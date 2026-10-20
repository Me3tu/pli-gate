#!/usr/bin/env python3
"""Acceptance-aware windowed-moment inversion for the external reference materials
(Steinberger disks and 7 MBq water cylinder), same test as acceptance_tm.py for Sim-1.

Re-runs the published clinical extractor (steinberger_extract_ap.clinical_run logic, same
sequential RNG, so the sharp-edge values must reproduce the published ones exactly) while
keeping per-shard histograms, then inverts the windowed mean
  (i) sharp edge (published),
  (ii) acceptance measured on the SAME dataset's randoms arm (upper-edge width and position),
  (iii) a-priori sigma_e = 359 ps.
"""
import json, os, sys
from pathlib import Path
import numpy as np, h5py
from scipy.optimize import curve_fit
from scipy.special import erfc

RES = Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[1])
sys.path.insert(0, str(RES / "2026-08-16_T6_nema"))
import nema_extract as X  # noqa
CFG = json.loads((RES / "2026-08-16_T6_nema/nema_voi_config.json").read_text())
Phi = lambda z: 0.5 * erfc(-z / np.sqrt(2.0))
bc = 0.5 * (X.TAU_BINS[:-1] + X.TAU_BINS[1:]); sel = (bc > 1.0) & (bc < 5.0); x = bc[sel]
tg = np.linspace(0.3, 12.0, 6000); t = np.linspace(1.0, 5.0, 8001)

def curve(sig=None, d=0.0):
    A = np.ones_like(t) if sig is None else Phi((5 + d - t) / sig) * Phi((t + 5 - d) / sig)
    E = np.exp(-t[None, :] / tg[:, None]) * A[None, :]
    return (E * t).sum(1) / E.sum(1)

def run(singles, phantom, rule):
    cfg = CFG[phantom]; vois, outline = cfg["vois"], cfg["phantom_outline"]
    f = h5py.File(singles, "r"); g = f["singles"]
    bounds = X.find_shard_boundaries(g["shard_id"][:])
    rng = np.random.default_rng(int(90.9 * 1000) or 1)
    HS = []
    for sh_id, s, e in bounds:
        out = X.clinical_per_shard(g["time_ns"][s:e], g["in_511"][s:e], g["in_prompt"][s:e],
                                   g["gx"][s:e], g["gy"][s:e], g["gz"][s:e], 0.0909, rng, vois, outline, matching=rule)
        Hs = np.zeros((len(vois) + 1, 2, len(bc)))
        if out is not None:
            for vi in range(-1, len(vois)):
                for a, arm in enumerate(out):
                    if arm[vi]:
                        Hs[vi + 1, a] = np.histogram(np.concatenate(arm[vi]), bins=X.TAU_BINS)[0]
        HS.append(Hs)
    f.close()
    return np.array(HS), ["bg"] + [v["name"] for v in vois], [np.nan] + [v.get("tau_cfg_ns", np.nan) for v in vois]

def invert(Hs, curves):
    hs = Hs[..., 0, :][..., sel] - Hs[..., 1, :][..., sel]
    mean = (hs * x).sum(-1) / hs.sum(-1)
    return {k: np.interp(mean, m, tg, left=np.nan, right=np.nan) for k, m in curves.items()}

out = {}
for tag, singles, phantom in [
        ("nema", Path(os.environ.get("PLI_RAW", "<missing>")) / "2026-08-16_T6_nema/singles_nema_all.h5", "nema_iq")]:
    if not Path(singles).exists() and not list(RES.glob(f"**/{Path(singles).name}")):
        raise SystemExit(f"{singles}: raw data not included in the release, see analysis/README.md")
    if not Path(singles).exists():
        cand = list(RES.glob(f"**/{Path(singles).name}"))
        singles = cand[0]
    for rule in ("closest", "all"):
        HS, names, cfgs = run(singles, phantom, rule)
        H = HS.sum(0)
        r = H[:, 1, :].sum(0); m = (bc > 2.5) & (bc < 7.5)
        p, c = curve_fit(lambda z, a, s, d: a * Phi((5 + d - z) / s), bc[m], r[m],
                         p0=[r[(bc > 1) & (bc < 4)].mean(), 0.3, 0.0], sigma=np.sqrt(np.clip(r[m], 1, None)))
        curves = {"sharp": curve(), "acc_measured": curve(p[1], p[2]), "acc359": curve(0.35897)}
        T0 = invert(H, curves)
        rng = np.random.default_rng(11); S = len(HS); flat = HS.reshape(S, -1); TB = []
        for _ in range(4):
            w = np.stack([np.bincount(rng.integers(0, S, S), minlength=S) for _ in range(250)]).astype(float)
            TB.append(invert((w @ flat).reshape((250,) + HS.shape[1:]), curves))
        res = dict(sigma_e_ps=1000 * p[1], sigma_e_err_ps=1000 * np.sqrt(c[1, 1]), edge_shift_ps=1000 * p[2], n_shards=S, vois={})
        print(f"\n{tag} {rule}: randoms-arm edge sigma_e {1000*p[1]:.1f} +- {1000*np.sqrt(c[1,1]):.1f} ps, shift {1000*p[2]:+.1f} ps, {S} shards")
        for iv, nm in enumerate(names):
            row = {}
            for k in curves:
                sd = float(np.nanstd(np.concatenate([tb[k][:, iv] for tb in TB]), ddof=1))
                row[k] = [float(T0[k][iv]), sd]
            res["vois"][nm] = dict(tau_cfg=cfgs[iv], **row)
            print(f"   {nm:14s} cfg {cfgs[iv]!s:>6}  " + "  ".join(f"{k} {row[k][0]:.4f}±{row[k][1]:.4f}" for k in curves))
        out[f"{tag}_{rule}"] = res
json.dump(out, open(RES / "2026-09-30_referee2/nema_acceptance.json", "w"), indent=1, default=float)
print("wrote nema_acceptance.json")
