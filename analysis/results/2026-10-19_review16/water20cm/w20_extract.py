#!/usr/bin/env python3
"""Review round 17 (2026-10-19): clinical extraction of the 20 cm water cylinder (Steinberger et al 2024 size).

Same extractor and readouts as 2026-09-30_referee2/steinberger_acceptance.py, which produced the Table 5 water row
of the 20 mm cylinder (VOI_15x25, sharp readout, all-prompts): steinberger_extract_ap.clinical_per_shard with the
same sequential RNG seed, 50 ps histograms, fit window 1 to 5 ns, windowed-moment inversion
  (i) sharp edge (production), (ii) acceptance measured on the same dataset's randoms arm, (iii) a priori 359 ps,
and a shard bootstrap (4 x 250 resamples, seed 11). Only the geometry differs: VOIs and phantom outline scaled from
the 20 mm design to the published 20 cm cylinder (VOI diameters 5, 10, 15 cm; lengths unchanged 15, 20, 25 cm).
Usage: w20_extract.py <singles.h5> [<singles.h5> ...]   (one file per batch; shards are taken per file)
Writes w20_extract.json next to this script."""
import json, os, sys
from pathlib import Path
import numpy as np, h5py
from scipy.optimize import curve_fit
from scipy.special import erfc

HERE = Path(__file__).resolve().parent
RES = Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[2])
sys.path.insert(0, str(RES / "2026-08-12_steinberger_allprompts"))
import steinberger_extract_ap as X  # noqa

OUTLINE = {"type": "cylinder", "R_mm": 102.0, "halfz_mm": 160.0, "buffer_mm": 2.0, "lor_tol_mm": 10.0,
           "_note": "20 cm diameter x 30 cm water cylinder + small buffer (same buffers as the 20 mm design)"}
VOIS = [dict(name=n, shape="cylinder", cx_mm=0.0, cy_mm=0.0, cz_mm=0.0, r_mm=r, half_z_mm=hz, tau_cfg_ns=1.84,
             tau_steinberger_ns=1.815) for n, r, hz in (("VOI_5x15", 25.0, 75.0), ("VOI_10x20", 50.0, 100.0),
                                                       ("VOI_15x25", 75.0, 125.0))]
Phi = lambda z: 0.5 * erfc(-z / np.sqrt(2.0))
bc = 0.5 * (X.TAU_BINS[:-1] + X.TAU_BINS[1:]); sel = (bc > 1.0) & (bc < 5.0); x = bc[sel]
tg = np.linspace(0.3, 12.0, 6000); t = np.linspace(1.0, 5.0, 8001)


def curve(sig=None, d=0.0):
    A = np.ones_like(t) if sig is None else Phi((5 + d - t) / sig) * Phi((t + 5 - d) / sig)
    E = np.exp(-t[None, :] / tg[:, None]) * A[None, :]
    return (E * t).sum(1) / E.sum(1)


def run(files, rule):
    rng = np.random.default_rng(int(90.9 * 1000) or 1)
    HS = []; ntrip = []
    for fn in files:
        f = h5py.File(fn, "r"); g = f["singles"]
        for sh_id, s, e in X.find_shard_boundaries(g["shard_id"][:]):
            out = X.clinical_per_shard(g["time_ns"][s:e], g["in_511"][s:e], g["in_prompt"][s:e],
                                       g["gx"][s:e], g["gy"][s:e], g["gz"][s:e], 0.0909, rng, VOIS, OUTLINE,
                                       matching=rule)
            Hs = np.zeros((len(VOIS) + 1, 2, len(bc)))
            if out is not None:
                for vi in range(-1, len(VOIS)):
                    for a, arm in enumerate(out):
                        if arm[vi]:
                            Hs[vi + 1, a] = np.histogram(np.concatenate(arm[vi]), bins=X.TAU_BINS)[0]
            HS.append(Hs)
        f.close()
    return np.array(HS), ["bg"] + [v["name"] for v in VOIS]


def invert(Hs, curves):
    hs = Hs[..., 0, :][..., sel] - Hs[..., 1, :][..., sel]
    mean = (hs * x).sum(-1) / hs.sum(-1)
    return {k: np.interp(mean, m, tg, left=np.nan, right=np.nan) for k, m in curves.items()}


files = sys.argv[1:]
out = {"files": files, "outline": OUTLINE, "vois": VOIS}
for rule in ("closest", "all"):
    HS, names = run(files, rule)
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
    res = dict(sigma_e_ps=1000 * p[1], sigma_e_err_ps=1000 * np.sqrt(c[1, 1]), edge_shift_ps=1000 * p[2],
               n_shards=S, vois={})
    print(f"\nwater20 {rule}: randoms-arm edge sigma_e {1000*p[1]:.1f} +- {1000*np.sqrt(c[1,1]):.1f} ps, "
          f"shift {1000*p[2]:+.1f} ps, {S} shards")
    for iv, nm in enumerate(names):
        n_tot = float(H[iv, 0, sel].sum()); n_rnd = float(H[iv, 1, sel].sum())
        row = dict(n_fitwin_total=n_tot, n_fitwin_random=n_rnd)
        for k in curves:
            sd = float(np.nanstd(np.concatenate([tb[k][:, iv] for tb in TB]), ddof=1))
            row[k] = [float(T0[k][iv]), sd]
        res["vois"][nm] = row
        print(f"   {nm:10s} N(1-5ns) {n_tot:9.0f} rnd {n_rnd:8.0f} ({100*n_rnd/max(n_tot-n_rnd,1):5.1f}% of signal)  "
              + "  ".join(f"{k} {row[k][0]:.4f}+-{row[k][1]:.4f}" for k in curves))
    out[rule] = res
json.dump(out, open(HERE / "w20_extract.json", "w"), indent=1, default=float)
print("wrote", HERE / "w20_extract.json")
