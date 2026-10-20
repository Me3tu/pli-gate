#!/usr/bin/env python3
"""Review round 17 (2026-10-19, second report Sec. 4.4): size of the scanner-radius effect.

The model's inner crystal-face radius is 392.83 mm (MUSIRE template) against 410 mm for the scanner. This script
(1) computes the geometric gain of the model over a 410 mm scanner per photon, per pair and per triple for the
    eight-tube source geometry: axial acceptance (|z| <= 530 mm on the mid-crystal cylinder) by Monte Carlo, times
    the transaxial sector coverage of a 38-sided polygon with 64.2 mm sectors (side = 2 r tan(pi/38)),
(2) recomputes the a priori acceptance width sigma_e of the eight-tube geometry with the frozen method of
    2026-08-26_A2_cond2/sigma_geo_sim1.py at both radii (the 392.83 mm value must reproduce 358.97 ps),
(3) re-reads the corrected-yield zero-delay eight-tube all-prompts arm (C0, 128 shards, 303 ps) with the
    acceptance-aware windowed moment at both widths, as acceptance_tm_c0.py does, to show what the width alone
    changes in the 1.9% mean absolute residual. The events stay those of the 392.83 mm model.
Writes radius410.json."""
import json, os, sys
from pathlib import Path
import numpy as np
from scipy.special import erfc

RES = Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[1])
C = 299.792458
DEPTH = 414.2571 - 392.8322
HZ = 530.0
SIG_T_PAIR = 0.0909 / np.sqrt(2.0)
RING_R, TUBE_R, TUBE_HZ = 80.0, 15.0, 25.0
ANG = np.arange(8) * (np.pi / 4.0)


def hit(p, u, rdet):
    a = u[:, 0] ** 2 + u[:, 1] ** 2
    b = 2 * (p[:, 0] * u[:, 0] + p[:, 1] * u[:, 1])
    c = p[:, 0] ** 2 + p[:, 1] ** 2 - rdet ** 2
    t = (-b + np.sqrt(np.maximum(b * b - 4 * a * c, 0.0))) / (2 * a)
    return t, np.abs(p[:, 2] + t * u[:, 2]) <= HZ


def iso(rng, m):
    c = rng.uniform(-1, 1, m); s = np.sqrt(1 - c * c); ph = rng.uniform(0, 2 * np.pi, m)
    return np.c_[s * np.cos(ph), s * np.sin(ph), c]


def tube_points(rng, k, n):
    cx, cy = RING_R * np.cos(ANG[k]), RING_R * np.sin(ANG[k])
    r = TUBE_R * np.sqrt(rng.random(n)); th = rng.uniform(0, 2 * np.pi, n)
    return np.c_[cx + r * np.cos(th), cy + r * np.sin(th), rng.uniform(-TUBE_HZ, TUBE_HZ, n)]


def geometry(r_in, n=2_000_000, seed=20260826):
    rdet = r_in + 0.5 * DEPTH
    sg, ax1, ax2, ax3 = [], [], [], []
    for k in range(8):
        rng = np.random.default_rng(seed + k)
        p = tube_points(rng, k, n); u = iso(rng, n)
        d1, o1 = hit(p, u, rdet); d2, o2 = hit(p, -u, rdet)
        v = iso(rng, n); dp, op = hit(p, v, rdet)
        ok = o1 & o2 & op
        sg.append(float(((dp[ok] - 0.5 * (d1[ok] + d2[ok])) / C).std(ddof=1)))
        ax1.append(o1.mean()); ax2.append((o1 & o2).mean()); ax3.append(ok.mean())
    side = 2 * r_in * np.tan(np.pi / 38); cov = 64.2 / side
    sgm = float(np.mean(sg))
    return dict(r_in=r_in, r_det=rdet, sector_side_mm=side, gap_mm=side - 64.2, coverage=cov,
                axial_single=float(np.mean(ax1)), axial_pair=float(np.mean(ax2)), axial_triple=float(np.mean(ax3)),
                photon=float(np.mean(ax1)) * cov, pair=float(np.mean(ax2)) * cov ** 2,
                triple=float(np.mean(ax3)) * cov ** 3, sigma_geo_ps=1000 * sgm,
                sigma_e_ps=1000 * float(np.hypot(sgm, SIG_T_PAIR)))


G = {"model": geometry(392.8322), "scanner": geometry(410.0)}
gain = {k: G["model"][k] / G["scanner"][k] - 1 for k in ("photon", "pair", "triple")}
print(json.dumps(G, indent=1)); print("model gain over 410 mm:", {k: f"{100*v:+.1f}%" for k, v in gain.items()})

# (3) C0 acceptance-aware readout at both widths
D = RES / "2026-09-30_feedback3/c0_sim1"
Z = np.load(D / "photon_conventions_c0/photon_conventions.npz")
H = Z["hist"][:, 0].astype(float); BINS, NAMES, TCFG = Z["bins"], list(Z["voi_names"]), Z["tau_cfg"]
bc = 0.5 * (BINS[:-1] + BINS[1:]); sel = (bc > 1.0) & (bc < 5.0); x = bc[sel]
Phi = lambda z: 0.5 * erfc(-z / np.sqrt(2.0)); tg = np.linspace(0.3, 12.0, 6000); t = np.linspace(1.0, 5.0, 8001)


def curve(se):
    A = Phi((5 - t) / se) * Phi((t + 5) / se); E = np.exp(-t[None, :] / tg[:, None]) * A[None, :]
    return (E * t).sum(1) / E.sum(1)


CUR = {"359": curve(0.35897), "model_recomputed": curve(G["model"]["sigma_e_ps"] / 1000),
       "410mm": curve(G["scanner"]["sigma_e_ps"] / 1000)}
TI = [i for i, nm in enumerate(NAMES) if nm != "bg"]


def readout(Hs):
    hs = Hs[..., 0, :][..., sel] - Hs[..., 1, :][..., sel]; m = (hs * x).sum(-1) / hs.sum(-1)
    return {k: np.interp(m, c, tg) for k, c in CUR.items()}


T0 = readout(H.sum(0)); rng = np.random.default_rng(4); S = H.shape[0]; flat = H.reshape(S, -1); TB = []
for _ in range(4):
    w = np.stack([np.bincount(rng.integers(0, S, S), minlength=S) for _ in range(250)]).astype(float)
    TB.append(readout((w @ flat).reshape((250,) + H.shape[1:])))
out = dict(geometry=G, gain_model_over_410=gain, c0_all_prompts={})
ir = 1  # all prompts
for k in CUR:
    b = 100 * (T0[k][ir, TI] - TCFG[TI]) / TCFG[TI]
    bb = np.concatenate([100 * (tb[k][:, ir][:, TI] - TCFG[TI]) / TCFG[TI] for tb in TB])
    d = np.concatenate([100 * (tb[k][:, ir][:, TI] - tb["359"][:, ir][:, TI]) / TCFG[TI] for tb in TB]).mean(1)
    out["c0_all_prompts"][k] = dict(mean_abs=float(np.abs(b).mean()), mean_abs_sd=float(np.abs(bb).mean(1).std(ddof=1)),
                                    mean_signed=float(b.mean()), mean_signed_sd=float(bb.mean(1).std(ddof=1)),
                                    change_vs_359_pp=float(np.mean(b - 100 * (T0["359"][ir, TI] - TCFG[TI]) / TCFG[TI])),
                                    change_sd=float(d.std(ddof=1)))
    print(k, {kk: round(v, 3) for kk, v in out["c0_all_prompts"][k].items()})
json.dump(out, open(Path(__file__).resolve().parent / "radius410.json", "w"), indent=1)
