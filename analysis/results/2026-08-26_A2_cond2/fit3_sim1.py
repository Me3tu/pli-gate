#!/usr/bin/env python3
"""Three-component + acceptance fits for the eight Sim-1 tubes, per matching
rule, on the existing C2 (388 ps) histograms. No simulation, fitting only.

The estimator is C1's, imported unchanged in structure: three EMG components
with the para-Ps and direct lifetimes FIXED (b2_estimator.TAU_PPS = 0.125 ns,
TAU_DIR = 0.388 ns), a Dirichlet prior on the intensities, and the analytic
acceptance folded in per bin.

ONE thing changes from C1, and it is the point of Condition 2: the acceptance
width. C1 used G.SIGMA_E = 312.4 ps, which is the MERCOLLI/water cylinder's.
Sim-1's own geometry gives 359.0 ps (sigma_geo_sim1.json). Using the shipped
value here would be the same reuse error the acceptance memo records.

Errors are PAIRED: one shard resample, both rules refitted from exactly those
shards, so the rule difference carries the correlation.
"""
from __future__ import annotations
import os, sys, json, time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize

RES = Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[1])
for p in ("2026-08-12_steinberger_allprompts", "2026-08-22_B2_estimator",
          "2026-08-23_acceptance"):
    sys.path.insert(0, str(RES / p))
import steinberger_extract_ap as X   # noqa
import b2_estimator as E             # noqa
import b2c_gate_acc as G             # noqa
from acceptance import A_analytic    # noqa

HERE = RES / "2026-08-26_A2_cond2"
SIGMA_E = json.load(open(HERE / "sigma_geo_sim1.json"))["sigma_e_ns"]   # Sim-1's own
BC = 0.5 * (X.TAU_BINS[:-1] + X.TAU_BINS[1:]); WB = BC[1] - BC[0]
M = (BC > -1.0) & (BC < 4.5)
EDGES = np.r_[BC[M] - 0.5 * WB, BC[M][-1] + 0.5 * WB]
ACCB = np.array([np.trapezoid(A_analytic(np.linspace(a, b, 60), SIGMA_E),
                              np.linspace(a, b, 60)) / (b - a)
                 for a, b in zip(EDGES[:-1], EDGES[1:])])
STARTS = (0.9, 1.4, 1.9, 2.4)


def _nll(th, counts):
    if not (0.3 < th[0] < 8.0 and 0.01 < th[2] < 1.0 and abs(th[1]) < 2.0):
        return 1e12
    w = np.exp(np.array([th[3], th[4], 0.0])); I = w / w.sum()
    p = np.zeros(len(EDGES) - 1)
    for Ik, tk in zip(I, (E.TAU_PPS, E.TAU_DIR, th[0])):
        p += Ik * np.diff(E.emg_cdf(EDGES, tk, th[1], th[2]))
    p = np.clip(p * ACCB, 1e-15, None); p /= p.sum()
    mu = counts.sum() * p
    return float(np.sum(mu - counts * np.log(mu))) - G.dirichlet_logpdf(I, E.ALPHA)


def fit3(counts, x0=None):
    """Full multi-start for point estimates, single warm start for replicates."""
    inits = ([np.array([x0, 0.0, .0909, np.log(.15 / .23), np.log(.62 / .23)])]
             if x0 is not None else
             [np.array([s, 0.0, .0909, np.log(.15 / .23), np.log(.62 / .23)])
              for s in STARTS])
    best = None
    for st in inits:
        r = minimize(_nll, st, args=(counts,), method="Nelder-Mead",
                     options=dict(maxiter=30000, xatol=1e-7, fatol=1e-7))
        if best is None or r.fun < best.fun:
            best = r
    t = float(best.x[0])
    return t if 0.35 < t < 7.5 else float("nan")


def spectrum(sig, rnd):
    return np.clip((sig - rnd)[M], 0, None)


if __name__ == "__main__":
    z = np.load(RES / "2026-08-26_B2_paired/hists_sim1_388.npz", allow_pickle=True)
    names = [str(s) for s in z["voi_names"]]
    i = names.index("T1_water_shibuya2020_n2")
    c = spectrum(z["closest_sig"][:, i, :].sum(0), z["closest_rnd"][:, i, :].sum(0))
    t0 = time.time(); v = fit3(c); t1 = time.time()
    print(f"  acceptance sigma_e used : {1000*SIGMA_E:.1f} ps  (Sim-1's own; "
          f"C1 used {1000*G.SIGMA_E:.1f} ps)")
    print(f"  T1 closest, 4-start fit : tau3 = {v:.4f} ns   in {t1-t0:.1f} s")
    t0 = time.time(); v2 = fit3(c, x0=v); t1 = time.time()
    print(f"  T1 closest, warm start  : tau3 = {v2:.4f} ns   in {t1-t0:.1f} s "
          f"(delta {1000*(v2-v):+.2f} ps)")
