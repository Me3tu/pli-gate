#!/usr/bin/env python3
"""B3: fit each Mercolli tube with BOTH estimators, acceptance-corrected.

Comparison axis: recovered o-Ps lifetime in ns against Mercolli's published
value for that sample. Two estimators on identical events: the production
truncated mean, and the three-component fit with the a priori acceptance.
Nothing is retuned per sample -- one acceptance, one model, four tubes.
"""
from __future__ import annotations
import json, os, sys
from pathlib import Path
import numpy as np
from scipy.optimize import minimize

RES = Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[1])
ACC = RES / "2026-08-23_acceptance"
sys.path.insert(0, str(ACC))
sys.path.insert(0, str(RES / "2026-08-22_B2_estimator"))
sys.path.insert(0, str(RES / "2026-08-12_steinberger_allprompts"))
import steinberger_extract_ap as X
import b2_estimator as E
import b2c_gate_acc as G
from acceptance import A_analytic

BC = 0.5 * (X.TAU_BINS[:-1] + X.TAU_BINS[1:]); WB = BC[1] - BC[0]
SEL = (BC > X.TAU_MIN_NS) & (BC < X.W_PROMPT_NS)
W = X.W_PROMPT_NS - X.TAU_MIN_NS
# Mercolli Table 2, published tau_oPs and BR3 per sample
PUB = {"T1_dry": (2.52, 0.211), "T2_wet": (2.37, 0.217),
       "T3_gelatin": (2.27, 0.222), "T4_water": (1.82, 0.268)}
FIT_LO, FIT_HI = -1.0, 4.5          # inside the flat-acceptance region


def trunc_mean(hs, hr):
    h = hs[SEL] - hr[SEL]; n = h.sum()
    if n <= 0: return float("nan")
    return X.truncation_correct(float((BC[SEL] * h).sum() / n) - X.TAU_MIN_NS, W)


def fit3(counts, edges):
    accb = np.array([np.trapezoid(A_analytic(np.linspace(a, b, 60), G.SIGMA_E),
                                  np.linspace(a, b, 60)) / (b - a)
                     for a, b in zip(edges[:-1], edges[1:])])

    def nll(th):
        if not (0.3 < th[0] < 8.0 and 0.01 < th[2] < 1.0 and abs(th[1]) < 2.0):
            return 1e12
        w = np.exp(np.array([th[3], th[4], 0.0])); I = w / w.sum()
        p = np.zeros(len(edges) - 1)
        for Ik, tk in zip(I, (E.TAU_PPS, E.TAU_DIR, th[0])):
            p += Ik * np.diff(E.emg_cdf(edges, tk, th[1], th[2]))
        p = np.clip(p * accb, 1e-15, None); p = p / p.sum()
        mu = counts.sum() * p
        return float(np.sum(mu - counts * np.log(mu))) - G.dirichlet_logpdf(I, E.ALPHA)

    best = None
    for st in (0.9, 1.4, 1.9, 2.4, 2.9):
        x0 = np.array([st, 0.0, 0.0909, np.log(0.15/0.23), np.log(0.62/0.23)])
        r = minimize(nll, x0, method="Nelder-Mead",
                     options=dict(maxiter=30000, xatol=1e-7, fatol=1e-7))
        if best is None or r.fun < best.fun: best = r
    t, f0 = float(best.x[0]), float(best.fun)
    if not (0.35 < t < 7.5): return float("nan"), float("nan"), None
    d = max(0.004, 0.004*t)
    pr = lambda tt: float(minimize(lambda q: nll(np.r_[tt, q]), best.x[1:],
                                   method="Nelder-Mead",
                                   options=dict(maxiter=8000, xatol=1e-7, fatol=1e-7)).fun)
    curv = (pr(t+d) - 2*f0 + pr(t-d)) / (d*d)
    return t, (float(1/np.sqrt(curv)) if curv > 0 else float("nan")), best.x


z = np.load("hists_mercolli.npz", allow_pickle=True)
names = [str(s) for s in z["voi_names"]]
M = (BC > FIT_LO) & (BC < FIT_HI)
edges = np.r_[BC[M] - 0.5*WB, BC[M][-1] + 0.5*WB]
out = {"fit_range_ns": [FIT_LO, FIT_HI], "sigma_e_ns": G.SIGMA_E, "tubes": {}}
print(f"fit range [{FIT_LO}, {FIT_HI}] ns; acceptance sigma_e = {1000*G.SIGMA_E:.1f} ps, "
      f"a priori, identical for all four\n")
print(f"{'tube':12s} {'rule':>8s} {'published':>10s} {'truncated':>10s} "
      f"{'3-comp':>18s} {'N':>10s}")
for tube, (pub, br3) in PUB.items():
    i = names.index(tube)
    for rule in ("closest", "all"):
        hs = z[f"{rule}_sig"][:, i, :].sum(0).astype(float)
        hr = z[f"{rule}_rnd"][:, i, :].sum(0).astype(float)
        tm = trunc_mean(hs, hr)
        sig = np.clip((hs - hr)[M], 0.0, None)
        t3, sd3, par = fit3(sig, edges)
        out["tubes"].setdefault(tube, {})[rule] = dict(
            published=pub, br3=br3, trunc=float(tm), tau3=float(t3), sd3=float(sd3),
            n=float(sig.sum()), spectrum=(hs - hr).tolist(),
            params=(None if par is None else [float(v) for v in par]))
        print(f"{tube:12s} {rule:>8s} {pub:10.2f} {tm:10.3f} "
              f"{t3:11.3f} +/- {sd3:.3f} {sig.sum():10,.0f}")
out["bin_centres_ns"] = BC.tolist()
json.dump(out, open("b3_fit.json", "w"))
print("\nwrote b3_fit.json")
