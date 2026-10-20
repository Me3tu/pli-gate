#!/usr/bin/env python3
"""B2 gate, re-run THROUGH the acceptance -- the failure mode it now must represent.

Synthetic spectra are generated WITH A(tau) applied and fitted by the
acceptance-corrected model. Same four conditions as the original gate.

A(tau) uses sigma_e = 312.4 ps: the quadrature sum of the a priori geometric
path-length spread (305.7 ps, from the scanner geometry alone) and the timing
term (64.3 ps, from the production CRT). Neither is fitted here.
"""
from __future__ import annotations
import json, os, sys
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from scipy.special import gammaln

sys.path.insert(0, ".")
sys.path.insert(0, str(Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[1]) / "2026-08-22_B2_estimator"))
import b2_estimator as E
from acceptance import A_analytic

SIGMA_E = 0.3124                  # ns, a priori; see acceptance_water.py
EDGES = E.EDGES
CEN = E.CEN
ACC_BIN = None


def acc_bins(edges, sigma_e=SIGMA_E, n=60):
    out = np.empty(len(edges) - 1)
    for i, (a, b) in enumerate(zip(edges[:-1], edges[1:])):
        x = np.linspace(a, b, n)
        out[i] = np.trapezoid(A_analytic(x, sigma_e), x) / (b - a)
    return out


ACC_BIN = acc_bins(EDGES)


def dirichlet_logpdf(x, a):
    x = np.clip(np.asarray(x, float), 1e-12, None); x = x / x.sum()
    return float(np.sum((a - 1.0) * np.log(x)) + gammaln(a.sum()) - np.sum(gammaln(a)))


def model_counts(theta, N, acc=True):
    tau3, t0, sig, w1, w2 = theta
    w = np.exp(np.array([w1, w2, 0.0])); I = w / w.sum()
    p = np.zeros(len(CEN))
    for Ik, tk in zip(I, (E.TAU_PPS, E.TAU_DIR, tau3)):
        p += Ik * np.diff(E.emg_cdf(EDGES, tk, t0, sig))
    if acc:
        p = p * ACC_BIN
    p = np.clip(p, 1e-15, None); p = p / p.sum()
    return N * p, I


def nll(theta, counts, use_prior=True, acc=True):
    if not (0.3 < theta[0] < 8.0 and 0.01 < theta[2] < 1.0 and abs(theta[1]) < 2.0):
        return 1e12
    mu, I = model_counts(theta, counts.sum(), acc)
    v = float(np.sum(mu - counts * np.log(mu)))
    if use_prior:
        v -= dirichlet_logpdf(I, E.ALPHA)
    return v


def fit(counts, acc=True, seed_tau=2.0):
    best = None
    for st in (seed_tau, 1.6, 2.6):
        x0 = np.array([st, 0.0, 0.0909, np.log(0.15/0.23), np.log(0.62/0.23)])
        r = minimize(nll, x0, args=(counts, True, acc), method="Nelder-Mead",
                     options=dict(maxiter=25000, xatol=1e-7, fatol=1e-7))
        if best is None or r.fun < best.fun:
            best = r
    t, f0 = float(best.x[0]), float(best.fun)
    d = max(0.004, 0.004*t)
    pr = lambda tt: float(minimize(lambda q: nll(np.r_[tt, q], counts, True, acc),
                                   best.x[1:], method="Nelder-Mead",
                                   options=dict(maxiter=8000, xatol=1e-7, fatol=1e-7)).fun)
    curv = (pr(t+d) - 2*f0 + pr(t-d)) / (d*d)
    sd = float(1/np.sqrt(curv)) if curv > 0 else float("nan")
    _, I = model_counts(best.x, counts.sum(), acc)
    return dict(tau=t, sd=sd, sigma=float(best.x[2]), I=[float(v) for v in I])


def synth_acc(tau3, N, rng):
    I = E.ALPHA / E.ALPHA.sum()
    th = np.r_[tau3, 0.0, 0.0909, np.log(I[0]/I[2]), np.log(I[1]/I[2])]
    mu, _ = model_counts(th, N, acc=True)      # generated WITH the acceptance
    return rng.poisson(mu).astype(float), I


if __name__ == "__main__":
    rng = np.random.default_rng(20260823)
    TRUTH = [1.50, 1.80, 2.04, 2.40, 2.80]
    out = {"sigma_e_ns": SIGMA_E, "runs": []}
    print(f"acceptance sigma_e = {1000*SIGMA_E:.1f} ps, a priori. Spectra are GENERATED")
    print(f"with the acceptance, then fitted with it.\n")
    for N in (10**6, 10**5):
        print(f"N = {N:,}")
        print(f"{'tau_true':>9s} {'tau_hat':>9s} {'sd':>8s} {'err %':>8s} "
              f"{'err/sd':>8s} {'I_oPs':>7s}")
        for tt in TRUTH:
            counts, Itrue = synth_acc(tt, N, rng)
            r = fit(counts, acc=True, seed_tau=tt*1.05)
            err = r["tau"] - tt
            out["runs"].append(dict(N=N, tau_true=tt, **r, err_ns=float(err),
                                    err_pct=float(100*err/tt),
                                    err_over_sd=float(abs(err)/r["sd"]),
                                    dI_oPs=float(r["I"][2]-Itrue[2])))
            print(f"{tt:9.2f} {r['tau']:9.4f} {r['sd']:8.4f} {100*err/tt:8.2f} "
                  f"{abs(err)/r['sd']:8.2f} {r['I'][2]:7.3f}")
        print()
    r6 = [x for x in out["runs"] if x["N"] == 10**6]
    r5 = [x for x in out["runs"] if x["N"] == 10**5]
    C = [("within 1 % at 1e6", max(abs(x["err_pct"]) for x in r6), 1.0),
         ("within 3 sigma at 1e6", max(x["err_over_sd"] for x in r6), 3.0),
         ("o-Ps share within 0.02", max(abs(x["dI_oPs"]) for x in r6), 0.02),
         ("within 3 sigma at 1e5", max(x["err_over_sd"] for x in r5), 3.0)]
    print("gate, re-run through the acceptance")
    allok = True
    for lab, worst, bar in C:
        ok = worst <= bar; allok &= ok
        print(f"  {lab:26s} worst {worst:8.3f}  bar {bar:5.2f}  {'PASS' if ok else 'FAIL'}")
    out["conditions"] = [dict(name=l, worst=float(w), bar=float(b), ok=bool(w <= b))
                         for l, w, b in C]
    out["all_pass"] = bool(allok)
    json.dump(out, open("b2c_gate_acc.json", "w"), indent=2)
    print(f"\nVERDICT: {'PASS' if allok else 'FAIL'}   wrote b2c_gate_acc.json")
