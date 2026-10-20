#!/usr/bin/env python3
"""Referee round 2, MC11 and statistical audit: synthetic validation of the lifetime
estimators under the exact randoms-subtraction process.

Generator (known truth). Delay spectrum = sum_k I_k Exp(tau_k) (x) N(0, sigma_d), times the
window acceptance A(tau; sigma_e) of the production extractor. tau_pPs = 0.125 ns and
tau_dir = 0.388 ns as in the component estimator. sigma_d = sqrt(1.5)*90.9 ps (pair-mean and
prompt single), sigma_e = 358.97 ps (Sim-1 geometry). Randoms: flat in delay, prompt arm
n_t ~ Poisson(s + r), delayed arm n_d ~ Poisson(r), independent. Randoms fraction f_R is
defined inside the component-fit range (-1, 4.5) ns. Each replicate is built from 32
independent shards (as the production data), so the published shard bootstrap can be run.

Estimators
  PUB   published component estimator: clip(n_t - n_d, 0), multinomial Poisson objective
        with Dirichlet prior (alpha = 0.75, 3.1, 1.15), five free parameters, four starts.
  NOCLIP  same objective without clipping negative bins.
  JOINT   count-level joint likelihood of both arms, per-bin randoms mean r_b profiled
        analytically (no flatness assumption, no clipping), same prior, six parameters.
  TM    windowed-moment (truncated-mean) estimator on (1, 5) ns with the analytic
        truncation inversion and the production delta-method sigma.
Intervals: shard bootstrap (B = 100, warm-started single-start refits) for PUB and JOINT,
percentile 68/95 %, and delta method for TM.
"""
from __future__ import annotations
import json, os, sys, time
from pathlib import Path
from multiprocessing import Pool
import numpy as np
from scipy.optimize import minimize
from scipy.special import erfc, erfcx, gammaln

OUT = Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[1]) / "2026-09-30_referee2/estimator_synth"
TAU_PPS, TAU_DIR = 0.125, 0.388
ALPHA = np.array([0.75, 3.1, 1.15])
SIG_D = np.sqrt(1.5) * 0.0909
SIG_E = 0.35896857566233875
W = 5.0
NSH = 32

def emg_cdf(t, tau, t0, sig):
    t = np.asarray(t, float); z = (t - t0) / sig; a = sig / tau; r2 = np.sqrt(2.0)
    Phi = 0.5 * erfc(-z / r2); u = (a - z) / r2
    lo = np.where(u >= 0.0, u, 0.0)
    term_lo = 0.5 * np.exp(-0.5 * z * z) * erfcx(lo)
    ex = np.clip(0.5 * a * a - a * z, -700.0, 700.0)
    term_hi = 0.5 * np.exp(ex) * erfc(np.where(u < 0.0, u, 0.0))
    return np.clip(Phi - np.where(u >= 0.0, term_lo, term_hi), 0.0, 1.0)

def Phi(x):
    return 0.5 * erfc(-x / np.sqrt(2.0))

def acc(tau):
    return Phi((W - tau) / SIG_E) * Phi((tau + W) / SIG_E)

def dir_logpdf(x, a):
    x = np.clip(x, 1e-12, None)
    return float(((a - 1) * np.log(x)).sum() + gammaln(a.sum()) - gammaln(a).sum())

class Grid:
    """Bins over (-1, 5) ns. Component fits use (-1, 4.5), the moment estimator (1, 5)."""
    def __init__(self, bw):
        self.bw = bw
        self.edges = np.round(np.arange(-1.0, 5.0 + 1e-9, bw), 6)
        self.bc = 0.5 * (self.edges[:-1] + self.edges[1:])
        self.M = self.bc < 4.5
        self.T = (self.bc > 1.0) & (self.bc < 5.0)
        e = self.edges[: int(self.M.sum()) + 1]
        self.eM = e
        self.accb = np.array([np.trapezoid(acc(np.linspace(a, b, 40)), np.linspace(a, b, 40)) / (b - a)
                              for a, b in zip(self.edges[:-1], self.edges[1:])])

    def shape(self, th, edges, accb):
        w = np.exp(np.array([th[3], th[4], 0.0])); I = w / w.sum()
        p = np.zeros(len(edges) - 1)
        for Ik, tk in zip(I, (TAU_PPS, TAU_DIR, th[0])):
            p += Ik * np.diff(emg_cdf(edges, tk, th[1], th[2]))
        return np.clip(p * accb, 1e-15, None), I

def truth_counts(g, tau3, I, Ns):
    th = [tau3, 0.0, SIG_D, np.log(I[0] / I[2]), np.log(I[1] / I[2])]
    p, _ = g.shape(th, g.edges, g.accb)
    pM = p[g.M].sum()
    return Ns * p / pM          # Ns signal triples expected inside (-1, 4.5)

# ---------------- estimators ----------------
def _bad(th):
    return not (0.3 < th[0] < 8.0 and 0.01 < th[2] < 1.0 and abs(th[1]) < 2.0)

def nll_pub(th, c, g):
    if _bad(th): return 1e12
    p, I = g.shape(th, g.eM, g.accb[g.M]); p = p / p.sum()
    mu = c.sum() * p
    return float(np.sum(mu - c * np.log(mu))) - dir_logpdf(I, ALPHA)

def nll_joint(th, nt, nd, g):
    if _bad(th) or not (-5 < th[5] < 25): return 1e12
    p, I = g.shape(th, g.eM, g.accb[g.M]); p = p / p.sum()
    s = np.exp(th[5]) * p
    b = 2 * s - nt - nd
    r = (-b + np.sqrt(b * b + 8 * nd * s)) / 4.0
    r = np.clip(r, 1e-12, None); m = s + r
    ll = (m - nt * np.log(m)) + (r - nd * np.log(r))
    return float(ll.sum()) - dir_logpdf(I, ALPHA)

STARTS = (0.9, 1.4, 1.9, 2.4)
def x0(t, extra=None):
    v = [t, 0.0, 0.0909, np.log(.15 / .23), np.log(.62 / .23)]
    return np.array(v + ([] if extra is None else [extra]))

def fit(obj, args, starts, extra=None):
    best = None
    for s in starts:
        st = s if isinstance(s, np.ndarray) else x0(s, extra)
        r = minimize(obj, st, args=args, method="Nelder-Mead",
                     options=dict(maxiter=30000 if len(starts) > 1 else 8000, xatol=1e-7, fatol=1e-7))
        if best is None or r.fun < best.fun: best = r
    t = float(best.x[0])
    return (t if 0.35 < t < 7.5 else np.nan), best.x

lam_g = np.geomspace(0.02, 50, 200001); u_ = 4.0 / lam_g
obs_g = lam_g - 4.0 * np.exp(-u_) / (-np.expm1(-u_))
def tm(nt, nd, g):
    ht, hr = nt[g.T], nd[g.T]; hs = ht - hr; ns = hs.sum()
    if ns <= 0: return np.nan, np.nan
    tobs = (g.bc[g.T] * hs).sum() / ns - 1.0
    lam = float(np.interp(tobs, obs_g, lam_g, left=np.nan, right=np.nan)) if tobs > 0 else np.nan
    sig_naive = abs(tobs) * np.sqrt(ht.sum() + hr.sum()) / max(ns, 1)
    lp = float(np.interp(tobs + 1e-3, obs_g, lam_g))
    return lam, abs((lp - lam) / 1e-3) * sig_naive

def one_rep(args):
    cond, rep, do_boot = args
    g = Grid(cond["bw"])
    rng = np.random.default_rng([20260930, cond["id"], rep])
    mu_s = truth_counts(g, cond["tau3"], cond["I"], cond["Ns"])
    fR = cond["fR"]
    r_per_bin = fR / (1 - fR) * cond["Ns"] / g.M.sum() if fR > 0 else 0.0
    sh_t = rng.poisson(np.broadcast_to((mu_s + r_per_bin) / NSH, (NSH, len(mu_s)))).astype(float)
    sh_d = rng.poisson(np.broadcast_to(np.full(len(mu_s), r_per_bin) / NSH, (NSH, len(mu_s)))).astype(float)
    nt, nd = sh_t.sum(0), sh_d.sum(0)
    res = dict(cond=cond["id"], rep=rep)
    diff = (nt - nd)[g.M]
    res["neg_frac"] = float((diff < 0).mean())
    ident = g.M.sum() > 6
    if cond.get("component", True) and ident:
        t_pub, x_pub = fit(nll_pub, (np.clip(diff, 0, None), g), STARTS)
        t_nc, _ = fit(nll_pub, (diff, g), STARTS)
        ext = np.log(max(diff.sum(), 10.0))
        t_j, x_j = fit(nll_joint, (nt[g.M], nd[g.M], g), STARTS, extra=ext)
        res.update(pub=t_pub, noclip=t_nc, joint=t_j)
        if do_boot and np.isfinite(t_pub) and np.isfinite(t_j):
            bp, bj = [], []
            for b in range(100):
                k = rng.integers(0, NSH, NSH)
                bt, bd = sh_t[k].sum(0), sh_d[k].sum(0)
                bp.append(fit(nll_pub, (np.clip((bt - bd)[g.M], 0, None), g), [x_pub])[0])
                bj.append(fit(nll_joint, (bt[g.M], bd[g.M], g), [x_j])[0])
            bp, bj = np.array(bp), np.array(bj)
            res["pub_boot"] = [float(np.nanstd(bp, ddof=1))] + [float(np.nanpercentile(bp, q)) for q in (2.5, 16, 84, 97.5)]
            res["joint_boot"] = [float(np.nanstd(bj, ddof=1))] + [float(np.nanpercentile(bj, q)) for q in (2.5, 16, 84, 97.5)]
    elif cond.get("component", True):
        res.update(pub=np.nan, noclip=np.nan, joint=np.nan, nonident=True)
    # moment estimator on a pure o-Ps + flat randoms spectrum (the zero-delay production model)
    if cond.get("tm", True):
        p1 = np.diff(emg_cdf(g.edges, cond["tau3"], 0.0, SIG_D)) * g.accb
        mu1 = cond["Ns"] * p1 / p1[g.T].sum()
        r1 = fR / (1 - fR) * cond["Ns"] / g.T.sum() if fR > 0 else 0.0
        a_t = rng.poisson(mu1 + r1).astype(float); a_d = rng.poisson(np.full(len(mu1), r1)).astype(float)
        res["tm"], res["tm_sig"] = tm(a_t, a_d, g)
    return res

def build_conditions():
    C = []; i = 0
    base = dict(tau3=2.0, I=[0.15, 0.62, 0.23])
    for fR in (0.0, 0.1, 0.2, 0.3, 0.5):
        for Ns in (2e4, 1e5, 5e5):
            C.append(dict(id=i, fR=fR, Ns=Ns, bw=0.05, boot=True, **base)); i += 1
    for bw in (0.1, 0.2, 0.5, 1.0, 2.0):
        for fR in (0.0, 0.2, 0.5):
            C.append(dict(id=i, fR=fR, Ns=1e5, bw=bw, boot=False, **base)); i += 1
    for t in (1.5, 2.5, 3.0):
        C.append(dict(id=i, fR=0.2, Ns=1e5, bw=0.05, boot=True, tau3=t, I=[0.15, 0.62, 0.23])); i += 1
    C.append(dict(id=i, fR=0.2, Ns=1e5, bw=0.05, boot=True, tau3=2.0, I=[0.10, 0.45, 0.45], note="prior-mismatched intensities")); i += 1
    return C

if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    NREP = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    conds = build_conditions()
    json.dump(conds, open(OUT / "conditions.json", "w"), indent=1)
    jobs = [(c, r, c["boot"]) for c in conds for r in range(NREP)]
    t0 = time.time(); rows = []
    with Pool(30) as pool:
        for i, res in enumerate(pool.imap_unordered(one_rep, jobs, chunksize=1)):
            rows.append(res)
            if i % 500 == 0:
                print(f"  {i}/{len(jobs)}  {time.time()-t0:.0f} s", flush=True)
                json.dump(rows, open(OUT / "reps_partial.json", "w"), default=float)
    json.dump(rows, open(OUT / "reps.json", "w"), default=float)
    print(f"done {time.time()-t0:.0f} s", flush=True)
