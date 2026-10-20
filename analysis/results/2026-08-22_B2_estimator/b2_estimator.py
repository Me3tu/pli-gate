#!/usr/bin/env python3
"""B2 - Mercolli's three-component estimator, and the synthetic-recovery gate.

Comparison axis: recovered o-Ps lifetime against the KNOWN value used to
generate the spectrum, in ns. Tolerances are fixed in REGISTERED_B2.md and are
not adjusted here.
"""
from __future__ import annotations
import json
import numpy as np
from scipy.special import erfcx, erfc
from scipy.optimize import minimize
from scipy.stats import dirichlet

TAU_PPS, TAU_DIR = 0.125, 0.388          # ns, FIXED per the block
ALPHA = np.array([0.75, 3.1, 1.15])      # Dirichlet, order (p-Ps, direct, o-Ps), see registration
EDGES = np.linspace(-2.0, 12.0, 561)     # 25 ps bins
CEN = 0.5 * (EDGES[:-1] + EDGES[1:])


def emg_cdf(t, tau, t0, sig):
    """CDF of exponential(tau) convolved with N(t0, sig).

    Two branches, because neither closed form is stable over the whole range:
      z <= a : erfcx form, whose argument is non-negative so erfcx cannot overflow
      z >  a : plain erfc form, whose exponent a^2/2 - a z is then negative
    """
    t = np.asarray(t, dtype=float)
    z = (t - t0) / sig
    a = sig / tau
    r2 = np.sqrt(2.0)
    Phi = 0.5 * erfc(-z / r2)
    u = (a - z) / r2
    lo = np.where(u >= 0.0, u, 0.0)              # safe arg for erfcx
    term_lo = 0.5 * np.exp(-0.5 * z * z) * erfcx(lo)
    ex = np.clip(0.5 * a * a - a * z, -700.0, 700.0)
    term_hi = 0.5 * np.exp(ex) * erfc(np.where(u < 0.0, u, 0.0))
    term = np.where(u >= 0.0, term_lo, term_hi)
    return np.clip(Phi - term, 0.0, 1.0)


def model_counts(theta, N):
    """Expected counts per bin. theta = [tau_oPs, t0, sig, w1, w2] (softmax weights)."""
    tau3, t0, sig, w1, w2 = theta
    w = np.exp(np.array([w1, w2, 0.0])); I = w / w.sum()
    taus = (TAU_PPS, TAU_DIR, tau3)
    p = np.zeros(len(CEN))
    for Ik, tk in zip(I, taus):
        c = emg_cdf(EDGES, tk, t0, sig)
        p += Ik * np.diff(c)
    # normalise within the fitted range: the model is a SHAPE fit conditional on
    # the event landing in the window, so a long-tau component whose tail leaves
    # the window is not penalised for it
    p = np.clip(p, 1e-15, None)
    p = p / p.sum()
    return N * p, I


def nll(theta, counts, use_prior):
    tau3, t0, sig = theta[0], theta[1], theta[2]
    if not (0.3 < tau3 < 8.0 and 0.01 < sig < 1.0 and abs(t0) < 2.0):
        return 1e12
    mu, I = model_counts(theta, counts.sum())
    val = float(np.sum(mu - counts * np.log(mu)))
    if use_prior:
        val -= float(dirichlet.logpdf(np.clip(I, 1e-9, None) / I.sum(), ALPHA))
    return val


def fit(counts, use_prior, seed_tau=2.0):
    best = None
    for st in (seed_tau, 1.6, 2.6):
        x0 = np.array([st, 0.0, 0.0909, np.log(0.15 / 0.23), np.log(0.62 / 0.23)])
        r = minimize(nll, x0, args=(counts, use_prior), method="Nelder-Mead",
                     options=dict(maxiter=20000, xatol=1e-6, fatol=1e-6))
        if best is None or r.fun < best.fun:
            best = r
    # error on tau3 from the curvature of the PROFILE likelihood: re-optimise the
    # other four parameters at tau3 +/- delta and fit a parabola. Two extra fits
    # instead of a full scan.
    t_hat = float(best.x[0]); f0 = float(best.fun)
    delta = max(0.004, 0.004 * t_hat)

    def prof(t):
        r = minimize(lambda q: nll(np.r_[t, q], counts, use_prior), best.x[1:],
                     method="Nelder-Mead",
                     options=dict(maxiter=6000, xatol=1e-7, fatol=1e-7))
        return float(r.fun)

    fm, fp = prof(t_hat - delta), prof(t_hat + delta)
    curv = (fp - 2.0 * f0 + fm) / (delta * delta)
    sd = float(1.0 / np.sqrt(curv)) if curv > 0 else float("nan")
    _, I = model_counts(best.x, counts.sum())
    return dict(tau=float(t_hat), sd=float(sd), t0=float(best.x[1]),
                sigma=float(best.x[2]), I=[float(v) for v in I], nll=float(f0))


def synth(tau3, N, rng, sigma=0.0909, t0=0.0, I=None):
    if I is None:
        I = ALPHA / ALPHA.sum()
    theta = np.r_[tau3, t0, sigma, np.log(I[0] / I[2]), np.log(I[1] / I[2])]
    mu, Itrue = model_counts(theta, N)
    return rng.poisson(mu).astype(float), Itrue


if __name__ == "__main__":
    rng = np.random.default_rng(20260822)
    TRUTH = [1.50, 1.80, 2.04, 2.40, 2.80]
    out = {"tau_pPs_ns": TAU_PPS, "tau_direct_ns": TAU_DIR,
           "dirichlet_alpha": ALPHA.tolist(), "bin_ns": float(EDGES[1] - EDGES[0]),
           "runs": []}
    for N in (10**6, 10**5):
        print(f"\nN = {N:,}")
        print(f"{'tau_true':>9s} {'prior?':>7s} {'tau_hat':>9s} {'sd':>8s} "
              f"{'err %':>8s} {'err/sd':>8s} {'I_oPs':>7s} {'dI':>7s}")
        for tt in TRUTH:
            counts, Itrue = synth(tt, N, rng)
            for up in (True, False):
                r = fit(counts, up, seed_tau=tt * 1.05)
                err = r["tau"] - tt
                row = dict(N=N, tau_true=tt, prior=up, **r,
                           err_ns=float(err), err_pct=float(100 * err / tt),
                           err_over_sd=float(abs(err) / r["sd"]),
                           I_oPs_true=float(Itrue[2]),
                           dI_oPs=float(r["I"][2] - Itrue[2]))
                out["runs"].append(row)
                print(f"{tt:9.2f} {str(up):>7s} {r['tau']:9.4f} {r['sd']:8.4f} "
                      f"{100*err/tt:8.2f} {abs(err)/r['sd']:8.2f} {r['I'][2]:7.3f} "
                      f"{r['I'][2]-Itrue[2]:+7.3f}")
    json.dump(out, open("b2_gate.json", "w"), indent=2)
    print("\nwrote b2_gate.json")
