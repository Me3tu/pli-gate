#!/usr/bin/env python3
"""N4 — derive the bias analytically for the estimator as implemented.

THE ESTIMATOR (subtract_and_mle, verbatim):
    fit window  W = [tau_min, T] = [1.0, 5.0] ns
    h_s = h_total - h_random                       (bin-wise subtraction)
    tau_obs = sum(bc * h_s)/sum(h_s) - tau_min     (truncated first moment)
    tau_hat = truncation_correct(tau_obs, T - tau_min)

THE PERTURBATION: the delayed window removes MORE than the true accidental
content. The excess is flat in delay -- accidentals are uniform in the window by
construction -- so the estimator sees

    signal:   N counts, density proportional to exp(-t/tau) on [tau_min, T]
    removed:  P = phi*N counts, density UNIFORM on [tau_min, T]

A uniform pedestal has mean (tau_min+T)/2, which lies ABOVE the truncated
exponential mean for every tau of interest here. Removing it therefore removes
proportionally more from the long-delay side, the first moment falls, and
tau_hat falls with it. That is the whole mechanism.

    m1(tau)   = truncated-exponential mean on [tau_min, T], minus tau_min
    m1_ped    = (T - tau_min)/2
    tau_obs(phi) = [m1(tau) - phi*m1_ped] / (1 - phi)
    tau_hat(phi) = truncation_correct(tau_obs(phi), T - tau_min)

    predicted fractional bias(phi) = tau_hat(phi)/tau_hat(0) - 1

NO FREE PARAMETERS. tau_min, T and the estimator come from the production code;
phi is measured per ROI in N2.

Tested exactly as eq. (2) was: prediction against measurement, identity line,
all phantoms, no tuning.
"""
from __future__ import annotations
import json
import os
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

S = Path(__file__).parent
RES = Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[1])
EXDIR = RES / "2026-08-12_steinberger_allprompts"
STYLE = RES / "2026-05-06_1715_phase3r2_d3_and_figure_baseline"
sys.path.insert(0, str(EXDIR))
import steinberger_extract_ap as EX  # noqa: E402

TAU_MIN = EX.TAU_MIN_NS      # 1.0
T_HI = EX.W_PROMPT_NS        # 5.0
WIDTH = T_HI - TAU_MIN       # 4.0


def trunc_mean(tau, lo=TAU_MIN, hi=T_HI):
    """E[t | lo < t < hi] - lo  for t ~ Exp(1/tau)."""
    a, b = lo / tau, hi / tau
    # Stable form: factor out exp(-a) so nothing underflows for small tau, where
    # a can reach 10^3 and exp(-a) flushes to zero, giving 0/0.
    e = np.exp(-(b - a))
    num = (1 + a) - (1 + b) * e
    den = 1.0 - e
    return tau * num / den - lo


def invert_trunc(m_target, lo=TAU_MIN, hi=T_HI):
    """Recover tau from the truncated first moment -- this is exactly what
    EX.truncation_correct does, solved here directly so the two can be
    cross-checked against each other."""
    if not np.isfinite(m_target) or m_target <= 0:
        return float("nan")
    f = lambda t: trunc_mean(t, lo, hi) - m_target
    try:
        return brentq(f, 1e-3, 500.0)
    except ValueError:
        return float("nan")


def predicted_bias(tau_true, phi):
    """Fractional bias in tau_hat caused by over-subtracting phi*N uniform counts."""
    m1 = trunc_mean(tau_true)
    m_ped = WIDTH / 2.0                      # uniform pedestal mean, minus tau_min
    if phi >= 1:
        return float("nan")
    m_obs = (m1 - phi * m_ped) / (1.0 - phi)
    tau_hat = invert_trunc(m_obs)
    tau_hat0 = invert_trunc(m1)              # = tau_true, by construction
    if not (np.isfinite(tau_hat) and np.isfinite(tau_hat0)):
        return float("nan")
    return tau_hat / tau_hat0 - 1.0


def main():
    # self-consistency: the analytic inverse must reproduce the production one
    print("cross-check of the inversion against the production truncation_correct:")
    for t in (1.56, 1.80, 2.04, 2.72):
        m = trunc_mean(t)
        print(f"  tau={t:.2f}  trunc_mean={m:.5f}  invert->{invert_trunc(m):.5f}  "
              f"EX.truncation_correct->{EX.truncation_correct(m, WIDTH):.5f}")

    # small-phi slope, the quantity N3b measured as 1.15
    print("\nanalytic slope d(bias)/d(phi) as phi -> 0:")
    for t in (1.56, 1.80, 2.04, 2.72):
        h = 1e-5
        sl = (predicted_bias(t, h) - predicted_bias(t, 0.0)) / h
        print(f"  tau={t:.2f} ns : {sl:+.3f}    (measured, all tau pooled: -1.15)")

    d = json.load(open(S / "delta_test.json"))["points"]
    rows = []
    for p in d:
        pb = predicted_bias(p["tau"], p["d_phi"])
        rows.append(dict(**p, pred_d_bias=-100.0 * pb))
    print(f"\n{'phantom':>18} {'ROI':>26} {'tau':>6} {'d_phi':>8} "
          f"{'pred':>8} {'meas':>8} {'+/-':>5} {'m/p':>6}")
    for r in sorted(rows, key=lambda x: x["d_phi"]):
        print(f"{r['phantom']:>18} {r['roi']:>26} {r['tau']:6.3f} {r['d_phi']:+8.4f} "
              f"{r['pred_d_bias']:8.2f} {r['d_bias']:+8.2f} {r['d_bias_sigma']:5.2f} "
              f"{r['d_bias']/r['pred_d_bias'] if r['pred_d_bias'] else float('nan'):6.2f}")

    w = [r for r in rows if r["d_bias_sigma"] < 3.0]
    m = np.array([r["d_bias"] for r in w]); q = np.array([r["pred_d_bias"] for r in w])
    from scipy.stats import pearsonr, spearmanr
    slope = float((q * m).sum() / (q * q).sum())
    r2 = 1 - float(((m - slope * q) ** 2).sum()) / float(((m - m.mean()) ** 2).sum())
    r2_id = 1 - float(((m - q) ** 2).sum()) / float(((m - m.mean()) ** 2).sum())
    print(f"\n=== well-measured, n={len(w)} ===")
    print(f"  Pearson r(pred, meas) = {pearsonr(q, m)[0]:+.4f}")
    print(f"  Spearman rho          = {spearmanr(q, m)[0]:+.4f}")
    print(f"  best-fit slope through origin = {slope:.3f}   R^2 = {r2:+.4f}")
    print(f"  R^2 about the IDENTITY LINE (slope fixed at 1, no fit) = {r2_id:+.4f}")
    print(f"  mean meas/pred = {np.mean(m/q):.3f}   spread {np.min(m/q):.3f}-{np.max(m/q):.3f}")

    json.dump(dict(points=rows, slope_vs_prediction=slope, r2_fitted=r2,
                   r2_identity=r2_id, tau_min=TAU_MIN, t_hi=T_HI),
              open(S / "n4_derivation.json", "w"), indent=1)

    sys.path.insert(0, str(STYLE))
    from style_baseline import setup_publication_rcparams, save_figure, figure_dims  # noqa
    setup_publication_rcparams()
    import matplotlib.pyplot as plt
    COL = {"Sim-1": "#1C5D86", "Steinberger disks": "#8F5510",
           "Steinberger water": "#2E6F5E", "NEMA IQ": "#8C2F3A"}
    MRK = {"Sim-1": "o", "Steinberger disks": "s", "Steinberger water": "^", "NEMA IQ": "D"}
    fig, ax = plt.subplots(figsize=figure_dims(120, 100))
    fin = np.isfinite(q) & np.isfinite(m)
    lim = [0, float(max(q[fin].max(), m[fin].max()) * 1.12)]
    ax.plot(lim, lim, color="0.35", lw=1.0, ls="--", zorder=1, label="identity (no free parameter)")
    for ph in COL:
        sel = [r for r in rows if r["phantom"] == ph]
        if not sel:
            continue
        ax.errorbar([r["pred_d_bias"] for r in sel], [r["d_bias"] for r in sel],
                    yerr=[r["d_bias_sigma"] for r in sel], fmt=MRK[ph], ms=5,
                    color=COL[ph], mfc=COL[ph], mec="white", mew=0.6,
                    ecolor=COL[ph], elinewidth=0.8, capsize=2, lw=0, zorder=3, label=ph)
    ax.set_xlim(*lim); ax.set_ylim(-8, lim[1])
    ax.set_xlabel(r"predicted $\Delta|$bias$|$ from uniform-pedestal over-subtraction (%)")
    ax.set_ylabel(r"measured $\Delta|$bias$|$  (%)")
    ax.set_title("Pedestal over-subtraction: derivation vs measurement", fontsize=9.5)
    ax.legend(fontsize=6.8, loc="upper left", framealpha=0.95)
    ax.grid(True, alpha=0.18, lw=0.4)
    (S / "figures").mkdir(exist_ok=True)
    save_figure(fig, S / "figures/n4_derivation")
    plt.close(fig)
    print(f"\nwrote {S}/figures/n4_derivation.{{pdf,png}} and n4_derivation.json")


if __name__ == "__main__":
    main()
