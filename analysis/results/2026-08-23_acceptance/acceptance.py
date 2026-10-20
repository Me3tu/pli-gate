#!/usr/bin/env python3
"""A(tau), the acceptance of the +/-5 ns prompt box in the histogrammed variable.

Two independent routes, per the 2026-08-23 decision:
  route 1  analytic, A(tau) = Phi((W-tau)/sigma_e) * Phi((tau+W)/sigma_e)
           with W = 5 ns exact by construction and sigma_e = sigma_t/sqrt(2)
           = 64.3 ps from the production timing model. Zero free parameters.
  route 2  empirical, the delayed/randoms arm normalised by its flat interior.
           Flat by truth, so its shape IS the acceptance.

Comparison axis: acceptance relative to the flat interior, against delay in ns.
"""
from __future__ import annotations
import json, os, sys
from pathlib import Path
import numpy as np
from scipy.special import erf
from scipy.optimize import curve_fit

RES = Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[1]) / "2026-08-22_residual_structure"
W = 5.0                       # ns, W_PROMPT_NS, from the extractor
SIGMA_T = 0.0909              # ns, production per-single timing sigma
SIGMA_E = SIGMA_T / np.sqrt(2.0)   # 64.3 ps, derived in REGISTERED_ACC.md


def Phi(x):
    return 0.5 * (1.0 + erf(x / np.sqrt(2.0)))


def A_analytic(tau, sigma_e=SIGMA_E, w=W):
    tau = np.asarray(tau, float)
    return Phi((w - tau) / sigma_e) * Phi((tau + w) / sigma_e)


def A_binned(edges, sigma_e=SIGMA_E, w=W, n=200):
    """Bin-averaged acceptance, so it can be compared with binned measurements."""
    edges = np.asarray(edges, float)
    out = np.empty(len(edges) - 1)
    for i, (a, b) in enumerate(zip(edges[:-1], edges[1:])):
        x = np.linspace(a, b, n)
        out[i] = np.trapezoid(A_analytic(x, sigma_e, w), x) / (b - a)
    return out


if __name__ == "__main__":
    E = json.loads((RES / "a1_edge_spectrum.json").read_text())
    cen = np.array(E["bin_centres_ns"])
    meas = np.array(E["ratio_to_core"])
    err = np.array(E["err_ratio"])
    half = 0.5 * (cen[1] - cen[0])
    edges = np.r_[cen - half, cen[-1] + half]

    print(f"W = {W} ns exact; sigma_t = {1000*SIGMA_T:.1f} ps -> "
          f"sigma_e = sigma_t/sqrt(2) = {1000*SIGMA_E:.1f} ps  (NOT fitted)\n")
    ana = A_binned(edges)
    print(f"{'bin (ns)':>14s} {'measured':>10s} {'+/-':>7s} {'route 1':>9s} {'diff':>8s}")
    for a, b, m, e, v in zip(edges[:-1], edges[1:], meas, err, ana):
        mark = "  <-- edge" if (a < -4.4 or b > 4.4) else ""
        print(f"[{a:+5.1f},{b:+5.1f}) {m:10.3f} {e:7.3f} {v:9.3f} {m-v:+8.3f}{mark}")

    # R1: fit sigma_e to route 2 purely as the agreement statistic
    def model(c, s):
        return A_binned(edges, sigma_e=s)
    p, cov = curve_fit(lambda c, s: model(c, s), cen, meas, p0=[0.07],
                       sigma=np.maximum(err, 1e-4), absolute_sigma=True)
    s_fit, s_err = float(p[0]), float(np.sqrt(cov[0, 0]))
    print(f"\nR1 agreement test")
    print(f"  route 1, derived : sigma_e = {1000*SIGMA_E:.1f} ps")
    print(f"  route 2, fitted  : sigma_e = {1000*s_fit:.1f} +/- {1000*s_err:.1f} ps")
    z = abs(s_fit - SIGMA_E) / s_err
    print(f"  difference       : {1000*(s_fit-SIGMA_E):+.1f} ps = {z:.1f} sigma")
    print(f"  VERDICT: {'AGREE' if z < 3 else 'DISAGREE'}")
    interior = (cen > -4.4) & (cen < 4.4)
    chi2 = float((((meas - ana) / np.maximum(err, 1e-4)) ** 2)[interior].sum()) / interior.sum()
    print(f"  chi2/bin over the flat interior, route 1 vs measured: {chi2:.2f}")
    json.dump(dict(W_ns=W, sigma_t_ns=SIGMA_T, sigma_e_analytic_ns=SIGMA_E,
                   sigma_e_fitted_ns=s_fit, sigma_e_fitted_err_ns=s_err,
                   agreement_sigma=float(z), bin_edges_ns=edges.tolist(),
                   measured=meas.tolist(), measured_err=err.tolist(),
                   route1_binned=ana.tolist()),
              open("acceptance.json", "w"), indent=2)
    print("\nwrote acceptance.json")
