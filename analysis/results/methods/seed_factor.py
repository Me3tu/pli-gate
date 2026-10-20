#!/usr/bin/env python3
"""Run-to-run factor of the eight-tube lifetimes (Methods), written for release v1.1.0.

The published number had no script on disk (inventory flag F2). This script recomputes it from the three
independent GATE seeds of 2026-08-12_seed_repeatability/seed_results.json (44Sc eight-tube phantom, seeds
2000000, 3000000, 4000000; per tube the event-identity tail lifetime tau_hat and its tail count n_tail).

  chi2    = sum over tubes k and seeds s of (tau_ks - mean_s tau_ks)^2 / sigma_ks^2,  sigma_ks = tau_ks / sqrt(n_ks)
            (counting uncertainty of each run, deviation from the unweighted seed mean), dof = 8 x (3 - 1) = 16
  factor  = sqrt(chi2 / dof)
  95 %    = sqrt(chi2 / q) for q the 97.5 % and 2.5 % quantiles of chi2(dof)
  spreads = per tube, seed standard deviation (ddof 1) over the seed mean, against the mean relative counting
            uncertainty 1/sqrt(n); both averaged over the eight tubes, and the ratio of the two means.
Note: the chi2 is 27.150. With the inverse-variance weighted seed mean instead of the plain mean it is 27.137, which
prints as 27.1, so the printed 27.2 corresponds to the plain seed mean. The factor (1.30) and its interval
(0.97 to 1.98) do not depend on that choice at the printed precision.
Writes seed_factor.json next to this script."""
import json, os
from pathlib import Path
import numpy as np
from scipy.stats import chi2

RES = Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[1])
d = json.load(open(RES / "2026-08-12_seed_repeatability/seed_results.json"))
seeds = list(d)
tubes = list(d[seeds[0]])
T = np.array([[d[s][k]["tau_hat"] for s in seeds] for k in tubes])
N = np.array([[d[s][k]["n_tail"] for s in seeds] for k in tubes], dtype=float)
sig = T / np.sqrt(N)
c2 = float((((T - T.mean(1, keepdims=True)) / sig) ** 2).sum())
dof = T.shape[0] * (T.shape[1] - 1)
f = float(np.sqrt(c2 / dof))
lo, hi = float(np.sqrt(c2 / chi2.ppf(0.975, dof))), float(np.sqrt(c2 / chi2.ppf(0.025, dof)))
spread = 100 * T.std(1, ddof=1) / T.mean(1)
count = 100 * (1 / np.sqrt(N)).mean(1)
out = dict(doc=__doc__, seeds=seeds, tubes=tubes, chi2=c2, dof=dof, p_value=float(chi2.sf(c2, dof)), factor=f, factor_95=[lo, hi],
           per_tube=[dict(tube=k, tau_mean=float(T[i].mean()), spread_pct=float(spread[i]), counting_pct=float(count[i]))
                     for i, k in enumerate(tubes)],
           mean_spread_pct=float(spread.mean()), mean_counting_pct=float(count.mean()),
           ratio_of_means=float(spread.mean() / count.mean()))
json.dump(out, open(Path(__file__).resolve().parent / "seed_factor.json", "w"), indent=1)
print(f"chi2 {c2:.1f} for {dof} dof (p {out['p_value']:.3f}); factor {f:.2f} [{lo:.2f}, {hi:.2f}]")
print(f"mean seed spread {spread.mean():.2f} % against counting {count.mean():.2f} %, ratio of means {spread.mean() / count.mean():.2f}")
