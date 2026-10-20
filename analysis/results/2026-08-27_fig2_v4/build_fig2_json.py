#!/usr/bin/env python3
"""Single source of truth for fig2, built from the FREEZE pipeline.

The 25 August fig2 plotted `2026-08-19_V3b_offset/v3b_vs_closure.json`, which
descends from the MAY phase3s_2 CSVs -- a different, pre-correction recipe. The
caption, the abstract and `numbers_freeze_ejmp.json:closure_12_tubes` all quote
the POST-CORRECTION O4 rerun. This rebuilds the figure on that rerun so the
panel, the caption and the freeze are the same numbers.

Inputs (the exact files the freeze entry names):
  2026-08-13_overnight/O4_rerun/sim1_calibration/analysis/post_taumin1p0.csv
  2026-08-13_overnight/O4_rerun/sim2_mercolli_te/analysis/post_taumin1p0.csv
  recipe: filter_triples_tof_pann --r-cut 75 ; analyze_per_tube_tau --tau-min 1.0
  binary 3522d567...  materials DB 2a26b68e...
"""
from __future__ import annotations
import csv, json, os
from pathlib import Path
import numpy as np
RES = Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[1])
O4 = RES / "2026-08-13_overnight/O4_rerun"
OUT = RES / "2026-08-27_fig2_v4"
SETS = [("Sc-44", "Sim-1", O4 / "sim1_calibration/analysis/post_taumin1p0.csv"),
        ("I-124", "Sim-2 TE", O4 / "sim2_mercolli_te/analysis/post_taumin1p0.csv")]

tubes = []
for iso, sim, p in SETS:
    for r in csv.DictReader(open(p)):
        if not r["tube"]:
            continue
        cfg = float(r["tau_TE_configured"]); rec = float(r["tau_recovered"])
        sd = float(r["sigma"]); res = float(r["residual_pct"])
        tubes.append(dict(isotope=iso, sim=sim, tube=r["tube"], material=r["material"],
                          tau_cfg=cfg, tau_rec=rec, sigma_ns=sd, residual_pct=res,
                          residual_ps=1000.0 * (rec - cfg),
                          sigma_ps=1000.0 * sd, n_used=int(r["n_used"])))

per_iso = {}
for iso in ("Sc-44", "I-124"):
    a = np.array([abs(t["residual_pct"]) for t in tubes if t["isotope"] == iso])
    per_iso[iso] = dict(n=len(a), mean_abs_pct=float(a.mean()), worst_pct=float(a.max()),
                        within_1pct=int((a <= 1).sum()))
allp = np.array([t["residual_pct"] for t in tubes])

# trend of the residual with the configured lifetime, weighted
x = np.array([t["tau_cfg"] for t in tubes])
y = np.array([t["residual_ps"] for t in tubes])
s = np.array([t["sigma_ps"] for t in tubes])
w = 1.0 / s ** 2
X = np.c_[np.ones_like(x), x]
C = np.linalg.inv(X.T @ (w[:, None] * X)); b = C @ (X.T @ (w * y))
r = y - X @ b
dof = len(x) - 2
chi2_dof = float((w * r ** 2).sum() / dof)
se = np.sqrt(np.diag(C) * max(chi2_dof, 1.0))

# systematic floor bringing chi2/dof to 1 about the fitted trend
def chi2(ss):
    return float((r ** 2 / (s ** 2 + ss ** 2)).sum() / dof)
if chi2(0.0) <= 1.0:
    ss = 0.0
else:
    lo, hi = 0.0, 200.0
    for _ in range(300):
        m = 0.5 * (lo + hi)
        if chi2(m) > 1.0: lo = m
        else: hi = m
    ss = 0.5 * (lo + hi)
share = ss ** 2 / (ss ** 2 + float(np.median(s)) ** 2)

out = dict(
    source="POST-CORRECTION O4 rerun (the freeze pipeline)",
    inputs=[str(p.relative_to(RES.parent)) for _, _, p in SETS],
    recipe="filter_triples_tof_pann --r-cut 75 ; analyze_per_tube_tau --tau-min 1.0",
    binary="3522d567...", materials_db="2a26b68e...",
    n_tubes=len(tubes), tubes=tubes, per_isotope=per_iso,
    pooled=dict(rms_pct=float(np.sqrt((allp ** 2).mean())),
                mean_abs_pct=float(np.abs(allp).mean()),
                within_1pct=int((np.abs(allp) <= 1).sum())),
    trend=dict(slope_ps_per_ns=float(b[1]), slope_err=float(se[1]),
               intercept_ps=float(b[0]), intercept_err=float(se[0]),
               chi2_dof_stat_only=chi2_dof, sigma_sys_ps=float(ss),
               median_sigma_ps=float(np.median(s)), systematic_share=float(share)),
    supersedes=dict(file="2026-08-25_fig2_regen/fig2_v3_summary.json",
                    plotted="v3b / May phase3s_2 pipeline",
                    numbers="Sc-44 0.387 %, I-124 0.446 %, pooled RMS 0.422 %",
                    why="a different, pre-correction recipe from the one the caption, the "
                        "abstract and the freeze all quote"))
json.dump(out, open(OUT / "fig2_recipe_validation.json", "w"), indent=1)
for iso, v in per_iso.items():
    print(f"  {iso:6}  mean |resid| {v['mean_abs_pct']:.4f} %   worst {v['worst_pct']:.3f} %   "
          f"within +-1 %: {v['within_1pct']}/{v['n']}")
print(f"  pooled  RMS {out['pooled']['rms_pct']:.4f} %   within +-1 %: "
      f"{out['pooled']['within_1pct']}/{out['n_tubes']}")
t = out["trend"]
print(f"  trend   slope {t['slope_ps_per_ns']:+.2f} +- {t['slope_err']:.2f} ps/ns, "
      f"intercept {t['intercept_ps']:+.1f} ps")
print(f"  floor   chi2/dof {t['chi2_dof_stat_only']:.3f} -> sigma_sys {t['sigma_sys_ps']:.2f} ps "
      f"({100*t['systematic_share']:.0f} % of variance)")
