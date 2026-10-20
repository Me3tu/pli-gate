#!/usr/bin/env python3
"""Review round 8 (2026-10-08): per-insert intervals and mean paired difference for the oxygen phantom
rerun with database v8 (plasma k_O2) and the production binary (corrected I3 convention).
Same estimators as 2026-09-30_feedback3/oxygen/oxygen_boot.py:
  per insert: tail mean above tau_min = 1 ns minus tau_min, sigma = tau_hat/sqrt(N); exact chi2(2N) interval;
  mean difference over the five inserts: parametric bootstrap of Gamma(N, tau_hat/N) per insert (B = 200000);
  figure source: re-extraction from the v8 triples (travel-time corrected), sphere of 16 mm about each insert
  centre on the source truth position (ASSIGN_R_MM of per_insert_hypoxia.py), event-level percentile bootstrap.
Writes oxygen_boot_v8.json next to this script."""
import os, sys, json, csv
from pathlib import Path
import numpy as np, h5py
from scipy.stats import chi2
RES = Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[1])
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from reconstruction.correct_travel_time import correct as travel_time_correct
HERE = Path(__file__).resolve().parent
OUT = Path(os.environ.get("PLI_RAW", "<missing>")) / "oxygen_v8"
if not Path(OUT / "triples_detector_hypoxia_te_v8_rcut075mm.h5").exists(): raise SystemExit(f"{OUT / "triples_detector_hypoxia_te_v8_rcut075mm.h5"}: raw data not included in the release, see analysis/README.md")
INS = [("TumorTissue_pO2_2", -60.0, 2.0), ("TumorTissue_pO2_6", -30.0, 6.0), ("TumorTissue_pO2_10", 0.0, 10.0),
       ("TumorTissue_pO2_20", 30.0, 20.0), ("TumorTissue_pO2_40", 60.0, 40.0)]
fr = {r["insert"]: r for r in csv.DictReader(open(HERE / "hypoxia_te_v8_per_insert.csv"))}
tau, _ = travel_time_correct(OUT / "triples_detector_hypoxia_te_v8_rcut075mm.h5")
with h5py.File(OUT / "triples_truth_hypoxia_te_v8_rcut075mm.h5", "r") as f:
    sx, sy, sz = (f["truth"][k][:] for k in ("source_x", "source_y", "source_z"))
cfg_by = {n: float(fr[n]["tau_TE_configured"]) for n, _, _ in INS}
contrast_ps = (cfg_by["TumorTissue_pO2_2"] - cfg_by["TumorTissue_pO2_40"]) * 1e3
rng = np.random.default_rng(20261008)
rows = []
for name, cx, po2 in INS:
    n = int(fr[name]["n_used"]); th = float(fr[name]["tau_recovered"]); cfg = cfg_by[name]
    sig = th / np.sqrt(n)
    elo, ehi = 2 * n * th / chi2.ppf(0.975, 2 * n), 2 * n * th / chi2.ppf(0.025, 2 * n)
    rows.append(dict(insert=name, pO2=po2, n_tail=n, n_total=int(fr[name]["n_total"]), I3=float(fr[name]["I3"]),
                     tau_cfg=cfg, tau_hat=th, sigma_analytic_ps=sig * 1e3, diff_ps=(th - cfg) * 1e3,
                     z=(th - cfg) / sig, exact95_ps=[(elo - cfg) * 1e3, (ehi - cfg) * 1e3]))
    r = rows[-1]
    print(f"{name:19s} N={n}/{r['n_total']} diff={r['diff_ps']:+7.1f} sig={r['sigma_analytic_ps']:.1f} "
          f"z={r['z']:+.2f} exact[{r['exact95_ps'][0]:+5.0f},{r['exact95_ps'][1]:+5.0f}]")
Bm = 200000
m = np.zeros(Bm)
for r in rows:
    n, th = r["n_tail"], r["tau_hat"]
    m += (rng.gamma(n, th / n, size=Bm) - r["tau_cfg"]) * 1e3
m /= len(rows)
mean_d = float(np.mean([r["diff_ps"] for r in rows]))
lo, hi = np.percentile(m, [2.5, 97.5])
sd = float(np.sqrt(sum(r["sigma_analytic_ps"] ** 2 for r in rows)) / len(rows))
chi2_5 = float(sum(r["z"] ** 2 for r in rows))
# slope of tau_hat against pO2 (weighted least squares), against the configured slope
x = np.array([r["pO2"] for r in rows]); y = np.array([r["tau_hat"] for r in rows]) * 1e3
w = 1 / np.array([r["sigma_analytic_ps"] for r in rows]) ** 2
X = np.vstack([np.ones_like(x), x]).T
cov = np.linalg.inv(X.T @ (X * w[:, None])); beta = cov @ (X.T @ (w * y))
ycfg = np.array([r["tau_cfg"] for r in rows]) * 1e3
bcfg = np.polyfit(x, ycfg, 1)[0]
res = dict(materials_sha12="23f607ea14c5", per_insert=rows, mean_diff_ps=mean_d,
           mean_diff_param95_ps=[float(lo), float(hi)], mean_diff_sd_ps=sd,
           mean_diff_normal95_ps=[mean_d - 1.96 * sd, mean_d + 1.96 * sd],
           chi2_vs_cfg=chi2_5, dof=len(rows),
           configured_contrast_ps=contrast_ps,
           mean_sigma_insert_ps=float(np.mean([r["sigma_analytic_ps"] for r in rows])),
           ratio_sigma_insert_to_contrast=float(np.mean([r["sigma_analytic_ps"] for r in rows]) / contrast_ps),
           slope_ps_per_mmHg=float(beta[1]), slope_sd_ps_per_mmHg=float(np.sqrt(cov[1, 1])),
           slope_cfg_ps_per_mmHg=float(bcfg))
print(f"mean diff {mean_d:+.1f} ps  parametric 95% [{lo:+.1f},{hi:+.1f}]  sd {sd:.1f}  chi2 {chi2_5:.2f}/5")
print(f"contrast {contrast_ps:.3f} ps  mean sigma {res['mean_sigma_insert_ps']:.1f}  ratio {res['ratio_sigma_insert_to_contrast']:.1f}")
print(f"slope {beta[1]:+.2f} +- {np.sqrt(cov[1,1]):.2f} ps/mmHg (configured {bcfg:+.4f})")

rng2 = np.random.default_rng(7)
B = 20000
allm = np.zeros(B)
rex = []
for name, cx, po2 in INS:
    cfg = cfg_by[name]
    sel = (sx - cx) ** 2 + sy ** 2 + sz ** 2 <= 16.0 ** 2
    s = tau[sel]; s = s[(s > 1.0) & np.isfinite(s)]
    th = s.mean() - 1.0
    bt = s[rng2.integers(0, s.size, size=(B, s.size))].mean(axis=1) - 1.0
    allm += (bt - cfg) * 1e3
    blo, bhi = np.percentile(bt, [2.5, 97.5])
    rex.append(dict(insert=name, pO2=po2, n_tail=int(s.size), n_insert=int(sel.sum()), tau_cfg=cfg, tau_hat=float(th),
                    diff_ps=float((th - cfg) * 1e3), sigma_analytic_ps=float(th / np.sqrt(s.size) * 1e3),
                    boot_sd_ps=float(bt.std(ddof=1) * 1e3), boot95_ps=[float((blo - cfg) * 1e3), float((bhi - cfg) * 1e3)]))
allm /= len(INS)
md = float(np.mean([r["diff_ps"] for r in rex]))
mlo, mhi = np.percentile(allm, [2.5, 97.5])
res["reextraction"] = dict(per_insert=rex, mean_diff_ps=md, mean_diff_boot95_ps=[float(mlo), float(mhi)],
                           mean_sigma_ps=float(np.mean([r["sigma_analytic_ps"] for r in rex])),
                           mean_bootsd_ps=float(np.mean([r["boot_sd_ps"] for r in rex])),
                           halfwidth_range_ps=[float(min((r["boot95_ps"][1] - r["boot95_ps"][0]) / 2 for r in rex)),
                                               float(max((r["boot95_ps"][1] - r["boot95_ps"][0]) / 2 for r in rex))],
                           B=B, note="event-level percentile bootstrap: tail events resampled within each insert; mean over inserts per replicate")
for r in rex:
    print(f"REX {r['insert']:19s} N={r['n_tail']} diff={r['diff_ps']:+7.1f} sig={r['sigma_analytic_ps']:.1f} "
          f"bootsd={r['boot_sd_ps']:.1f} boot95=[{r['boot95_ps'][0]:+.0f},{r['boot95_ps'][1]:+.0f}]")
print(f"REX mean diff {md:+.1f}  event-boot 95% [{mlo:+.1f},{mhi:+.1f}]  mean sigma {res['reextraction']['mean_sigma_ps']:.1f}")
old = json.load(open(RES / "2026-09-30_feedback3/oxygen/oxygen_boot.json"))
res["previous_v2"] = dict(mean_diff_ps=old["mean_diff_ps"], mean_diff_param95_ps=old["mean_diff_param95_ps"],
                          n_tail=[r["n_tail"] for r in old["per_insert"]], configured_contrast_ps=old["configured_contrast_ps"],
                          ratio_sigma_insert_to_contrast=old["ratio_sigma_insert_to_contrast"])
res["tail_count_ratio_v8_over_v2"] = float(sum(r["n_tail"] for r in rows) / sum(res["previous_v2"]["n_tail"]))
print("tail count ratio v8/v2", res["tail_count_ratio_v8_over_v2"])
json.dump(res, open(HERE / "oxygen_boot_v8.json", "w"), indent=1)
