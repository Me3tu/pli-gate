#!/usr/bin/env python3
"""Review round 6 (2026-10-08): the Eq. 3 consistency regression on the SIGNED rule
difference.

The published regression (2026-08-27_fig9_v2/build_fig9_json.py) compares
  meas = |b_closest| - |b_all|            (folded, test_delta.py:62)
against
  pred = -100 * predicted_bias(tau, d_phi) (derive_pedestal_bias.py, frozen),
which is the lifetime reduction that the differential pedestal d_phi causes, i.e. a
prediction of the SIGNED difference b_all - b_closest. Where b_all <= 0 and
b_closest <= 0 the two measured statistics coincide; where a bias is positive or
straddles zero the folded statistic differs from the quantity predicted, and its
bootstrap distribution folds.

Here the measured statistic is the signed difference d_b = b_all - b_closest, with
its paired shard-bootstrap uncertainty (same histograms, same B = 4000, same seed
20260826 and same resampling as b2c_final.py), and the frozen prediction. The
statistics of the published regression are recomputed unchanged in form.
Nothing is refitted except the single scatter constant, as before.
Post hoc (review), not registered.
"""
from __future__ import annotations
import os, sys, re, json
from pathlib import Path
import numpy as np
RES = Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[1])
S = RES / "2026-08-17_derivation_test"
sys.path.insert(0, str(S)); sys.path.insert(0, str(RES / "2026-08-12_steinberger_allprompts"))
import steinberger_extract_ap as X  # noqa
from derive_pedestal_bias import predicted_bias  # noqa
exec(compile(re.search(r"^def tau_from_hists\(.*?(?=^def |\Z)",
     (RES / "2026-08-19_V4_stats/v4d_paired.py").read_text(), re.S | re.M).group(0), "f", "exec"))
OUT = RES / "2026-10-08_signed_eq3"
ARMS = {"Sim-1": RES / "2026-08-22_G3_consequences/hists_sim1.npz",
        "NEMA IQ": RES / "2026-08-22_G3_consequences/hists_nema.npz",
        "Steinberger water": RES / "2026-08-19_V4_stats/hists_water.npz",
        "Steinberger disks": RES / "2026-08-19_V4_stats/hists_disks.npz"}
B = 4000
pts = json.loads((S / "delta_test.json").read_text())["points"]
fig9 = {p["roi"]: p for p in json.load(open(RES / "2026-08-27_fig9_v2/fig9_pedestal_derivation.json"))["points"]}
cache = {}
rows = []
for p in pts:
    ph, roi = p["phantom"], p["roi"]
    if ph not in cache:
        cache[ph] = np.load(ARMS[ph], allow_pickle=True)
    d = cache[ph]
    names = [str(s) for s in d["voi_names"]]
    if roi not in names:
        print("MISSING", ph, roi); continue
    ki = names.index(roi); tau = float(d["tau_cfg"][ki])
    cs, cr, as_, ar = d["closest_sig"], d["closest_rnd"], d["all_sig"], d["all_rnd"]
    ns = cs.shape[0]
    t0c = tau_from_hists(cs[:, ki].sum(0), cr[:, ki].sum(0))
    t0a = tau_from_hists(as_[:, ki].sum(0), ar[:, ki].sum(0))
    b0c = 100 * (t0c - tau) / tau; b0a = 100 * (t0a - tau) / tau
    rng = np.random.default_rng(20260826)
    tc = np.empty(B); ta = np.empty(B)
    for b in range(B):
        i = rng.integers(0, ns, ns)
        tc[b] = tau_from_hists(cs[i, ki].sum(0), cr[i, ki].sum(0))
        ta[b] = tau_from_hists(as_[i, ki].sum(0), ar[i, ki].sum(0))
    ok = np.isfinite(tc) & np.isfinite(ta)
    bc = 100 * (tc[ok] - tau) / tau; ba = 100 * (ta[ok] - tau) / tau
    sgn = ba - bc; fold = np.abs(bc) - np.abs(ba)
    q = np.percentile(sgn, [16, 84]); qf = np.percentile(fold, [16, 84])
    pred = -100.0 * predicted_bias(p["tau"], p["d_phi"])
    rows.append(dict(phantom=ph, roi=roi, tau=tau, b_closest=float(b0c), b_all=float(b0a),
                     meas_signed=float(b0a - b0c), meas_folded=float(abs(b0c) - abs(b0a)),
                     meas_folded_published=fig9[roi]["meas"], pred=float(pred),
                     pred_published=fig9[roi]["pred"],
                     sig_signed=float(0.5 * (q[1] - q[0])), sig_folded=float(0.5 * (qf[1] - qf[0])),
                     sig_folded_published=fig9[roi]["sigma_paired"],
                     straddle_all=float(min((ba > 0).mean(), (ba < 0).mean())),
                     straddle_closest=float(min((bc > 0).mean(), (bc < 0).mean())),
                     n_nan=int((~ok).sum())))
    r = rows[-1]
    print(f"{ph[:12]:12} {roi[:22]:22} b_c {b0c:+7.2f} b_a {b0a:+7.2f} signed {r['meas_signed']:7.2f} "
          f"folded {r['meas_folded']:7.2f} (pub {r['meas_folded_published']:7.2f}) pred {pred:6.2f} "
          f"sig_s {r['sig_signed']:5.2f} sig_f {r['sig_folded']:5.2f} (pub {r['sig_folded_published']:5.2f})", flush=True)


def stats(sel, key_m, key_s):
    m = np.array([r[key_m] for r in sel]); qv = np.array([r["pred"] for r in sel])
    s = np.array([r[key_s] for r in sel]); res = m - qv; n = len(sel)
    c0 = float(np.sum((res / s) ** 2) / n)
    w = 1 / s ** 2; sl = float((w * qv * m).sum() / (w * qv * qv).sum())
    c_sl = float(np.sum(((m - sl * qv) / s) ** 2) / (n - 1))
    def chi(ss): return float(np.sum(res ** 2 / (s ** 2 + ss ** 2)) / n)
    if chi(0.0) <= 1.0: ss = 0.0
    else:
        lo, hi = 0.0, 200.0
        for _ in range(300):
            mid = 0.5 * (lo + hi)
            if chi(mid) > 1.0: lo = mid
            else: hi = mid
        ss = 0.5 * (lo + hi)
    tot = np.sqrt(s ** 2 + ss ** 2); wt = 1 / tot ** 2
    sl_f = float((wt * qv * m).sum() / (wt * qv * qv).sum())
    med = float(np.median(s)); share = ss ** 2 / (ss ** 2 + med ** 2)
    pulls = res / tot
    r2 = 1 - float((res ** 2).sum()) / float(((m - m.mean()) ** 2).sum())
    return dict(n=n, chi2_per_region=c0, slope=sl, chi2_dof_slope=c_sl, scatter_pp=ss,
                slope_with_scatter=sl_f, median_sigma=med, stat_share=1 - share, scatter_share=share,
                max_abs_pull=float(np.abs(pulls).max()),
                max_pull_roi=sel[int(np.argmax(np.abs(pulls)))]["roi"],
                n_beyond_2=int((np.abs(pulls) > 2).sum()), r2_identity=r2,
                sign_agree=int(np.sum(np.sign(m) == np.sign(qv))),
                pulls={r["roi"]: float(x) for r, x in zip(sel, pulls)})


resolved = [r for r in rows if not r["roi"].startswith("sphere_")]
out = dict(note=__doc__, B=B, rows=rows,
           signed=dict(all20=stats(rows, "meas_signed", "sig_signed"),
                       resolved14=stats(resolved, "meas_signed", "sig_signed")),
           folded_recomputed=dict(all20=stats(rows, "meas_folded", "sig_folded"),
                                  resolved14=stats(resolved, "meas_folded", "sig_folded")))
for k in ("signed", "folded_recomputed"):
    for s in ("all20", "resolved14"):
        x = out[k][s]
        print(f"{k:18} {s:10} n={x['n']} chi2/reg {x['chi2_per_region']:.2f} slope {x['slope']:.3f} "
              f"chi2/dof(slope) {x['chi2_dof_slope']:.2f} scatter {x['scatter_pp']:.3f} pp slope|scatter "
              f"{x['slope_with_scatter']:.3f} median sig {x['median_sigma']:.3f} stat share {100*x['stat_share']:.0f}% "
              f"max pull {x['max_abs_pull']:.2f} ({x['max_pull_roi']}) >2: {x['n_beyond_2']} sign agree {x['sign_agree']}/{x['n']} R2 {x['r2_identity']:.3f}")
json.dump(out, open(OUT / "signed_eq3.json", "w"), indent=1)
print("wrote", OUT / "signed_eq3.json")
