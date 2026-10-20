#!/usr/bin/env python3
"""Review round 13 (2026-10-14), referee M1: the signed pedestal regression of Fig. 9c (2026-10-08_signed_eq3/signed_eq3.py)
with the eight tubes of the corrected-yield zero-delay arm C0 in place of the eight tubes of the discovery stream.

Everything else is unchanged: the twelve regions of the separately simulated configurations (disks, water cylinder,
NEMA) are taken from signed_eq3.json as published; for the C0 tubes the measured signed rule difference and its paired
shard-bootstrap uncertainty are computed from hists_c0_0ps.npz exactly as signed_eq3.py does (tau_from_hists of
v4d_paired.py, B = 4000, seed 20260826), and the prediction is the frozen pedestal relation (derive_pedestal_bias.
predicted_bias) at the over-subtraction measured on C0 (c0_compare.json, d_phi). The regression statistics are those of
signed_eq3.py (function copied). Post hoc (review)."""
from __future__ import annotations
import os, sys, re, json
from pathlib import Path
import numpy as np
RES = Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[1])
S = RES / "2026-08-17_derivation_test"
sys.path.insert(0, str(S)); sys.path.insert(0, str(RES / "2026-08-12_steinberger_allprompts"))
import steinberger_extract_ap as X  # noqa  (tau_from_hists uses X)
from derive_pedestal_bias import predicted_bias  # noqa
exec(compile(re.search(r"^def tau_from_hists\(.*?(?=^def |\Z)",
     (RES / "2026-08-19_V4_stats/v4d_paired.py").read_text(), re.S | re.M).group(0), "f", "exec"))
src = (RES / "2026-10-08_signed_eq3/signed_eq3.py").read_text()
exec(compile(re.search(r"^def stats\(.*?(?=^resolved = )", src, re.S | re.M).group(0), "stats", "exec"))
OUT = RES / "2026-10-14_known_randoms"
B = 4000
PUB = json.load(open(RES / "2026-10-08_signed_eq3/signed_eq3.json"))
C0C = {r["region"]: r["c0"] for r in json.load(open(RES / "2026-09-30_feedback3/c0_sim1/c0_compare.json"))["rows"]}
d = np.load(RES / "2026-09-30_feedback3/c0_sim1/hists_c0_0ps.npz", allow_pickle=True)
names = [str(s) for s in d["voi_names"]]
cs, cr, as_, ar = d["closest_sig"], d["closest_rnd"], d["all_sig"], d["all_rnd"]; ns = cs.shape[0]
rows_c0 = []
for roi in [n for n in names if n.startswith("T")]:
    ki = names.index(roi); tau = float(d["tau_cfg"][ki])
    t0c = tau_from_hists(cs[:, ki].sum(0), cr[:, ki].sum(0)); t0a = tau_from_hists(as_[:, ki].sum(0), ar[:, ki].sum(0))
    b0c = 100 * (t0c - tau) / tau; b0a = 100 * (t0a - tau) / tau
    rng = np.random.default_rng(20260826); tc = np.empty(B); ta = np.empty(B)
    for b in range(B):
        i = rng.integers(0, ns, ns)
        tc[b] = tau_from_hists(cs[i, ki].sum(0), cr[i, ki].sum(0)); ta[b] = tau_from_hists(as_[i, ki].sum(0), ar[i, ki].sum(0))
    ok = np.isfinite(tc) & np.isfinite(ta)
    sgn = (100 * (ta[ok] - tau) / tau) - (100 * (tc[ok] - tau) / tau); q = np.percentile(sgn, [16, 84])
    pred = -100.0 * predicted_bias(tau, C0C[roi]["d_phi"])
    rows_c0.append(dict(phantom="Sim-1 (C0)", roi=roi, tau=tau, b_closest=float(b0c), b_all=float(b0a), meas_signed=float(b0a - b0c),
                        pred=float(pred), pred_c0_compare=float(C0C[roi]["pedestal_pred"]), meas_c0_compare=float(C0C[roi]["did"]),
                        sig_signed=float(0.5 * (q[1] - q[0])), sig_c0_compare=float(C0C[roi]["sd_did"])))
    r = rows_c0[-1]
    print(f"{roi[:24]:24} signed {r['meas_signed']:6.2f} (c0_compare {r['meas_c0_compare']:6.2f}) pred {pred:6.2f} "
          f"(c0_compare {r['pred_c0_compare']:6.2f}) sig {r['sig_signed']:.3f} (c0_compare {r['sig_c0_compare']:.3f})", flush=True)
others = [r for r in PUB["rows"] if r["phantom"] != "Sim-1"]
rows = rows_c0 + others
resolved = [r for r in rows if not r["roi"].startswith("sphere_")]
out = dict(note=__doc__, B=B, rows_c0=rows_c0, signed=dict(all20=stats(rows, "meas_signed", "sig_signed"),
                                                       resolved14=stats(resolved, "meas_signed", "sig_signed")),
           published_discovery=dict(all20={k: PUB["signed"]["all20"][k] for k in ("slope", "scatter_pp", "slope_with_scatter", "chi2_per_region", "max_abs_pull", "r2_identity")}))
for s in ("all20", "resolved14"):
    x = out["signed"][s]
    print(f"C0 tubes {s:10} n={x['n']} chi2/reg {x['chi2_per_region']:.2f} slope {x['slope']:.3f} scatter {x['scatter_pp']:.3f} pp "
          f"slope|scatter {x['slope_with_scatter']:.3f} max pull {x['max_abs_pull']:.2f} ({x['max_pull_roi']}) >2: {x['n_beyond_2']} "
          f"sign agree {x['sign_agree']}/{x['n']} R2 {x['r2_identity']:.3f}")
json.dump(out, open(OUT / "signed_eq3_c0.json", "w"), indent=1)
