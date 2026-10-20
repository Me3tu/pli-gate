#!/usr/bin/env python3
"""Referee round 2, MC11: summary of the synthetic estimator study (estimator_synth.py).

For every condition and estimator: valid fits, failure fraction, bias, empirical SD and RMSE
(all in % of the true o-Ps lifetime, with the Monte Carlo standard error of the bias), and
the coverage of the nominal 68 and 95 % intervals: shard-bootstrap percentile intervals for
the component estimators (PUB, JOINT), delta-method +-1 and +-1.96 sigma for the moment
estimator (TM). Also the mean fraction of negative randoms-subtracted bins, and the ratio of
the mean bootstrap SD to the empirical SD.
Writes estimator_synth/summary.json and estimator_synth/tableS10_rows.tex."""
import json, os, sys
from pathlib import Path
import numpy as np

D = Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[1]) / "2026-09-30_referee2/estimator_synth"
src = D / ("reps.json" if (D / "reps.json").exists() and "--partial" not in sys.argv else "reps_partial.json")
rows = json.load(open(src)); C = {c["id"]: c for c in json.load(open(D / "conditions.json"))}
print("source", src.name, "replicates", len(rows))


def f(x):
    return np.nan if x is None else float(x)


out = []
for cid, c in C.items():
    R = [r for r in rows if r["cond"] == cid]
    if not R: continue
    t = c["tau3"]; e = dict(id=cid, fR=c["fR"], Ns=c["Ns"], bw=c["bw"], tau3=t, I=c["I"], n=len(R),
                            neg_frac=float(np.mean([r["neg_frac"] for r in R])), nonident=any(r.get("nonident") for r in R))
    for est in ("pub", "noclip", "joint", "tm"):
        v = np.array([f(r.get(est)) for r in R])
        ok = np.isfinite(v)
        if ok.sum() == 0:
            e[est] = None; continue
        d = 100 * (v[ok] - t) / t
        s = dict(nvalid=int(ok.sum()), fail=float(1 - ok.mean()), bias=float(d.mean()),
                 mcse=float(d.std(ddof=1) / np.sqrt(ok.sum())) if ok.sum() > 1 else np.nan,
                 sd=float(d.std(ddof=1)) if ok.sum() > 1 else np.nan, rmse=float(np.sqrt(np.mean(d ** 2))))
        if est in ("pub", "joint"):
            B = [r.get(est + "_boot") for r in R if r.get(est + "_boot") is not None and np.isfinite(f(r.get(est)))]
            if B:
                B = np.array(B, float)
                s["cov68"] = float(np.mean((B[:, 2] <= t) & (t <= B[:, 3])))
                s["cov95"] = float(np.mean((B[:, 1] <= t) & (t <= B[:, 4])))
                s["sd_ratio"] = float(np.nanmean(B[:, 0]) / (s["sd"] * t / 100))
                s["nboot"] = int(len(B))
        if est == "tm":
            sg = np.array([f(r.get("tm_sig")) for r in R])[ok]
            s["cov68"] = float(np.mean(np.abs(v[ok] - t) <= sg))
            s["cov95"] = float(np.mean(np.abs(v[ok] - t) <= 1.96 * sg))
            s["sd_ratio"] = float(np.nanmean(sg) / (s["sd"] * t / 100))
        e[est] = s
    out.append(e)
json.dump(out, open(D / "summary.json", "w"), indent=1)


def fmt(s, key, p=1):
    if s is None or key not in s or s[key] is None or not np.isfinite(s[key]): return "--"
    return f"{s[key]:+.{p}f}" if key == "bias" else f"{s[key]:.{p}f}"


def cov(s):
    if s is None or "cov95" not in s: return "--"
    return f"{s['cov68']:.2f}/{s['cov95']:.2f}"


L = []
for e in out:
    lab = []
    if e["bw"] != 0.05: lab.append(f"bin {1000 * e['bw']:.0f} ps")
    if e["tau3"] != 2.0: lab.append(f"$\\tau_3={e['tau3']:.1f}$ ns")
    if e["I"] != [0.15, 0.62, 0.23]: lab.append("prior-mismatched $I$")
    ns = f"{e['Ns']:.0e}".replace("e+0", "\\times10^{") + "}"
    row = [f"{e['fR']:.1f}", f"${ns}$", "; ".join(lab) if lab else "base",
           f"{e['neg_frac']:.2f}"]
    for est in ("pub", "noclip", "joint", "tm"):
        s = e.get(est)
        row.append("n.i." if e["nonident"] and est != "tm" else (fmt(s, "bias", 2) if s else "--"))
    for est in ("pub", "joint"):
        s = e.get(est)
        row.append("--" if e["nonident"] else fmt(s, "sd", 2))
    for est in ("pub", "joint", "tm"):
        row.append(cov(e.get(est)))
    L.append(" & ".join(row) + " \\\\")
open(D / "tableS10_rows.tex", "w").write("\n".join(L) + "\n")
for e in out:
    print(e["id"], e["fR"], int(e["Ns"]), e["bw"], e["tau3"], "neg", round(e["neg_frac"], 3),
          *[(k, None if e.get(k) is None else {kk: (round(vv, 3) if isinstance(vv, float) else vv) for kk, vv in e[k].items()}) for k in ("pub", "noclip", "joint", "tm")], sep=" | ")
