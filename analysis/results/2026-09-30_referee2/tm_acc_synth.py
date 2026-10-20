#!/usr/bin/env python3
"""Referee round 2, MC11 companion: the windowed-moment (truncated-mean) estimator on the
synthetic spectra of estimator_synth.py (pure o-Ps EMG times the window acceptance plus flat
randoms, delayed-arm subtraction), read out two ways on the same replicates:
  sharp   the production inversion (sharp truncation on (1, 5) ns)
  accept  the acceptance-aware inversion: the binned expectation of the window mean,
          computed with the same EMG, bins and acceptance (sigma_e = 359 ps), inverted
          numerically.
200 replicates per condition. Bias, SD and delta-method 68/95 % coverage.
Writes 2026-09-30_referee2/estimator_synth/tm_acc.json."""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import estimator_synth as S  # noqa

OUT = S.OUT / "tm_acc.json"
NREP = 200


def expected_mean(tau, g):
    p = np.diff(S.emg_cdf(g.edges, tau, 0.0, S.SIG_D)) * g.accb
    w = p[g.T]
    return float((g.bc[g.T] * w).sum() / w.sum())


def run(cond):
    g = S.Grid(cond["bw"]); t3 = cond["tau3"]; fR = cond["fR"]; Ns = cond["Ns"]
    taus = np.linspace(0.8, 5.0, 841); mA = np.array([expected_mean(t, g) for t in taus])
    rng = np.random.default_rng([20260930, 77, cond["id"]])
    p1 = np.diff(S.emg_cdf(g.edges, t3, 0.0, S.SIG_D)) * g.accb
    mu1 = Ns * p1 / p1[g.T].sum()
    r1 = fR / (1 - fR) * Ns / g.T.sum() if fR > 0 else 0.0
    sharp, acc, sg_s, sg_a = [], [], [], []
    for _ in range(NREP):
        a_t = rng.poisson(mu1 + r1).astype(float); a_d = rng.poisson(np.full(len(mu1), r1)).astype(float)
        l, s = S.tm(a_t, a_d, g); sharp.append(l); sg_s.append(s)
        ht, hr = a_t[g.T], a_d[g.T]; hs = ht - hr; ns = hs.sum()
        m = (g.bc[g.T] * hs).sum() / ns
        la = float(np.interp(m, mA, taus, left=np.nan, right=np.nan))
        # delta method: sigma of the window mean propagated through d tau / d m
        sig_m = np.sqrt(((g.bc[g.T] - m) ** 2 * (ht + hr)).sum()) / ns
        dt = float(np.interp(m + 1e-3, mA, taus) - la) / 1e-3
        acc.append(la); sg_a.append(abs(dt) * sig_m)
    res = dict(id=cond["id"], fR=fR, Ns=Ns, bw=cond["bw"], tau3=t3)
    for k, v, sg in (("sharp", sharp, sg_s), ("accept", acc, sg_a)):
        v = np.array(v); sg = np.array(sg); ok = np.isfinite(v); d = 100 * (v[ok] - t3) / t3
        res[k] = dict(nvalid=int(ok.sum()), bias=float(d.mean()), mcse=float(d.std(ddof=1) / np.sqrt(ok.sum())),
                      sd=float(d.std(ddof=1)), rmse=float(np.sqrt(np.mean(d ** 2))),
                      cov68=float(np.mean(np.abs(v[ok] - t3) <= sg[ok])), cov95=float(np.mean(np.abs(v[ok] - t3) <= 1.96 * sg[ok])))
    return res


if __name__ == "__main__":
    conds = [c for c in S.build_conditions() if c["bw"] <= 0.5]
    rows = [run(c) for c in conds]
    json.dump(rows, open(OUT, "w"), indent=1)
    for r in rows:
        print(r["id"], r["fR"], int(r["Ns"]), r["bw"], r["tau3"],
              "sharp %+.2f (sd %.2f, cov %.2f/%.2f)" % (r["sharp"]["bias"], r["sharp"]["sd"], r["sharp"]["cov68"], r["sharp"]["cov95"]),
              "accept %+.2f+-%.2f (sd %.2f, cov %.2f/%.2f)" % (r["accept"]["bias"], r["accept"]["mcse"], r["accept"]["sd"], r["accept"]["cov68"], r["accept"]["cov95"]))
