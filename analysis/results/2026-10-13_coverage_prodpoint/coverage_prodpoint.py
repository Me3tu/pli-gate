#!/usr/bin/env python3
"""Interval coverage of the lifetime estimators at the production operating point.

Extends the synthetic study of Note S11 (2026-09-30_referee2/estimator_synth.py, imported unchanged as S),
whose randoms-fraction grid starts at f_R = 0.1, to the operating point of the production component fits
(C2 all-prompts arm, comp_joint_prod.json: f_R = 0.074 to 0.077 in (-1, 4.5) ns, N_s = 3.0 to 3.2e5 signal
triples per tube, configured lifetimes 1.77 to 2.72 ns), and adds shard-bootstrap intervals for the
windowed-moment readout at its own production operating point (C0 all-prompts arm, fit window (1, 5) ns,
measured per tube by moment_window_point.py).

Conditions (same generator, binning, acceptance and idealized subtraction as Note S11):
  100  component: f_R 0.075, N_s 3.2e5, tau3 2.0 ns
  101  component: f_R 0.075, N_s 3.2e5, tau3 2.7 ns
  102  moment:    f_R and N_s of the moment window from moment_window_point.json, tau3 2.0 ns
Per replicate: PUB (clipped, prespecified), NOCLIP and JOINT point fits with the four starts of Note S11,
shard bootstrap (32 shards, B = 200, warm-started refits) for PUB and JOINT, moment readout sharp (TM) and
acceptance-aware (TM_A) with delta-method sigma and with a shard bootstrap (B = 200).
Interval types: bootstrap percentile (as Note S11) and normal intervals +-z * bootstrap SD (as the quoted
uncertainties of the paper). 400 replicates per condition.
Writes reps.json and summary.json next to this script."""
import json, os, sys, time
from pathlib import Path
from multiprocessing import Pool
import numpy as np
sys.path.insert(0, str(Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[1]) / "2026-09-30_referee2"))
import estimator_synth as S  # noqa: E402

OUT = Path(__file__).resolve().parent
NREP = int(sys.argv[1]) if len(sys.argv) > 1 else 400
NB = 200
NSH = S.NSH
BASE_I = [0.15, 0.62, 0.23]


def conditions():
    mw = json.load(open(OUT / "moment_window_point.json"))
    return [dict(id=100, kind="component", fR=0.075, Ns=3.2e5, bw=0.05, tau3=2.0, I=BASE_I),
            dict(id=101, kind="component", fR=0.075, Ns=3.2e5, bw=0.05, tau3=2.7, I=BASE_I),
            dict(id=102, kind="moment", fR=round(mw["fR_median"], 3), Ns=float(round(mw["Ns_median"], -3)), bw=0.05, tau3=2.0, I=BASE_I)]


def mean_table(g):
    taus = np.linspace(0.8, 5.0, 841)
    mA = []
    for t in taus:
        p = np.diff(S.emg_cdf(g.edges, t, 0.0, S.SIG_D)) * g.accb
        w = p[g.T]; mA.append(float((g.bc[g.T] * w).sum() / w.sum()))
    return taus, np.array(mA)


def tm_both(at, ad, g, taus, mA):
    lam, sig = S.tm(at, ad, g)
    ht, hr = at[g.T], ad[g.T]; hs = ht - hr; ns = hs.sum()
    m = (g.bc[g.T] * hs).sum() / ns
    la = float(np.interp(m, mA, taus, left=np.nan, right=np.nan))
    sig_m = np.sqrt(((g.bc[g.T] - m) ** 2 * (ht + hr)).sum()) / ns
    dt = float(np.interp(m + 1e-3, mA, taus) - la) / 1e-3
    return lam, sig, la, abs(dt) * sig_m


def one_rep(args):
    cond, rep = args
    g = S.Grid(cond["bw"]); fR, Ns, t3 = cond["fR"], cond["Ns"], cond["tau3"]
    rng = np.random.default_rng([20261013, cond["id"], rep])
    res = dict(cond=cond["id"], rep=rep)
    if cond["kind"] == "component":
        mu_s = S.truth_counts(g, t3, cond["I"], Ns)
        r_b = fR / (1 - fR) * Ns / g.M.sum()
        sh_t = rng.poisson(np.broadcast_to((mu_s + r_b) / NSH, (NSH, len(mu_s)))).astype(float)
        sh_d = rng.poisson(np.broadcast_to(np.full(len(mu_s), r_b) / NSH, (NSH, len(mu_s)))).astype(float)
        nt, nd = sh_t.sum(0), sh_d.sum(0); diff = (nt - nd)[g.M]
        res["neg_frac"] = float((diff < 0).mean())
        t_pub, x_pub = S.fit(S.nll_pub, (np.clip(diff, 0, None), g), S.STARTS)
        t_nc, _ = S.fit(S.nll_pub, (diff, g), S.STARTS)
        t_j, x_j = S.fit(S.nll_joint, (nt[g.M], nd[g.M], g), S.STARTS, extra=np.log(max(diff.sum(), 10.0)))
        res.update(pub=t_pub, noclip=t_nc, joint=t_j)
        bp, bj = [], []
        for _ in range(NB):
            k = rng.integers(0, NSH, NSH); bt, bd = sh_t[k].sum(0), sh_d[k].sum(0)
            bp.append(S.fit(S.nll_pub, (np.clip((bt - bd)[g.M], 0, None), g), [x_pub])[0])
            bj.append(S.fit(S.nll_joint, (bt[g.M], bd[g.M], g), [x_j])[0])
        res["pub_boot"] = bp; res["joint_boot"] = bj
    else:
        taus, mA = mean_table(g)
        p1 = np.diff(S.emg_cdf(g.edges, t3, 0.0, S.SIG_D)) * g.accb
        mu1 = Ns * p1 / p1[g.T].sum(); r1 = fR / (1 - fR) * Ns / g.T.sum()
        sh_t = rng.poisson(np.broadcast_to((mu1 + r1) / NSH, (NSH, len(mu1)))).astype(float)
        sh_d = rng.poisson(np.broadcast_to(np.full(len(mu1), r1) / NSH, (NSH, len(mu1)))).astype(float)
        l, s, la, sa = tm_both(sh_t.sum(0), sh_d.sum(0), g, taus, mA)
        res.update(tm=l, tm_sig=s, tma=la, tma_sig=sa)
        bs, ba = [], []
        for _ in range(NB):
            k = rng.integers(0, NSH, NSH)
            l2, _, la2, _ = tm_both(sh_t[k].sum(0), sh_d[k].sum(0), g, taus, mA)
            bs.append(l2); ba.append(la2)
        res["tm_boot"] = bs; res["tma_boot"] = ba
    return res


def summarize(rows, C):
    out = []
    for c in C:
        R = [r for r in rows if r["cond"] == c["id"]]; t = c["tau3"]
        e = dict(c, n=len(R))
        if c["kind"] == "component":
            e["neg_frac"] = float(np.mean([r["neg_frac"] for r in R]))
            ests = (("pub", "pub_boot"), ("noclip", None), ("joint", "joint_boot"))
        else:
            ests = (("tm", "tm_boot"), ("tma", "tma_boot"))
        for est, bk in ests:
            v = np.array([r[est] for r in R], float); ok = np.isfinite(v); d = 100 * (v[ok] - t) / t
            s = dict(nvalid=int(ok.sum()), bias=float(d.mean()), mcse=float(d.std(ddof=1) / np.sqrt(ok.sum())),
                     sd=float(d.std(ddof=1)))
            if bk:
                B = np.array([r[bk] for r in R], float)[ok]; vv = v[ok]
                bsd = np.nanstd(B, axis=1, ddof=1)
                q = np.nanpercentile(B, [2.5, 16, 84, 97.5], axis=1)
                s["cov68_pct"] = float(np.mean((q[1] <= t) & (t <= q[2])))
                s["cov95_pct"] = float(np.mean((q[0] <= t) & (t <= q[3])))
                s["cov68_norm"] = float(np.mean(np.abs(vv - t) <= bsd))
                s["cov95_norm"] = float(np.mean(np.abs(vv - t) <= 1.96 * bsd))
                s["sd_ratio"] = float(np.nanmean(bsd) / (s["sd"] * t / 100))
                # factor by which the bootstrap SD must be scaled for 95% normal coverage
                z = np.abs(vv - t) / bsd
                s["k95"] = float(np.quantile(z, 0.95) / 1.96)
                s["mcse_cov95"] = float(np.sqrt(0.95 * 0.05 / ok.sum()))
            if est in ("tm", "tma"):
                sg = np.array([r[est + "_sig"] for r in R], float)[ok]
                s["cov68_delta"] = float(np.mean(np.abs(v[ok] - t) <= sg))
                s["cov95_delta"] = float(np.mean(np.abs(v[ok] - t) <= 1.96 * sg))
            e[est] = s
        out.append(e)
    return out


if __name__ == "__main__":
    C = conditions()
    if "--summarize" in sys.argv:
        rows = json.load(open(OUT / "reps.json"))
    else:
        json.dump(C, open(OUT / "conditions.json", "w"), indent=1)
        jobs = [(c, r) for c in C for r in range(NREP)]
        t0 = time.time(); rows = []
        with Pool(30) as pool:
            for i, res in enumerate(pool.imap_unordered(one_rep, jobs, chunksize=1)):
                rows.append(res)
                if i % 100 == 0:
                    print(f"  {i}/{len(jobs)}  {time.time() - t0:.0f} s", flush=True)
                    json.dump(rows, open(OUT / "reps_partial.json", "w"), default=float)
        json.dump(rows, open(OUT / "reps.json", "w"), default=float)
        print(f"done {time.time() - t0:.0f} s", flush=True)
    S_ = summarize(rows, C)
    json.dump(S_, open(OUT / "summary.json", "w"), indent=1)
    for e in S_:
        print(json.dumps({k: (round(v, 4) if isinstance(v, float) else v) for k, v in e.items() if not isinstance(v, dict)}))
        for k, v in e.items():
            if isinstance(v, dict): print("  ", k, {kk: round(vv, 4) if isinstance(vv, float) else vv for kk, vv in v.items()})
