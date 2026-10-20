#!/usr/bin/env python3
"""Referee round 2, MC11 follow-up: the count-level joint likelihood applied to the
production eight-tube 388 ps histograms (2026-08-26_B2_paired/hists_sim1_388.npz).

The synthetic study (estimator_synth.py) showed that clipping the negative bins of the
randoms-subtracted histogram biases the published three-component estimator low once
randoms are present, and that the joint likelihood of the prompt and delayed arms (per-bin
randoms mean profiled analytically, no clipping) removes that bias. Here the same three
objectives are applied to the production data, per tube and per matching rule, with the
component shape, fixed lifetimes, prior and Sim-1 acceptance of fit3_sim1.py unchanged:
  PUB     clip(n_t - n_d, 0), multinomial Poisson objective (the published estimator)
  NOCLIP  the same objective without clipping
  JOINT   n_t ~ Poisson(s_b + r_b), n_d ~ Poisson(r_b), r_b profiled, s_b = S p_b(theta)
Paired shard bootstrap (B draws, both rules and both estimators refitted on the same
resample) gives the uncertainty of every bias and of the rule shift.
Writes 2026-09-30_referee2/comp_joint_prod.json."""
from __future__ import annotations
import os, sys, json, time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy.optimize import minimize

RES = Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[1])
sys.path.insert(0, str(RES / "2026-08-26_A2_cond2"))
import fit3_sim1 as F  # noqa  (EDGES, ACCB, M, fit3, _nll, E, G)

NPZ = RES / "2026-08-26_B2_paired/hists_sim1_388.npz"
OUT = RES / "2026-09-30_referee2/comp_joint_prod.json"
B = int(sys.argv[1]) if len(sys.argv) > 1 else 200


def shape(th):
    w = np.exp(np.array([th[3], th[4], 0.0])); I = w / w.sum()
    p = np.zeros(len(F.EDGES) - 1)
    for Ik, tk in zip(I, (F.E.TAU_PPS, F.E.TAU_DIR, th[0])):
        p += Ik * np.diff(F.E.emg_cdf(F.EDGES, tk, th[1], th[2]))
    p = np.clip(p * F.ACCB, 1e-15, None)
    return p / p.sum(), I


def nll_joint(th, nt, nd):
    if not (0.3 < th[0] < 8.0 and 0.01 < th[2] < 1.0 and abs(th[1]) < 2.0 and -5 < th[5] < 25):
        return 1e12
    p, I = shape(th)
    s = np.exp(th[5]) * p
    b = 2 * s - nt - nd
    r = np.clip((-b + np.sqrt(b * b + 8 * nd * s)) / 4.0, 1e-12, None)
    m = s + r
    return float(np.sum(m - nt * np.log(m)) + np.sum(r - nd * np.log(r))) - F.G.dirichlet_logpdf(I, F.E.ALPHA)


def fit_joint(nt, nd, x0=None):
    ext = np.log(max((nt - nd).sum(), 10.0))
    inits = [x0] if x0 is not None else [np.array([s, 0.0, .0909, np.log(.15 / .23), np.log(.62 / .23), ext]) for s in F.STARTS]
    best = None
    for st in inits:
        r = minimize(nll_joint, st, args=(nt, nd), method="Nelder-Mead",
                     options=dict(maxiter=40000, xatol=1e-7, fatol=1e-7))
        if best is None or r.fun < best.fun: best = r
    t = float(best.x[0])
    return (t if 0.35 < t < 7.5 else float("nan")), best.x


def fit_noclip(nt, nd):
    c = (nt - nd)[F.M]
    best = None
    for s in F.STARTS:
        st = np.array([s, 0.0, .0909, np.log(.15 / .23), np.log(.62 / .23)])
        r = minimize(F._nll, st, args=(c,), method="Nelder-Mead", options=dict(maxiter=30000, xatol=1e-7, fatol=1e-7))
        if best is None or r.fun < best.fun: best = r
    t = float(best.x[0]); return t if 0.35 < t < 7.5 else float("nan")


def one_tube(args):
    ki, name, tau = args
    z = np.load(NPZ, allow_pickle=True)
    arms = {"closest": (z["closest_sig"][:, ki, :], z["closest_rnd"][:, ki, :]),
            "all": (z["all_sig"][:, ki, :], z["all_rnd"][:, ki, :])}
    out = dict(region=name, tau_cfg=tau)
    pts, xj = {}, {}
    for rule, (S, R) in arms.items():
        nt, nd = S.sum(0), R.sum(0)
        ntM, ndM = nt[F.M], nd[F.M]
        out[rule] = dict(integer_arms=bool(np.allclose(S, np.round(S)) and np.allclose(R, np.round(R))),
                         n_prompt=float(ntM.sum()), n_delayed=float(ndM.sum()),
                         f_R=float(ndM.sum() / ntM.sum()), N_s=float((ntM - ndM).sum()),
                         neg_frac=float(((ntM - ndM) < 0).mean()))
        t_pub = F.fit3(F.spectrum(nt, nd)); t_nc = fit_noclip(nt, nd); t_j, x = fit_joint(ntM, ndM)
        pts[rule] = dict(pub=t_pub, noclip=t_nc, joint=t_j); xj[rule] = x
        out[rule].update(tau_pub=t_pub, tau_noclip=t_nc, tau_joint=t_j,
                         bias_pub_pct=100 * (t_pub - tau) / tau, bias_noclip_pct=100 * (t_nc - tau) / tau,
                         bias_joint_pct=100 * (t_j - tau) / tau)
    ns = arms["all"][0].shape[0]
    rng = np.random.default_rng(20260930 + ki)
    bt = {(r, e): np.full(B, np.nan) for r in arms for e in ("pub", "joint")}
    for b in range(B):
        k = rng.integers(0, ns, ns)
        for rule, (S, R) in arms.items():
            nt, nd = S[k].sum(0), R[k].sum(0)
            bt[(rule, "pub")][b] = F.fit3(F.spectrum(nt, nd), x0=pts[rule]["pub"])
            bt[(rule, "joint")][b] = fit_joint(nt[F.M], nd[F.M], x0=xj[rule])[0]
    for e in ("pub", "joint"):
        for rule in arms:
            v = bt[(rule, e)]; ok = np.isfinite(v); q = np.percentile(100 * (v[ok] - tau) / tau, [16, 84])
            out[rule][f"sig_{e}_pct"] = float(0.5 * (q[1] - q[0]))
        d = 100 * (bt[("all", e)] - bt[("closest", e)]) / tau; ok = np.isfinite(d); q = np.percentile(d[ok], [16, 84])
        out[f"shift_{e}_pp"] = 100 * (pts["all"][e] - pts["closest"][e]) / tau
        out[f"shift_{e}_sig_pp"] = float(0.5 * (q[1] - q[0])); out[f"n_nan_{e}"] = int((~ok).sum())
    out["shift_noclip_pp"] = 100 * (pts["all"]["noclip"] - pts["closest"]["noclip"]) / tau
    return out


if __name__ == "__main__":
    z = np.load(NPZ, allow_pickle=True)
    names = [str(s) for s in z["voi_names"]]; cfg = z["tau_cfg"]
    jobs = [(i, n, float(cfg[i])) for i, n in enumerate(names) if n != "bg" and np.isfinite(cfg[i])]
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=8) as ex:
        rows = list(ex.map(one_tube, jobs))
    rows.sort(key=lambda r: r["tau_cfg"])
    for r in rows:
        print(f"{r['region']:28} cfg {r['tau_cfg']:.3f}  fR(all) {r['all']['f_R']:.3f}  Ns {r['all']['N_s']:.0f}  neg {r['all']['neg_frac']:.2f}"
              f" | all: pub {r['all']['bias_pub_pct']:+.2f} noclip {r['all']['bias_noclip_pct']:+.2f} joint {r['all']['bias_joint_pct']:+.2f}+-{r['all']['sig_joint_pct']:.2f}"
              f" | closest: pub {r['closest']['bias_pub_pct']:+.2f} joint {r['closest']['bias_joint_pct']:+.2f}"
              f" | shift pub {r['shift_pub_pp']:.2f}+-{r['shift_pub_sig_pp']:.2f} joint {r['shift_joint_pp']:.2f}+-{r['shift_joint_sig_pp']:.2f}")
    res = dict(B=B, rows=rows)
    for e in ("pub", "joint"):
        s = np.array([r[f"shift_{e}_pp"] for r in rows]); w = 1 / np.array([r[f"shift_{e}_sig_pp"] for r in rows]) ** 2
        res[f"ivw_shift_{e}_pp"] = float((w * s).sum() / w.sum()); res[f"ivw_shift_{e}_err_pp"] = float(np.sqrt(1 / w.sum()))
        a = np.array([r["all"][f"bias_{e}_pct"] for r in rows])
        res[f"all_bias_{e}_mean_pct"] = float(a.mean()); res[f"all_bias_{e}_n_within_2p3"] = int((np.abs(a) <= 2.3).sum())
        print(e, "IVW shift", round(res[f"ivw_shift_{e}_pp"], 3), "+-", round(res[f"ivw_shift_{e}_err_pp"], 3),
              "all-prompts mean bias", round(a.mean(), 2), "within 2.3%:", int((np.abs(a) <= 2.3).sum()))
    json.dump(res, open(OUT, "w"), indent=1)
    print(f"done {time.time() - t0:.0f} s")
