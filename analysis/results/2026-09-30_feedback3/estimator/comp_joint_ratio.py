#!/usr/bin/env python3
"""Feedback round 3: the count-level joint likelihood as the primary component estimator.
For the eight-tube 388 ps histograms (2026-08-26_B2_paired/hists_sim1_388.npz), per tube and on the
SAME shard resamples: rule shift 100 (tau_all - tau_closest)/tau_cfg under the truncated mean (TM,
tau_from_hists of v4d_paired), the published pseudo-likelihood (PUB, fit3_sim1.fit3) and the joint
likelihood (JOINT, comp_joint_prod.fit_joint). Pooled fixed-effect (IVW) shifts, the JOINT/TM and
PUB/TM ratios with paired-bootstrap percentile intervals, and for JOINT the heterogeneity statistics
(Cochran Q, I^2, DerSimonian-Laird tau^2 and mean, Knapp-Hartung interval). B is a constant here
(Python 3.14 forkserver workers do not inherit sys.argv). Writes comp_joint_ratio.json."""
import os, sys, re, json, time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy.stats import t as tdist
RES = Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[2])
sys.path.insert(0, str(RES / "2026-09-30_referee2"))
sys.path.insert(0, str(RES / "2026-08-12_steinberger_allprompts"))
import comp_joint_prod as CJ  # noqa  (fit_joint, F = fit3_sim1)
import steinberger_extract_ap as X  # noqa
exec(compile(re.search(r"^def tau_from_hists\(.*?(?=^def |\Z)",
     (RES / "2026-08-19_V4_stats/v4d_paired.py").read_text(), re.S | re.M).group(0), "f", "exec"))
F = CJ.F
NPZ = RES / "2026-08-26_B2_paired/hists_sim1_388.npz"
OUT = Path(__file__).with_name("comp_joint_ratio.json")
B = 400


def one_tube(args):
    ki, name, tau = args
    z = np.load(NPZ, allow_pickle=True)
    arms = {r: (z[f"{r}_sig"][:, ki, :], z[f"{r}_rnd"][:, ki, :]) for r in ("closest", "all")}
    def est(nt, nd, x0p=None, x0j=None):
        tm = tau_from_hists(nt, nd)
        pub = F.fit3(F.spectrum(nt, nd), x0=x0p) if x0p is not None else F.fit3(F.spectrum(nt, nd))
        tj, xj = CJ.fit_joint(nt[F.M], nd[F.M], x0=x0j)
        return tm, pub, tj, xj
    pt, x0 = {}, {}
    for r, (S, Rr) in arms.items():
        tm, pub, tj, xj = est(S.sum(0), Rr.sum(0)); pt[r] = (tm, pub, tj); x0[r] = (pub, xj)
    shift = lambda a, c: 100 * (a - c) / tau
    out = dict(region=name, tau_cfg=tau,
               point={e: shift(pt["all"][i], pt["closest"][i]) for i, e in enumerate(("tm", "pub", "joint"))},
               bias_all={e: 100 * (pt["all"][i] - tau) / tau for i, e in enumerate(("tm", "pub", "joint"))},
               bias_closest={e: 100 * (pt["closest"][i] - tau) / tau for i, e in enumerate(("tm", "pub", "joint"))})
    ns = arms["all"][0].shape[0]
    rng = np.random.default_rng(20261001 + ki)
    bs = {e: np.full(B, np.nan) for e in ("tm", "pub", "joint")}
    ba = {e: np.full(B, np.nan) for e in ("tm", "pub", "joint")}
    for b in range(B):
        k = rng.integers(0, ns, ns)
        v = {}
        for r, (S, Rr) in arms.items():
            tm, pub, tj, _ = est(S[k].sum(0), Rr[k].sum(0), x0p=x0[r][0], x0j=x0[r][1]); v[r] = (tm, pub, tj)
        for i, e in enumerate(("tm", "pub", "joint")):
            bs[e][b] = shift(v["all"][i], v["closest"][i]); ba[e][b] = 100 * (v["all"][i] - tau) / tau
    out["boot_shift"] = {e: bs[e].tolist() for e in bs}
    out["boot_bias_all"] = {e: ba[e].tolist() for e in ba}
    out["sig_shift"] = {e: float(0.5 * np.subtract(*np.nanpercentile(bs[e], [84, 16]))) for e in bs}
    out["sig_bias_all"] = {e: float(0.5 * np.subtract(*np.nanpercentile(ba[e], [84, 16]))) for e in ba}
    return out


def hetero(y, s):
    w = 1 / s ** 2; k = len(y); mu = (w * y).sum() / w.sum(); Q = float((w * (y - mu) ** 2).sum())
    I2 = max(0.0, (Q - (k - 1)) / Q) if Q > 0 else 0.0
    C = w.sum() - (w ** 2).sum() / w.sum(); t2 = max(0.0, (Q - (k - 1)) / C)
    ws = 1 / (s ** 2 + t2); mre = (ws * y).sum() / ws.sum(); se_re = np.sqrt(1 / ws.sum())
    qkh = (ws * (y - mre) ** 2).sum() / (k - 1); se_kh = np.sqrt(qkh / ws.sum()); tq = tdist.ppf(0.975, k - 1)
    return dict(fe=float(mu), fe_se=float(np.sqrt(1 / w.sum())), Q=Q, df=k - 1, I2=float(I2), tau2=float(t2),
                re=float(mre), re_se=float(se_re), kh95=[float(mre - tq * se_kh), float(mre + tq * se_kh)])


if __name__ == "__main__":
    z = np.load(NPZ, allow_pickle=True)
    names = [str(s) for s in z["voi_names"]]; cfg = z["tau_cfg"]
    jobs = [(i, n, float(cfg[i])) for i, n in enumerate(names) if n != "bg" and np.isfinite(cfg[i])]
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=8) as ex:
        rows = list(ex.map(one_tube, jobs))
    rows.sort(key=lambda r: r["tau_cfg"])
    res = dict(B=B, rows=[{k: v for k, v in r.items() if not k.startswith("boot")} for r in rows])
    for e in ("tm", "pub", "joint"):
        y = np.array([r["point"][e] for r in rows]); s = np.array([r["sig_shift"][e] for r in rows])
        res[f"hetero_{e}"] = hetero(y, s)
        res[f"range_{e}"] = [float(y.min()), float(y.max())]
        a = np.array([r["bias_all"][e] for r in rows]); sa = np.array([r["sig_bias_all"][e] for r in rows])
        res[f"all_bias_{e}"] = dict(mean=float(a.mean()), median=float(np.median(a)), n_within_2p3=int((np.abs(a) <= 2.3).sum()),
                                   ivw=float((a / sa ** 2).sum() / (1 / sa ** 2).sum()), ivw_se=float(np.sqrt(1 / (1 / sa ** 2).sum())),
                                   per_tube={r["region"]: float(r["bias_all"][e]) for r in rows})
    # paired ratios of the pooled shifts: fixed IVW weights from the per-tube bootstrap sigmas
    for e in ("pub", "joint"):
        wt = np.array([1 / r["sig_shift"]["tm"] ** 2 for r in rows]); we = np.array([1 / r["sig_shift"][e] ** 2 for r in rows])
        bt = np.array([r["boot_shift"]["tm"] for r in rows]); be = np.array([r["boot_shift"][e] for r in rows])
        num = (we[:, None] * be).sum(0) / we.sum(); den = (wt[:, None] * bt).sum(0) / wt.sum()
        ratio = num / den; ok = np.isfinite(ratio)
        res[f"ratio_{e}_over_tm"] = dict(point=res[f"hetero_{e}"]["fe"] / res["hetero_tm"]["fe"],
                                         ci95=[float(x) for x in np.percentile(ratio[ok], [2.5, 97.5])], n_ok=int(ok.sum()))
    for e in ("tm", "pub", "joint"):
        h = res[f"hetero_{e}"]
        print(f"{e:5s} FE {h['fe']:.2f}+-{h['fe_se']:.2f}  range {res[f'range_{e}'][0]:.1f}-{res[f'range_{e}'][1]:.1f}  Q {h['Q']:.1f}/{h['df']} I2 {100*h['I2']:.1f}%  "
              f"RE {h['re']:.2f}+-{h['re_se']:.2f} KH [{h['kh95'][0]:.1f},{h['kh95'][1]:.1f}]  all-prompts bias mean {res[f'all_bias_{e}']['mean']:+.2f} "
              f"median {res[f'all_bias_{e}']['median']:+.2f} within2.3 {res[f'all_bias_{e}']['n_within_2p3']}")
    for e in ("pub", "joint"):
        print(e, "/TM ratio", res[f"ratio_{e}_over_tm"])
    json.dump(res, open(OUT, "w"), indent=1)
    print(f"done {time.time()-t0:.0f} s")
