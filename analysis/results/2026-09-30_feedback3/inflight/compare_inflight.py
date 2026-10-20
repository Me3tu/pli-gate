#!/usr/bin/env python3
"""In-flight sensitivity: hybrid arm (stock in-flight annihilation kept) against C2 (shipped
physics), identical otherwise. Scores the predictions of REGISTERED_INFLIGHT.md.
Per tube and pooled (inverse-variance), hybrid minus C2, with independent-arm errors:
  truth-grouped ('oracle') lifetime of sim1_extract, clinical truncated-mean lifetimes under both
  rules, the rule artifact Delta|b| = |b_closest| - |b_all| (truncated mean), the three-component
  rule shift and all-prompts bias under the joint likelihood and the clipped pseudo-likelihood,
  genuine triples, randoms fraction, phi, 511-window and prompt-window singles, and mu.
Writes compare_inflight.json."""
import os, sys, re, json, csv
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np, h5py
RES = Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[2])
D = RES / "2026-09-30_feedback3/inflight"; C2 = RES / "2026-08-24_C2_sim1_388"
sys.path.insert(0, str(RES / "2026-09-30_referee2")); sys.path.insert(0, str(RES / "2026-08-12_steinberger_allprompts"))
import comp_joint_prod as CJ  # noqa
import steinberger_extract_ap as X  # noqa
exec(compile(re.search(r"^def tau_from_hists\(.*?(?=^def |\Z)",
     (RES / "2026-08-19_V4_stats/v4d_paired.py").read_text(), re.S | re.M).group(0), "f", "exec"))
F = CJ.F
ARMS = {"c2": dict(npz=RES / "2026-08-26_B2_paired/hists_sim1_388.npz", out=C2, ov=C2 / "oversub_sim1_{}.json",
                   singles=Path(os.environ.get("PLI_RAW", "<missing>")) / "c2_sim1_388/singles_c2_sim1_388.h5"),
        "hybrid": dict(npz=D / "hists_hybrid_388.npz", out=D, ov=D / "oversub_{}.json",
                       singles=Path(os.environ.get("PLI_RAW", "<missing>")) / "feedback3_inflight/c2_hybrid/singles_c2_hybrid.h5")}
if not all(Path(ARMS[a]["singles"]).exists() for a in ARMS): raise SystemExit("singles of " + ", ".join(ARMS) + ": raw data not included in the release, see analysis/README.md")
B = 200


def est_tube(args):
    arm, ki, name, tau = args
    z = np.load(ARMS[arm]["npz"], allow_pickle=True)
    H = {r: (z[f"{r}_sig"][:, ki, :], z[f"{r}_rnd"][:, ki, :]) for r in ("closest", "all")}
    def one(S, Rr, x0=None):
        nt, nd = S.sum(0), Rr.sum(0)
        tm = tau_from_hists(nt, nd)
        pub = F.fit3(F.spectrum(nt, nd)) if x0 is None else F.fit3(F.spectrum(nt, nd), x0=x0[0])
        tj, xj = CJ.fit_joint(nt[F.M], nd[F.M], x0=None if x0 is None else x0[1])
        return np.array([tm, pub, tj]), (pub, xj)
    pt, x0 = {}, {}
    for r in H:
        pt[r], x0[r] = one(*H[r])
    b = lambda v: 100 * (v - tau) / tau
    def quantities(p):
        bc, ba = b(p["closest"]), b(p["all"])
        return dict(tm_all=ba[0], tm_closest=bc[0], did_tm=abs(bc[0]) - abs(ba[0]),
                    pub_all=ba[1], joint_all=ba[2], shift_pub=ba[1] - bc[1], shift_joint=ba[2] - bc[2])
    q0 = quantities(pt)
    ns = H["all"][0].shape[0]; rng = np.random.default_rng(777 + ki)
    bq = {k: np.full(B, np.nan) for k in q0}
    for i in range(B):
        k = rng.integers(0, ns, ns)
        p = {r: one(H[r][0][k], H[r][1][k], x0=x0[r])[0] for r in H}
        for kk, v in quantities(p).items(): bq[kk][i] = v
    sd = {kk: float(0.5 * np.subtract(*np.nanpercentile(v, [84, 16]))) for kk, v in bq.items()}
    return dict(arm=arm, region=name, tau_cfg=tau, q={k: float(v) for k, v in q0.items()}, sd=sd)


def per_voi(arm, rule):
    return {r["voi"]: r for r in csv.DictReader(open(ARMS[arm]["out"] / f"out_{rule}" / "per_voi.csv"))}


def singles_stats(path):
    with h5py.File(path, "r") as h:
        g = h["singles"]; sh = g["shard_id"][:]; t = g["time_ns"][:]
        n5 = int(g["in_511"][:].sum()); npr = int(g["in_prompt"][:].sum()); n = int(sh.size)
        span = 0.0
        for sid in np.unique(sh):
            m = sh == sid; tt = t[m]; span += float(tt.max() - tt.min())
    return dict(n_singles=n, n_511=n5, n_prompt=npr, span_ns=span, mu=npr / span * 10.0)


if __name__ == "__main__":
    z = np.load(ARMS["c2"]["npz"], allow_pickle=True)
    names = [str(s) for s in z["voi_names"]]; cfg = z["tau_cfg"]
    tubes = [(i, n, float(cfg[i])) for i, n in enumerate(names) if n != "bg" and np.isfinite(cfg[i])]
    jobs = [(arm, i, n, t) for arm in ARMS for (i, n, t) in tubes]
    with ProcessPoolExecutor(max_workers=16) as ex:
        est = list(ex.map(est_tube, jobs))
    E = {(e["arm"], e["region"]): e for e in est}
    out = dict(tubes=[], pooled={}, singles={})
    for arm in ARMS:
        out["singles"][arm] = singles_stats(ARMS[arm]["singles"])
    rows = []
    for i, n, tau in tubes:
        r = dict(region=n, tau_cfg=tau)
        for arm in ARMS:
            pa, pc = per_voi(arm, "all")[n], per_voi(arm, "closest")[n]
            ova = {x["voi"]: x for x in json.load(open(str(ARMS[arm]["ov"]).format("all")))["rows"]}[n]
            ovc = {x["voi"]: x for x in json.load(open(str(ARMS[arm]["ov"]).format("closest")))["rows"]}[n]
            e = E[(arm, n)]
            r[arm] = dict(oracle_tau=float(pa["oracle_tau_hat_ns"]), oracle_sig=float(pa["oracle_sigma_ns"]),
                          oracle_n=float(pa["oracle_n_signal"]), n_total_all=float(pa["clinical_n_total"]),
                          n_random_all=float(pa["clinical_n_random"]), n_signal_all=float(pa["clinical_n_signal"]),
                          fR_all=float(pa["clinical_n_random"]) / float(pa["clinical_n_total"]),
                          n_true_sig=float(ova["n_true_sig"]), phi_all=float(ova["phi"]), phi_closest=float(ovc["phi"]),
                          **{k: v for k, v in e["q"].items()}, **{f"sd_{k}": v for k, v in e["sd"].items()})
        h, c = r["hybrid"], r["c2"]
        r["d"] = dict(oracle_pct=100 * (h["oracle_tau"] - c["oracle_tau"]) / tau,
                      oracle_pct_se=100 * np.hypot(h["oracle_sig"], c["oracle_sig"]) / tau,
                      oracle_n_rel=h["oracle_n"] / c["oracle_n"] - 1, true_sig_rel=h["n_true_sig"] / c["n_true_sig"] - 1,
                      fR_all_diff=h["fR_all"] - c["fR_all"])
        for k in ("tm_all", "tm_closest", "did_tm", "pub_all", "joint_all", "shift_pub", "shift_joint"):
            r["d"][k] = h[k] - c[k]; r["d"][k + "_se"] = float(np.hypot(h["sd_" + k], c["sd_" + k]))
        rows.append(r)
        print(f"{n:28s} oracle {r['d']['oracle_pct']:+.2f}+-{r['d']['oracle_pct_se']:.2f}%  TM all {r['d']['tm_all']:+.2f}+-{r['d']['tm_all_se']:.2f}"
              f"  closest {r['d']['tm_closest']:+.2f}+-{r['d']['tm_closest_se']:.2f}  did {r['d']['did_tm']:+.2f}+-{r['d']['did_tm_se']:.2f}"
              f"  joint shift {r['d']['shift_joint']:+.2f}+-{r['d']['shift_joint_se']:.2f}  joint all {r['d']['joint_all']:+.2f}+-{r['d']['joint_all_se']:.2f}"
              f"  genuine {100*r['d']['true_sig_rel']:+.2f}%  oracle_n {100*r['d']['oracle_n_rel']:+.2f}%  fR {r['d']['fR_all_diff']:+.4f}")
    out["tubes"] = rows
    for k in ("oracle_pct", "tm_all", "tm_closest", "did_tm", "pub_all", "joint_all", "shift_pub", "shift_joint"):
        v = np.array([r["d"][k] for r in rows]); s = np.array([r["d"][(k + "_se") if k != "oracle_pct" else "oracle_pct_se"] for r in rows])
        w = 1 / s ** 2; m = float((w * v).sum() / w.sum()); se = float(np.sqrt(1 / w.sum()))
        chi2 = float((w * (v - m) ** 2).sum())
        out["pooled"][k] = dict(mean=m, se=se, z=m / se, chi2=chi2, df=len(v) - 1)
        print(f"POOLED {k:12s} {m:+.3f} +- {se:.3f}  z {m/se:+.2f}  chi2 {chi2:.1f}/{len(v)-1}")
    for k in ("true_sig_rel", "oracle_n_rel"):
        v = np.array([r["d"][k] for r in rows]); out["pooled"][k] = dict(mean=float(v.mean()), min=float(v.min()), max=float(v.max()))
        print(f"{k}: mean {100*v.mean():+.2f}% range {100*v.min():+.2f} to {100*v.max():+.2f}%")
    sh, sc = out["singles"]["hybrid"], out["singles"]["c2"]
    out["singles"]["rel"] = {k: sh[k] / sc[k] - 1 for k in ("n_singles", "n_511", "n_prompt", "mu")}
    print("singles rel (hybrid/c2 - 1):", {k: round(100 * v, 3) for k, v in out["singles"]["rel"].items()}, "%")
    P = out["pooled"]; S = out["singles"]["rel"]
    out["scores"] = {
        "P1 511 singles fall 1-3%": bool(-0.03 <= S["n_511"] <= -0.01),
        "P2 prompt singles |change| < 0.5%": bool(abs(S["n_prompt"]) < 0.005),
        "P3 genuine triples fall 1.5-3.5%": bool(-0.035 <= P["true_sig_rel"]["mean"] <= -0.015),
        "P4 lifetimes pooled |z| < 2 (oracle, TM all, TM closest)": bool(all(abs(P[k]["z"]) < 2 for k in ("oracle_pct", "tm_all", "tm_closest"))),
        "P5 Delta|b| pooled |z| < 2": bool(abs(P["did_tm"]["z"]) < 2),
        "P6 mu |change| < 0.5%": bool(abs(S["mu"]) < 0.005)}
    for k, v in out["scores"].items(): print(("PASS " if v else "FAIL ") + k)
    json.dump(out, open(D / "compare_inflight.json", "w"), indent=1)
