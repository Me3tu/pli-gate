#!/usr/bin/env python3
"""C0 (corrected-physics zero-delay eight-tube arm) against the frozen zero-delay stream (June,
superseded I3 convention) and against C2 (388 ps, corrected physics). Scores REGISTERED_C0.md.
Per tube: truncated-mean biases under both rules and the rule artifact Delta|b| = |b_closest| - |b_all|
with paired shard-bootstrap sigma (both rules from the same shard draw, as in B2), truth-signal counts
in [1,5] ns and over-subtraction phi (measure_oversubtraction), truth-grouped (oracle) lifetimes,
the pedestal-law prediction from each arm's own phi (derive_pedestal_bias.predicted_bias), and the
retention Delta|b|(388, C2) / Delta|b|(0, arm). Prompt competition mu from the singles of each arm.
Writes c0_compare.json."""
import os, sys, re, csv, json
from pathlib import Path
import numpy as np, h5py
RES = Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[2])
L = RES / "2026-09-30_feedback3/c0_sim1"
sys.path.insert(0, str(RES / "2026-08-12_steinberger_allprompts")); sys.path.insert(0, str(RES / "2026-08-17_derivation_test"))
import steinberger_extract_ap as X  # noqa
import derive_pedestal_bias as PB  # noqa
exec(compile(re.search(r"^def tau_from_hists\(.*?(?=^def |\Z)", (RES / "2026-08-19_V4_stats/v4d_paired.py").read_text(), re.S | re.M).group(0), "f", "exec"))
ARMS = {
    "june0": dict(npz=RES / "2026-08-22_G3_consequences/hists_sim1.npz",
                  ov=RES / "2026-08-17_derivation_test/oversub_sim1_{}.json", pervoi=L / "june_out_all/per_voi.csv",
                  singles=Path(os.environ.get("PLI_RAW", "<missing>")) / "2026-06-10_brain_tumor_1B_truetime/sim1_singles_truetime.h5"),
    "c0": dict(npz=L / "hists_c0_0ps.npz", ov=L / "oversub_{}.json", pervoi=L / "out_all/per_voi.csv",
               singles=Path(os.environ.get("PLI_RAW", "<missing>")) / "feedback3_c0/c0_sim1_0ps/singles_c0_sim1_0ps.h5"),
    "c2": dict(npz=RES / "2026-08-26_B2_paired/hists_sim1_388.npz", ov=RES / "2026-08-24_C2_sim1_388/oversub_sim1_{}.json",
               pervoi=RES / "2026-08-24_C2_sim1_388/out_all/per_voi.csv", singles=Path(os.environ.get("PLI_RAW", "<missing>")) / "c2_sim1_388/singles_c2_sim1_388.h5"),
}
if not all(Path(ARMS[a]["singles"]).exists() for a in ARMS): raise SystemExit("singles of " + ", ".join(ARMS) + ": raw data not included in the release, see analysis/README.md")
B = 2000


def did_arm(npz):
    d = np.load(npz, allow_pickle=True); names = [str(s) for s in d["voi_names"]]; cfg = d["tau_cfg"]
    cs, cr, as_, ar = d["closest_sig"], d["closest_rnd"], d["all_sig"], d["all_rnd"]; ns = cs.shape[0]
    rng = np.random.default_rng(20261001); out = {}
    for ki, nm in enumerate(names):
        if nm == "bg" or not np.isfinite(cfg[ki]):
            continue
        tau = float(cfg[ki])
        t0c = tau_from_hists(cs[:, ki].sum(0), cr[:, ki].sum(0)); t0a = tau_from_hists(as_[:, ki].sum(0), ar[:, ki].sum(0))
        bc0, ba0 = 100 * (t0c - tau) / tau, 100 * (t0a - tau) / tau
        dd = np.empty(B); bcs = np.empty(B); bas = np.empty(B)
        for b in range(B):
            i = rng.integers(0, ns, ns)
            tc = tau_from_hists(cs[i, ki].sum(0), cr[i, ki].sum(0)); ta = tau_from_hists(as_[i, ki].sum(0), ar[i, ki].sum(0))
            bcs[b], bas[b] = 100 * (tc - tau) / tau, 100 * (ta - tau) / tau; dd[b] = abs(bcs[b]) - abs(bas[b])
        ok = np.isfinite(dd)
        h = lambda v: float(0.5 * np.subtract(*np.nanpercentile(v, [84, 16])))  # manuscript convention: half the 16-84 width
        out[nm] = dict(tau_cfg=tau, b_closest=float(bc0), b_all=float(ba0), did=float(abs(bc0) - abs(ba0)),
                       sd_did=h(dd[ok]), sd_did_std=float(np.std(dd[ok], ddof=1)), sd_b_closest=h(bcs),
                       sd_b_all=h(bas), n_ok=int(ok.sum()))
    return out


def ov(arm):
    r = {}
    for rule in ("closest", "all"):
        for x in json.load(open(str(ARMS[arm]["ov"]).format(rule)))["rows"]:
            r.setdefault(x["voi"], {})[rule] = x
    return r


def pervoi(arm):
    p = ARMS[arm]["pervoi"]
    return {r["voi"]: r for r in csv.DictReader(open(p))} if p.exists() else {}


def mu(path):
    with h5py.File(path, "r") as h:
        g = h["singles"]; sh = g["shard_id"][:]; t = g["time_ns"][:]; npr = int(g["in_prompt"][:].sum())
    order = np.argsort(sh, kind="stable"); sh = sh[order]; t = t[order]
    b = np.flatnonzero(np.diff(sh)) + 1; st = np.r_[0, b]; en = np.r_[b, sh.size]
    span = float(sum(t[s:e].max() - t[s:e].min() for s, e in zip(st, en)))
    return dict(n_prompt=npr, span_ns=span, mu=npr / span * 10.0)


def ivw(v, s):
    v, s = np.asarray(v, float), np.asarray(s, float); w = 1 / s ** 2
    m = float((w * v).sum() / w.sum()); return m, float(np.sqrt(1 / w.sum()))


D = {a: did_arm(ARMS[a]["npz"]) for a in ARMS}
O = {a: ov(a) for a in ARMS}; P = {a: pervoi(a) for a in ARMS}
M = {a: mu(ARMS[a]["singles"]) for a in ARMS}
tubes = sorted(D["c0"])
rows = []
for t in tubes:
    r = dict(region=t, tau_cfg=D["c0"][t]["tau_cfg"])
    for a in ARMS:
        o = O[a][t]; dphi = o["closest"]["phi"] - o["all"]["phi"]
        pred = 100 * abs(PB.predicted_bias(D[a][t]["tau_cfg"], dphi))
        r[a] = dict(**D[a][t], n_true_sig=o["all"]["n_true_sig"], phi_closest=o["closest"]["phi"], phi_all=o["all"]["phi"],
                    d_phi=dphi, pedestal_pred=pred,
                    oracle_tau=float(P[a][t]["oracle_tau_hat_ns"]) if t in P[a] else None,
                    oracle_sig=float(P[a][t]["oracle_sigma_ns"]) if t in P[a] else None)
    r["tail_count_ratio_c0_over_june"] = r["c0"]["n_true_sig"] / r["june0"]["n_true_sig"]
    if r["c0"]["oracle_tau"] is not None and r["june0"]["oracle_tau"] is not None:
        r["oracle_diff_ps"] = 1e3 * (r["c0"]["oracle_tau"] - r["june0"]["oracle_tau"])
        r["oracle_diff_se_ps"] = 1e3 * float(np.hypot(r["c0"]["oracle_sig"], r["june0"]["oracle_sig"]))
    r["did_ratio_c0_over_june"] = r["c0"]["did"] / r["june0"]["did"]
    r["did_ratio_se"] = abs(r["did_ratio_c0_over_june"]) * float(np.hypot(r["c0"]["sd_did"] / r["c0"]["did"], r["june0"]["sd_did"] / r["june0"]["did"]))
    for base in ("c0", "june0"):
        q = r["c2"]["did"] / r[base]["did"]
        r[f"retention_vs_{base}"] = q
        r[f"retention_vs_{base}_se"] = abs(q) * float(np.hypot(r["c2"]["sd_did"] / r["c2"]["did"], r[base]["sd_did"] / r[base]["did"]))
    sig_pp = float(np.hypot(r["c0"]["sd_did"], 0.37))
    r["c0_pedestal_pull"] = (r["c0"]["did"] - r["c0"]["pedestal_pred"]) / sig_pp
    rows.append(r)
    print(f"{t:28s} tail C0/June {r['tail_count_ratio_c0_over_june']:.3f}  DiD June {r['june0']['did']:.2f}+-{r['june0']['sd_did']:.2f}"
          f"  C0 {r['c0']['did']:.2f}+-{r['c0']['sd_did']:.2f}  C2 {r['c2']['did']:.2f}+-{r['c2']['sd_did']:.2f}"
          f"  ratio {r['did_ratio_c0_over_june']:.3f}  ret(C0) {r['retention_vs_c0']:.3f}+-{r['retention_vs_c0_se']:.3f}"
          f"  ret(June) {r['retention_vs_june0']:.3f}  ped pull C0 {r['c0_pedestal_pull']:+.2f}"
          + (f"  oracle d {r['oracle_diff_ps']:+.1f}+-{r['oracle_diff_se_ps']:.1f} ps" if 'oracle_diff_ps' in r else ""))
S = {}
S["retention_ivw_c0"] = ivw([r["retention_vs_c0"] for r in rows], [r["retention_vs_c0_se"] for r in rows])
S["retention_ivw_june0"] = ivw([r["retention_vs_june0"] for r in rows], [r["retention_vs_june0_se"] for r in rows])
S["did_ratio_ivw"] = ivw([r["did_ratio_c0_over_june"] for r in rows], [r["did_ratio_se"] for r in rows])
S["did_ivw"] = {a: ivw([r[a]["did"] for r in rows], [r[a]["sd_did"] for r in rows]) for a in ARMS}
S["b_closest_mean_abs"] = {a: float(np.mean([abs(r[a]["b_closest"]) for r in rows])) for a in ARMS}
S["b_all_mean_abs"] = {a: float(np.mean([abs(r[a]["b_all"]) for r in rows])) for a in ARMS}
S["b_all_ivw"] = {a: ivw([r[a]["b_all"] for r in rows], [r[a]["sd_b_all"] for r in rows]) for a in ARMS}
S["tail_ratio_range"] = [min(r["tail_count_ratio_c0_over_june"] for r in rows), max(r["tail_count_ratio_c0_over_june"] for r in rows)]
if all("oracle_diff_ps" in r for r in rows):
    S["oracle_diff_ivw_ps"] = ivw([r["oracle_diff_ps"] for r in rows], [r["oracle_diff_se_ps"] for r in rows])
    S["oracle_max_abs_z"] = max(abs(r["oracle_diff_ps"] / r["oracle_diff_se_ps"]) for r in rows)
S["mu"] = M; S["mu_rel_c0_vs_c2"] = M["c0"]["mu"] / M["c2"]["mu"] - 1; S["mu_rel_c0_vs_june0"] = M["c0"]["mu"] / M["june0"]["mu"] - 1
S["pedestal_max_abs_pull_c0"] = max(abs(r["c0_pedestal_pull"]) for r in rows)
print(json.dumps(S, indent=1))
sc = {}
sc["Q1 tail count C0/June 1.25-1.42 per tube"] = bool(all(1.25 <= r["tail_count_ratio_c0_over_june"] <= 1.42 for r in rows))
if "oracle_diff_ivw_ps" in S:
    m, s = S["oracle_diff_ivw_ps"]
    sc["Q2 oracle lifetimes unchanged (per tube |z|<3, pooled |z|<2)"] = bool(S["oracle_max_abs_z"] < 3 and abs(m / s) < 2)
sc["Q3 mu C0 within 0.5% of C2"] = bool(abs(S["mu_rel_c0_vs_c2"]) < 0.005)
sc["Q4 Delta|b| smaller in C0 in every tube, pooled ratio 0.65-0.85"] = bool(all(r["c0"]["did"] < r["june0"]["did"] for r in rows) and 0.65 <= S["did_ratio_ivw"][0] <= 0.85)
sc["Q5 pedestal law predicts C0 within 3 sigma per tube"] = bool(S["pedestal_max_abs_pull_c0"] < 3)
sc["Q6 retention vs C0 pooled 0.42-0.54"] = bool(0.42 <= S["retention_ivw_c0"][0] <= 0.54)
for k, v in sc.items(): print(("PASS " if v else "FAIL ") + k)
json.dump(dict(rows=rows, summary=S, scores=sc, B=B, sigma_convention='half the 16-84 percentile width of the paired shard bootstrap'), open(L / "c0_compare.json", "w"), indent=1)
