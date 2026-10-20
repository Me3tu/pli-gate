#!/usr/bin/env python3
"""Review round 5: replication (Fig. 3) with the acceptance width of the replication geometry.
Re-runs, at a given sigma_e, everything the paper quotes from the replication fit:
 b3 main fits (all tubes, both rules), the b4 endpoint sweep (stat, half-spread syst, total, z
 against the published values with their uncertainties), the Fig. 3a spectrum of dry XAD-4,
 and the lower-edge sensitivity. Identical code paths (b3_fit defs exec'd unchanged). Only
 G.SIGMA_E is set. Inputs read-only from 2026-08-23_B3, outputs written here.
Usage: run_repl.py SIGMA_E_NS"""
import json, os, sys
from pathlib import Path
import numpy as np
SIG = float(sys.argv[1])
OUT = Path(__file__).resolve().parent
B3 = Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[1]) / "2026-08-23_B3"
os.chdir(B3); sys.argv = ["b3"]
src = open(B3 / "b3_fit.py").read().split('z = np.load("hists_mercolli.npz"')[0]
exec(compile(src, "b3_fit_defs", "exec"))          # defines X, E, G, BC, WB, fit3, PUB, trunc_mean
G.SIGMA_E = SIG
from acceptance import A_analytic
PUBERR = {"T1_dry": 0.03, "T2_wet": 0.03, "T3_gelatin": 0.04, "T4_water": 0.02}   # fig12_replication.json forest
z = np.load("hists_mercolli.npz", allow_pickle=True)
names = [str(s) for s in z["voi_names"]]
res = dict(sigma_e_ns=SIG, fits={}, sweep={}, forest=[], edge={})
HIS = [3.5, 4.0, 4.5, 5.0]
for tube in PUB:
    i = names.index(tube)
    for rule in ("closest", "all"):
        hs = z[f"{rule}_sig"][:, i, :].sum(0).astype(float); hr = z[f"{rule}_rnd"][:, i, :].sum(0).astype(float)
        M = (BC > -1.0) & (BC < 4.5)
        edges = np.r_[BC[M] - 0.5 * WB, BC[M][-1] + 0.5 * WB]
        t3, sd3, par = fit3(np.clip((hs - hr)[M], 0.0, None), edges)
        res["fits"].setdefault(tube, {})[rule] = dict(tau3=float(t3), sd3=float(sd3), params=[float(p) for p in par])
        vals, sds = [], []
        for hi in HIS:
            m = (BC > -1.0) & (BC < hi)
            e2 = np.r_[BC[m] - 0.5 * WB, BC[m][-1] + 0.5 * WB]
            t, sd, _ = fit3(np.clip((hs - hr)[m], 0.0, None), e2)
            vals.append(t); sds.append(sd)
        v = np.array(vals); s = np.array(sds); g = np.isfinite(v)
        syst = float(v[g].max() - v[g].min()); stat = float(np.median(s[g])); tot = float(np.hypot(stat, syst / 2))
        res["sweep"].setdefault(tube, {})[rule] = dict(tau=[float(x) for x in v], sd=[float(x) for x in s],
                                                       syst_full=syst, syst_half=syst / 2, stat=stat, total=tot)
        pub = PUB[tube][0]; pe = PUBERR[tube]
        res["forest"].append(dict(tube=tube, rule=rule, pub=pub, pub_err=pe, ours=float(t3), stat=stat,
                                  syst_half=syst / 2, syst_full=syst, total=tot, diff_ps=1000 * (t3 - pub),
                                  z=float((t3 - pub) / np.hypot(tot, pe))))
        print(tube, rule, f"tau3={t3:.4f}+/-{sd3:.4f} sweep={1000*syst:.1f}ps stat={1000*stat:.1f} z={res['forest'][-1]['z']:+.2f}", flush=True)
    # lower-edge sensitivity (all-prompts), as fit_range_sensitivity.py
    hs = z["all_sig"][:, i, :].sum(0).astype(float); hr = z["all_rnd"][:, i, :].sum(0).astype(float)
    res["edge"][tube] = {}
    for lo in (-1.0, -0.5, -0.3, 0.0, 0.5):
        M = (BC > lo) & (BC < 4.5)
        edges = np.r_[BC[M] - 0.5 * WB, BC[M][-1] + 0.5 * WB]
        t3, sd3, par = fit3(np.clip((hs - hr)[M], 0.0, None), edges)
        res["edge"][tube][str(lo)] = dict(tau3=float(t3), sd3=float(sd3))
# Fig. 3a spectrum, dry XAD-4 all-prompts (as fig3_spectrum.py)
i = names.index("T1_dry")
hs = z["all_sig"][:, i, :].sum(0).astype(float); hr = z["all_rnd"][:, i, :].sum(0).astype(float)
M = (BC > -1.0) & (BC < 4.5)
edges = np.r_[BC[M] - 0.5 * WB, BC[M][-1] + 0.5 * WB]
sig = np.clip((hs - hr)[M], 0.0, None)
t3, sd3, par = fit3(sig, edges)
accb = np.array([np.trapezoid(A_analytic(np.linspace(a, b, 60), G.SIGMA_E), np.linspace(a, b, 60)) / (b - a)
                 for a, b in zip(edges[:-1], edges[1:])])
w = np.exp(np.array([par[3], par[4], 0.0])); I = w / w.sum()
comps = [Ik * np.diff(E.emg_cdf(edges, tk, par[1], par[2])) * accb for Ik, tk in zip(I, (E.TAU_PPS, E.TAU_DIR, par[0]))]
tot = np.clip(np.sum(comps, axis=0), 1e-15, None); norm = sig.sum() / tot.sum()
y = (hs - hr)[M]; var = hs[M] + hr[M]
pull = (y - norm * tot) / np.sqrt(np.clip(var, 1, None))
res["spectrum"] = dict(tube="T1_dry", rule="all", fit_range_ns=[-1.0, 4.5], tau3=float(t3), sd3=float(sd3),
                       params=[float(v) for v in par], bin_centres_ns=BC[M].tolist(), data=y.tolist(), var=var.tolist(),
                       model_total=(norm * tot).tolist(), model_ops=(norm * comps[2]).tolist(),
                       model_short=(norm * (comps[0] + comps[1])).tolist(), pull=pull.tolist(),
                       chi2=float((pull ** 2).sum()), ndf=int(M.sum() - 5), max_abs_pull=float(np.abs(pull).max()),
                       t_max_abs_pull=float(BC[M][np.argmax(np.abs(pull))]))
print("spectrum chi2", res["spectrum"]["chi2"], res["spectrum"]["ndf"], "max|pull|", res["spectrum"]["max_abs_pull"])
json.dump(res, open(OUT / f"repl_sig{SIG:.4f}.json", "w"), indent=1)
print("wrote", OUT / f"repl_sig{SIG:.4f}.json")
