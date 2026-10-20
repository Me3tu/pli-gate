#!/usr/bin/env python3
"""V2 support — the three context numbers the registered DiD does not carry.

1. Prediction 1 asked whether the ABSOLUTE bias worsens. Bias has to be
   measured against a comparison axis, and the right one is each arm's OWN
   oracle (truth-tagged genuine triples), not the configured 1.84 ns: the
   early-time model changes what a single-exponential estimator targets, so
   comparing arm B to 1.84 would score the model change as if it were bias.
2. The artefact is quoted in ns, but both arms' tau differ by ~0.55 ns, so an
   absolute shrinkage is not the same statement as a fractional one. Both are
   reported; the ABSOLUTE one is what was registered.
3. Prediction 2's mechanism said the artefact should fall by the fractional
   growth of the fit-window genuine count. That is a testable number.
"""
from __future__ import annotations
import csv, json, os
from pathlib import Path

RES = Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[1])
CSVS = {
    ("A_0ps", "closest"): RES / "2026-08-17_derivation_test/water_closest_o4/per_voi.csv",
    ("A_0ps", "all"):     RES / "2026-08-17_derivation_test/water_all_o4/per_voi.csv",
    ("B_388ps", "closest"): RES / "2026-08-19_V2_earlytime/out_388_closest/per_voi.csv",
    ("B_388ps", "all"):     RES / "2026-08-19_V2_earlytime/out_388_all/per_voi.csv",
}
D = {k: {r["voi"]: r for r in csv.DictReader(open(p))} for k, p in CSVS.items()}
did = json.load(open("v2_did.json"))["results"]
names = ["bg", "VOI_5x15", "VOI_10x20", "VOI_15x25"]
out = {}

print("=== 1. bias vs each arm's own oracle (ns) ===")
print(f"{'VOI':10s} {'oracle_0':>9s} {'oracle_388':>11s} | "
      f"{'bias_c_0':>9s} {'bias_c_388':>11s} | {'bias_a_0':>9s} {'bias_a_388':>11s}")
for n in names:
    o0 = float(D[("A_0ps", "closest")][n]["oracle_tau_hat_ns"])
    o1 = float(D[("B_388ps", "closest")][n]["oracle_tau_hat_ns"])
    b = {}
    for arm, o in (("A_0ps", o0), ("B_388ps", o1)):
        for m in ("closest", "all"):
            b[(arm, m)] = float(D[(arm, m)][n]["clinical_tau_trunc_ns"]) - o
    print(f"{n:10s} {o0:9.4f} {o1:11.4f} | {b[('A_0ps','closest')]:9.4f} "
          f"{b[('B_388ps','closest')]:11.4f} | {b[('A_0ps','all')]:9.4f} "
          f"{b[('B_388ps','all')]:11.4f}")
    out.setdefault(n, {})["oracle_0ps"] = o0
    out[n]["oracle_388ps"] = o1
    out[n]["bias_vs_oracle"] = {f"{a}_{m}": v for (a, m), v in b.items()}

print("\n=== 2. artefact, absolute vs fractional ===")
print(f"{'VOI':10s} {'abs_0':>9s} {'abs_388':>9s} {'abs_ret':>8s} | "
      f"{'frac_0':>8s} {'frac_388':>9s} {'frac_ret':>9s}")
for n in names:
    r = did[n]
    f0 = r["artefact_0ps"] / r["tau_all_0ps"]
    f1 = r["artefact_388ps"] / r["tau_all_388ps"]
    print(f"{n:10s} {r['artefact_0ps']:9.4f} {r['artefact_388ps']:9.4f} "
          f"{r['retained_fraction']*100:7.1f}% | {f0*100:7.2f}% {f1*100:8.2f}% "
          f"{f1/f0*100:8.1f}%")
    out[n]["artefact_frac_0ps"] = f0
    out[n]["artefact_frac_388ps"] = f1
    out[n]["retained_fractional"] = f1 / f0

print("\n=== 3. prediction-2 mechanism check: retention vs count growth ===")
print("   predicted retention = 1 / (growth of the genuine count)")
print(f"{'VOI':10s} {'clin_n_sig x':>13s} {'pred_ret':>9s} | "
      f"{'oracle_n x':>11s} {'pred_ret':>9s} | {'measured':>9s}")
for n in names:
    g_cl = (float(D[("B_388ps", "closest")][n]["clinical_n_signal"])
            / float(D[("A_0ps", "closest")][n]["clinical_n_signal"]))
    g_or = (float(D[("B_388ps", "closest")][n]["oracle_n_signal"])
            / float(D[("A_0ps", "closest")][n]["oracle_n_signal"]))
    m = did[n]["retained_fraction"]
    print(f"{n:10s} {g_cl:13.3f} {1/g_cl*100:8.1f}% | {g_or:11.3f} "
          f"{1/g_or*100:8.1f}% | {m*100:8.1f}%")
    out[n]["growth_clinical_n_signal"] = g_cl
    out[n]["growth_oracle_n_signal"] = g_or

Path("v2_context.json").write_text(json.dumps(out, indent=2))
print("\nwrote v2_context.json")

# ---------------------------------------------------------------------------
# 4. Does the oracle shift close with no free parameters?
#
# The oracle estimator has NO upper cut (steinberger_extract_ap.py:252-258):
# it is mean(tau | tau > 1 ns) - 1 ns. For an exponential that is exactly the
# lifetime, which is why arm A returns the configured 1.84 ns.
#
# In arm B the same population gains leaked direct annihilations. Above the
# 1 ns censor an exponential is memoryless, so each component contributes its
# own lifetime and the estimator returns the COUNT-WEIGHTED MEAN:
#
#     tau_oracle(388) = f_d * 0.388 + (1 - f_d) * tau_oracle(0)
#
# f_d, the leaked-direct share of the fit-window genuine population, is not
# fitted: arm A's direct events sit at exactly 0 delay and cannot reach 1 ns
# (that is 11 sigma at the 90.9 ps timing smear), so every extra genuine
# triple in arm B is a leaked direct one, and f_d = 1 - N_A/N_B.
TAU_DIRECT = 0.388
print("\n=== 4. oracle shift, closed with no free parameters ===")
print(f"{'VOI':10s} {'f_d':>7s} {'predicted':>10s} {'measured':>9s} {'rel err':>9s}")
for n in names:
    nA = float(D[("A_0ps", "closest")][n]["oracle_n_signal"])
    nB = float(D[("B_388ps", "closest")][n]["oracle_n_signal"])
    f_d = 1.0 - nA / nB
    pred = f_d * TAU_DIRECT + (1.0 - f_d) * out[n]["oracle_0ps"]
    meas = out[n]["oracle_388ps"]
    print(f"{n:10s} {f_d:7.3f} {pred:10.4f} {meas:9.4f} {(pred-meas)/meas*100:8.2f}%")
    out[n]["f_direct_leaked"] = f_d
    out[n]["oracle_388ps_predicted"] = pred
    out[n]["oracle_388ps_rel_err"] = (pred - meas) / meas

Path("v2_context.json").write_text(json.dumps(out, indent=2))
print("\nrewrote v2_context.json with the closure check")
