#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-only
# Copyright (C) 2026 Dhritman Dhar and German Cancer Research Center (DKFZ)
"""PLI package test suite. Nine checks, one command, no network, no big data.

Every reference value lives in reference_values.json beside this file. Nothing
reaches outside the package.

    python3 test_pli.py            run everything
    python3 test_pli.py -k oxygen  run the tests whose name contains "oxygen"
"""
from __future__ import annotations
import json, math, os, subprocess, sys
from pathlib import Path

PKG = Path(__file__).resolve().parent.parent
# tao_eldrup.py and pedestal.py belong to the analysis part of the release.
# Search the layouts it has used, first match wins; PLI_ANALYSIS_PY overrides.
for _d in reversed([os.environ.get("PLI_ANALYSIS_PY", ""),
                    PKG / "analysis" / "pli_analysis",
                    PKG / "analysis",
                    PKG / "analysis" / "legacy_v1.0.1"]):
    if _d and Path(_d).is_dir():
        sys.path.insert(0, str(_d))
REF = json.loads((Path(__file__).parent / "reference_values.json").read_text())
TESTS: list[tuple[str, callable]] = []


def test(name):
    def deco(fn):
        TESTS.append((name, fn))
        return fn
    return deco


# ------------------------------------------------------------------------- 1
@test("tao_eldrup_numeric")
def t1():
    """R -> tau -> R round-trips, and water sits at its literature value."""
    from tao_eldrup import tau_eff, R_from_tau
    for R in (0.20, 0.27, 0.35, 0.50):
        assert abs(R_from_tau(tau_eff(R, 0.0), 0.0) - R) < 1e-6, f"round trip failed at R={R}"
    R_water = R_from_tau(1.80, 0.0)
    assert abs(tau_eff(R_water, 0.0) - 1.80) < 1e-9
    return f"round-trip < 1e-6 nm; water anchor {tau_eff(R_water, 0.0):.4f} ns"


# ------------------------------------------------------------------------- 2
@test("o2_quenching_conversion")
def t2():
    """Oxygen shortens the lifetime, monotonically, by the documented amount."""
    from tao_eldrup import tau_eff, R_from_tau, K_O2_NS_MMHG
    r = REF["oxygen_quenching"]
    db_k = json.loads((PKG / "materials/pli_materials.json").read_text())["constants"]["k_O2_ns_mmHg"]
    assert abs(r["kappa_ns_mmHg"] - db_k) < 1e-12, \
        f"reference kappa {r['kappa_ns_mmHg']} differs from the shipped database {db_k}"
    assert abs(K_O2_NS_MMHG - db_k) < 1e-12, \
        f"tao_eldrup kappa {K_O2_NS_MMHG} differs from the shipped database {db_k}"
    R = R_from_tau(r["tau_at_pO2_0_ns"], 0.0)
    prev = 1e9
    for p in (0, 10, 20, 40, 100, 200):
        t = tau_eff(R, p)
        assert t < prev, f"tau not decreasing at pO2={p}"
        prev = t
    t100 = tau_eff(R, 100.0)
    assert abs(t100 - r["tau_at_pO2_100_ns"]) < 1e-3, f"tau(100 mmHg) = {t100:.4f}"
    return (f"monotone; tau({r['tau_at_pO2_0_ns']} ns at 0) -> {t100:.4f} ns at "
            f"100 mmHg, a shift of {(r['tau_at_pO2_0_ns']-t100)*1000:.1f} ps")


# ------------------------------------------------------------------------- 3
@test("i3_convention_is_declared")
def t3():
    """Every database entry declares its I3 convention. Absence is a hard error
    in the loader because the two conventions differ by 4/3."""
    n = 0
    for f in ("materials/pli_materials.json", "materials/pli_materials_vipman.json"):
        d = json.loads((PKG / f).read_text())
        for m in d["materials"]:
            c = m.get("I3_convention")
            assert c in ("ops_intensity", "total_formation"), \
                f"{f}: {m['name']} declares {c!r}"
            i3 = m["I3"] * (4 / 3 if c == "ops_intensity" else 1.0)
            assert 0.0 <= i3 <= 1.0, f"{m['name']}: implied formation prob {i3:.3f}"
            n += 1
    return f"{n} entries all declare a valid convention with P(form) in [0,1]"


# ------------------------------------------------------------------------- 4
@test("three_gamma_branching [recorded result]")
def t4():
    """RECORDED RESULT, not a live check of your build: it validates the shipped
    constant and the consistency of a tally measured on the reference run. It
    does not execute package code. Re-deriving it on your own build means
    counting 3-gamma finals in a run log."""
    r = REF["three_gamma"]
    T3 = r["tau_3gamma_ns"]
    # the ortho-positronium vacuum lifetime, the constant the model is built on
    assert abs(T3 - 142.05) < 1e-9, f"shipped tau_3gamma {T3} != 142.05 ns"
    implied = r["production_P3g_given_oPs_percent"] / 100.0 * T3
    assert 2.0 < implied < 2.5, f"implied mean tau_eff {implied:.3f} ns out of range"
    return (f"P(3g)=tau/{T3}; production "
            f"{r['production_P3g_given_oPs_percent']}% -> mean tau_eff {implied:.3f} ns")


# ------------------------------------------------------------------------- 5
@test("ore_powell_spectrum")
def t5():
    """The three-photon final state, checked by the compiled C++ test."""
    for exe in (PKG / "custom_physics/build/test_orepowell",
                PKG / "build/test_orepowell"):
        if exe.exists():
            p = subprocess.run([str(exe)], capture_output=True, text=True, timeout=600)
            assert p.returncode == 0, f"test_orepowell rc={p.returncode}\n{p.stdout[-500:]}"
            return f"{exe.name} passed (energy conservation + Dalitz envelope)"
    return "SKIP: test_orepowell not built (build custom_physics to enable)"


# ------------------------------------------------------------------------- 6
@test("null_transport_vs_source [recorded result]")
def t6():
    """RECORDED RESULT, not a live check of your build: it compares a stored
    measurement against a stored tolerance, both from this file. It does not
    execute package code. Reproducing it means re-running the null arm."""
    r = REF["null_transport_vs_source"]
    assert abs(r["difference_percent"]) < r["tolerance_percent"], \
        f"null arm {r['difference_percent']}% exceeds {r['tolerance_percent']}%"
    return (f"measured null difference {r['difference_ns']:+.5f} ns "
            f"({r['difference_percent']:+.2f} %), within {r['tolerance_percent']} %")


# ------------------------------------------------------------------------- 7
@test("boundary_case_materials")
def t7():
    """Air forms no positronium; every declared lifetime is finite and positive."""
    d = json.loads((PKG / "materials/pli_materials_vipman.json").read_text())
    air = [m for m in d["materials"] if m["name"] == "Air"]
    assert air, "no Air entry in the atlas database"
    assert air[0]["I3"] == REF["boundary_case"]["air_I3"], "air I3 changed"
    n = 0
    for m in d["materials"]:
        t = m.get("tau_TE_ns")
        assert t is not None and math.isfinite(t) and t > 0, f"{m['name']} tau={t}"
        n += 1
    return f"Air I3=0; all {n} tau_eff finite and positive"


# ------------------------------------------------------------------------- 8
@test("pedestal_analytic_identity")
def t8():
    """The truncated-mean inversion is exact, phi=0 is a no-op, and
    over-subtraction biases the lifetime low by the documented amount."""
    import pedestal as P
    r = REF["pedestal_law"]
    for tau in (1.56, 1.80, 2.04, 2.72):
        m1 = P.trunc_mean(tau)
        assert abs(P.invert_trunc(m1) - tau) < 1e-4, f"inversion failed at {tau}"
        assert abs(P.predicted_bias(tau, 0.0)) < 1e-9, "phi=0 must give zero bias"
    b = 100.0 * P.predicted_bias(1.80, 0.10)
    assert b < 0, "over-subtraction must bias tau LOW"
    assert abs(b - r["bias_tau1p8_phi0p1_percent"]) < r["tolerance_percent"], \
        f"bias {b:.2f}% != reference {r['bias_tau1p8_phi0p1_percent']}%"
    return f"inversion exact to 1e-4; phi=0 no-op; bias(1.8 ns, phi=0.1) = {b:.2f}%"


# ------------------------------------------------------------------------- 9
@test("fixed_seed_regression")
def t9():
    """The quickstart, if it has been run, must reproduce its committed value."""
    log = PKG / "examples/quickstart/quickstart_run.log"
    exp = json.loads((PKG / "examples/quickstart/expected_output.json").read_text())
    if not log.exists():
        return "SKIP: run examples/quickstart/run_quickstart.sh first"
    tau = [float(l.split(":")[1]) for l in log.read_text(errors="ignore").splitlines()
           if l.startswith("PLI_LIFETIME:")]
    above = [t for t in tau if t > 1.0]
    assert above, "no annihilations above the censor in the quickstart log"
    hat = sum(above) / len(above) - 1.0
    d = hat - exp["tau_oPs_ns"]
    assert abs(d) <= exp["tolerance_ns"], \
        f"tau_oPs {hat:.4f} differs from {exp['tau_oPs_ns']} by {d:+.4f} ns"
    return (f"{len(tau):,} annihilations, {len(above):,} above censor, "
            f"tau_oPs = {hat:.4f} ns ({d:+.4f} vs configured)")


def main() -> int:
    only = sys.argv[sys.argv.index("-k") + 1] if "-k" in sys.argv else None
    sel = [(n, f) for n, f in TESTS if not only or only in n]
    w = max(len(n) for n, _ in sel) if sel else 10
    npass = nfail = nskip = 0
    print(f"{'test':<{w}}  result")
    print("-" * (w + 62))
    for name, fn in sel:
        try:
            msg = fn()
            if str(msg).startswith("SKIP"):
                nskip += 1; status = "SKIP"
            else:
                npass += 1; status = "PASS"
        except Exception as e:
            nfail += 1; status = "FAIL"; msg = f"{type(e).__name__}: {e}"
        print(f"{name:<{w}}  {status}  {msg}")
    print("-" * (w + 62))
    print(f"{npass} passed, {nfail} failed, {nskip} skipped")
    return 1 if nfail else 0


if __name__ == "__main__":
    sys.exit(main())
