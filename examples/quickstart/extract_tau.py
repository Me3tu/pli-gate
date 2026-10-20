#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-only
# Copyright (C) 2026 Dhritman Dhar and German Cancer Research Center (DKFZ)
"""Recover the o-Ps lifetime from a quickstart run's stdout.

The physics writes one record per annihilation when PLI_LIFETIME_LOG=1:

    PLI_LIFETIME:<tau_ns>:<material>:<x>:<y>:<z>:<global_time_ns>

`tau_ns` is the delay the physics assigned to that annihilation: 0 for direct
annihilation, a draw from the 125 ps exponential for p-Ps, and a draw from the
material's tau_eff for o-Ps.

Estimator: mean(tau | tau > 1 ns) - 1 ns. For an exponential this is exactly
the lifetime, because the exponential is memoryless. The 1 ns censor removes
direct annihilation entirely and leaves a p-Ps contamination of exp(-1/0.125)
= 3.4e-4, which is far below the statistical error of a quickstart-sized run.

Usage:  python3 extract_tau.py run.log [--material NAME]
"""
from __future__ import annotations
import argparse, json, math, sys
from pathlib import Path

TAU_MIN_NS = 1.0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("logfile", type=Path)
    ap.add_argument("--material", default=None,
                    help="restrict to one material (default: the commonest)")
    ap.add_argument("--expected", type=Path, default=Path(__file__).with_name("expected_output.json"))
    a = ap.parse_args()

    if not a.logfile.is_file():
        print(f"ERROR: no such file: {a.logfile}", file=sys.stderr)
        return 2

    per_mat: dict[str, list[float]] = {}
    for line in a.logfile.read_text(errors="ignore").splitlines():
        if not line.startswith("PLI_LIFETIME:"):
            continue
        f = line.split(":")
        if len(f) < 3:
            continue
        try:
            per_mat.setdefault(f[2], []).append(float(f[1]))
        except ValueError:
            pass

    if not per_mat:
        print("ERROR: no PLI_LIFETIME records found.", file=sys.stderr)
        print("       Did you run with PLI_LIFETIME_LOG=1 ?", file=sys.stderr)
        return 1

    mat = a.material or max(per_mat, key=lambda m: len(per_mat[m]))
    tau = per_mat[mat]
    above = [t for t in tau if t > TAU_MIN_NS]
    n = len(above)
    if n < 100:
        print(f"ERROR: only {n} annihilations above {TAU_MIN_NS} ns; run longer.",
              file=sys.stderr)
        return 1

    tau_hat = sum(above) / n - TAU_MIN_NS
    sigma = tau_hat / math.sqrt(n)

    print(f"material            : {mat}")
    print(f"annihilations       : {len(tau):,}")
    print(f"above {TAU_MIN_NS} ns censor  : {n:,}")
    print()
    print(f"  tau_oPs = {tau_hat:.4f} +/- {sigma:.4f} ns")
    print()

    if a.expected.is_file():
        exp = json.loads(a.expected.read_text())
        ref, tol = exp["tau_oPs_ns"], exp["tolerance_ns"]
        d = tau_hat - ref
        ok = abs(d) <= tol
        print(f"expected            : {ref:.4f} ns (configured value)")
        print(f"difference          : {d:+.4f} ns   tolerance +/-{tol:.4f}")
        print()
        print("RESULT: PASS" if ok else "RESULT: FAIL")
        return 0 if ok else 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
