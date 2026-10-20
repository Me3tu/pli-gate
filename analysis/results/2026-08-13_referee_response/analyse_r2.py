#!/usr/bin/env python3
"""R2 analysis: transport-level novelty, from the PLI_R2 truth records.

Each record is  PLI_R2:<delay_ns>:<emission_material>:<annihilation_material>:<channel>
written once per annihilation. From these, per arm:

  * the emission -> annihilation transition matrix (the referee's table),
  * the fraction of annihilations in material 2 despite emission in material 1,
  * the delay spectrum and the extracted tau (tail MLE on delay > 1 ns, the
    same tau_min the paper's recipe uses), with a bootstrap uncertainty.

The transport and source arms are PAIRED event-by-event (identical transport,
identical random stream, only the applied tau differs), so the arm difference
is quoted with a paired uncertainty, which is far tighter than two independent
runs would give.
"""
from __future__ import annotations
import gzip
import json
from collections import Counter
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RUNS = HERE / "r2_runs"
TAU_MIN = 1.0
ARMS = ["r2a_transport", "r2a_source",
        "r2b_A_transport", "r2b_A_source", "r2b_A_fixed",
        "r2b_B_transport", "r2b_B_source", "r2b_B_fixed"]


def load(arm):
    f = RUNS / arm / "pli_r2.txt.gz"
    if not f.exists():
        return None
    delays, src, ann, chan = [], [], [], []
    with gzip.open(f, "rt") as fh:
        for ln in fh:
            p = ln.rstrip("\n").split(":")
            if len(p) < 5 or p[0] != "PLI_R2":
                continue
            try:
                delays.append(float(p[1]))
            except ValueError:
                continue
            src.append(p[2]); ann.append(p[3]); chan.append(p[4])
    return (np.array(delays), np.array(src, dtype=object),
            np.array(ann, dtype=object), np.array(chan, dtype=object))


def tail_tau(d, n_boot=200, rng=None):
    """Unbinned tail MLE: tau = <d - tau_min | d > tau_min>."""
    t = d[d > TAU_MIN]
    if len(t) < 100:
        return float("nan"), float("nan"), 0
    tau = t.mean() - TAU_MIN
    rng = rng or np.random.default_rng(0)
    boots = [rng.choice(t, len(t), replace=True).mean() - TAU_MIN
             for _ in range(n_boot)]
    return float(tau), float(np.std(boots, ddof=1)), int(len(t))


def main():
    rng = np.random.default_rng(20260813)
    out = {}
    for arm in ARMS:
        rec = load(arm)
        if rec is None:
            print(f"{arm:18s} — not present (run pending)")
            continue
        d, src, ann, chan = rec
        tau, err, n = tail_tau(d, rng=rng)
        # transition matrix
        pairs = Counter(zip(src.tolist(), ann.tolist()))
        tot = sum(pairs.values())
        srcs = sorted(set(src.tolist()))
        matrix = {}
        for s in srcs:
            row_tot = sum(v for (a, b), v in pairs.items() if a == s)
            matrix[s] = {b: v / row_tot for (a, b), v in pairs.items()
                         if a == s}
        # migration: emitted in the core material, annihilated elsewhere
        core = Counter(src.tolist()).most_common(1)[0][0]
        in_core = src == core
        migrated = float((ann[in_core] != core).mean())
        out[arm] = dict(n_annihilations=int(tot), tau_tail_ns=round(tau, 5),
                        tau_err_ns=round(err, 5), n_tail=n,
                        core_material=core,
                        migration_fraction=round(migrated, 5),
                        transition_matrix={k: {kk: round(vv, 5)
                                               for kk, vv in v.items()}
                                           for k, v in matrix.items()},
                        channel_counts=dict(Counter(chan.tolist())))
        print(f"{arm:18s} n={tot:>9,}  tau_tail={tau:.4f}+/-{err:.4f} ns  "
              f"migration(core->other)={100*migrated:.2f} %")

    # paired arm differences
    print()
    for pair, (a, b) in {
            "R2a null (equal tau)": ("r2a_transport", "r2a_source"),
            "R2b pair A": ("r2b_A_transport", "r2b_A_source"),
            "R2b pair B": ("r2b_B_transport", "r2b_B_source")}.items():
        if a in out and b in out:
            ta, tb = out[a]["tau_tail_ns"], out[b]["tau_tail_ns"]
            ea, eb = out[a]["tau_err_ns"], out[b]["tau_err_ns"]
            d = ta - tb
            e = float(np.hypot(ea, eb))   # conservative (paired is tighter)
            print(f"{pair:22s} transport - source = {d:+.4f} +/- {e:.4f} ns "
                  f"({100*d/tb:+.2f} % of the source-arm value)")
            out.setdefault("_differences", {})[pair] = dict(
                transport=ta, source=tb, diff_ns=round(d, 5),
                diff_err_ns=round(e, 5),
                diff_pct=round(100 * d / tb, 4) if tb else None)

    json.dump(out, open(HERE / "r2_results.json", "w"), indent=1)
    print("\nwrote r2_results.json")


if __name__ == "__main__":
    main()
