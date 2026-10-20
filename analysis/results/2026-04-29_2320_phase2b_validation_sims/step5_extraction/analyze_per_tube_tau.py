"""Phase 2b Step 5: per-tube oracle MLE tau on eventID-tagged triples.

Reads detector + truth H5 from extract_triples_chunked / extract_triples_i124,
applies the project's existing travel-time correction (importing
reconstruction.correct_travel_time.correct), tags each triple to a tube
via the truth-side source position, and computes per-tube MLE tau.

MLE formula (matches dhar_pli_sensitivity_2026/pipeline/xad4_complete_pipeline.py
lines 565-583, the project's standard exponential MLE):
    tau_o = mean(tau_corrected[tau_corrected > tau_min]) - tau_min
    sigma = tau_o / sqrt(N)
where tau_min is a low cutoff to exclude prompt-coincidence pile-up.

Usage:
    python analyze_per_tube_tau.py --detector <det.h5> --truth <tru.h5>
        --tube-config <sim1|sim2> --tau-min 0.3
        --out-md <validation_table.md> --out-csv <table.csv>
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import h5py
import numpy as np

REPO = Path(os.environ.get("PLI_REPO") or Path(__file__).resolve().parents[4])
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

# Import the project's existing travel-time correction (read-only, just call
# the correct() function on our paths). The module's main() hardcodes a path,
# but correct() is path-parameterised.
from reconstruction.correct_travel_time import correct as travel_time_correct  # noqa: E402

# Tube layouts: must mirror the GATE macros' /placement/setTranslation values.
TUBE_LAYOUTS = {
    "sim1": [
        # name in materials json,                centre (x_mm, y_mm), half-extent
        ("water_shibuya2020_n2",       80.000,   0.000),
        ("adipose_avachat2024",        56.569,  56.569),
        ("adipose_moskal2023",          0.000,  80.000),
        ("hepatic_avachat2024",       -56.569,  56.569),
        ("muscle_avachat2024",        -80.000,   0.000),
        ("brain_healthy_moskal2024",  -56.569, -56.569),
        ("brain_tumor_moskal2024",      0.000, -80.000),
        ("myxoma_moskal2023",          56.569, -56.569),
    ],
    "sim2": [
        # Mercolli 4 tubes (see gate/geometry/phantom_mercolli_xad4.mac)
        ("XAD4_dry",      -45.0, 0.0),
        ("XAD4_wet",      -15.0, 0.0),
        ("XAD4_gelatin",  +15.0, 0.0),
        ("XAD4_water",    +45.0, 0.0),
    ],
}

# Tube assignment radius: per config, slightly larger than the geometric tube
# radius to allow for positron-range smearing of source position. Sim 1 tubes
# are r=15 mm. Sim 2 tubes are r=6 mm with adjacent centres 30 mm apart so
# 9 mm assignment leaves a 12 mm guard band.
TUBE_R_MM_BY_CONFIG = {"sim1": 18.0, "sim2": 9.0}


def assign_tube(src_x: np.ndarray, src_y: np.ndarray,
                tubes: list[tuple[str, float, float]],
                tube_r_mm: float) -> np.ndarray:
    """Assign each triple to its origin tube (-1 if no tube within tube_r_mm)."""
    assignment = np.full(len(src_x), -1, dtype=np.int8)
    for i, (_name, cx, cy) in enumerate(tubes):
        d2 = (src_x - cx) ** 2 + (src_y - cy) ** 2
        mask = d2 <= tube_r_mm ** 2
        assignment[mask] = i
    return assignment


def mle_tau(tau_arr: np.ndarray, tau_min: float) -> tuple[float, float, int]:
    """Standard exponential MLE on a censored sample.
    tau_arr: corrected tau in ns (one per triple).
    Returns (tau_estimate, sigma, n_used).
    """
    sample = tau_arr[(tau_arr > tau_min) & np.isfinite(tau_arr)]
    n = sample.size
    if n < 20:
        return float("nan"), float("nan"), n
    tau_o = float(sample.mean()) - tau_min
    sigma = tau_o / np.sqrt(n)
    return tau_o, sigma, n


def load_materials_json(path: Path) -> dict:
    with open(path) as f:
        d = json.load(f)
    return {m["name"]: m for m in d["materials"]}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--detector", required=True, type=Path)
    ap.add_argument("--truth",    required=True, type=Path)
    ap.add_argument("--tube-config", required=True, choices=["sim1", "sim2"])
    ap.add_argument("--tau-min", type=float, default=0.3)
    ap.add_argument("--materials-json", type=Path,
                    default=REPO / "materials" / "pli_materials.json")
    ap.add_argument("--out-md",  required=True, type=Path)
    ap.add_argument("--out-csv", required=True, type=Path)
    ap.add_argument("--label",   default="")
    args = ap.parse_args()

    print(f"[step5] detector: {args.detector}")
    print(f"[step5] truth:    {args.truth}")
    print(f"[step5] tube_config: {args.tube_config}  tau_min={args.tau_min}")

    # Travel-time correction: call existing implementation.
    print("[step5] applying TOF travel-time correction (correct_travel_time.correct)")
    tau, P_annih = travel_time_correct(args.detector)
    print(f"[step5]   N triples = {len(tau):,}")
    print(f"[step5]   tau median = {np.median(tau):.3f} ns; "
          f"frac >0 = {(tau > 0).mean():.3f}")

    # Truth source positions
    with h5py.File(args.truth, "r") as ft:
        src_x = ft["truth"]["source_x"][:]
        src_y = ft["truth"]["source_y"][:]

    tubes = TUBE_LAYOUTS[args.tube_config]
    tube_r_mm = TUBE_R_MM_BY_CONFIG[args.tube_config]
    assignment = assign_tube(src_x, src_y, tubes, tube_r_mm)
    print(f"[step5]   tube assignment counts (tube_r={tube_r_mm} mm):")
    for i, (name, cx, cy) in enumerate(tubes):
        n = int((assignment == i).sum())
        print(f"     T{i+1} {name:30s} ({cx:+.1f},{cy:+.1f}): {n:,}")
    n_unassigned = int((assignment == -1).sum())
    print(f"     unassigned: {n_unassigned:,}")

    # Per-tube MLE
    materials = load_materials_json(args.materials_json)
    print(f"[step5] per-tube MLE (tau_min={args.tau_min} ns):")

    rows = []
    for i, (name, cx, cy) in enumerate(tubes):
        m = materials.get(name, {})
        tau_pub = m.get("tau_pub_ns", float("nan"))
        tau_TE  = m.get("tau_TE_ns",  float("nan"))
        tau_unc = m.get("tau_pub_uncertainty_ns",
                        m.get("tau_err_ns") or float("nan"))
        I3       = m.get("I3", float("nan"))
        tau_src  = m.get("tau_source", "")
        idx_mask = (assignment == i)
        tau_recovered, sigma, n_used = mle_tau(tau[idx_mask], args.tau_min)
        residual_pct = (tau_recovered - tau_TE) / tau_TE * 100 if np.isfinite(tau_recovered) else float("nan")
        n_total = int(idx_mask.sum())
        print(f"     T{i+1:1d} {name[:25]:25s}  "
              f"cfg={tau_TE:.3f}  rec={tau_recovered:.3f}+/-{sigma:.4f}  "
              f"resid={residual_pct:+.2f}%  N={n_used:,}/{n_total:,}")
        rows.append(dict(
            tube=f"T{i+1}", material=name, tau_pub_lit=tau_pub,
            tau_pub_unc=tau_unc, I3=I3,
            tau_TE_configured=tau_TE,
            tau_recovered=tau_recovered, sigma=sigma,
            residual_pct=residual_pct,
            n_used=n_used, n_total=n_total, tau_source=tau_src,
        ))

    # CSV
    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    keys = ["tube", "material", "tau_pub_lit", "tau_pub_unc", "I3",
            "tau_TE_configured", "tau_recovered", "sigma", "residual_pct",
            "n_used", "n_total", "tau_source"]
    with open(args.out_csv, "w") as f:
        f.write(",".join(keys) + "\n")
        for row in rows:
            f.write(",".join(str(row[k]).replace(",", ";") for k in keys) + "\n")
    print(f"[step5] wrote {args.out_csv}")

    # Markdown
    with open(args.out_md, "w") as f:
        f.write(f"# Step 5 validation table — {args.label}\n\n")
        f.write(f"Source: `{args.detector.name}` + `{args.truth.name}`  ")
        f.write(f"`tau_min={args.tau_min}` ns  `tube_R={tube_r_mm}` mm\n\n")
        f.write("| tube | material | tau_pub (ns) | unc | I3 | tau_TE_cfg | "
                "tau_recovered | sigma | resid % | N | source |\n")
        f.write("|---|---|---|---|---|---|---|---|---|---|---|\n")
        for r in rows:
            unc = f"±{r['tau_pub_unc']}" if np.isfinite(r['tau_pub_unc']) else "—"
            f.write(f"| {r['tube']} | `{r['material']}` | "
                    f"{r['tau_pub_lit']:.4f} | {unc} | {r['I3']:.3f} | "
                    f"{r['tau_TE_configured']:.4f} | "
                    f"{r['tau_recovered']:.4f} | "
                    f"{r['sigma']:.4f} | "
                    f"{r['residual_pct']:+.2f} | "
                    f"{r['n_used']:,} | {r['tau_source']} |\n")
        # Summary stat
        valid = [r for r in rows if np.isfinite(r["residual_pct"])]
        if valid:
            mean_abs = np.mean([abs(r["residual_pct"]) for r in valid])
            f.write(f"\n**mean |residual| across {len(valid)} tubes: "
                    f"{mean_abs:.3f} %**\n")
            print(f"[step5] mean |residual| = {mean_abs:.3f} %")
    print(f"[step5] wrote {args.out_md}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
