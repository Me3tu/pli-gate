"""Phase 3n PART 2: TOF P_ann positron-escape filter for phantom triples.

Removes triples whose TOF-derived annihilation position is too far from
the per-event truth source position. Intended for phantom geometries
(Sim 1 multi-tissue Sc-44, Sim 2 I-124 XAD-4) where each event has a
discrete source position recorded in the truth H5. Phase 3m showed
this filter closes the Sim 2 T4 +13.65 % τ bias to within ±5 % of cfg
by removing positron-escape contamination, while preserving Sim 1
(low-escape) within ~2 % of baseline.

SCOPE LIMITATION
================
This filter is **phantom-specific**. It requires per-event discrete
source positions, which exist in `triples_truth_*.h5` for Sim 1 / Sim 2
but are meaningless for VIPman recons (voxelized activity distribution,
no per-event "source"). For VIPman, positron escape is sufficiently
small (< 1 % per Phase 3k cross-isotope check) that the filter is
unnecessary. If needed, a different criterion such as
`P_ann ∈ patient body` via µmap intersection would be required, and
that is Phase 3p territory.

USAGE
=====
    python reconstruction/filter_triples_tof_pann.py \\
        --triples-detector  PATH_IN_DETECTOR.h5 \\
        --triples-truth     PATH_IN_TRUTH.h5 \\
        --r-cut             75.0 \\
        --output-detector   PATH_OUT_DETECTOR.h5 \\
        --output-truth      PATH_OUT_TRUTH.h5

    # report only (no write), useful for R_cut sweeps
    python reconstruction/filter_triples_tof_pann.py \\
        --triples-detector PATH --triples-truth PATH \\
        --r-cut 75.0 --report-only

OUTPUT NAMING CONVENTION
========================
Suggested: ``<input_stem>_rcut075mm.h5`` alongside the input.

R_cut DEFAULT: 75 mm
====================
Brackets a 7-mm tube radius + ~4σ TOF blur (σ_LOR ≈ 19 mm from CRT
247 ps FWHM). Loose enough to keep TOF-noise-broadened in-tube events,
tight enough to cut LSO-ring escapes (radius ≈ 393 mm) and most
air-gap escapes. See Phase 3m findings.md for the empirical motivation.

SIZE LIMIT
==========
Single-batch in-memory read is acceptable up to ~50 M triples per H5
(~few GB peak). Beyond that, chunked I/O would be needed. Not yet
implemented because no phantom currently exceeds that scale (Sim 1 =
3.1 M, Sim 2 = 2.1 M, Phase 2e VIPman time-grouped = 24.6 M).
"""
from __future__ import annotations
import argparse
import h5py
import numpy as np
from pathlib import Path

C_MM_NS = 299.792458
GATE_TIME_SCALE_NS = 1e9   # GATE singles "time" is in seconds, × 1e9 → ns


def detect_time_units(t1_raw, t2_raw):
    """Sanity-check time units. Print median |t1 - t2| in raw H5 units
    and decide whether GATE-style seconds-to-ns conversion is needed.

    For CRT ≈ 250 ps FWHM, expected median |t1-t2| is ~150-400 ps. In
    seconds: ~1.5e-10 to 4e-10. In ns: ~0.15 to 0.4. Three orders of
    magnitude apart, so the discrimination is robust.
    """
    raw_dt = np.abs(t1_raw - t2_raw)
    raw_med = float(np.median(raw_dt))
    print(f"  startup sanity: median |t1 - t2| (raw H5 units) = {raw_med:.4e}")
    if raw_med < 1e-7:
        # Far below 100 ns, consistent with seconds-storage (typical
        # 1.5e-10 to 4e-10 for the phantoms we run). Apply × 1e9
        # conversion.
        print(f"  → interpreted as SECONDS; applying × {GATE_TIME_SCALE_NS:.0e} → ns "
               f"(expected median |Δt| ≈ {raw_med * GATE_TIME_SCALE_NS:.3f} ns)")
        return "seconds"
    if 1e-3 < raw_med < 1e2:
        # Plausible ns-storage range (e.g., median 0.247 ns).
        print(f"  → interpreted as NANOSECONDS; NO conversion applied "
               f"(would have been × {GATE_TIME_SCALE_NS:.0e} too large)")
        return "nanoseconds"
    raise ValueError(
        f"Time-units sanity check failed: median |t1 - t2| = {raw_med:.4e} "
        "is in neither the 'seconds' regime (< 1e-7) nor the 'ns' regime "
        "(1e-3 to 1e2). Inspect the input H5; aborting before silent "
        "unit-conversion error.")


def compute_pann_mm(p1, p2, t1_ns, t2_ns):
    mp = 0.5 * (p1 + p2)
    lor = p2 - p1
    ll = np.linalg.norm(lor, axis=-1, keepdims=True)
    lu = lor / np.maximum(ll, 1e-9)
    off = 0.5 * C_MM_NS * (t1_ns - t2_ns)[:, None]
    return mp + off * lu


def filter_h5(det_in: Path, tru_in: Path,
                det_out: Path | None, tru_out: Path | None,
                r_cut: float, report_only: bool) -> dict:
    print(f"Filter TOF P_ann at R_cut = {r_cut:.2f} mm")
    print(f"  detector in : {det_in}")
    print(f"  truth in    : {tru_in}")

    with h5py.File(det_in, "r") as fd, h5py.File(tru_in, "r") as ft:
        d = fd["detector"]
        t = ft["truth"]
        n = d["a1_x"].shape[0]
        print(f"  N input     : {n:,}")

        # Time units sanity (uses a sample to be cheap)
        sample_size = min(50_000, n)
        units = detect_time_units(d["a1_time"][:sample_size],
                                    d["a2_time"][:sample_size])
        scale = GATE_TIME_SCALE_NS if units == "seconds" else 1.0

        # Load full arrays
        p1 = np.stack([d["a1_x"][:], d["a1_y"][:], d["a1_z"][:]], axis=-1)
        p2 = np.stack([d["a2_x"][:], d["a2_y"][:], d["a2_z"][:]], axis=-1)
        t1_ns = d["a1_time"][:] * scale
        t2_ns = d["a2_time"][:] * scale
        src = np.stack([t["source_x"][:], t["source_y"][:], t["source_z"][:]], axis=-1)

        pann = compute_pann_mm(p1, p2, t1_ns, t2_ns)
        dist = np.linalg.norm(pann - src, axis=-1)
        keep = dist <= r_cut
        n_kept = int(keep.sum())
        retention = n_kept / n if n else 0.0

        print(f"  N kept      : {n_kept:,}  ({retention*100:.2f} %)")
        print(f"  N dropped   : {n - n_kept:,}")
        if n_kept > 0:
            print(f"  median |P_ann−src| (kept)   : {np.median(dist[keep]):.2f} mm")
        if n_kept < n:
            print(f"  median |P_ann−src| (dropped): {np.median(dist[~keep]):.2f} mm")

        report = dict(n_input=n, n_kept=n_kept, retention=float(retention),
                       r_cut_mm=float(r_cut),
                       median_pann_src_kept_mm=(float(np.median(dist[keep])) if n_kept > 0 else None),
                       median_pann_src_dropped_mm=(float(np.median(dist[~keep])) if n_kept < n else None),
                       time_units_detected=units)

        if report_only:
            return report

        if det_out is None or tru_out is None:
            raise ValueError("output paths required when not --report-only")

        det_out.parent.mkdir(parents=True, exist_ok=True)
        tru_out.parent.mkdir(parents=True, exist_ok=True)

        # Detector output
        with h5py.File(det_out, "w") as fdo:
            g = fdo.create_group("detector")
            for k in d.keys():
                g.create_dataset(k, data=d[k][:][keep])
            for k, v in d.attrs.items():
                g.attrs[k] = v
            if "derived" in fd.keys():
                gd = fdo.create_group("derived")
                for k in fd["derived"].keys():
                    arr = fd["derived"][k][:]
                    gd.create_dataset(k, data=arr[keep])
                for k, v in fd["derived"].attrs.items():
                    gd.attrs[k] = v
            mg = fdo.create_group("filter_metadata")
            mg.attrs["filter"] = "tof_pann"
            mg.attrs["r_cut_mm"] = float(r_cut)
            mg.attrs["n_input"] = n
            mg.attrs["n_kept"] = n_kept
            mg.attrs["retention"] = float(retention)
            mg.attrs["source_h5_in"] = str(det_in)
            mg.attrs["filter_phase"] = "3n"
            mg.attrs["time_units_detected"] = units

        # Truth output
        with h5py.File(tru_out, "w") as fto:
            g = fto.create_group("truth")
            for k in t.keys():
                g.create_dataset(k, data=t[k][:][keep])
            for k, v in t.attrs.items():
                g.attrs[k] = v
            mg = fto.create_group("filter_metadata")
            mg.attrs["filter"] = "tof_pann"
            mg.attrs["r_cut_mm"] = float(r_cut)
            mg.attrs["n_input"] = n
            mg.attrs["n_kept"] = n_kept
            mg.attrs["retention"] = float(retention)
            mg.attrs["source_h5_in"] = str(tru_in)
            mg.attrs["filter_phase"] = "3n"
            mg.attrs["time_units_detected"] = units

        print(f"  wrote {det_out}")
        print(f"  wrote {tru_out}")

    return report


def _cli():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--triples-detector", required=True, type=Path)
    ap.add_argument("--triples-truth",    required=True, type=Path)
    ap.add_argument("--r-cut",            type=float, default=75.0,
                     help="Distance cut |P_ann_TOF - source| in mm. Default 75 mm.")
    ap.add_argument("--output-detector",  type=Path, default=None)
    ap.add_argument("--output-truth",     type=Path, default=None)
    ap.add_argument("--report-only",      action="store_true",
                     help="Print retention stats; do not write filtered H5.")
    args = ap.parse_args()
    if not args.report_only and (args.output_detector is None or args.output_truth is None):
        ap.error("--output-detector and --output-truth required unless --report-only")
    filter_h5(args.triples_detector, args.triples_truth,
               args.output_detector, args.output_truth,
               args.r_cut, args.report_only)


if __name__ == "__main__":
    _cli()
