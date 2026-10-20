"""Task B v2: Per-voxel-then-per-tissue aggregation.

The reviewer-relevant measurement is the per-voxel τ̂ map (computed at each
grid size), aggregated by tissue label via per-voxel modal assignment from
the 4 mm atlas. This isolates partial-volume-at-voxel-level: a coarser grid
mixes more atlas cells per voxel → more PVE per voxel → more bias.
"""
from __future__ import annotations
import os
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

RES = Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[2])
P3Q = Path(os.environ.get("PLI_RAW", "<missing>")) / "2026-08-13_overnight/O4_rerun/vipman"
RUN_B = RES / "simple/run_B_4mm_10iter_geosens_norm"
OUT = RES / "2026-08-13_overnight/V1_vipman_postcorrection"

DET_H5 = P3Q / "triples_detector_phase3q.h5"
TRU_H5 = P3Q / "triples_truth_phase3q.h5"
TL_4MM_PATH = RUN_B / "tissue_label_3d.npy"

EXT_X = (-280.0, 280.0)
EXT_Y = (-140.0, 140.0)
EXT_Z = (-40.0, 40.0)
VOX_4 = 4.0

CFG_TAU_NS = {
    "Brain": 2.70, "Liver": 2.04, "Body": 2.00, "Blood": 1.80,
    "Lung": 1.85, "Heart": 2.00, "Kidney": 2.00, "Pancreas": 2.00,
    "Spleen": 2.00, "Intestine": 2.00, "Testis": 2.00,
    "Cartilage": 2.10, "SpineBone": 1.60, "RibBone": 1.60,
    "Muscle": 2.03, "Adipose": 2.55, "Air": 1.00,
}

TAU_MIN_NS = 1.0
R_CUT_MM = 75.0
N_MIN_PER_VOXEL = 30
MIN_VOXELS_PER_TISSUE = 5    # report tissues with at least this many valid voxels

VOXEL_SIZES = [4.0, 6.0, 8.0]


def voxel_tau_map(p_ann, tau, voxel_mm):
    """Compute per-voxel cumulant MLE map at the given grid size."""
    NX = int(np.floor((EXT_X[1] - EXT_X[0]) / voxel_mm))
    NY = int(np.floor((EXT_Y[1] - EXT_Y[0]) / voxel_mm))
    NZ = int(np.floor((EXT_Z[1] - EXT_Z[0]) / voxel_mm))
    i = np.floor((p_ann[:, 0] - EXT_X[0]) / voxel_mm).astype(np.int32)
    j = np.floor((p_ann[:, 1] - EXT_Y[0]) / voxel_mm).astype(np.int32)
    k = np.floor((p_ann[:, 2] - EXT_Z[0]) / voxel_mm).astype(np.int32)
    in_grid = (i >= 0) & (i < NX) & (j >= 0) & (j < NY) & (k >= 0) & (k < NZ)
    in_range = (tau > TAU_MIN_NS) & np.isfinite(tau)
    keep = in_grid & in_range
    flat = i[keep] * NY * NZ + j[keep] * NZ + k[keep]
    tau_sh = tau[keep] - TAU_MIN_NS
    SIZE = NX * NY * NZ
    sum_tau = np.bincount(flat, weights=tau_sh, minlength=SIZE).reshape(NX, NY, NZ)
    n_eff = np.bincount(flat, minlength=SIZE).reshape(NX, NY, NZ).astype(np.int32)
    mask_valid = n_eff >= N_MIN_PER_VOXEL
    with np.errstate(divide="ignore", invalid="ignore"):
        tau_map = np.where(mask_valid, sum_tau / n_eff, np.nan)
    return tau_map, n_eff, mask_valid, (NX, NY, NZ)


def voxel_to_tissue_map(tl_4mm, voxel_mm, dims):
    """For each voxel at given grid size, compute the modal atlas tissue label.

    Each voxel of size voxel_mm covers (voxel_mm / 4)³ atlas cells in
    principle (with possible fractional overlap at grid boundaries). We
    use a sampling approach: for each voxel center, look up the 4mm-atlas
    cell containing it.
    """
    NX, NY, NZ = dims
    out = np.empty((NX, NY, NZ), dtype=object)
    NX_4 = tl_4mm.shape[0]
    NY_4 = tl_4mm.shape[1]
    NZ_4 = tl_4mm.shape[2]
    for i in range(NX):
        for j in range(NY):
            for k in range(NZ):
                # Find all 4mm-atlas cells overlapping this voxel
                x_lo = EXT_X[0] + i * voxel_mm
                x_hi = EXT_X[0] + (i + 1) * voxel_mm
                y_lo = EXT_Y[0] + j * voxel_mm
                y_hi = EXT_Y[0] + (j + 1) * voxel_mm
                z_lo = EXT_Z[0] + k * voxel_mm
                z_hi = EXT_Z[0] + (k + 1) * voxel_mm
                ix_lo = max(0, int((x_lo - EXT_X[0]) / VOX_4))
                ix_hi = min(NX_4, int(np.ceil((x_hi - EXT_X[0]) / VOX_4)))
                iy_lo = max(0, int((y_lo - EXT_Y[0]) / VOX_4))
                iy_hi = min(NY_4, int(np.ceil((y_hi - EXT_Y[0]) / VOX_4)))
                iz_lo = max(0, int((z_lo - EXT_Z[0]) / VOX_4))
                iz_hi = min(NZ_4, int(np.ceil((z_hi - EXT_Z[0]) / VOX_4)))
                block = tl_4mm[ix_lo:ix_hi, iy_lo:iy_hi, iz_lo:iz_hi]
                if block.size == 0:
                    out[i, j, k] = "OutOfFOV"
                    continue
                vals, counts = np.unique(block, return_counts=True)
                out[i, j, k] = vals[counts.argmax()]
    return out


def per_tissue_summary(tau_map, n_eff, mask_valid, tissue_labels, vox_mm):
    """For each unique tissue label, summarize the per-voxel τ̂ values."""
    rows = []
    unique_t = sorted(set(tissue_labels.ravel()) - {"OutOfFOV"})
    for t in unique_t:
        sel = (tissue_labels == t) & mask_valid
        n_vox = int(sel.sum())
        if n_vox < MIN_VOXELS_PER_TISSUE:
            continue
        tau_vals = tau_map[sel]
        median_n_per_vox = int(np.median(n_eff[sel]))
        # Per-voxel median + mean
        median_tau = float(np.median(tau_vals))
        mean_tau = float(tau_vals.mean())
        std_tau = float(tau_vals.std())
        # Combined cumulant MLE on all triples in tissue voxels: weight by n_eff
        weighted_tau_hat = float((tau_vals * n_eff[sel]).sum() / n_eff[sel].sum())
        cfg = CFG_TAU_NS.get(t, float("nan"))
        residual_pct = 100 * (weighted_tau_hat - cfg) / cfg if cfg > 0 else float("nan")
        rows.append(dict(
            voxel_mm=vox_mm,
            tissue=t,
            cfg_tau_eff_ns=cfg,
            n_voxels=n_vox,
            median_n_per_voxel=median_n_per_vox,
            median_tau_hat_ns=median_tau,
            mean_tau_hat_ns=mean_tau,
            std_tau_hat_ns=std_tau,
            weighted_tau_hat_ns=weighted_tau_hat,
            residual_pct=residual_pct,
        ))
    return rows


def main():
    if not Path(DET_H5).exists(): raise SystemExit(f"{DET_H5}: raw data not included in the release, see analysis/README.md")
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    print(f"[{time.strftime('%H:%M:%S')}] loading P_annih_tof + tau_corrected")
    with h5py.File(DET_H5, "r") as f:
        P = f["derived"]["P_annih_tof_mm"][:].astype(np.float32)
        tau = f["derived"]["tau_corrected_ns"][:].astype(np.float32)
    print(f"  N = {len(tau):,}; elapsed {time.time()-t0:.1f}s")

    print(f"[{time.strftime('%H:%M:%S')}] loading truth source for R_cut")
    with h5py.File(TRU_H5, "r") as f:
        sx = f["truth"]["source_x"][:].astype(np.float32)
        sy = f["truth"]["source_y"][:].astype(np.float32)
        sz = f["truth"]["source_z"][:].astype(np.float32)

    dist2 = (P[:, 0] - sx)**2 + (P[:, 1] - sy)**2 + (P[:, 2] - sz)**2
    keep = dist2 <= R_CUT_MM**2
    n_kept = int(keep.sum())
    print(f"  R_cut={R_CUT_MM} mm: {n_kept:,} / {len(tau):,} = {100*n_kept/len(tau):.2f} %")

    del sx, sy, sz, dist2
    P = P[keep]
    tau = tau[keep]

    print(f"[{time.strftime('%H:%M:%S')}] loading 4 mm tissue label atlas")
    tl_4mm = np.load(TL_4MM_PATH, allow_pickle=True)
    print(f"  shape={tl_4mm.shape} dtype={tl_4mm.dtype} unique={np.unique(tl_4mm)}")

    all_rows = []
    for vox_mm in VOXEL_SIZES:
        print(f"\n[{time.strftime('%H:%M:%S')}] === VOXEL {vox_mm} mm ===")
        tau_map, n_eff, mask_valid, dims = voxel_tau_map(P, tau, vox_mm)
        n_v = int(mask_valid.sum())
        med = float(np.nanmedian(tau_map))
        print(f"  grid {dims}: valid voxels {n_v:,}/{np.prod(dims):,}  median τ̂={med:.3f} ns")

        print(f"  computing tissue labels at {vox_mm} mm grid...")
        tissue_labels = voxel_to_tissue_map(tl_4mm, vox_mm, dims)
        rows = per_tissue_summary(tau_map, n_eff, mask_valid, tissue_labels, vox_mm)
        for r in rows:
            print(f"  {r['tissue']:<10}  n_vox={r['n_voxels']:>5}  med_N={r['median_n_per_voxel']:>5}  "
                  f"τ̂_w={r['weighted_tau_hat_ns']:.3f}  cfg={r['cfg_tau_eff_ns']:.2f}  "
                  f"resid={r['residual_pct']:+6.2f}%")
        all_rows.extend(rows)

    df = pd.DataFrame(all_rows)
    df.to_csv(OUT / "B_voxel_sweep.csv", index=False)
    print(f"\nwrote {OUT}/B_voxel_sweep.csv ({len(df)} rows) in {(time.time()-t0)/60:.1f} min")


if __name__ == "__main__":
    main()
