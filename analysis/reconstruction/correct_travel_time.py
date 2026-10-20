"""
Apply TOF-based travel-time correction to the extracted triples.

For each triple:
  1. Estimate annihilation position P from TOF of the two 511s:
     P = midpoint(p1, p2) + (c/2) * (t1 - t2) * unit(p1 - p2)
  2. Compute t_annih = t1 - |p1 - P| / c  (equivalent via p2)
  3. Compute t_decay = t_prompt - |p_prompt - P| / c
     (prompt emitted at approximately P since positron range is sub-mm for Sc-44)
  4. Corrected τ = t_annih - t_decay

Writes corrected_tau back into triples_detector_vipman_500M.h5 as a new dataset.
No truth fields are read.
"""
from __future__ import annotations
import h5py
import numpy as np
from pathlib import Path

C_MM_NS = 299.792458  # speed of light in mm/ns
GATE_TIME_SCALE_NS = 1e9  # GATE 9.4 stores time in seconds; convert to ns

DET = Path("data/triples_detector_vipman_500M.h5")


def correct(det_path: Path) -> tuple[np.ndarray, np.ndarray]:
    with h5py.File(det_path, "r") as f:
        g = f["detector"]
        p1 = np.stack([g["a1_x"][:], g["a1_y"][:], g["a1_z"][:]], axis=-1)
        p2 = np.stack([g["a2_x"][:], g["a2_y"][:], g["a2_z"][:]], axis=-1)
        pp = np.stack([g["p_x"][:],  g["p_y"][:],  g["p_z"][:]],  axis=-1)
        t1 = g["a1_time"][:] * GATE_TIME_SCALE_NS
        t2 = g["a2_time"][:] * GATE_TIME_SCALE_NS
        tp = g["p_time"][:]  * GATE_TIME_SCALE_NS

    # 1. Annihilation position P from TOF (along the LOR p1->p2)
    midpoint = 0.5 * (p1 + p2)
    lor = p2 - p1
    lor_len = np.linalg.norm(lor, axis=-1, keepdims=True)
    lor_unit = lor / np.maximum(lor_len, 1e-9)
    # TOF offset along LOR, from midpoint toward p2 if t1 > t2
    offset = 0.5 * C_MM_NS * (t1 - t2)[:, None]
    P = midpoint + offset * lor_unit

    # 2. Annihilation time, using a1 (symmetrically equivalent with a2)
    d1 = np.linalg.norm(p1 - P, axis=-1)
    t_annih = t1 - d1 / C_MM_NS

    # 3. Prompt decay time, assuming prompt is emitted at P
    dp = np.linalg.norm(pp - P, axis=-1)
    t_decay = tp - dp / C_MM_NS

    # 4. Corrected τ
    tau = t_annih - t_decay
    return tau, P


def main():
    tau, P_annih = correct(DET)

    with h5py.File(DET, "a") as f:
        if "derived" in f:
            del f["derived"]
        d = f.create_group("derived")
        d.create_dataset("tau_corrected_ns", data=tau,
                         compression="gzip", compression_opts=4)
        d.create_dataset("P_annih_tof_mm", data=P_annih,
                         compression="gzip", compression_opts=4)
        d.attrs["note"] = ("TOF-corrected τ and TOF-derived annihilation "
                           "position. Detector-only.")

    print(f"[correct] wrote {len(tau):,} corrected τ values")
    print(f"[correct] τ median = {np.median(tau):.3f} ns")
    print(f"[correct] τ mean   = {np.mean(tau):.3f} ns")
    print(f"[correct] τ > 0 fraction: {(tau > 0).mean():.3f}")
    print(f"[correct] τ in [0.5, 10]: {((tau > 0.5) & (tau < 10)).mean():.3f}")
    print(f"[correct] τ percentiles  5/50/95: "
          f"{np.percentile(tau, 5):.3f} / "
          f"{np.percentile(tau, 50):.3f} / "
          f"{np.percentile(tau, 95):.3f} ns")


if __name__ == "__main__":
    main()
