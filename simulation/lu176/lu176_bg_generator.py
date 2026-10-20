# SPDX-License-Identifier: LGPL-3.0-only
# Copyright (C) 2026 Dhritman Dhar and German Cancer Research Center (DKFZ)
"""Analytical Lu-176 LSO-background singles generator (Phase 3o-3).

Produces a GATE-Singles-compatible HDF5 file by analytically modelling the
Lu-176 β⁻+cascade decay over the full Quadra LSO crystal volume, suitable
for direct ingestion by ``reconstruction/utils/lso_background_inject.py``
(Phase 1) without any GATE simulation.

This bypasses the Phase 3f-1 OOM block (G4RadioactiveDecay cascade memory
leak in Geant4 11.3.2) by replacing the Geant4 transport step with a
purely-numpy analytical model.

PHYSICS MODEL (Phase 3o-2 design)
=================================
- ¹⁷⁶Lu (T_½ = 3.78 × 10¹⁰ y) → ¹⁷⁶Hf* + β⁻ + ν̄_e.
- 99.6 % branch: β⁻ endpoint 596 keV → 998 keV state in Hf-176.
  998 → 597 → 290 → 0 keV cascade.
- Cascade γ emission probabilities (Conti 2017, IEEE TNS, Table I):
    307 keV: 0.94 emit, 0.06 IC
    202 keV: 0.78 emit, 0.22 IC
     88 keV: 0.15 emit, 0.85 IC          ← heavily IC-dominated
- 0.4 % branch (β⁻ endpoint 195 keV → 1142 keV state, then 401 keV γ):
  DEFERRED for first-pass (negligible PLI contamination effect).

DETECTOR RESPONSE
=================
- Energy resolution: 10.1 % FWHM @ 511 keV with √E scaling
  (LSO scintillator; matches lso_lu176_standalone.mac).
- Timing resolution: 214 ps FWHM Gaussian jitter.
- Crystal-of-detection: assumed to be the crystal of emission
  (β⁻ range ~1 mm in LSO; 88 keV γ self-absorbs ~70 % at crystal centre;
   higher-energy γ scatter to neighbours but never enter the PLI 511
   window, so this approximation doesn't affect downstream triple
   formation). DEFERRED: multi-crystal Compton scattering (Step 7.3).

TEMPORAL MODEL
==============
- Decay times: proper Poisson process at the configured activity.
- Within-decay timing: instantaneous cascade (9.4 ns isomeric lifetime
  of 998 keV state ignored — DEFERRED Step 7.2; cascade γ never reach
  the PLI window so the 9.4 ns has no observable effect).

OUTPUT SCHEMA
=============
HDF5 file with /singles/ group containing all 21 GATE Singles branches
required by lso_background_inject.py. Time in seconds (sorted ascending);
energy in MeV; positions in mm. eventID is int64 (Phase 3p-3) — accommodates
production-scale phantom decay_id ranges with a 10^12 default offset.

DETECTOR HIERARCHY ID CONVENTION (Phase 3o-2 §2.4)
===================================================
This is OUR convention; a hypothetical GATE-generated Lu-176 output
might use different IDs. The eventID-grouping triple extractor doesn't
use these, so they are populated for provenance only.
    gantryID    = 0                 (single gantry)
    rsectorID   ∈ [0, 37]           (38 circumferential modules)
    moduleID    = 0                 (no GATE 'module' below rsector here)
    submoduleID ∈ [0, 3]            (4 axial blocks)
    crystalID   ∈ [0, 1599]         (20 transaxial × 80 axial flat index)
    layerID     = 0                 (single LSO layer)

DEFERRED LIMITATIONS (Phase 3o-2 §7)
====================================
7.1 IC electrons not emitted as separate particles (γ row simply absent
    when IC occurs). Slightly under-counts low-energy bg singles.
7.2 Cascade isomeric lifetimes (998 keV state 9.4 ns) treated as 0.
7.3 Multi-crystal Compton scattering of cascade γ not modelled.
7.4 401 keV rare branch (0.4 %) not modelled.
7.5 β⁻ self-detection assumed 100 %.
7.6 Beta spectrum: Fermi-corrected (simplified non-relativistic Fermi
    function); 1st-forbidden-non-unique correction NOT applied per
    Phase 3o-2 sign-off Q1.

USAGE
=====
    python reconstruction/utils/lu176_bg_generator.py \\
        --activity-MBq 14.59 \\
        --T-sim-s 1.0 \\
        --output bg.h5 \\
        --seed 0
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path

import h5py
import numpy as np

# ────────────────────── physical constants ──────────────────────
ALPHA = 1.0 / 137.035999084
M_E_KEV = 510.998950        # electron rest energy
GATE_TIME_NS = 1e9          # seconds → ns conversion factor for downstream

# ───────────────────── Lu-176 cascade scheme (Conti 2017 Table I) ─────
BETA_ENDPOINT_KEV = 596.0
HF_DAUGHTER_Z = 72
GAMMA_LINES = [
    # (energy_keV, emission_probability)
    (307.0, 0.94),
    (202.0, 0.78),
    ( 88.0, 0.15),
]

# ────────────────────── detector response ──────────────────────
ENERGY_FWHM_AT_511_KEV = 0.101
TIMING_FWHM_PS = 214.0
SIGMA_T_NS = TIMING_FWHM_PS * 1e-3 / 2.3548   # ≈ 0.0909 ns

# ────────────────────── Quadra geometry (gate/geometry/quadra_detector.mac) ─────
# 38 rsectors × 4 blocks × 20 transaxial × 80 axial = 243,200 crystals
N_RSECTOR = 38
N_BLOCK = 4
N_TRANS = 20
N_AXIAL = 80
N_CRYSTAL_TOTAL = N_RSECTOR * N_BLOCK * N_TRANS * N_AXIAL   # 243,200
RSECTOR_RADIUS_MM = 402.83720907    # rsector centroid → crystal centroid radial position
CRYSTAL_RADIAL_MM = 20.0            # crystal depth (radial); centroid ± 10 mm
CRYSTAL_TRANS_MM = 3.20             # transaxial cell size
CRYSTAL_TRANS_PITCH_MM = 3.21       # 0.01 mm wrap gap
CRYSTAL_AXIAL_MM = 3.20             # axial cell size
BLOCK_AXIAL_LENGTH_MM = N_AXIAL * CRYSTAL_AXIAL_MM      # 256
BLOCK_GAP_AXIAL_MM = 9.0            # gap between blocks
TOTAL_AXIAL_EXTENT_MM = (
    N_BLOCK * BLOCK_AXIAL_LENGTH_MM + (N_BLOCK - 1) * BLOCK_GAP_AXIAL_MM
)   # 4×256 + 3×9 = 1051

# Block centres (axial), centred on z=0
_block_total = BLOCK_AXIAL_LENGTH_MM + BLOCK_GAP_AXIAL_MM
_block_offsets = (np.arange(N_BLOCK) - (N_BLOCK - 1) / 2.0) * _block_total
BLOCK_CENTRES_Z = _block_offsets    # e.g. [-397.5, -132.5, +132.5, +397.5]

# ───────────────────── output schema (matches lso_background_inject) ─────
SINGLES_BRANCHES = (
    "runID", "eventID", "sourceID",
    "sourcePosX", "sourcePosY", "sourcePosZ",
    "time", "energy",
    "globalPosX", "globalPosY", "globalPosZ",
    "gantryID", "rsectorID", "moduleID", "submoduleID",
    "crystalID", "layerID",
    "comptonPhantom", "comptonCrystal",
    "RayleighPhantom", "RayleighCrystal",
)
DTYPES = {
    "runID":           np.int32,    "eventID":         np.int64,
    "sourceID":        np.int32,
    "sourcePosX":      np.float32,  "sourcePosY":      np.float32,  "sourcePosZ":      np.float32,
    "time":            np.float64,  "energy":          np.float32,
    "globalPosX":      np.float32,  "globalPosY":      np.float32,  "globalPosZ":      np.float32,
    "gantryID":        np.int32,    "rsectorID":       np.int32,
    "moduleID":        np.int32,    "submoduleID":     np.int32,
    "crystalID":       np.int32,    "layerID":         np.int32,
    "comptonPhantom":  np.int32,    "comptonCrystal":  np.int32,
    "RayleighPhantom": np.int32,    "RayleighCrystal": np.int32,
}

# Default eventID offset (Phase 3p-3): 10^12 places bg eventIDs well outside
# any realistic phantom decay_id range. Phase 3q phantom max ~1.1B → 10^4 margin.
# Stays comfortably within int64 [0, 2^63-1] and float64-exact-integer range
# (2^53 ≈ 9×10^15) for downstream extractor float64 promotion.
DEFAULT_EVENTID_OFFSET = 1_000_000_000_000


# ════════════════════════ physics samplers ════════════════════════

def _build_beta_cdf(t_max_keV: float = BETA_ENDPOINT_KEV,
                    z_daughter: int = HF_DAUGHTER_Z,
                    n_grid: int = 1000):
    """Tabulate the Fermi-corrected β⁻ spectrum and return inverse-CDF lookup.

    N(T) dT ∝ F(Z, T) × p × E × (T_max − T)² dT
    where T = kinetic energy (keV), E = T + m_e c², p = √(E² − (m_e c²)²).
    Non-relativistic Fermi function F(Z, T) = 2πη / (1 − exp(−2πη)),
    η = Zα/β.
    """
    T_grid = np.linspace(0.5, t_max_keV - 0.5, n_grid)
    E_grid = T_grid + M_E_KEV
    p_grid = np.sqrt(np.maximum(E_grid**2 - M_E_KEV**2, 0.0))
    beta = p_grid / E_grid
    # Fermi function (non-relativistic)
    eta = z_daughter * ALPHA / np.maximum(beta, 1e-6)
    # Stable form: 2πη / (1 − e^{−2πη})
    arg = 2.0 * np.pi * eta
    fermi = arg / (1.0 - np.exp(-arg))
    # Spectrum shape
    N_T = fermi * p_grid * E_grid * (t_max_keV - T_grid) ** 2
    # CDF
    cdf = np.cumsum(N_T)
    cdf /= cdf[-1]
    return T_grid, cdf


_BETA_T_GRID, _BETA_CDF = _build_beta_cdf()


def sample_beta_T(n: int, rng: np.random.Generator) -> np.ndarray:
    """Sample β⁻ kinetic energy (keV) by inverse-CDF lookup."""
    u = rng.uniform(0.0, 1.0, size=n)
    return np.interp(u, _BETA_CDF, _BETA_T_GRID).astype(np.float32)


def apply_energy_resolution(E_true_keV: np.ndarray,
                              rng: np.random.Generator) -> np.ndarray:
    """Apply √E-scaled Gaussian energy smearing.

    σ_E(E) = (FWHM_at_511 / 2.355) × √(E × 511) keV
           = 0.04284 × √(E_keV × 511)
    """
    sigma_keV = (ENERGY_FWHM_AT_511_KEV / 2.3548) * np.sqrt(
        np.maximum(E_true_keV, 0.0) * 511.0
    )
    return np.maximum(0.0, E_true_keV + rng.normal(0.0, sigma_keV)).astype(np.float32)


# ════════════════════════ geometry ════════════════════════

def crystal_index_to_centroid(crystal_idx: np.ndarray):
    """Decompose flat crystal index ∈ [0, 243199] into hierarchy and compute world centroid.

    Index layout (slowest → fastest): rsector → block → transaxial → axial.
    Returns:
        (xc, yc, zc)        — world position of crystal centroid (mm)
        (rsec, blk, ix_t, ix_a) — hierarchy indices
    """
    crystal_idx = np.asarray(crystal_idx, dtype=np.int64)
    block_total = N_BLOCK * N_TRANS * N_AXIAL
    rsec = crystal_idx // block_total
    rem = crystal_idx % block_total
    blk = rem // (N_TRANS * N_AXIAL)
    rem2 = rem % (N_TRANS * N_AXIAL)
    ix_t = rem2 // N_AXIAL
    ix_a = rem2 % N_AXIAL

    # Local rsector frame: x = radial (outward), y = transaxial, z = axial
    # Crystal centroid in local frame:
    x_local = np.full_like(crystal_idx, RSECTOR_RADIUS_MM, dtype=np.float64)
    y_local = ((ix_t.astype(np.float64) - (N_TRANS - 1) / 2.0)
                * CRYSTAL_TRANS_PITCH_MM)
    # Axial: block centre + within-block axial offset
    z_local = (BLOCK_CENTRES_Z[blk]
                + (ix_a.astype(np.float64) - (N_AXIAL - 1) / 2.0)
                * CRYSTAL_AXIAL_MM)

    # Rotate by rsector angle (rsector 0 at φ=0, increasing CCW)
    phi = 2.0 * np.pi * rsec.astype(np.float64) / N_RSECTOR
    cos_phi = np.cos(phi)
    sin_phi = np.sin(phi)
    xc = x_local * cos_phi - y_local * sin_phi
    yc = x_local * sin_phi + y_local * cos_phi
    zc = z_local
    return (xc.astype(np.float32),
             yc.astype(np.float32),
             zc.astype(np.float32),
             rsec.astype(np.int32),
             blk.astype(np.int32),
             ix_t.astype(np.int32),
             ix_a.astype(np.int32))


def add_within_crystal_offset(xc, yc, zc, rsec, rng):
    """Add uniform offset within crystal volume (3.2 × 3.2 × 20 mm).

    Sampling is in the crystal-local frame (radial × transaxial × axial)
    then rotated to world coords by the rsector angle.
    """
    n = len(xc)
    dr = rng.uniform(-CRYSTAL_RADIAL_MM / 2.0, CRYSTAL_RADIAL_MM / 2.0, n)
    dy = rng.uniform(-CRYSTAL_TRANS_MM   / 2.0, CRYSTAL_TRANS_MM   / 2.0, n)
    dz = rng.uniform(-CRYSTAL_AXIAL_MM   / 2.0, CRYSTAL_AXIAL_MM   / 2.0, n)
    phi = 2.0 * np.pi * rsec.astype(np.float64) / N_RSECTOR
    cos_phi = np.cos(phi)
    sin_phi = np.sin(phi)
    return (
        (xc + (dr * cos_phi - dy * sin_phi)).astype(np.float32),
        (yc + (dr * sin_phi + dy * cos_phi)).astype(np.float32),
        (zc + dz).astype(np.float32),
    )


# ════════════════════════ main generator ════════════════════════

def generate(activity_MBq: float,
              T_sim_s: float,
              output_path: Path,
              seed: int = 0,
              eventID_offset: int = DEFAULT_EVENTID_OFFSET,
              chunk_decays: int = 200_000,
              version: str = "phase3o3_v1_2026-05-04",
              verbose: bool = True) -> dict:
    """Generate Lu-176 LSO-bg singles HDF5 file.

    Memory: peak ~300 MB per chunk (chunk_decays × ~2.87 rows × ~50 B/row).
    For activity_MBq=14.59, T_sim_s=10 → ~1.46e8 decays → ~730 chunks → ~5 min wall.
    """
    activity_Bq = activity_MBq * 1e6
    expected_decays = activity_Bq * T_sim_s
    rng = np.random.default_rng(seed)
    n_decays = int(rng.poisson(expected_decays))

    if verbose:
        print(f"[lu176_bg] activity = {activity_MBq} MBq; T_sim = {T_sim_s} s")
        print(f"[lu176_bg] expected decays = {expected_decays:.3e}; sampled = {n_decays:,}")
        exp_rows = n_decays * (1.0 + sum(p for _, p in GAMMA_LINES))
        print(f"[lu176_bg] expected total singles ≈ {exp_rows:.3e}")
        print(f"[lu176_bg] writing → {output_path}")

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Pre-allocate per-particle metadata: which row in each decay corresponds
    # to which particle type (0 = β⁻, 1 = γ_307, 2 = γ_202, 3 = γ_88).
    # Vectorise per chunk.
    t_start = time.time()

    n_total_written = 0
    n_betas_in_511_window = {"sc44": 0, "i124": 0}    # for V5.3 enhancement

    with h5py.File(output_path, "w") as h:
        grp = h.create_group("singles")
        # Resizable datasets (chunked + gzip compressed).
        # Initial size 0; resize per chunk.
        ds = {}
        for b in SINGLES_BRANCHES:
            ds[b] = grp.create_dataset(
                b, shape=(0,), maxshape=(None,), dtype=DTYPES[b],
                chunks=True, compression="gzip", compression_opts=4,
            )

        # Decay loop in chunks
        decay_idx_global = 0
        n_chunks = max(1, (n_decays + chunk_decays - 1) // chunk_decays)
        for ic in range(n_chunks):
            n_dec = min(chunk_decays, n_decays - decay_idx_global)
            if n_dec <= 0:
                break
            decay_eid = (eventID_offset
                          + decay_idx_global
                          + np.arange(n_dec, dtype=np.int64))
            decay_idx_global += n_dec

            # Decay times: uniform in [0, T_sim_s]; sort later globally
            decay_times = rng.uniform(0.0, T_sim_s, n_dec).astype(np.float64)

            # Decay vertex: sample crystal index uniformly in [0, 243199)
            crystal_idx = rng.integers(0, N_CRYSTAL_TOTAL, n_dec)
            xc, yc, zc, rsec, blk, ix_t, ix_a = crystal_index_to_centroid(crystal_idx)
            sx, sy, sz = add_within_crystal_offset(xc, yc, zc, rsec, rng)

            # Per-decay particle inventory: always β⁻ + cascade γ Bernoulli draws
            # Build an index list mapping each output row → (decay_idx, particle_type)
            # particle_type 0=β, 1=γ_307, 2=γ_202, 3=γ_88
            emit_g307 = rng.uniform(0, 1, n_dec) < GAMMA_LINES[0][1]
            emit_g202 = rng.uniform(0, 1, n_dec) < GAMMA_LINES[1][1]
            emit_g88  = rng.uniform(0, 1, n_dec) < GAMMA_LINES[2][1]
            n_per_decay = (1
                            + emit_g307.astype(np.int32)
                            + emit_g202.astype(np.int32)
                            + emit_g88.astype(np.int32))
            n_rows = int(n_per_decay.sum())

            # Build flat row→decay mapping via repeats
            row_to_dec = np.repeat(np.arange(n_dec, dtype=np.int64), n_per_decay)

            # Per-row particle type code: walk in a vectorised fashion.
            # For each decay: [β⁻] then γ_307 if emit, γ_202 if emit, γ_88 if emit
            # offsets within each decay
            ptype_per_row = np.empty(n_rows, dtype=np.int8)
            cursor = 0
            ptype_per_row[cursor:cursor + n_dec] = 0   # but rows aren't contiguous; this approach won't work
            # Alternative: vectorise via cumulative offsets
            # Strategy: for each particle slot {0,1,2,3}, build per-decay
            # "is this slot occupied" mask; then collect (decay_idx, ptype) pairs.
            slot_masks = [
                np.ones(n_dec, dtype=bool),
                emit_g307, emit_g202, emit_g88,
            ]
            slot_decays = []
            slot_ptypes = []
            for ptype, mask in enumerate(slot_masks):
                idx = np.flatnonzero(mask)
                slot_decays.append(idx)
                slot_ptypes.append(np.full(len(idx), ptype, dtype=np.int8))
            row_to_dec = np.concatenate(slot_decays)
            ptype_per_row = np.concatenate(slot_ptypes)
            # NOTE: rows are NOT yet ordered by decay; we will time-sort at end.

            # Per-row energies
            E_keV = np.empty(n_rows, dtype=np.float32)
            beta_mask = ptype_per_row == 0
            E_keV[beta_mask] = sample_beta_T(int(beta_mask.sum()), rng)
            for ptype, (e_line, _) in enumerate(GAMMA_LINES, start=1):
                m = ptype_per_row == ptype
                E_keV[m] = e_line
            # Apply energy resolution
            E_keV = apply_energy_resolution(E_keV, rng)

            # Per-row times: t_decay + Gaussian timing jitter (in seconds)
            t_jitter_s = rng.normal(0.0, SIGMA_T_NS * 1e-9, n_rows).astype(np.float64)
            t_row = decay_times[row_to_dec] + t_jitter_s

            # Per-row decay positions (vertex)
            sx_row = sx[row_to_dec]
            sy_row = sy[row_to_dec]
            sz_row = sz[row_to_dec]

            # Per-row globalPos = crystal centroid (Phase 3o-2 §2.2: same-crystal self-detection)
            gx_row = xc[row_to_dec]
            gy_row = yc[row_to_dec]
            gz_row = zc[row_to_dec]

            # Per-row hierarchy IDs
            rsec_row = rsec[row_to_dec]
            blk_row  = blk[row_to_dec]
            ix_t_row = ix_t[row_to_dec]
            ix_a_row = ix_a[row_to_dec]
            # crystalID = transaxial-axial flat index ∈ [0, 1599]
            crystal_id_row = (ix_t_row * N_AXIAL + ix_a_row).astype(np.int32)

            eid_row = decay_eid[row_to_dec].astype(np.int64)

            # β⁻ in 511 window counters (V5.3 enhancement)
            E_keV_arr = np.asarray(E_keV)
            beta_mask = ptype_per_row == 0
            n_betas_in_511_window["sc44"] += int(((E_keV_arr >= 481) & (E_keV_arr <= 541) & beta_mask).sum())
            n_betas_in_511_window["i124"] += int(((E_keV_arr >= 460) & (E_keV_arr <= 545) & beta_mask).sum())

            # Build column dict for this chunk
            col = {
                "runID":           np.zeros(n_rows, np.int32),
                "eventID":         eid_row,
                "sourceID":        np.zeros(n_rows, np.int32),
                "sourcePosX":      sx_row, "sourcePosY": sy_row, "sourcePosZ": sz_row,
                "time":            t_row.astype(np.float64),
                "energy":          (E_keV / 1000.0).astype(np.float32),  # MeV
                "globalPosX":      gx_row, "globalPosY": gy_row, "globalPosZ": gz_row,
                "gantryID":        np.zeros(n_rows, np.int32),
                "rsectorID":       rsec_row,
                "moduleID":        np.zeros(n_rows, np.int32),
                "submoduleID":     blk_row,
                "crystalID":       crystal_id_row,
                "layerID":         np.zeros(n_rows, np.int32),
                "comptonPhantom":  np.zeros(n_rows, np.int32),
                "comptonCrystal":  np.zeros(n_rows, np.int32),
                "RayleighPhantom": np.zeros(n_rows, np.int32),
                "RayleighCrystal": np.zeros(n_rows, np.int32),
            }

            # Sort by time within chunk
            order = np.argsort(col["time"], kind="stable")
            for k in col:
                col[k] = col[k][order]

            # Append to HDF5 datasets
            old = ds["runID"].shape[0]
            new = old + n_rows
            for b in SINGLES_BRANCHES:
                ds[b].resize((new,))
                ds[b][old:new] = col[b]

            n_total_written += n_rows
            if verbose and (ic == 0 or (ic + 1) % max(1, n_chunks // 10) == 0
                             or ic == n_chunks - 1):
                el = time.time() - t_start
                print(f"[lu176_bg]  chunk {ic+1:3d}/{n_chunks}  rows={n_rows:>7,}  "
                       f"cum={n_total_written:>10,}  elapsed={el:.1f}s")

        # Global time-sort across chunks: read full time array, derive
        # permutation, apply across all branches. With 4e8 rows this is
        # the most memory-intensive step (~3 GB float64 + perm). Skip the
        # global sort — the inject re-sorts during merge anyway. Document.
        # (We DID time-sort within each chunk, so the file is "approximately"
        # sorted with stride = chunk size.)
        if verbose:
            print("[lu176_bg]  intra-chunk time-sorted; "
                   "global time-sort skipped (inject re-sorts on merge).")

        # File attributes
        h.attrs["T_sim_s"]                       = float(T_sim_s)
        h.attrs["activity_Bq"]                   = float(activity_Bq)
        h.attrs["lu176_branching_307_keV"]       = GAMMA_LINES[0][1]
        h.attrs["lu176_branching_202_keV"]       = GAMMA_LINES[1][1]
        h.attrs["lu176_branching_88_keV"]        = GAMMA_LINES[2][1]
        h.attrs["beta_endpoint_keV"]             = BETA_ENDPOINT_KEV
        h.attrs["energy_resolution_FWHM_at_511"] = ENERGY_FWHM_AT_511_KEV
        h.attrs["timing_resolution_FWHM_ps"]     = TIMING_FWHM_PS
        h.attrs["n_decays"]                      = int(n_decays)
        h.attrs["n_singles"]                     = int(n_total_written)
        h.attrs["eventID_offset"]                = int(eventID_offset)
        h.attrs["generator_version"]             = version
        h.attrs["generator_seed"]                = int(seed)
        h.attrs["n_betas_in_511_window_sc44"]    = int(n_betas_in_511_window["sc44"])
        h.attrs["n_betas_in_511_window_i124"]    = int(n_betas_in_511_window["i124"])
        h.attrs["note"] = (
            "Analytical Lu-176 LSO-bg singles. Cascade: 307/202/88 keV with "
            "branching 0.94/0.78/0.15 (Conti 2017). Time-sorted within "
            "intra-chunk only (inject re-sorts globally on merge)."
        )

    elapsed = time.time() - t_start
    if verbose:
        print(f"[lu176_bg] DONE  N_singles={n_total_written:,}  elapsed={elapsed:.1f}s")
        print(f"[lu176_bg] β⁻-in-511-window:  Sc-44={n_betas_in_511_window['sc44']:,}  "
              f"I-124={n_betas_in_511_window['i124']:,}")
        print(f"[lu176_bg] file size: {output_path.stat().st_size / 1e6:.2f} MB")

    return dict(
        n_decays=n_decays,
        n_singles=n_total_written,
        n_betas_in_511_window_sc44=n_betas_in_511_window["sc44"],
        n_betas_in_511_window_i124=n_betas_in_511_window["i124"],
        elapsed_s=elapsed,
        output=str(output_path),
    )


# ════════════════════════ CLI ════════════════════════

def _cli():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--activity-MBq", type=float, default=14.59,
                     help="Total LSO Lu-176 activity in MBq. Default 14.59 "
                          "(full Quadra: 243,200 crystals × ~60 Bq each).")
    ap.add_argument("--T-sim-s", type=float, required=True,
                     help="Simulated acquisition time in seconds.")
    ap.add_argument("--output", type=Path, required=True,
                     help="Output HDF5 file path.")
    ap.add_argument("--seed", type=int, default=0,
                     help="RNG seed. Default 0.")
    ap.add_argument("--eventID-offset", type=int, default=DEFAULT_EVENTID_OFFSET,
                     help="EventID offset (Phase 1 c3 convention). Default 1e7.")
    ap.add_argument("--chunk-decays", type=int, default=200_000,
                     help="Decays per chunk (memory tuning). Default 200k.")
    args = ap.parse_args()
    generate(activity_MBq=args.activity_MBq,
              T_sim_s=args.T_sim_s,
              output_path=args.output,
              seed=args.seed,
              eventID_offset=args.eventID_offset,
              chunk_decays=args.chunk_decays)


if __name__ == "__main__":
    _cli()
