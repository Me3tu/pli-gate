# SPDX-License-Identifier: LGPL-3.0-only
# Copyright (C) 2026 Dhritman Dhar and German Cancer Research Center (DKFZ)
"""Lu-176 LSO background injection — post-processing merge into phantom singles.

Reads a phantom GATE singles ROOT (or HDF5 dump) and a Lu-176 standalone
singles ROOT, time-aligns them by treating the Lu-176 stream as a continuous
background, and writes a merged HDF5 singles file with a per-event provenance
flag so downstream consumers (extract_triples_chunked.py) can separate
"phantom-only" from "phantom+bg" triple counts.

This module is purely additive: it does not modify any production-current
script (extract_triples_chunked.py, simple_recon.py, ...). It produces a
NEW HDF5 singles file with the same schema as a phantom GATE Singles tree
plus one extra column `source_flag` (uint8: 0 = phantom signal, 1 = LSO bg).

Schema of output H5:
    /singles/eventID                                   int64     ← Phase 3p-3
    /singles/runID, sourceID                           int32
    /singles/sourcePosX/Y/Z                            float32
    /singles/time                                      float64
    /singles/energy                                    float32
    /singles/globalPosX/Y/Z                            float32
    /singles/gantryID, rsectorID, moduleID,
              submoduleID, crystalID, layerID          int32
    /singles/comptonPhantom, comptonCrystal,
              RayleighPhantom, RayleighCrystal         int32
    /singles/source_flag                               uint8     ← NEW
        0 = phantom origin
        1 = LSO Lu-176 bg origin
    attrs:
      n_phantom, n_bg, n_total
      phantom_T_sim_s, bg_T_sim_s, scale_bg_to_phantom_T
      bg_activity_Bq, phantom_activity_Bq

Time alignment:
    The Lu-176 standalone sim ran for T_bg simulated seconds at A_bg activity,
    producing N_bg singles. The phantom sim ran for T_phantom simulated
    seconds. To inject a physically-correct fraction of bg singles into the
    phantom timeline, we sample
        N_bg_inject = N_bg × (T_phantom / T_bg)
    from the bg singles uniformly with replacement (Poisson-sampling with
    unit-mean rate; if N_bg >> N_bg_inject we sample without replacement
    instead). Times of selected bg singles are uniformly redistributed over
    [0, T_phantom].

The module exposes:
    inject(phantom_path, bg_path, output_path, ...) — main entry
    _PhantomReader / _BgReader                       — small helpers
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

# Branches we expect on a GATE Singles tree (matches what
# extract_triples_chunked.py reads).
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

# dtype map for HDF5 output. Matches GATE Singles tree types.
DTYPES = {
    "runID":           np.int32,
    "eventID":         np.int64,
    "sourceID":        np.int32,
    "sourcePosX":      np.float32,
    "sourcePosY":      np.float32,
    "sourcePosZ":      np.float32,
    "time":            np.float64,
    "energy":          np.float32,
    "globalPosX":      np.float32,
    "globalPosY":      np.float32,
    "globalPosZ":      np.float32,
    "gantryID":        np.int32,
    "rsectorID":       np.int32,
    "moduleID":        np.int32,
    "submoduleID":     np.int32,
    "crystalID":       np.int32,
    "layerID":         np.int32,
    "comptonPhantom":  np.int32,
    "comptonCrystal":  np.int32,
    "RayleighPhantom": np.int32,
    "RayleighCrystal": np.int32,
    "source_flag":     np.uint8,
}


@dataclass
class _Stream:
    """In-memory view of a singles stream (one of: phantom, bg)."""
    arrays: dict[str, np.ndarray]
    T_sim_s: float       # simulated acquisition time in seconds
    activity_Bq: float   # source activity (or None if not known)
    label: str           # "phantom" or "bg" (for logs)

    def __len__(self) -> int:
        return len(self.arrays["time"])

    def select(self, idx: np.ndarray) -> dict[str, np.ndarray]:
        """Return arrays subsetted by `idx`."""
        return {k: v[idx] for k, v in self.arrays.items()}


def _load_h5_singles(path: Path,
                     T_sim_s: float | None = None,
                     activity_Bq: float | None = None,
                     label: str = "h5") -> _Stream:
    """Load a singles HDF5 file (test/synthetic format)."""
    import h5py
    arrs = {}
    with h5py.File(path, "r") as h:
        grp = h["singles"] if "singles" in h else h
        for b in SINGLES_BRANCHES:
            if b in grp:
                arrs[b] = grp[b][:]
            elif b in h:
                arrs[b] = h[b][:]
            else:
                # missing branches default to zeros (synthetic test inputs
                # may not populate every detector field)
                n = len(arrs.get("time", []))
                if n == 0 and "time" in grp:
                    n = len(grp["time"])
                arrs[b] = np.zeros(n, dtype=DTYPES[b])
        if T_sim_s is None:
            T_sim_s = float(h.attrs.get("T_sim_s", 1.0))
        if activity_Bq is None:
            activity_Bq = float(h.attrs.get("activity_Bq", 0.0))
    return _Stream(arrays=arrs, T_sim_s=T_sim_s,
                   activity_Bq=activity_Bq, label=label)


def _load_root_singles(path: Path,
                       T_sim_s: float,
                       activity_Bq: float,
                       label: str = "root") -> _Stream:
    """Load a GATE Singles tree from a ROOT file via uproot."""
    import uproot
    f = uproot.open(path)
    s = f["Singles"]
    # All branches loaded eagerly. For very large bg streams (>100 M
    # singles), prefer chunked iteration — but the LSO bg target is
    # ≤250 M and ~12.5 GB raw which fits on this machine; a future
    # optimisation can chunk if needed.
    arrs = {b: s[b].array(library="np") for b in SINGLES_BRANCHES if b in s.keys()}
    # Defensive default-fill for missing branches (older GATE configs may
    # omit some).
    n = len(arrs["time"])
    for b in SINGLES_BRANCHES:
        if b not in arrs:
            arrs[b] = np.zeros(n, dtype=DTYPES[b])
    return _Stream(arrays=arrs, T_sim_s=T_sim_s,
                   activity_Bq=activity_Bq, label=label)


def _redistribute_bg_times(n_inject: int, T_phantom_s: float,
                           rng: np.random.Generator) -> np.ndarray:
    """Sample n_inject new bg-event times uniform on [0, T_phantom]."""
    return rng.uniform(0.0, T_phantom_s, size=n_inject)


def _select_bg_indices(bg: _Stream, T_phantom_s: float,
                       rng: np.random.Generator) -> np.ndarray:
    """Sample a Poisson-mean number of bg events for the phantom acquisition.

    N_target = N_bg × (T_phantom / T_bg)        (deterministic mean)
    N_inject ~ Poisson(N_target)                (Poisson realisation)
    Indices sampled without replacement when N_inject ≤ N_bg, else
    with-replacement to satisfy the Poisson rate at higher injection
    fractions.
    """
    if bg.T_sim_s <= 0:
        raise ValueError(f"bg.T_sim_s must be > 0, got {bg.T_sim_s}")
    N_bg = len(bg)
    mean = N_bg * (T_phantom_s / bg.T_sim_s)
    N_inject = int(rng.poisson(mean))
    replace = N_inject > N_bg
    return rng.choice(N_bg, size=N_inject, replace=replace)


def inject(phantom_path: str | Path,
           bg_path: str | Path,
           output_path: str | Path,
           phantom_T_sim_s: float | None = None,
           bg_T_sim_s: float | None = None,
           bg_activity_Bq: float | None = None,
           phantom_activity_Bq: float | None = None,
           seed: int = 0,
           verbose: bool = True) -> dict:
    """Merge a phantom singles file with an Lu-176 bg singles file.

    Parameters
    ----------
    phantom_path : path to phantom singles, .h5 (synthetic test) or
        .root (GATE Singles tree). If .root, must supply phantom_T_sim_s.
    bg_path : path to LSO bg singles. Same format options.
    output_path : .h5 file to write, with /singles group + source_flag column.
    phantom_T_sim_s, bg_T_sim_s : simulated acquisition times (seconds).
        For .h5 inputs these are read from file attrs if not given.
    bg_activity_Bq, phantom_activity_Bq : source activities, recorded as
        attrs on the output but not used in the merge math (only the
        time ratio matters once activities are baked into the singles
        counts).
    seed : RNG seed for reproducible bg sampling.
    verbose : print progress.

    Returns
    -------
    summary dict with n_phantom, n_bg, n_bg_injected, n_total, scale.
    """
    import h5py
    phantom_path = Path(phantom_path)
    bg_path      = Path(bg_path)
    output_path  = Path(output_path)
    rng = np.random.default_rng(seed)

    # Load phantom
    if phantom_path.suffix in (".h5", ".hdf5"):
        phantom = _load_h5_singles(phantom_path, phantom_T_sim_s,
                                    phantom_activity_Bq, "phantom")
    elif phantom_path.suffix == ".root":
        if phantom_T_sim_s is None:
            raise ValueError("phantom_T_sim_s required for .root input")
        phantom = _load_root_singles(phantom_path, phantom_T_sim_s,
                                      phantom_activity_Bq or 0.0, "phantom")
    else:
        raise ValueError(f"unsupported phantom format: {phantom_path.suffix}")

    # Load bg
    if bg_path.suffix in (".h5", ".hdf5"):
        bg = _load_h5_singles(bg_path, bg_T_sim_s, bg_activity_Bq, "bg")
    elif bg_path.suffix == ".root":
        if bg_T_sim_s is None:
            raise ValueError("bg_T_sim_s required for .root input")
        bg = _load_root_singles(bg_path, bg_T_sim_s,
                                 bg_activity_Bq or 0.0, "bg")
    else:
        raise ValueError(f"unsupported bg format: {bg_path.suffix}")

    if verbose:
        print(f"[lso_inject] phantom: {len(phantom):,} singles "
              f"over T = {phantom.T_sim_s} s")
        print(f"[lso_inject] bg:      {len(bg):,} singles over T = {bg.T_sim_s} s")

    # Sample bg subset
    bg_idx = _select_bg_indices(bg, phantom.T_sim_s, rng)
    bg_sub = bg.select(bg_idx)
    n_bg_inj = len(bg_idx)
    bg_sub["time"] = _redistribute_bg_times(n_bg_inj, phantom.T_sim_s, rng)
    if verbose:
        scale = phantom.T_sim_s / bg.T_sim_s
        print(f"[lso_inject] injecting {n_bg_inj:,} bg singles "
              f"(scale = T_p/T_bg = {scale:.4g})")

    # Build merged stream by concatenation, then sort by time
    merged: dict[str, np.ndarray] = {}
    for b in SINGLES_BRANCHES:
        merged[b] = np.concatenate([
            phantom.arrays[b].astype(DTYPES[b], copy=False),
            bg_sub[b].astype(DTYPES[b], copy=False),
        ])
    flag = np.concatenate([
        np.zeros(len(phantom), dtype=np.uint8),
        np.ones(n_bg_inj, dtype=np.uint8),
    ])
    merged["source_flag"] = flag

    # Sort by time
    order = np.argsort(merged["time"], kind="stable")
    for k in list(merged):
        merged[k] = merged[k][order]

    # Write
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(output_path, "w") as h:
        grp = h.create_group("singles")
        for b in SINGLES_BRANCHES + ("source_flag",):
            grp.create_dataset(b, data=merged[b],
                                compression="gzip", compression_opts=4,
                                chunks=True)
        h.attrs["n_phantom"]            = len(phantom)
        h.attrs["n_bg"]                 = n_bg_inj
        h.attrs["n_total"]              = len(phantom) + n_bg_inj
        h.attrs["phantom_T_sim_s"]      = phantom.T_sim_s
        h.attrs["bg_T_sim_s"]           = bg.T_sim_s
        h.attrs["scale_bg_to_phantom_T"]= phantom.T_sim_s / bg.T_sim_s
        h.attrs["bg_activity_Bq"]       = bg.activity_Bq
        h.attrs["phantom_activity_Bq"]  = phantom.activity_Bq
        h.attrs["seed"]                 = seed
        h.attrs["note"] = ("Merged phantom + LSO Lu-176 bg singles. "
                            "source_flag: 0=phantom, 1=bg.")

    if verbose:
        print(f"[lso_inject] wrote {output_path} "
              f"({output_path.stat().st_size/1e6:.2f} MB)")

    return dict(
        n_phantom=len(phantom),
        n_bg=n_bg_inj,
        n_total=len(phantom) + n_bg_inj,
        scale=phantom.T_sim_s / bg.T_sim_s,
        output=str(output_path),
    )


# ----------------------------------------------------------------------------
# Unit tests / synthetic dry run helpers
# ----------------------------------------------------------------------------

def _make_synthetic_singles(path: Path, n: int, T_sim_s: float,
                             activity_Bq: float,
                             energy_mean: float = 0.511, energy_sigma: float = 0.05,
                             seed: int = 0) -> None:
    """Write a tiny synthetic singles H5 (used by dry-run tests).

    Energy is in MeV (matching GATE convention).
    Time is uniform on [0, T_sim_s].
    eventID is monotonic increasing (one decay per integer).
    Position fields are zero — we only care about time/energy/eventID for
    dry-run validation of the inject algebra and downstream extraction.
    """
    import h5py
    rng = np.random.default_rng(seed)
    arrs = {b: np.zeros(n, dtype=DTYPES[b]) for b in SINGLES_BRANCHES}
    arrs["time"]    = np.sort(rng.uniform(0.0, T_sim_s, size=n)).astype(np.float64)
    arrs["energy"]  = rng.normal(energy_mean, energy_sigma, size=n).astype(np.float32)
    arrs["eventID"] = np.arange(n, dtype=np.int64)
    arrs["runID"]   = np.zeros(n, dtype=np.int32)
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as h:
        grp = h.create_group("singles")
        for b in SINGLES_BRANCHES:
            grp.create_dataset(b, data=arrs[b])
        h.attrs["T_sim_s"]      = T_sim_s
        h.attrs["activity_Bq"]  = activity_Bq


def _self_test() -> None:
    """Run a small in-process unit test on synthetic data."""
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        # Phantom: 10 k singles in 1 s at 7 MBq, energy peaked near 511 keV
        _make_synthetic_singles(td / "phantom.h5",
                                 n=10_000, T_sim_s=1.0, activity_Bq=7e6,
                                 energy_mean=0.511, energy_sigma=0.030, seed=1)
        # Bg: 14 590 singles in 1 s at 14.59 MBq, energy peaked near 0.3 MeV
        _make_synthetic_singles(td / "bg.h5",
                                 n=14_590, T_sim_s=1.0, activity_Bq=14.59e6,
                                 energy_mean=0.300, energy_sigma=0.080, seed=2)
        out = td / "merged.h5"
        summary = inject(td / "phantom.h5", td / "bg.h5", out,
                         seed=42, verbose=True)

        # Validate
        import h5py
        with h5py.File(out, "r") as h:
            n_total = h["singles/time"].shape[0]
            n_phantom_ds = int((h["singles/source_flag"][:] == 0).sum())
            n_bg_ds = int((h["singles/source_flag"][:] == 1).sum())
            t = h["singles/time"][:]
            assert np.all(np.diff(t) >= -1e-12), "output is not time-sorted"
            assert n_total == n_phantom_ds + n_bg_ds
            assert n_phantom_ds == 10_000
            assert abs(n_bg_ds - 14_590) < 5 * np.sqrt(14_590)  # Poisson 5σ
            assert h.attrs["n_phantom"] == 10_000
            print(f"[self_test] PASS  n_phantom={n_phantom_ds:,}  "
                  f"n_bg={n_bg_ds:,}  total={n_total:,}")

        return summary


if __name__ == "__main__":
    _self_test()
