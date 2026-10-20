"""
Chunked, fixed-RAM rewrite of extract_triples.py.

Stream-reads ROOT files in CHUNK_ENTRIES-sized blocks, maintains a
single-event tail buffer across chunk boundaries, and appends to
resizable HDF5 datasets so peak RSS stays bounded regardless of input
size.

Output H5 schema is bit-identical to extract_triples.py:
  triples_detector<tag>.h5  group "detector": 15 datasets (float64)
                                              + same attrs as original
  triples_truth<tag>.h5     group "truth":    11 datasets
                                              (eventID, source_file_index int64,
                                               others float64)

CLI:
  --input-format {root, h5}      default 'root'.
  --root-dir DIR                  required for root mode.
  --input PATH                    required for h5 mode (single .h5 file
                                  with a /singles/ group, e.g., the
                                  output of lso_background_inject).
  --out-dir DIR
  --tag TAG

Phase 3p-1 added the h5 input mode for the lso_background_inject
round-trip (Phase 3o-4). H5 input is pre-sorted by (eventID, time)
in-memory before chunking. Warning logged if input exceeds 50 M rows
(pre-sort peak ~ N × 100 B). For ROOT, behaviour is bit-identical to
the pre-Phase-3p-1 implementation (regression PASS verified against
NEMA 88M reference 2026-05-04).

Design notes are in
reconstruction/results/2026-04-26_0137_extract_triples_chunked/extract_design.md
and Phase 3p-1 design at
reconstruction/results/2026-05-04_2220_phase3p1_h5_extractor/findings.md.

KNOWN BEHAVIOUR (Phase 3p-1 R3 finding, NOT FIXED, preserved exactly):
  T_COINC_NS = 5.0 is compared against `time` field values that are in
  SECONDS (GATE Singles convention). The threshold is therefore
  effectively 5 SECONDS, not 5 ns. The threshold is currently a no-op
  because eventID grouping already constrains all candidates to a
  single decay (within ~ns naturally). Methodology description:
  "coincidence selection via eventID grouping. The T_COINC_NS variable
  is currently effectively a no-op due to time-units mismatch and is
  retained pending a future fix." Fix candidate: scale T_COINC_NS by
  1e-9 in the comparisons. Out of scope for Phase 3p-1 (would change
  all historical results).
"""

from __future__ import annotations

import argparse
import os
import resource
import time
from pathlib import Path

import h5py
import numpy as np
import ROOT

ROOT.gErrorIgnoreLevel = ROOT.kError

E_511_LO, E_511_HI = 481.0, 541.0
E_PROMPT_LO, E_PROMPT_HI = 1100.0, 1200.0
T_COINC_NS = 5.0

DETECTOR_FIELDS = [
    "a1_time", "a1_x", "a1_y", "a1_z", "a1_energy",
    "a2_time", "a2_x", "a2_y", "a2_z", "a2_energy",
    "p_time",  "p_x",  "p_y",  "p_z",  "p_energy",
]
TRUTH_FIELDS = [
    "eventID",
    "source_x", "source_y", "source_z",
    "a1_compton_phantom", "a2_compton_phantom", "p_compton_phantom",
    "a1_rayleigh_phantom", "a2_rayleigh_phantom", "p_rayleigh_phantom",
    "source_file_index",
]
INT64_TRUTH = {"eventID", "source_file_index"}

NEEDED = [
    "eventID",
    "time", "energy",
    "globalPosX", "globalPosY", "globalPosZ",
    "sourcePosX", "sourcePosY", "sourcePosZ",
    "comptonPhantom", "RayleighPhantom",
]

CHUNK_ENTRIES = 5_000_000
H5_FLUSH_BATCH = 50_000


def rss_gb() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 / 1024


def iter_root_files(root_dir: Path) -> list[Path]:
    files = sorted(root_dir.glob("sim_*.root"))
    if not files:
        raise FileNotFoundError(f"No sim_*.root in {root_dir}")
    return files


def _read_h5_singles_sorted(path: Path, log) -> tuple[dict, bool]:
    """Phase 3p-1: load /singles/ branches from an HDF5 file, then sort
    by (eventID, time) so the existing chunk-pending logic can apply
    unchanged.

    Inject's merged H5 is time-sorted (`lso_background_inject` line 276),
    which interleaves rows from different decays. The pending-buffer
    chunk logic assumes consecutive rows with the same eventID belong
    to the same decay, true after sorting by (eventID, time).

    Returns (arrs, has_source_flag). All NEEDED branches are stored as
    float64 to match the ROOT-mode storage convention.
    """
    log(f"  [h5 mode] reading {path}")
    with h5py.File(path, "r") as f:
        grp = f["singles"] if "singles" in f else f
        arrs = {}
        for b in NEEDED:
            if b not in grp:
                raise KeyError(f"Required branch '{b}' missing in {path}")
            arrs[b] = grp[b][:]
        has_source_flag = "source_flag" in grp
        if has_source_flag:
            arrs["source_flag"] = grp["source_flag"][:]

    n = arrs["time"].size
    if n > 50_000_000:
        log(f"  [h5 mode] WARNING: {n:,} rows exceeds 50 M soft limit; "
            f"in-memory pre-sort peak RAM ≈ {n * 100 / 1e9:.1f} GB.")
    log(f"  [h5 mode] loaded {n:,} rows; "
        f"source_flag {'present' if has_source_flag else 'absent'}")

    # Promote NEEDED to float64 to match ROOT-mode storage (read_chunk:107).
    for b in NEEDED:
        arrs[b] = arrs[b].astype(np.float64, copy=False)

    # Sort by (eventID, time). lexsort: last key is primary.
    order = np.lexsort((arrs["time"], arrs["eventID"]))
    for k in list(arrs.keys()):
        arrs[k] = arrs[k][order]
    log(f"  [h5 mode] pre-sorted by (eventID, time)")
    return arrs, has_source_flag


def open_singles(path: Path):
    f = ROOT.TFile.Open(str(path), "READ")
    if f.IsZombie():
        raise IOError(f"Cannot open {path}")
    t = f.Get("Singles")
    if not t or int(t.GetEntries()) == 0:
        f.Close()
        return None, None
    t.SetBranchStatus("*", 0)
    for b in NEEDED:
        t.SetBranchStatus(b, 1)
    bufs = {}
    for b in NEEDED:
        if b == "eventID" or b in ("comptonPhantom", "RayleighPhantom"):
            bufs[b] = np.zeros(1, dtype=np.int32)
        elif b == "time":
            bufs[b] = np.zeros(1, dtype=np.float64)
        else:
            bufs[b] = np.zeros(1, dtype=np.float32)
        t.SetBranchAddress(b, bufs[b])
    return f, (t, bufs)


def read_chunk(t, bufs, start: int, end: int) -> dict:
    """Read entries [start, end) into a dict of float64 arrays.

    All branches are stored as float64 to match the original's storage
    layout (extract_triples.py line 93 promotes everything to float64).
    """
    n = end - start
    out = {b: np.empty(n, dtype=np.float64) for b in NEEDED}
    for i in range(n):
        t.GetEntry(start + i)
        for b in NEEDED:
            out[b][i] = bufs[b][0]
    return out


class TripleAccumulator:
    """In-RAM batch of triples. Flushes to resizable H5 datasets."""

    def __init__(self, det_grp: h5py.Group, tru_grp: h5py.Group,
                 batch: int = H5_FLUSH_BATCH):
        self.batch_size = batch
        self.batch_det = {k: [] for k in DETECTOR_FIELDS}
        self.batch_tru = {k: [] for k in TRUTH_FIELDS}
        self.dset_det = {}
        self.dset_tru = {}
        self.det_grp = det_grp
        self.tru_grp = tru_grp
        self._n = 0
        self._n_total = 0

    def _ensure_dsets(self):
        if self.dset_det:
            return
        for k in DETECTOR_FIELDS:
            self.dset_det[k] = self.det_grp.create_dataset(
                k, shape=(0,), maxshape=(None,), dtype=np.float64,
                chunks=(10_000,), compression="gzip", compression_opts=4,
            )
        for k in TRUTH_FIELDS:
            dt = np.int64 if k in INT64_TRUTH else np.float64
            self.dset_tru[k] = self.tru_grp.create_dataset(
                k, shape=(0,), maxshape=(None,), dtype=dt,
                chunks=(10_000,), compression="gzip", compression_opts=4,
            )

    def add(self, det_row: dict, tru_row: dict):
        for k in DETECTOR_FIELDS:
            self.batch_det[k].append(det_row[k])
        for k in TRUTH_FIELDS:
            self.batch_tru[k].append(tru_row[k])
        self._n += 1
        self._n_total += 1
        if self._n >= self.batch_size:
            self.flush()

    def flush(self):
        if self._n == 0:
            return
        self._ensure_dsets()
        for k, dset in self.dset_det.items():
            arr = np.asarray(self.batch_det[k], dtype=np.float64)
            old = dset.shape[0]
            dset.resize((old + arr.size,))
            dset[old:old + arr.size] = arr
            self.batch_det[k] = []
        for k, dset in self.dset_tru.items():
            dt = np.int64 if k in INT64_TRUTH else np.float64
            arr = np.asarray(self.batch_tru[k], dtype=dt)
            old = dset.shape[0]
            dset.resize((old + arr.size,))
            dset[old:old + arr.size] = arr
            self.batch_tru[k] = []
        self._n = 0


def find_triples_in_event(s: int, e: int, arrs: dict, file_idx: int,
                          acc: TripleAccumulator) -> int:
    """Apply the original's per-event triple finder to arrs[s:e].

    Mirrors lines 146-193 of extract_triples.py exactly.
    Returns number of triples produced.
    """
    if e - s < 3:
        return 0
    energy_keV = arrs["energy"][s:e] * 1000.0
    in_511 = (energy_keV >= E_511_LO) & (energy_keV <= E_511_HI)
    in_prompt = (energy_keV >= E_PROMPT_LO) & (energy_keV <= E_PROMPT_HI)
    n511 = int(in_511.sum())
    nprompt = int(in_prompt.sum())
    if n511 < 2 or nprompt < 1:
        return 0

    time_e = arrs["time"][s:e]
    x_e = arrs["globalPosX"][s:e]
    y_e = arrs["globalPosY"][s:e]
    z_e = arrs["globalPosZ"][s:e]
    sx_e = arrs["sourcePosX"][s:e]
    sy_e = arrs["sourcePosY"][s:e]
    sz_e = arrs["sourcePosZ"][s:e]
    cp_e = arrs["comptonPhantom"][s:e]
    rp_e = arrs["RayleighPhantom"][s:e]
    eid = arrs["eventID"][s]   # constant within an event

    idx_511 = np.where(in_511)[0]
    idx_prompt = np.where(in_prompt)[0]

    n_added = 0
    for i in range(n511):
        for j in range(i + 1, n511):
            a1 = int(idx_511[i]); a2 = int(idx_511[j])
            if abs(time_e[a1] - time_e[a2]) > T_COINC_NS:
                continue
            for pi in range(nprompt):
                p = int(idx_prompt[pi])
                t_mid = 0.5 * (time_e[a1] + time_e[a2])
                if abs(time_e[p] - t_mid) > T_COINC_NS * 2:
                    continue
                det = {
                    "a1_time": time_e[a1], "a2_time": time_e[a2], "p_time": time_e[p],
                    "a1_x": x_e[a1], "a2_x": x_e[a2], "p_x": x_e[p],
                    "a1_y": y_e[a1], "a2_y": y_e[a2], "p_y": y_e[p],
                    "a1_z": z_e[a1], "a2_z": z_e[a2], "p_z": z_e[p],
                    "a1_energy": energy_keV[a1],
                    "a2_energy": energy_keV[a2],
                    "p_energy":  energy_keV[p],
                }
                tru = {
                    "eventID": int(eid),
                    "source_x": sx_e[a1], "source_y": sy_e[a1], "source_z": sz_e[a1],
                    "a1_compton_phantom": cp_e[a1],
                    "a2_compton_phantom": cp_e[a2],
                    "p_compton_phantom":  cp_e[p],
                    "a1_rayleigh_phantom": rp_e[a1],
                    "a2_rayleigh_phantom": rp_e[a2],
                    "p_rayleigh_phantom":  rp_e[p],
                    "source_file_index": int(file_idx),
                }
                acc.add(det, tru)
                n_added += 1
    return n_added


def find_triples_in_pending(pending_arrs: dict, file_idx: int,
                            acc: TripleAccumulator) -> int:
    """Same as find_triples_in_event but on a standalone single-event dict."""
    n = pending_arrs["eventID"].size
    if n < 3:
        return 0
    return find_triples_in_event(0, n, pending_arrs, file_idx, acc)


def process_one_file(path: Path, file_idx: int, acc: TripleAccumulator,
                     log) -> tuple[int, int]:
    f, ts = open_singles(path)
    if ts is None:
        return 0, 0
    t, bufs = ts
    n = int(t.GetEntries())
    n_chunks = (n + CHUNK_ENTRIES - 1) // CHUNK_ENTRIES
    log(f"  [{path.name}] {n:,} entries in {n_chunks} chunks of {CHUNK_ENTRIES:,}")

    pending_eid = None       # python int or None
    pending_arrs = None      # dict of np.float64 arrays
    triples_this_file = 0
    max_pending = 0

    t0 = time.time()
    for chunk_idx, chunk_start in enumerate(range(0, n, CHUNK_ENTRIES)):
        chunk_end = min(chunk_start + CHUNK_ENTRIES, n)
        ck_t0 = time.time()
        arrs = read_chunk(t, bufs, chunk_start, chunk_end)
        eids = arrs["eventID"]
        n_chunk = len(eids)

        chunk_offset = 0

        # 1. absorb prefix matching pending
        if pending_eid is not None:
            if n_chunk > 0 and int(eids[0]) == pending_eid:
                mismatch = np.where(eids != float(pending_eid))[0]
                prefix_end = int(mismatch[0]) if mismatch.size else n_chunk
                for k in pending_arrs:
                    pending_arrs[k] = np.concatenate(
                        [pending_arrs[k], arrs[k][:prefix_end]])
                chunk_offset = prefix_end
            if chunk_offset == n_chunk:
                # entire chunk absorbed, still inside same event
                continue
            # event closes at chunk_offset
            triples_this_file += find_triples_in_pending(
                pending_arrs, file_idx, acc)
            pending_eid, pending_arrs = None, None

        # 2. find event boundaries in arrs[chunk_offset:]
        if chunk_offset == n_chunk:
            continue
        sub_eids = eids[chunk_offset:]
        diffs = np.diff(sub_eids)
        breaks = np.where(diffs != 0)[0] + 1
        starts = np.concatenate(([0], breaks))
        ends = np.concatenate((breaks, [len(sub_eids)]))

        # 3. process all events except the last
        for k_seg in range(len(starts) - 1):
            s = int(starts[k_seg]) + chunk_offset
            e = int(ends[k_seg]) + chunk_offset
            triples_this_file += find_triples_in_event(
                s, e, arrs, file_idx, acc)

        # 4. stash last segment as new pending
        s_last = int(starts[-1]) + chunk_offset
        pending_eid = int(eids[s_last])
        pending_arrs = {k: arrs[k][s_last:].copy() for k in arrs}
        max_pending = max(max_pending, pending_arrs["eventID"].size)

        log(f"    chunk {chunk_idx+1}/{n_chunks} "
            f"[{chunk_start:>10,}–{chunk_end:>10,}]  "
            f"{time.time()-ck_t0:5.1f}s  "
            f"triples={triples_this_file:>8,}  "
            f"pending_size={pending_arrs['eventID'].size:>5}  "
            f"rss={rss_gb():.2f}GB")

    if pending_eid is not None:
        triples_this_file += find_triples_in_pending(
            pending_arrs, file_idx, acc)

    f.Close()
    log(f"  [{path.name}] done. {triples_this_file:,} triples, "
        f"max_pending_size={max_pending}, wall={time.time()-t0:.1f}s")
    return n, triples_this_file


def process_one_h5_file(path: Path, file_idx: int, acc: TripleAccumulator,
                         log) -> tuple[int, int]:
    """Phase 3p-1: HDF5 equivalent of process_one_file.

    Loads, pre-sorts, then chunks through the (eventID, time) sorted
    arrays applying the existing pending-buffer logic. Triple-formation
    semantics are bit-identical to ROOT mode (same find_triples_in_event
    core).
    """
    arrs_full, _has_sflag = _read_h5_singles_sorted(path, log)
    n = int(arrs_full["eventID"].size)
    if n == 0:
        return 0, 0
    n_chunks = (n + CHUNK_ENTRIES - 1) // CHUNK_ENTRIES
    log(f"  [{path.name}] {n:,} entries (sorted) in {n_chunks} chunks "
        f"of {CHUNK_ENTRIES:,}")

    pending_eid = None
    pending_arrs = None
    triples_this_file = 0
    max_pending = 0
    t0 = time.time()

    for chunk_idx, chunk_start in enumerate(range(0, n, CHUNK_ENTRIES)):
        chunk_end = min(chunk_start + CHUNK_ENTRIES, n)
        ck_t0 = time.time()
        # View slice (no copy needed)
        arrs = {k: arrs_full[k][chunk_start:chunk_end] for k in arrs_full
                if k in NEEDED}
        eids = arrs["eventID"]
        n_chunk = len(eids)
        chunk_offset = 0

        # 1. absorb prefix matching pending
        if pending_eid is not None:
            if n_chunk > 0 and int(eids[0]) == pending_eid:
                mismatch = np.where(eids != float(pending_eid))[0]
                prefix_end = int(mismatch[0]) if mismatch.size else n_chunk
                for k in pending_arrs:
                    pending_arrs[k] = np.concatenate(
                        [pending_arrs[k], arrs[k][:prefix_end]])
                chunk_offset = prefix_end
            if chunk_offset == n_chunk:
                continue
            triples_this_file += find_triples_in_pending(
                pending_arrs, file_idx, acc)
            pending_eid, pending_arrs = None, None

        # 2. find event boundaries
        if chunk_offset == n_chunk:
            continue
        sub_eids = eids[chunk_offset:]
        diffs = np.diff(sub_eids)
        breaks = np.where(diffs != 0)[0] + 1
        starts = np.concatenate(([0], breaks))
        ends = np.concatenate((breaks, [len(sub_eids)]))

        # 3. process all events except the last
        for k_seg in range(len(starts) - 1):
            s = int(starts[k_seg]) + chunk_offset
            e = int(ends[k_seg]) + chunk_offset
            triples_this_file += find_triples_in_event(
                s, e, arrs, file_idx, acc)

        # 4. stash last segment as pending
        s_last = int(starts[-1]) + chunk_offset
        pending_eid = int(eids[s_last])
        pending_arrs = {k: arrs[k][s_last:].copy() for k in arrs}
        max_pending = max(max_pending, pending_arrs["eventID"].size)

        log(f"    chunk {chunk_idx+1}/{n_chunks} "
            f"[{chunk_start:>10,}–{chunk_end:>10,}]  "
            f"{time.time()-ck_t0:5.1f}s  "
            f"triples={triples_this_file:>8,}  "
            f"pending_size={pending_arrs['eventID'].size:>5}  "
            f"rss={rss_gb():.2f}GB")

    if pending_eid is not None:
        triples_this_file += find_triples_in_pending(
            pending_arrs, file_idx, acc)

    log(f"  [{path.name}] done. {triples_this_file:,} triples, "
        f"max_pending_size={max_pending}, wall={time.time()-t0:.1f}s")
    return n, triples_this_file


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-format", choices=["root", "h5"], default="root",
                    help="Input format. Default 'root' (backward-compat).")
    ap.add_argument("--root-dir", default=None,
                    help="ROOT mode: directory containing sim_*.root files.")
    ap.add_argument("--input", type=Path, default=None,
                    help="H5 mode: single .h5 file with /singles/ group "
                         "(e.g., output of lso_background_inject).")
    ap.add_argument("--out-dir", default="data")
    ap.add_argument("--tag", default="_chunked")
    args = ap.parse_args()

    if args.input_format == "root":
        if args.root_dir is None:
            ap.error("--root-dir is required when --input-format=root")
        root_dir = Path(args.root_dir).resolve()
        files = iter_root_files(root_dir)
        process_fn = process_one_file
        print(f"[extract_triples_chunked] {len(files)} ROOT file(s) "
              f"in {root_dir}")
    elif args.input_format == "h5":
        if args.input is None:
            ap.error("--input is required when --input-format=h5")
        h5_path = Path(args.input).resolve()
        if not h5_path.exists():
            ap.error(f"--input file not found: {h5_path}")
        files = [h5_path]
        process_fn = process_one_h5_file
        print(f"[extract_triples_chunked] 1 H5 file: {h5_path}")
    else:
        ap.error(f"unknown --input-format: {args.input_format}")

    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    for p in files:
        print(f"  {p.name}  ({p.stat().st_size/1e9:.2f} GB)")

    det_path = out_dir / f"triples_detector{args.tag}.h5"
    tru_path = out_dir / f"triples_truth{args.tag}.h5"
    print(f"[extract_triples_chunked] writing to:")
    print(f"  {det_path}")
    print(f"  {tru_path}")

    t_start = time.time()

    with h5py.File(det_path, "w") as fd, h5py.File(tru_path, "w") as ft:
        det_grp = fd.create_group("detector")
        tru_grp = ft.create_group("truth")
        acc = TripleAccumulator(det_grp, tru_grp)

        total_singles = 0
        total_triples = 0
        for idx, p in enumerate(files):
            print(f"\n[{idx+1}/{len(files)}] {p.name}")
            n, nt = process_fn(p, file_idx=idx, acc=acc, log=print)
            total_singles += n
            total_triples += nt

        acc.flush()

        fd.attrs["energy_window_511_keV"] = (E_511_LO, E_511_HI)
        fd.attrs["energy_window_prompt_keV"] = (E_PROMPT_LO, E_PROMPT_HI)
        fd.attrs["time_coinc_window_ns"] = T_COINC_NS
        fd.attrs["n_triples"] = total_triples
        fd.attrs["note"] = "Detector fields only. Usable in reconstruction."
        ft.attrs["n_triples"] = total_triples
        ft.attrs["note"] = (
            "TRUTH-ONLY fields. Must NEVER be used in reconstruction. "
            "Validation code may read these with caution."
        )

    wall = time.time() - t_start
    print(f"\n[extract_triples_chunked] {total_singles:,} singles → "
          f"{total_triples:,} triples in {wall:.1f}s "
          f"({wall/60:.1f} min)")
    print(f"[extract_triples_chunked] wrote {det_path}  "
          f"({det_path.stat().st_size/1e6:.1f} MB)")
    print(f"[extract_triples_chunked] wrote {tru_path}  "
          f"({tru_path.stat().st_size/1e6:.1f} MB)")
    print(f"[extract_triples_chunked] peak RSS = {rss_gb():.2f} GB")


if __name__ == "__main__":
    main()
