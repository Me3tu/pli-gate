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

CLI matches extract_triples.py: --root-dir, --out-dir, --tag.

Design notes are in
reconstruction/results/2026-04-26_0137_extract_triples_chunked/extract_design.md


[Phase 2b Step 5b variant: i124 windows substituted, rest is identical to extract_triples_chunked.py.]
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

E_511_LO, E_511_HI = 460.0, 545.0   # Mercolli I-124 LEW (was 481/541 for Sc-44)
E_PROMPT_LO, E_PROMPT_HI = 568.0, 639.0   # Mercolli I-124 HEW (was 1100/1200 for Sc-44)
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root-dir", required=True)
    ap.add_argument("--out-dir", default="data")
    ap.add_argument("--tag", default="_chunked")
    args = ap.parse_args()

    root_dir = Path(args.root_dir).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    files = iter_root_files(root_dir)
    print(f"[extract_triples_chunked] {len(files)} ROOT file(s) in {root_dir}")
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
            n, nt = process_one_file(p, file_idx=idx, acc=acc, log=print)
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
