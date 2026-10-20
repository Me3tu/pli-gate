"""Stream a GATE Singles ROOT (one per-batch sim-merged.root) into a compact
Singles HDF5, applying the energy pre-filter and detecting shard boundaries
via eventID resets.

Schema (as specified on 2026-06-10):
    eventID          int32     shard-local; resets per shard inside the ROOT
    shard_id         uint8     shard index within this batch (0..15 typically)
    time_ns          float64   TRUE single-photon time (UNSMEARED), within-shard
    energy_keV       float32   post energy-smear, pre-time-smear
    gx, gy, gz       float32   crystal hit position (mm)
    sx, sy, sz       float32   source truth position (mm)
    compton_phantom  uint8     count of in-phantom Compton scatters per single
    in_511           uint8     1 if energy_keV in [481, 541]
    in_prompt        uint8     1 if energy_keV in [1100, 1200]

Energy pre-filter (wider buffer, decided 2026-06-10):
    keep if energy_keV in [460, 562] or [1080, 1220]

Output HDF5 layout:
    /singles/{eventID, shard_id, time_ns, energy_keV, gx, gy, gz, sx, sy, sz,
              compton_phantom, in_511, in_prompt}
    /attrs: n_singles_in_root, n_singles_kept, n_shards_in_batch, batch_id,
            E_keep_511_window, E_keep_prompt_window, gate_root_path,
            gate_root_size_bytes

Streams in chunks of 5M Singles to keep RAM bounded.
"""
from __future__ import annotations
import argparse
import json
import os
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import uproot

# DEFAULT energy windows: Sc-44 (brain-tumor / Sim 1).
# For I-124 (Mercolli, Steinberger), override via the CLI flags below.
# Wider energy pre-filter (analysis windows + ±20 keV buffer)
E_KEEP_511 = (460.0, 562.0)
E_KEEP_PROMPT = (1080.0, 1220.0)
# Tight analysis windows (for the in_511 / in_prompt flag)
E_511_TIGHT = (481.0, 541.0)
E_PROMPT_TIGHT = (1100.0, 1200.0)

NEEDED = ["eventID", "time", "energy",
          "globalPosX", "globalPosY", "globalPosZ",
          "sourcePosX", "sourcePosY", "sourcePosZ",
          "comptonPhantom"]

CHUNK_ENTRIES = 5_000_000


def detect_shards_from_eventid(eid_global_sorted_by_load_order: np.ndarray) -> np.ndarray:
    """Given the eventID array IN LOAD ORDER (the order Singles were written),
    detect shard boundaries via eventID resets. Returns shard_id (uint8) per
    Single.

    A shard boundary is where eventID decreases (or stays at zero) relative
    to the previous Single. We use diff < 0 to detect resets and cumsum the
    boundary mask.
    """
    n = len(eid_global_sorted_by_load_order)
    if n == 0:
        return np.zeros(0, dtype=np.uint8)
    resets = np.where(np.diff(eid_global_sorted_by_load_order) < 0)[0] + 1
    shard_id = np.zeros(n, dtype=np.uint8)
    if len(resets) > 0:
        shard_id[resets] = 1
        np.cumsum(shard_id, out=shard_id)
    return shard_id


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", required=True, type=Path,
                    help="per-batch sim-merged.root (hadd of N shards)")
    ap.add_argument("--out", required=True, type=Path,
                    help="output HDF5 path")
    ap.add_argument("--batch-id", required=True, type=int,
                    help="batch index (recorded in attrs)")
    ap.add_argument("--gzip-level", type=int, default=4,
                    help="HDF5 gzip compression level (1-9; 4 is a good speed/ratio)")
    ap.add_argument("--keep-511", nargs=2, type=float, default=[460.0, 562.0],
                    metavar=("LO", "HI"),
                    help="keV pre-filter for 511-band (default 460..562, Sc-44)")
    ap.add_argument("--keep-prompt", nargs=2, type=float, default=[1080.0, 1220.0],
                    metavar=("LO", "HI"),
                    help="keV pre-filter for prompt-band (default 1080..1220, Sc-44)")
    ap.add_argument("--tight-511", nargs=2, type=float, default=[481.0, 541.0],
                    metavar=("LO", "HI"),
                    help="keV tight window for in_511 flag (default 481..541)")
    ap.add_argument("--tight-prompt", nargs=2, type=float, default=[1100.0, 1200.0],
                    metavar=("LO", "HI"),
                    help="keV tight window for in_prompt flag (default 1100..1200)")
    args = ap.parse_args()

    global E_KEEP_511, E_KEEP_PROMPT, E_511_TIGHT, E_PROMPT_TIGHT
    E_KEEP_511     = tuple(args.keep_511)
    E_KEEP_PROMPT  = tuple(args.keep_prompt)
    E_511_TIGHT    = tuple(args.tight_511)
    E_PROMPT_TIGHT = tuple(args.tight_prompt)
    print(f"[singles-h5] energy windows: keep_511={E_KEEP_511} keep_prompt={E_KEEP_PROMPT}")
    print(f"[singles-h5]                 tight_511={E_511_TIGHT} tight_prompt={E_PROMPT_TIGHT}")

    t0 = time.time()
    print(f"[singles-h5] ROOT: {args.root}")
    print(f"[singles-h5] OUT:  {args.out}")
    if not args.root.exists():
        print(f"ERROR: ROOT not found", file=sys.stderr); return 2
    root_size = args.root.stat().st_size
    print(f"[singles-h5] ROOT size: {root_size/1e9:.2f} GB")

    args.out.parent.mkdir(parents=True, exist_ok=True)

    f_root = uproot.open(args.root)
    t = f_root["Singles"]
    n_total = t.num_entries
    print(f"[singles-h5] N Singles in ROOT: {n_total:,}")

    # Pass 1: stream all entries, build the keep mask per-chunk, count.
    # We need eventID across the WHOLE file to detect shard resets, so
    # eventID is loaded fully (8 bytes/entry × ~3e9 = 24 GB worst case, too big).
    # Instead, we load eventID-only first (int32 → 4 bytes; ~12 GB for 3e9
    # would still be big). Better: stream and detect resets on the fly,
    # writing chunks as we go.

    # Approach: single streaming pass.
    #   - keep a running shard_id counter
    #   - keep last eventID seen
    #   - per chunk: apply energy filter, compute shard_id, append to HDF5

    with h5py.File(args.out, "w") as h5:
        g = h5.create_group("singles")
        # Resizable datasets
        ds_eid = g.create_dataset("eventID", shape=(0,), maxshape=(None,),
                                   dtype="<i4", chunks=(1_000_000,),
                                   compression="gzip", compression_opts=args.gzip_level,
                                   shuffle=True)
        ds_sh = g.create_dataset("shard_id", shape=(0,), maxshape=(None,),
                                  dtype="<u1", chunks=(1_000_000,),
                                  compression="gzip", compression_opts=args.gzip_level)
        ds_t = g.create_dataset("time_ns", shape=(0,), maxshape=(None,),
                                 dtype="<f8", chunks=(1_000_000,),
                                 compression="gzip", compression_opts=args.gzip_level,
                                 shuffle=True)
        ds_e = g.create_dataset("energy_keV", shape=(0,), maxshape=(None,),
                                 dtype="<f4", chunks=(1_000_000,),
                                 compression="gzip", compression_opts=args.gzip_level,
                                 shuffle=True)
        for name in ("gx", "gy", "gz", "sx", "sy", "sz"):
            g.create_dataset(name, shape=(0,), maxshape=(None,),
                             dtype="<f4", chunks=(1_000_000,),
                             compression="gzip", compression_opts=args.gzip_level,
                             shuffle=True)
        ds_cp = g.create_dataset("compton_phantom", shape=(0,), maxshape=(None,),
                                  dtype="<u1", chunks=(1_000_000,),
                                  compression="gzip", compression_opts=args.gzip_level)
        ds_i5 = g.create_dataset("in_511", shape=(0,), maxshape=(None,),
                                  dtype="<u1", chunks=(1_000_000,),
                                  compression="gzip", compression_opts=args.gzip_level)
        ds_ip = g.create_dataset("in_prompt", shape=(0,), maxshape=(None,),
                                  dtype="<u1", chunks=(1_000_000,),
                                  compression="gzip", compression_opts=args.gzip_level)

        n_kept = 0
        running_shard = 0
        last_eid = -1
        n_streamed = 0
        for arr in t.iterate(NEEDED, library="np", step_size=CHUNK_ENTRIES):
            eid = arr["eventID"]
            n_streamed += len(eid)

            # Detect shard resets WITHIN this chunk + at chunk boundary
            if len(eid) == 0:
                continue
            chunk_shard = np.empty(len(eid), dtype=np.uint8)
            # First entry: shard reset if eid[0] < last_eid (or unconditionally
            # increment if last_eid < 0 and we are starting fresh)
            if last_eid < 0:
                chunk_shard[0] = running_shard   # 0 at the very start
            else:
                if eid[0] < last_eid:
                    running_shard += 1
                chunk_shard[0] = running_shard
            # Subsequent entries
            for i in range(1, len(eid)):
                if eid[i] < eid[i - 1]:
                    running_shard += 1
                chunk_shard[i] = running_shard
            last_eid = int(eid[-1])

            # Convert units
            e_keV = arr["energy"] * 1000.0
            t_ns = arr["time"] * 1e9   # GATE time is seconds

            # Energy pre-filter (wider buffer)
            keep = (((e_keV >= E_KEEP_511[0]) & (e_keV <= E_KEEP_511[1])) |
                    ((e_keV >= E_KEEP_PROMPT[0]) & (e_keV <= E_KEEP_PROMPT[1])))
            if not keep.any():
                continue

            eid_k = eid[keep].astype(np.int32, copy=False)
            sh_k = chunk_shard[keep]
            t_k = t_ns[keep].astype(np.float64, copy=False)
            e_k = e_keV[keep].astype(np.float32, copy=False)
            gx_k = arr["globalPosX"][keep].astype(np.float32, copy=False)
            gy_k = arr["globalPosY"][keep].astype(np.float32, copy=False)
            gz_k = arr["globalPosZ"][keep].astype(np.float32, copy=False)
            sx_k = arr["sourcePosX"][keep].astype(np.float32, copy=False)
            sy_k = arr["sourcePosY"][keep].astype(np.float32, copy=False)
            sz_k = arr["sourcePosZ"][keep].astype(np.float32, copy=False)
            cp_k = arr["comptonPhantom"][keep].astype(np.uint8, copy=False)
            i5_k = (((e_keV[keep] >= E_511_TIGHT[0]) & (e_keV[keep] <= E_511_TIGHT[1]))
                    .astype(np.uint8))
            ip_k = (((e_keV[keep] >= E_PROMPT_TIGHT[0]) & (e_keV[keep] <= E_PROMPT_TIGHT[1]))
                    .astype(np.uint8))

            # Append
            n_new = len(eid_k)
            for ds, arr_k in [
                (ds_eid, eid_k), (ds_sh, sh_k), (ds_t, t_k),
                (ds_e, e_k),
                (g["gx"], gx_k), (g["gy"], gy_k), (g["gz"], gz_k),
                (g["sx"], sx_k), (g["sy"], sy_k), (g["sz"], sz_k),
                (ds_cp, cp_k), (ds_i5, i5_k), (ds_ip, ip_k)
            ]:
                old = ds.shape[0]
                ds.resize((old + n_new,))
                ds[old:] = arr_k
            n_kept += n_new

            if n_streamed % 50_000_000 == 0 or n_streamed == n_total:
                elapsed = time.time() - t0
                print(f"  streamed {n_streamed:>11,} / {n_total:,}  "
                      f"({100*n_streamed/n_total:5.1f} %)  "
                      f"kept {n_kept:>10,}  elapsed {elapsed:.0f}s")

        n_shards_in_batch = running_shard + 1
        # attrs
        h5["singles"].attrs["n_singles_in_root"] = n_total
        h5["singles"].attrs["n_singles_kept"] = n_kept
        h5["singles"].attrs["n_shards_in_batch"] = n_shards_in_batch
        h5["singles"].attrs["batch_id"] = args.batch_id
        h5["singles"].attrs["E_keep_511_window"] = list(E_KEEP_511)
        h5["singles"].attrs["E_keep_prompt_window"] = list(E_KEEP_PROMPT)
        h5["singles"].attrs["E_511_tight_window"] = list(E_511_TIGHT)
        h5["singles"].attrs["E_prompt_tight_window"] = list(E_PROMPT_TIGHT)
        h5["singles"].attrs["gate_root_path"] = str(args.root)
        h5["singles"].attrs["gate_root_size_bytes"] = root_size
        h5["singles"].attrs["schema_version"] = 1
        h5["singles"].attrs["note"] = (
            "TRUE (unsmeared) single-photon time. Energy is post energy-smear. "
            "Wider energy pre-filter for scatter studies; tight analysis "
            "windows in_511 / in_prompt flags.")

    elapsed = time.time() - t0
    out_size = args.out.stat().st_size
    print(f"\n[singles-h5] DONE  n_total={n_total:,}  n_kept={n_kept:,}  "
          f"({100*n_kept/max(n_total,1):.2f} %)")
    print(f"[singles-h5]   shards detected: {n_shards_in_batch}")
    print(f"[singles-h5]   ROOT  {root_size/1e9:.2f} GB  →  HDF5 {out_size/1e9:.2f} GB  "
          f"(ratio {out_size/root_size:.2f})")
    print(f"[singles-h5]   wall: {elapsed:.0f}s = {elapsed/60:.1f} min")

    # Per-batch summary JSON for the run-aggregator + M7 random-rate analysis
    summary = dict(
        batch_id=args.batch_id,
        n_singles_in_root=int(n_total),
        n_singles_kept=int(n_kept),
        n_shards_in_batch=int(n_shards_in_batch),
        h5_path=str(args.out),
        h5_size_bytes=int(out_size),
        root_size_bytes=int(root_size),
        E_keep_511_window=list(E_KEEP_511),
        E_keep_prompt_window=list(E_KEEP_PROMPT),
        wall_clock_s=float(elapsed),
    )
    args.out.with_suffix(".json").write_text(json.dumps(summary, indent=2))
    print(f"[singles-h5]   wrote {args.out.with_suffix('.json')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
