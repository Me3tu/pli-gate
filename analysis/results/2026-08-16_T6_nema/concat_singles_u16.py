#!/usr/bin/env python3
"""Concatenate the 16 NEMA per-batch singles files, re-keying shard_id globally.

WHY THIS EXISTS RATHER THAN 2026-08-13_overnight/concat_singles.py
That script stores the global shard id as uint8 and refuses above 15 batches of
16 shards. NEMA has 16 x 16 = 256 blocks, which is one past uint8. It stores
uint16 instead. `find_shard_boundaries` in the extractor only does
`np.diff(shard_id) != 0`, so the wider dtype is transparent to it.

WHY THE RE-KEY MATTERS AT ALL
eventID restarts at 0 in every shard. The extractor groups by eventID but only
WITHIN a contiguous run of equal shard_id, so without a global re-key the 16
batches' shard ids would collide and singles from different batches would be
fused into coincidences that never happened. Verified on batch_00 before any of
this ran: the merged ROOT has exactly 15 eventID resets in file order (16
blocks of ~991k), and within one block the mean eventID group size is 1.45 with
a p99 time spread of 4.9 ns -- physical. Across the whole merged file, ignoring
blocks, the mean group size is 16.28, i.e. all 16 shards fused.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

import h5py
import numpy as np

FIELDS = ["eventID", "time_ns", "energy_keV", "gx", "gy", "gz",
          "sx", "sy", "sz", "compton_phantom", "in_511", "in_prompt"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-dir", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    parts = sorted(args.in_dir.glob("singles_batch_*.h5"))
    if not parts:
        raise SystemExit(f"no singles_batch_*.h5 in {args.in_dir}")
    print(f"[concat] {len(parts)} per-batch files")

    counts, nshards = [], []
    for p in parts:
        with h5py.File(p, "r") as f:
            counts.append(f["singles"]["eventID"].shape[0])
            nshards.append(int(f["singles"]["shard_id"][:].max()) + 1)
    total = sum(counts)
    if len(set(nshards)) != 1:
        print(f"[concat] WARNING: shards per batch vary: {nshards}")
    print(f"[concat] {total:,} singles, shards/batch {nshards}")

    global_blocks = sum(nshards)
    if global_blocks > np.iinfo(np.uint16).max:
        raise SystemExit(f"{global_blocks} blocks exceeds uint16")
    print(f"[concat] {global_blocks} global shard blocks -> uint16")

    dtypes = {"eventID": np.int32, "time_ns": np.float64, "energy_keV": np.float32,
              "gx": np.float32, "gy": np.float32, "gz": np.float32,
              "sx": np.float32, "sy": np.float32, "sz": np.float32,
              "compton_phantom": np.uint8, "in_511": np.uint8, "in_prompt": np.uint8}

    with h5py.File(args.out, "w") as fo:
        g = fo.create_group("singles")
        for k in FIELDS:
            g.create_dataset(k, shape=(total,), dtype=dtypes[k],
                             compression="gzip", compression_opts=4)
        g.create_dataset("shard_id", shape=(total,), dtype=np.uint16,
                         compression="gzip", compression_opts=4)
        cur, offset = 0, 0
        for bi, (p, n) in enumerate(zip(parts, counts)):
            with h5py.File(p, "r") as fi:
                gi = fi["singles"]
                for k in FIELDS:
                    g[k][cur:cur + n] = gi[k][:]
                loc = gi["shard_id"][:].astype(np.uint16)
                g["shard_id"][cur:cur + n] = loc + offset
                offset += nshards[bi]
            cur += n
            print(f"[concat]   {p.name}: +{n:,} -> {cur:,}", flush=True)
        # Carry the provenance attrs forward. Without these the merged file --
        # the one every analysis actually opens -- has no record of which
        # energy bands it contains, and a question about a band the converter
        # discarded returns data instead of failing. That has produced one
        # bogus result already; see energy_window_guard.py.
        CARRY = ("E_keep_511_window", "E_keep_prompt_window",
                 "E_511_tight_window", "E_prompt_tight_window")
        with h5py.File(parts[0], "r") as f0:
            first = {k: f0["singles"].attrs[k] for k in CARRY
                     if k in f0["singles"].attrs}
        for p_ in parts[1:]:
            with h5py.File(p_, "r") as fi:
                for k, v in first.items():
                    other = fi["singles"].attrs.get(k)
                    if other is None or not np.array_equal(np.asarray(v),
                                                           np.asarray(other)):
                        raise SystemExit(
                            f"[concat] ABORT: {k} differs between {parts[0].name} "
                            f"({v}) and {p_.name} ({other}). Refusing to merge "
                            f"batches with inconsistent energy selections.")
        for k, v in first.items():
            g.attrs[k] = v
        if not first:
            raise SystemExit(
                "[concat] ABORT: no E_keep_* attrs on the input batches, so the "
                "merged file's energy content would be unrecorded.")

        g.attrs["n_batches"] = len(parts)
        g.attrs["n_singles"] = total
        g.attrs["shard_id_dtype"] = "uint16"
        g.attrs["note"] = ("T6 NEMA concat; shard_id re-keyed globally as "
                           "sum(shards in preceding batches) + local")

    # verify the blocks came out contiguous, which is what the extractor assumes
    with h5py.File(args.out, "r") as f:
        sh = f["singles"]["shard_id"][:]
    nblocks = int((np.diff(sh) != 0).sum()) + 1
    print(f"[concat] contiguous blocks: {nblocks} (expect {global_blocks})")
    json.dump(dict(n_singles=int(total), n_batches=len(parts),
                   shards_per_batch=nshards, global_blocks=int(global_blocks),
                   contiguous_blocks=int(nblocks), ok=bool(nblocks == global_blocks)),
              open(args.out.with_suffix(".concat.json"), "w"), indent=1)
    if nblocks != global_blocks:
        raise SystemExit("block count mismatch -- shard re-key is wrong, refusing")
    print(f"[concat] wrote {args.out} ({args.out.stat().st_size/1e9:.2f} GB)")


if __name__ == "__main__":
    main()
