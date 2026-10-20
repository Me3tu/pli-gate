# SPDX-License-Identifier: LGPL-3.0-only
# Copyright (C) 2026 Dhritman Dhar and German Cancer Research Center (DKFZ)
"""Phase 3q production driver — VIPman 1.4B per-batch sequential pipeline.

Pattern (validated by Phase 3p-2 smoketest at 50M scale):
  for batch_id in range(32):
    GATE sim (16 threads × 4 GB cap)
      → ROOT → HDF5 (int64 decay_id re-key)
      → bg gen (default offset 10^12)
      → inject merge
      → extract triples (h5 mode)
      → cleanup intermediates (keep only triples + logs + JSON)
  concat 32 per-batch triples
  correct_travel_time
  simple_recon (4 mm × 10 iter, --mumap, --sensitivity 4mm-fullfid)
  per-tissue analysis + figures

Resume: detects existing triples_detector_batch_NN.h5 with valid n_triples
attribute and skips that batch. Per-batch is the failure-recovery granule.

Per-batch sanity range: [500_000, 1_500_000] triples (smoketest gave
22.6 triples/k-primary; per-batch 43.75M × 22.6 ≈ 990K).

Phase 3q seed scheme:
  GATE shard:  1_000_000 + batch_id * 1000 + thread_id
  bg gen:      12_345 + batch_id
  inject:      42 + batch_id

NO production code changes. Calls existing utilities:
  - scripts/run_gate_parallel_musire.sh  (GATE sim runner)
  - reconstruction/utils/lu176_bg_generator.py
  - reconstruction/utils/lso_background_inject.py  (in-process)
  - reconstruction/extract_triples_chunked.py  (subprocess)
  - reconstruction/correct_travel_time.py  (in-process; module-level DET path
    overridden via monkey-patch — see _correct_travel_time_phase3q below)
  - reconstruction/simple_recon.py single  (subprocess)

Phase 3q sha (must be unchanged):
  lu176_bg_generator.py    03f9953d   (Phase 3p-3)
  lso_background_inject.py a1458848   (Phase 3p-3)
  extract_triples_chunked.py 634ab75a (Phase 3p-1)
  simple_recon.py            df11628e (Phase 3n PART 1)
  correct_travel_time.py     269d39c8
"""
from __future__ import annotations
import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import h5py
import numpy as np

# v1.1.0: paths relative to the release and from the environment (simulation/README.md)
REPO = Path(__file__).resolve().parents[1]          # simulation/
ANALYSIS = Path(os.environ.get("PLI_ANALYSIS", REPO.parent / "analysis"))

# ── Production parameters ────────────────────────────────────────────────────
N_BATCHES = 32
THREADS_PER_BATCH = 16
T_TOTAL_S = 26.77            # = 1.4e9 / 52.29e6 (atlas activity)
T_PER_BATCH_S = T_TOTAL_S / N_BATCHES  # 0.836562 s
ATLAS_ACTIVITY_BQ = 52.29e6
BG_ACTIVITY_MBQ = 14.59       # full Quadra LSO

MACRO = REPO / "macros/vipman_slab.mac"
RUNNER = REPO / "runner/run_gate_parallel.sh"
SENSITIVITY = Path(os.environ.get("PLI_SENSITIVITY_NPY", ANALYSIS / "results/sensitivity_4mm_fullfidelity_norm.npy"))
MUMAP = Path(os.environ.get("PLI_VIPMAN_DIR", "PLI_VIPMAN_DIR-not-set")) / "vipman-atlas-lungs-sim-mumap-harmonized.mhd"

PER_BATCH_TRIPLE_MIN = 500_000
PER_BATCH_TRIPLE_MAX = 1_500_000
EVENTID_OFFSET_BG = 1_000_000_000_000

VENV_PY = Path(os.environ.get("PLI_PYTHON", sys.executable))
PY_PATH = os.environ.get("PLI_EXTRA_PYTHONPATH", "")


# ── Per-batch primitives ─────────────────────────────────────────────────────

def gate_sim(batch_id: int, batch_dir: Path, log) -> Path:
    """Run GATE sim for one batch (16 threads × 4 GB cap, sequential)."""
    out = batch_dir / "sim"
    sim_root = out / "sim-merged.root"
    if sim_root.exists() and sim_root.stat().st_size > 1_000_000_000:
        log(f"  [batch {batch_id}] sim ROOT already present ({sim_root.stat().st_size/1e9:.2f} GB), skipping GATE")
        return sim_root
    out.mkdir(parents=True, exist_ok=True)
    cmd = [str(RUNNER),
            "--macro", str(MACRO),
            "--total-time-sec", f"{T_PER_BATCH_S}",
            "--threads", str(THREADS_PER_BATCH),
            "--bounded-max-gb", "4",
            "--batch-id", str(batch_id),
            "--label", "vipman_1400M_poisson",
            "--output-dir", str(out)]
    log(f"  [batch {batch_id}] GATE: " + " ".join(cmd))
    t0 = time.time()
    subprocess.check_call(cmd, cwd=REPO)
    log(f"  [batch {batch_id}] GATE wall: {(time.time() - t0)/60:.1f} min, "
         f"output: {sim_root.stat().st_size/1e9:.2f} GB")
    return sim_root


def root_to_h5(sim_root: Path, phantom_h5: Path, t_sim_s: float, log):
    """STEP 5 equivalent. ROOT → HDF5 with int64 decay_id re-key (Phase 3p-3 pattern)."""
    if phantom_h5.exists() and phantom_h5.stat().st_size > 100_000_000:
        log(f"  phantom h5 already present, skipping STEP 5")
        return
    BRANCHES = (
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
        "runID": np.int32, "eventID": np.int64, "sourceID": np.int32,  # Phase 3p-3 int64
        "sourcePosX": np.float32, "sourcePosY": np.float32, "sourcePosZ": np.float32,
        "time": np.float64, "energy": np.float32,
        "globalPosX": np.float32, "globalPosY": np.float32, "globalPosZ": np.float32,
        "gantryID": np.int32, "rsectorID": np.int32, "moduleID": np.int32,
        "submoduleID": np.int32, "crystalID": np.int32, "layerID": np.int32,
        "comptonPhantom": np.int32, "comptonCrystal": np.int32,
        "RayleighPhantom": np.int32, "RayleighCrystal": np.int32,
    }
    import uproot
    t0 = time.time()
    f = uproot.open(sim_root)
    n_total = int(f["Singles"].num_entries)
    arrs = {b: np.empty(n_total, dtype=DTYPES[b]) for b in BRANCHES}
    pos = 0
    for chunk in f["Singles"].iterate(filter_name=list(BRANCHES), library="np",
                                            step_size=2_000_000):
        n = len(chunk["time"])
        for b in BRANCHES:
            arrs[b][pos:pos + n] = chunk[b].astype(DTYPES[b], copy=False) if b in chunk else 0
        pos += n
    # int64 re-key (Phase 3p-3)
    eid_orig = arrs["eventID"]
    times = arrs["time"]
    is_new = np.concatenate([[True], np.diff(eid_orig) != 0]) | \
                 np.concatenate([[False], np.diff(times) < -0.5])
    arrs["eventID"] = (np.cumsum(is_new.astype(np.int64)) - 1).astype(np.int64)
    n_decays = int(arrs["eventID"].max() + 1)
    with h5py.File(phantom_h5, "w") as h:
        g = h.create_group("singles")
        for b in BRANCHES:
            g.create_dataset(b, data=arrs[b], compression="gzip",
                                compression_opts=4, chunks=True)
        h.attrs["T_sim_s"] = t_sim_s
        h.attrs["activity_Bq"] = ATLAS_ACTIVITY_BQ
        h.attrs["n_singles"] = n_total
        h.attrs["n_decays"] = n_decays
        h.attrs["source_root"] = str(sim_root)
    log(f"  STEP 5 wall: {time.time()-t0:.1f}s, n_singles={n_total:,}, n_decays={n_decays:,}, "
         f"output: {phantom_h5.stat().st_size/1e6:.1f} MB")


def bg_gen(bg_h5: Path, t_sim_s: float, batch_id: int, log):
    """STEP 6 equivalent. Lu-176 bg gen with default offset 10^12 (Phase 3p-3)."""
    if bg_h5.exists() and bg_h5.stat().st_size > 100_000_000:
        log(f"  bg h5 already present, skipping STEP 6")
        return
    cmd = [str(VENV_PY),
            str(REPO / "lu176/lu176_bg_generator.py"),
            "--activity-MBq", str(BG_ACTIVITY_MBQ),
            "--T-sim-s", f"{t_sim_s}",
            "--output", str(bg_h5),
            "--seed", str(12_345 + batch_id)]
    env = os.environ.copy()
    env["PYTHONPATH"] = PY_PATH + ":" + env.get("PYTHONPATH", "")
    t0 = time.time()
    out = subprocess.check_output(cmd, env=env, cwd=REPO,
                                       stderr=subprocess.STDOUT, timeout=300)
    last = out.decode(errors="ignore").splitlines()[-3:]
    log("  " + "\n  ".join(last))
    log(f"  STEP 6 wall: {time.time()-t0:.1f}s, output: {bg_h5.stat().st_size/1e6:.1f} MB")


def inject_step(phantom_h5: Path, bg_h5: Path, merged_h5: Path,
                  batch_id: int, log):
    """STEP 7. In-process inject (peak ~25 GB)."""
    if merged_h5.exists() and merged_h5.stat().st_size > 1_000_000_000:
        log(f"  merged h5 already present, skipping STEP 7")
        return
    sys.path.insert(0, str(REPO / "lu176"))
    from lso_background_inject import inject  # noqa
    t0 = time.time()
    inject(phantom_h5, bg_h5, merged_h5, seed=42 + batch_id, verbose=False)
    log(f"  STEP 7 wall: {time.time()-t0:.1f}s, output: {merged_h5.stat().st_size/1e6:.1f} MB")


def extract_triples(merged_h5: Path, batch_dir: Path, batch_id: int, log) -> tuple[Path, Path]:
    """STEP 8. h5-mode triple extraction via subprocess."""
    det = batch_dir / "triples_detector_batch.h5"
    tru = batch_dir / "triples_truth_batch.h5"
    if det.exists() and tru.exists() and det.stat().st_size > 10_000_000:
        with h5py.File(det, "r") as f:
            if int(f.attrs.get("n_triples", 0)) > 0:
                log(f"  triples already present, skipping STEP 8")
                return det, tru
    cmd = [str(VENV_PY),
            str(ANALYSIS / "reconstruction/extract_triples_chunked.py"),
            "--input-format", "h5",
            "--input", str(merged_h5),
            "--out-dir", str(batch_dir),
            "--tag", "_batch"]
    env = os.environ.copy()
    env["PYTHONPATH"] = PY_PATH + ":" + env.get("PYTHONPATH", "")
    t0 = time.time()
    subprocess.check_call(cmd, env=env, cwd=REPO)
    log(f"  STEP 8 wall: {time.time()-t0:.1f}s, "
         f"detector: {det.stat().st_size/1e6:.1f} MB, "
         f"truth: {tru.stat().st_size/1e6:.1f} MB")
    return det, tru


def validate_batch(det: Path, batch_id: int, log) -> int:
    """Per-batch sanity check on triple count."""
    with h5py.File(det, "r") as f:
        n_triples = int(f.attrs["n_triples"])
    if n_triples < PER_BATCH_TRIPLE_MIN or n_triples > PER_BATCH_TRIPLE_MAX:
        raise ValueError(f"batch {batch_id} triple count {n_triples:,} "
                            f"outside sanity range [{PER_BATCH_TRIPLE_MIN:,}, "
                            f"{PER_BATCH_TRIPLE_MAX:,}]")
    log(f"  batch {batch_id} triple count: {n_triples:,} ✓ "
         f"(range [{PER_BATCH_TRIPLE_MIN/1e6:.1f}, {PER_BATCH_TRIPLE_MAX/1e6:.1f}]M)")
    return n_triples


def cleanup_batch(batch_dir: Path, keep_merged: bool, log):
    """Delete sim ROOTs + phantom h5 + bg h5 + merged h5 + per-batch shards.
    Keep: triples_*.h5, *.log, *.json."""
    candidates = [
        batch_dir / "phantom.h5",
        batch_dir / "bg.h5",
        batch_dir / "sim",  # contains all GATE outputs
    ]
    if not keep_merged:
        candidates.append(batch_dir / "merged.h5")
    freed = 0
    for c in candidates:
        if c.is_dir():
            sz = sum(p.stat().st_size for p in c.rglob("*") if p.is_file())
            shutil.rmtree(c)
            freed += sz
        elif c.exists():
            sz = c.stat().st_size
            c.unlink()
            freed += sz
    log(f"  cleanup freed: {freed/1e9:.2f} GB")


def write_batch_summary(batch_dir: Path, batch_id: int,
                          n_triples: int, walls: dict):
    """Per-batch JSON summary."""
    summary = dict(
        batch_id=batch_id,
        n_triples=n_triples,
        T_sim_s=T_PER_BATCH_S,
        threads=THREADS_PER_BATCH,
        seed_gate_base=1_000_000 + batch_id * 1_000,
        seed_bg=12_345 + batch_id,
        seed_inject=42 + batch_id,
        walls=walls,
        timestamp=time.strftime("%Y-%m-%d %H:%M:%S"),
    )
    (batch_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, default=str))


# ── Per-batch orchestrator ───────────────────────────────────────────────────

def run_batch(batch_id: int, sess: Path, keep_merged: bool, force: bool):
    batch_dir = sess / f"batch_{batch_id:02d}"
    batch_dir.mkdir(parents=True, exist_ok=True)
    log_path = batch_dir / "batch.log"

    def log(msg):
        ts = time.strftime("%H:%M:%S")
        line = f"[{ts}] {msg}"
        print(line)
        with log_path.open("a") as f:
            f.write(line + "\n")

    # Resume detection
    det = batch_dir / "triples_detector_batch.h5"
    if not force and det.exists() and det.stat().st_size > 10_000_000:
        try:
            with h5py.File(det, "r") as f:
                n = int(f.attrs.get("n_triples", 0))
            if PER_BATCH_TRIPLE_MIN <= n <= PER_BATCH_TRIPLE_MAX:
                log(f"batch {batch_id}: triples already valid (n={n:,}), SKIP")
                return n
        except Exception:
            pass

    log(f"=== batch {batch_id}/{N_BATCHES} START ===")
    walls = {}

    t = time.time(); sim_root = gate_sim(batch_id, batch_dir, log); walls["gate_min"] = (time.time()-t)/60
    phantom_h5 = batch_dir / "phantom.h5"
    t = time.time(); root_to_h5(sim_root, phantom_h5, T_PER_BATCH_S, log); walls["root2h5_s"] = time.time()-t
    bg_h5 = batch_dir / "bg.h5"
    t = time.time(); bg_gen(bg_h5, T_PER_BATCH_S, batch_id, log); walls["bg_s"] = time.time()-t
    merged_h5 = batch_dir / "merged.h5"
    t = time.time(); inject_step(phantom_h5, bg_h5, merged_h5, batch_id, log); walls["inject_s"] = time.time()-t
    t = time.time(); det_h5, tru_h5 = extract_triples(merged_h5, batch_dir, batch_id, log); walls["extract_s"] = time.time()-t
    n_triples = validate_batch(det_h5, batch_id, log)
    cleanup_batch(batch_dir, keep_merged, log)
    write_batch_summary(batch_dir, batch_id, n_triples, walls)
    log(f"=== batch {batch_id} DONE n_triples={n_triples:,} ===")
    return n_triples


# ── Concat across batches + downstream ───────────────────────────────────────

def concat_triples(sess: Path, batch_ids: list[int], log):
    """Concat 32 per-batch detector + truth H5s into final triples files."""
    out_det = sess / "triples_detector_phase3q.h5"
    out_tru = sess / "triples_truth_phase3q.h5"
    if out_det.exists() and out_tru.exists():
        with h5py.File(out_det, "r") as f:
            log(f"  concat already present (n_triples={int(f.attrs.get('n_triples', 0)):,}), skipping")
        return out_det, out_tru
    t0 = time.time()
    det_files = [sess / f"batch_{i:02d}/triples_detector_batch.h5" for i in batch_ids]
    tru_files = [sess / f"batch_{i:02d}/triples_truth_batch.h5" for i in batch_ids]
    DETECTOR_KEYS = sorted(h5py.File(det_files[0], "r")["detector"].keys())
    TRUTH_KEYS = sorted(h5py.File(tru_files[0], "r")["truth"].keys())
    INT64_TRUTH = {"eventID", "source_file_index"}

    n_total = 0
    for f in det_files:
        with h5py.File(f, "r") as h:
            n_total += int(h.attrs["n_triples"])

    with h5py.File(out_det, "w") as fd, h5py.File(out_tru, "w") as ft:
        det_grp = fd.create_group("detector")
        tru_grp = ft.create_group("truth")
        for k in DETECTOR_KEYS:
            det_grp.create_dataset(k, shape=(n_total,), dtype=np.float64,
                                       compression="gzip", compression_opts=4)
        for k in TRUTH_KEYS:
            dt = np.int64 if k in INT64_TRUTH else np.float64
            tru_grp.create_dataset(k, shape=(n_total,), dtype=dt,
                                       compression="gzip", compression_opts=4)
        cursor = 0
        for fd_in, ft_in, batch_id in zip(det_files, tru_files, batch_ids):
            with h5py.File(fd_in, "r") as h:
                n = int(h.attrs["n_triples"])
                for k in DETECTOR_KEYS:
                    det_grp[k][cursor:cursor+n] = h["detector"][k][:]
            with h5py.File(ft_in, "r") as h:
                for k in TRUTH_KEYS:
                    if k == "source_file_index":
                        # remap to batch_id (was 0 within each batch's extract)
                        tru_grp[k][cursor:cursor+n] = batch_id
                    else:
                        tru_grp[k][cursor:cursor+n] = h["truth"][k][:]
            cursor += n
        fd.attrs["n_triples"] = n_total
        fd.attrs["n_batches"] = len(batch_ids)
        fd.attrs["energy_window_511_keV"] = (481.0, 541.0)
        fd.attrs["energy_window_prompt_keV"] = (1100.0, 1200.0)
        fd.attrs["time_coinc_window_ns"] = 5.0
        fd.attrs["note"] = "Phase 3q concat across 32 batches"
        ft.attrs["n_triples"] = n_total
        ft.attrs["n_batches"] = len(batch_ids)
        ft.attrs["note"] = "Phase 3q truth; source_file_index = batch_id"
    log(f"  concat wall: {time.time()-t0:.1f}s, n_triples={n_total:,}, "
         f"detector: {out_det.stat().st_size/1e9:.2f} GB, "
         f"truth: {out_tru.stat().st_size/1e9:.2f} GB")
    return out_det, out_tru


def correct_travel_time(det: Path, log):
    """In-process travel-time correction on the concat detector h5."""
    sys.path.insert(0, str(ANALYSIS / "reconstruction"))
    from correct_travel_time import correct  # noqa
    t0 = time.time()
    tau, P_annih = correct(det)
    n_nan = int(np.isnan(tau).sum() + np.isinf(tau).sum() +
                  np.isnan(P_annih).any(axis=1).sum() +
                  np.isinf(P_annih).any(axis=1).sum())
    if n_nan:
        raise RuntimeError(f"correct_travel_time produced {n_nan} NaN/Inf")
    with h5py.File(det, "a") as f:
        if "derived" in f:
            del f["derived"]
        d = f.create_group("derived")
        d.create_dataset("tau_corrected_ns", data=tau,
                            compression="gzip", compression_opts=4)
        d.create_dataset("P_annih_tof_mm", data=P_annih,
                            compression="gzip", compression_opts=4)
    log(f"  correct_travel_time wall: {time.time()-t0:.1f}s, "
         f"τ median {np.median(tau):.3f} ns, "
         f"τ in [0.5,10]: {((tau>0.5)&(tau<10)).mean():.3f}")


def simple_recon(det: Path, sess: Path, log):
    """4 mm × 10 iter, --mumap, --sensitivity 4mm-fullfid (Run B baseline recipe)."""
    out = sess / "recon_phase3q"
    cmd = [str(VENV_PY),
            str(ANALYSIS / "reconstruction/simple_recon.py"), "single",
            "--label", "phase3q",
            "--voxel", "4.0",
            "--iters", "10",
            "--tau-range", "-1.0", "15.0",
            "--outdir", str(out),
            "--sensitivity", str(SENSITIVITY),
            "--mumap", str(MUMAP),
            "--triples-override", str(det)]
    env = os.environ.copy()
    env["PARALLELPROJ_C_LIB"] = os.environ.get("PARALLELPROJ_C_LIB", "/usr/local/lib/libparallelproj_c.so")
    env["PARALLELPROJ_CUDA_LIB"] = os.environ.get("PARALLELPROJ_CUDA_LIB", "/usr/local/lib/libparallelproj_cuda.so")
    env["PYTHONPATH"] = PY_PATH + ":" + env.get("PYTHONPATH", "")
    t0 = time.time()
    subprocess.check_call(cmd, env=env, cwd=REPO)
    log(f"  simple_recon wall: {time.time()-t0:.1f}s, output: {out}/")


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--session", required=True, help="output session dir")
    ap.add_argument("--start-batch", type=int, default=0)
    ap.add_argument("--end-batch", type=int, default=N_BATCHES)
    ap.add_argument("--force-batch", type=int, action="append", default=[],
                     help="re-run this batch even if triples exist (repeatable)")
    ap.add_argument("--keep-merged", action="store_true",
                     help="DEBUG: keep per-batch merged H5 (109 GB across 32)")
    ap.add_argument("--skip-downstream", action="store_true",
                     help="run batches only, skip concat + recon")
    ap.add_argument("--dry-run", action="store_true",
                     help="validate inputs and exit; no GATE")
    args = ap.parse_args()

    sess = Path(args.session).resolve()
    sess.mkdir(parents=True, exist_ok=True)

    print(f"Phase 3q production driver")
    print(f"  session:   {sess}")
    print(f"  batches:   [{args.start_batch}, {args.end_batch})")
    print(f"  T_total:   {T_TOTAL_S:.3f} s = 1.4B at {ATLAS_ACTIVITY_BQ/1e6} MBq")
    print(f"  T/batch:   {T_PER_BATCH_S:.4f} s ({T_PER_BATCH_S/THREADS_PER_BATCH:.4f} s/shard)")

    if args.dry_run:
        for p in (MACRO, RUNNER, SENSITIVITY, MUMAP):
            print(f"  {'✓' if p.exists() else '✗'} {p}")
        print("dry-run complete")
        return

    t_session_start = time.time()
    n_triples_per_batch = []
    for bid in range(args.start_batch, args.end_batch):
        force = bid in args.force_batch
        n = run_batch(bid, sess, args.keep_merged, force)
        n_triples_per_batch.append(n)

    print(f"\nAll batches done. Total triples: {sum(n_triples_per_batch):,}")
    print(f"Total batch-loop wall: {(time.time()-t_session_start)/3600:.2f} h")

    if args.skip_downstream:
        return

    print("\n=== concat + downstream ===")
    batch_ids = list(range(args.start_batch, args.end_batch))
    det, tru = concat_triples(sess, batch_ids, print)
    correct_travel_time(det, print)
    simple_recon(det, sess, print)

    print(f"\nPhase 3q complete. Total session wall: {(time.time()-t_session_start)/3600:.2f} h")


if __name__ == "__main__":
    main()
