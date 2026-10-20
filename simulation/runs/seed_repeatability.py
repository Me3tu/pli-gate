#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-only
# Copyright (C) 2026 Dhritman Dhar and German Cancer Research Center (DKFZ)
"""S2: seed repeatability of the 12-tube closure (Sim-1 half).

3 additional seeds (bases 2e6/3e6/4e6, production base 1e6), reduced
from the planned 4 because 4 x (32 min GATE + convert/extract)
exceeds the ~2 h budget. Stated here rather than silently
truncated. Sim-2 (I-124) is not re-run: same recipe, costlier isotope, and
the Sim-1 half carries the closure claim's statistics (8 of 12 tubes).

Per seed: 4 sequential batches x 16 threads x 4.425 s (= 250 M decays
at 14.13 MBq, the production configuration), per-batch extraction
(eventID-grouped production extractor) with immediate ROOT deletion
(disk: 20 GB free vs 24 GB full-run ROOT), field-wise triples concat,
locked recipe (filter_triples_tof_pann --r-cut 75, V3 cumulant
tau_min = 1.0), per-tube tau_hat vs configured.
"""
from __future__ import annotations
import json
import shutil
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np

# v1.1.0: paths relative to the release and from the environment (simulation/README.md)
import os
ROOT = Path(__file__).resolve().parents[1]          # simulation/
ANALYSIS = Path(os.environ.get("PLI_ANALYSIS", ROOT.parent / "analysis"))
HERE = Path(os.environ.get("PLI_OUT") or sys.exit("set PLI_OUT (simulation/README.md)")) / "results/2026-08-12_seed_repeatability"
MACRO = ROOT / "macros/closure_sim1_sc44.mac"
RUNNER = ROOT / "runner/run_gate_parallel.sh"
SEEDS = [2_000_000, 3_000_000, 4_000_000]
TUBES = {"T1_water_shibuya2020_n2": 1.9215, "T2_adipose_avachat2024": 2.54,
         "T3_adipose_moskal2023": 2.72, "T4_hepatic_avachat2024": 2.04,
         "T5_muscle_avachat2024": 2.03, "T6_brain_healthy_moskal2024": 2.72,
         "T7_brain_tumor_moskal2024": 1.77, "T8_myxoma_moskal2023": 1.92}
TUBE_POS = {"T1_water_shibuya2020_n2": (80.0, 0.0),
            "T2_adipose_avachat2024": (56.569, 56.569),
            "T3_adipose_moskal2023": (0.0, 80.0),
            "T4_hepatic_avachat2024": (-56.569, 56.569),
            "T5_muscle_avachat2024": (-80.0, 0.0),
            "T6_brain_healthy_moskal2024": (-56.569, -56.569),
            "T7_brain_tumor_moskal2024": (0.0, -80.0),
            "T8_myxoma_moskal2023": (56.569, -56.569)}
TAU_MIN = 1.0


def sh(cmd, **kw):
    print("+", " ".join(str(c) for c in cmd), flush=True)
    subprocess.run([str(c) for c in cmd], check=True, **kw)


def make_runner(base_seed: int) -> Path:
    dst = ROOT / "runner" / f"runner_seed{base_seed}.sh"
    txt = RUNNER.read_text().replace("1000000 + BATCH_ID * 1000",
                                     f"{base_seed} + BATCH_ID * 1000")
    dst.write_text(txt)
    dst.chmod(0o755)
    return dst


def concat_h5(parts, out):
    with h5py.File(parts[0]) as h0:
        keys = [(g, k) for g in h0 for k in h0[g]]
    data = {gk: [] for gk in keys}
    for p in parts:
        with h5py.File(p) as h:
            for g, k in keys:
                data[(g, k)].append(h[g][k][:])
    with h5py.File(out, "w") as ho:
        for (g, k), arrs in data.items():
            ho.create_dataset(f"{g}/{k}", data=np.concatenate(arrs))


def v3_per_tube(det_h5):
    C = 299.792458
    with h5py.File(det_h5) as h:
        g = h["detector"]
        p1 = np.stack([g["a1_x"][:], g["a1_y"][:], g["a1_z"][:]], -1)
        p2 = np.stack([g["a2_x"][:], g["a2_y"][:], g["a2_z"][:]], -1)
        pp = np.stack([g["p_x"][:], g["p_y"][:], g["p_z"][:]], -1)
        t1 = g["a1_time"][:] * 1e9
        t2 = g["a2_time"][:] * 1e9
        tp = g["p_time"][:] * 1e9
    lor = p2 - p1
    dirn = lor / np.clip(np.linalg.norm(lor, axis=1, keepdims=True), 1e-9, None)
    P = 0.5 * (p1 + p2) + (0.5 * C * (t1 - t2))[:, None] * dirn
    d1 = np.linalg.norm(p1 - P, axis=1)
    dp = np.linalg.norm(pp - P, axis=1)
    tau = (t1 - d1 / C) - (tp - dp / C)
    out = {}
    for name, cfg in TUBES.items():
        cx, cy = TUBE_POS[name]
        sel = ((P[:, 0] - cx) ** 2 + (P[:, 1] - cy) ** 2 <= 19.0 ** 2) & \
              (np.abs(P[:, 2]) <= 29.0)
        t = tau[sel]
        tail = t[t > TAU_MIN]
        out[name] = dict(tau_hat=float(tail.mean() - TAU_MIN),
                         n_tail=int(len(tail)), tau_cfg=cfg)
    return out


def main():
    results = {}
    for base in SEEDS:
        sdir = HERE / f"seed_{base}"
        sdir.mkdir(parents=True, exist_ok=True)
        runner = make_runner(base)
        fdet = sdir / "triples_detector_filtered.h5"
        if fdet.exists():
            print(f"[seedrep] seed {base}: filtered h5 exists, analysis only",
                  flush=True)
            results[str(base)] = v3_per_tube(fdet)
            json.dump(results, open(HERE / "seed_results_partial.json", "w"),
                      indent=1)
            print(f"[seedrep] seed {base} done", flush=True)
            continue
        parts_det, parts_tru = [], []
        for b in range(8):
            bdir = sdir / f"batch_{b}"
            for attempt in (1, 2):
                try:
                    sh(["bash", runner, "--macro", MACRO,
                        "--total-primaries", "12000000", "--threads", "16",
                        "--bounded-max-gb", "5", "--batch-id", str(b),
                        "--label", f"seedrep{base}", "--output-dir", bdir])
                    break
                except subprocess.CalledProcessError:
                    if attempt == 2:
                        raise
                    print(f"[seedrep] batch {b} failed, retrying once",
                          flush=True)
                    shutil.rmtree(bdir, ignore_errors=True)
            merged = bdir / "sim-merged.root"
            (bdir / "sim_merged.root").symlink_to(merged)
            sh([sys.executable, ANALYSIS / "reconstruction/extract_triples_chunked.py",
                "--root-dir", bdir, "--out-dir", bdir, "--tag", f"_b{b}"])
            for f in bdir.glob("sim-00*.root"):
                f.unlink()
            merged.unlink(missing_ok=True)
            (bdir / "sim_merged.root").unlink(missing_ok=True)
            parts_det.append(next(bdir.glob("triples_detector*_b%d.h5" % b)))
            parts_tru.append(next(bdir.glob("triples_truth*_b%d.h5" % b)))
        det = sdir / "triples_detector_all.h5"
        tru = sdir / "triples_truth_all.h5"
        concat_h5(parts_det, det)
        concat_h5(parts_tru, tru)
        fdet = sdir / "triples_detector_filtered.h5"
        sh([sys.executable, ANALYSIS / "reconstruction/filter_triples_tof_pann.py",
            "--triples-detector", det, "--triples-truth", tru,
            "--r-cut", "75.0", "--output-detector", fdet,
            "--output-truth", sdir / "triples_truth_filtered.h5"])
        results[str(base)] = v3_per_tube(fdet)
        json.dump(results, open(HERE / "seed_results_partial.json", "w"),
                  indent=1)
        print(f"[seedrep] seed {base} done", flush=True)
    json.dump(results, open(HERE / "seed_results.json", "w"), indent=1)
    print("SEED_REPEATABILITY_DONE", flush=True)


if __name__ == "__main__":
    main()

# RETRACTION 2026-08-12: an earlier note here claimed a ~2.6x plugin
# slowdown vs April. FALSE: the primaries-mode patch had silently failed
# to apply, so all failed attempts ran time-mode 3.9M-primary shards
# (2.6x April's 1.5M). Runtime and 4GB-cap deaths fully explained by
# shard size. Config now verified: 8 batches x 12M primaries, cap 5GB.
