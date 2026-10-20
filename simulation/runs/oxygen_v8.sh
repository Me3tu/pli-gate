#!/bin/bash
# SPDX-License-Identifier: LGPL-3.0-only
# Copyright (C) 2026 Dhritman Dhar and German Cancer Research Center (DKFZ)
# Review round 8 (2026-10-08): rerun of the oxygen-gradient phantom with database v8
# (plasma-based k_O2 = 2.5704e-5 ns^-1 mmHg^-1, materials/pli_materials.json sha 23f607ea...) and the production
# binary (corrected I3 convention). Same recipe as Phase 3s-2 Stage 3 TE (archived
# 2026-05-06_2258_phase3s_2_paired_reruns/phase3s2_finish_orchestrator.py, stage_3_arm): Sc-44 uniform 1 MBq in the
# water cylinder, 250 M primaries = 250 s, 4 sequential batches x 16 threads (3.90625 s per shard), runner
# scripts/run_gate_parallel_musire.sh (seeds 1e6 + 1000 b + t), hadd, chunked extractor (eventID grouping,
# E_PROMPT [1100,1200], E_511 [481,541]), R_cut 75 mm filter, per-insert tail mean above tau_min = 1 ns.
# One GATE runner on the box at a time.
set -u
source "$(dirname "${BASH_SOURCE[0]}")/../runner/paths.sh"
L=$PLI_OUT/results/2026-10-08_oxygen_v8
OUT=$PLI_OUT/oxygen_v8
A=$PLI_ANALYSIS/pli-quadra-archive-2026-08-12/2026-05-06_2258_phase3s_2_paired_reruns
say(){ echo "[$(date +%F' '%T)] $*" >> "$L/run.log"; }
mkdir -p "$L" "$OUT/sim" "$OUT/extracted"
cd "$PLI_SIM"
source "$PLI_ENV_SCRIPT" >/dev/null 2>&1
[ -z "${PLI_VENV:-}" ] || source "$PLI_VENV/bin/activate"
export PLI_MATERIALS_JSON=${PLI_MATERIALS_JSON:-$PLI_PKG/materials/pli_materials.json}
unset PLI_DIRECT_TAU_PS PLI_INFLIGHT_STOCK GATE_BIN_DIR
export PYTHONPATH=${PLI_EXTRA_PYTHONPATH:+$PLI_EXTRA_PYTHONPATH:}${PYTHONPATH:-}
export PARALLELPROJ_C_LIB=${PARALLELPROJ_C_LIB:-/usr/local/lib/libparallelproj_c.so} PARALLELPROJ_CUDA_LIB=${PARALLELPROJ_CUDA_LIB:-/usr/local/lib/libparallelproj_cuda.so}
say "start; materials $(sha256sum $PLI_MATERIALS_JSON | cut -c1-12); Gate $(command -v Gate) $(sha256sum $(command -v Gate) | cut -c1-12); free=$(df -BG --output=avail "$PLI_OUT" | tail -1 | tr -dc '0-9')GB"
# v1.1.0: the macro this sed generated from the archived Phase 3s-2 macro is shipped as macros/oxygen_gradient.mac
MAC=$PLI_SIM/macros/oxygen_gradient.mac
#[v1.1.0 disabled] sed -e 's|^/gate/output/root/setFileName .*|/gate/output/root/setFileName hypoxia_te_v8|' \
#[v1.1.0 disabled]     -e 's|^# Phase 3s-2 — Hypoxia te arm Poisson master mac (auto-gen)|# Review round 8 oxygen rerun (database v8); copy of the Phase 3s-2 hypoxia te master mac|' \
#[v1.1.0 disabled]     "$A/hypoxia_te_v2/hypoxia_te_poisson.mac" > "$MAC"
grep -q "setTimeStop  3.906250 s" "$MAC" || { say "ABORT: unexpected setTimeStop"; echo OXY_FAILED >> "$L/run.log"; exit 1; }
for b in 0 1 2 3; do
  bd=$OUT/sim/batch_0$b
  [ -f "$bd/sim-merged.root" ] && { say "batch $b present, skip"; continue; }
  mkdir -p "$bd"
  runner/run_gate_parallel.sh --macro "$MAC" --total-time-sec 62.5 --threads 16 --bounded-max-gb 4 \
      --batch-id $b --label hypoxia_te_v8 --output-dir "$bd" >> "$L/gate.log" 2>&1
  rc=$?; [ $rc -ne 0 ] || [ ! -f "$bd/sim-merged.root" ] && { say "ABORT batch $b rc=$rc"; echo OXY_FAILED >> "$L/run.log"; exit 1; }
  grep -q "${GATE_INFLIGHT_BIN_DIR:-gate-inflight}" "$L/gate.log" && { say "ABORT: runner used the hybrid binary"; echo OXY_FAILED >> "$L/run.log"; exit 1; }
  say "batch $b done ($(du -h $bd/sim-merged.root | cut -f1))"
done
${ROOT_INSTALL:-$ROOTSYS}/bin/hadd -f "$OUT/sim/sim-merged.root" $OUT/sim/batch_0*/sim-merged.root >> "$L/hadd.log" 2>&1 \
  || { say "ABORT hadd"; echo OXY_FAILED >> "$L/run.log"; exit 1; }
rm -f $OUT/sim/batch_0*/sim-merged.root
ln -sf sim-merged.root "$OUT/sim/sim_merged.root"
say "GATE done ($(du -h $OUT/sim/sim-merged.root | cut -f1)); extracting"
python "$PLI_ANALYSIS/reconstruction/extract_triples_chunked.py" --input-format root --root-dir "$OUT/sim" --out-dir "$OUT/extracted" \
    --tag _hypoxia_te_v8 >> "$L/extract.log" 2>&1 || { say "ABORT extract"; echo OXY_FAILED >> "$L/run.log"; exit 1; }
python "$PLI_ANALYSIS/reconstruction/filter_triples_tof_pann.py" --triples-detector "$OUT/extracted/triples_detector_hypoxia_te_v8.h5" \
    --triples-truth "$OUT/extracted/triples_truth_hypoxia_te_v8.h5" \
    --output-detector "$OUT/triples_detector_hypoxia_te_v8_rcut075mm.h5" \
    --output-truth "$OUT/triples_truth_hypoxia_te_v8_rcut075mm.h5" --r-cut 75.0 >> "$L/filter.log" 2>&1 \
    || { say "ABORT filter"; echo OXY_FAILED >> "$L/run.log"; exit 1; }
python "$A/analyzers/per_insert_hypoxia.py" --detector "$OUT/triples_detector_hypoxia_te_v8_rcut075mm.h5" \
    --truth "$OUT/triples_truth_hypoxia_te_v8_rcut075mm.h5" --tau-min 1.0 --materials-json "$PLI_MATERIALS_JSON" \
    --out-md "$L/hypoxia_te_v8_per_insert.md" --out-csv "$L/hypoxia_te_v8_per_insert.csv" \
    --label "review round 8 oxygen rerun (database v8, 250M)" >> "$L/per_insert.log" 2>&1 \
    || { say "ABORT per-insert"; echo OXY_FAILED >> "$L/run.log"; exit 1; }
say "done"; echo OXY_DONE >> "$L/run.log"
