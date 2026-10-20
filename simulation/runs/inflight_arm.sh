#!/bin/bash
# SPDX-License-Identifier: LGPL-3.0-only
# Copyright (C) 2026 Dhritman Dhar and German Cancer Research Center (DKFZ)
# In-flight sensitivity arm: C2 replicated with the hybrid binary. See REGISTERED_INFLIGHT.md.
set -u
source "$(dirname "${BASH_SOURCE[0]}")/../runner/paths.sh"
L=$PLI_OUT/results/2026-09-30_feedback3/inflight
OUT=$PLI_OUT/feedback3_inflight/c2_hybrid
MACRO=$PLI_SIM/macros/eight_tube_sc44.mac
PER_SHARD_SEC=0.137933
cd "$PLI_SIM"; source "$PLI_ENV_SCRIPT" >/dev/null 2>&1; [ -z "${PLI_VENV:-}" ] || source "$PLI_VENV/bin/activate"
export PLI_MATERIALS_JSON=${PLI_MATERIALS_JSON:-$PLI_PKG/materials/pli_materials.json}
export PLI_DIRECT_TAU_PS=388
export PLI_INFLIGHT_STOCK=1
export GATE_BIN_DIR=${GATE_INFLIGHT_BIN_DIR:?set GATE_INFLIGHT_BIN_DIR to the bin/ directory of the in-flight variant build (integration/PATCH_NOTES.md)}
say(){ echo "[$(date +%H:%M:%S)] $*" >> "$L/run.log"; }
mkdir -p "$OUT" "$L"
say "=== HYBRID C2 START: in-flight kept, 388 ps, 10x12 + 1x8 shards x ${PER_SHARD_SEC}s ==="
for b in $(seq 0 10); do
  bd=$(printf %02d $b)
  THREADS=12; [ $b -eq 10 ] && THREADS=8
  TSEC=$(python3 -c "print($PER_SHARD_SEC*$THREADS)")
  [ -f "$OUT/singles_batch_$bd.h5" ] && [ -f "$OUT/singles_batch_$bd.json" ] && { say "  batch $bd present, skip"; continue; }
  rm -rf "$OUT/sim"; mkdir -p "$OUT/sim"
  "$PLI_SIM/runner/run_gate_parallel_bin.sh"       --macro "$MACRO"       --total-time-sec $TSEC --threads $THREADS --bounded-max-gb 5       --batch-id $b --label c2_hybrid --output-dir "$OUT/sim" >> "$L/gate.log" 2>&1
  rc=$?; [ $rc -ne 0 ] && { say "ABORT batch $bd gate rc=$rc"; echo HYB_FAILED >> "$L/run.log"; exit 1; }
  grep -q "$GATE_BIN_DIR" "$L/gate.log" || { say "ABORT: runner did not use the hybrid binary"; echo HYB_FAILED >> "$L/run.log"; exit 1; }
  python "$PLI_ANALYSIS/results/2026-06-10_brain_tumor_1B_truetime/scripts/gate_root_to_singles_h5.py"       --root "$OUT/sim/sim-merged.root" --out "$OUT/singles_batch_$bd.h5" --batch-id $b       --keep-511 435 870 --keep-prompt 870 1300       --tight-511 481 541 --tight-prompt 1100 1200 >> "$L/convert.log" 2>&1
  rc=$?; [ $rc -ne 0 ] && { say "ABORT batch $bd convert rc=$rc"; echo HYB_FAILED >> "$L/run.log"; exit 1; }
  rm -rf "$OUT/sim"
  say "  batch $bd done -> $(du -h "$OUT/singles_batch_$bd.h5" | cut -f1) free=$(df -BG --output=avail "$PLI_OUT" | tail -1 | tr -dc '0-9')GB"
done
say "concat"
python "$PLI_ANALYSIS/results/2026-08-16_T6_nema/concat_singles_u16.py"     --in-dir "$OUT" --out "$OUT/singles_c2_hybrid.h5" >> "$L/concat.log" 2>&1
rc=$?; say "concat rc=$rc"
[ $rc -ne 0 ] && { echo HYB_FAILED >> "$L/run.log"; exit 1; }
echo HYB_DONE >> "$L/run.log"
