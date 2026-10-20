#!/bin/bash
# SPDX-License-Identifier: LGPL-3.0-only
# Copyright (C) 2026 Dhritman Dhar and German Cancer Research Center (DKFZ)
# Sequential batched runner for the 20 cm water cylinder. Usage: driver.sh <batch ids...>
# Each batch: 16 threads x 0.183125 s (2.93 s), seeds 1,000,000 + 1000*batch + thread (batch ids 200+).
# After a verified conversion the merged ROOT is deleted (recovery: re-run that batch id, same seeds).
source "$(dirname "${BASH_SOURCE[0]}")/../runner/paths.sh"
cd "$PLI_SIM"
source "$PLI_ENV_SCRIPT" >/dev/null 2>&1
export PLI_MATERIALS_JSON=${PLI_MATERIALS_JSON:-$PLI_PKG/materials/pli_materials.json}
OUT=$PLI_OUT/results/2026-10-19_review16/water20cm
CONV=$PLI_ANALYSIS/results/2026-06-10_brain_tumor_1B_truetime/scripts/gate_root_to_singles_h5.py
mkdir -p "$OUT"
for b in "$@"; do
  t0=$(date +%s)
  runner/run_gate_parallel.sh --macro $PLI_SIM/macros/water_cylinder_20cm.mac \
    --total-time-sec 2.93 --threads 16 --bounded-max-gb 5 \
    --batch-id $b --label water20 --output-dir $OUT/sim/batch_$b >> $OUT/driver.log 2>&1
  rc=$?
  echo "[w20] batch $b gate rc=$rc wall=$(( $(date +%s)-t0 ))s $(date -Is)" >> $OUT/driver.log
  [ $rc -ne 0 ] && { echo "[w20] STOP at batch $b" >> $OUT/driver.log; exit 1; }
  ( [ -z "${PLI_VENV:-}" ] || source "$PLI_VENV/bin/activate"; python3 $CONV --root $OUT/sim/batch_$b/sim-merged.root \
      --out $OUT/water20_batch_$b.h5 --batch-id $b > $OUT/water20_batch_$b.convert.log 2>&1 )
  crc=$?
  echo "[w20] batch $b convert rc=$crc $(date -Is)" >> $OUT/driver.log
  if [ $crc -eq 0 ] && [ -s $OUT/water20_batch_$b.h5 ]; then rm -f "$OUT/sim/batch_$b/sim-merged.root"; else echo "[w20] STOP convert $b" >> $OUT/driver.log; exit 1; fi
done
echo "[w20] DONE $* $(date -Is)" >> $OUT/driver.log
