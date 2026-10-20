#!/bin/bash
# SPDX-License-Identifier: LGPL-3.0-only
# Copyright (C) 2026 Dhritman Dhar and German Cancer Research Center (DKFZ)
# C1: NEMA IQ under the realistic early-time model, PLI_DIRECT_TAU_PS=388.
# Same macro, activity, batching and thread count as the frozen 0 ps run
# (2026-08-11_nema_iec_250M/launch_nema_250M.sh) so the two arms differ ONLY in
# the direct-annihilation delay.
#
# ISOTOPE CHECK (the lesson from the I-124 loss): NEMA IQ is Sc-44, prompt
# 1157 keV, so the converter's DEFAULT windows are the correct ones here. They
# are passed explicitly anyway, so the record shows what was used rather than
# what was defaulted.
set -u
source "$(dirname "${BASH_SOURCE[0]}")/../runner/paths.sh"
L=$PLI_OUT/results/2026-08-24_C1_nema388
OUT=$PLI_OUT/c1_nema388
cd "$PLI_SIM"; source "$PLI_ENV_SCRIPT" >/dev/null 2>&1; [ -z "${PLI_VENV:-}" ] || source "$PLI_VENV/bin/activate"
export PLI_MATERIALS_JSON=${PLI_MATERIALS_JSON:-$PLI_PKG/materials/pli_materials.json}
export PLI_DIRECT_TAU_PS=388
say(){ echo "[$(date +%H:%M:%S)] $*" >> "$L/run.log"; }
mkdir -p "$OUT" "$L"
say "=== C1 START: NEMA IQ, Sc-44, direct at 388 ps, 16 batches x 2.2 s x 16 threads ==="
for b in $(seq 0 15); do
  bd=$(printf %02d $b)
  [ -f "$OUT/singles_batch_$bd.h5" ] && [ -f "$OUT/singles_batch_$bd.json" ] && { say "  batch $bd present, skip"; continue; }
  rm -rf "$OUT/sim"; mkdir -p "$OUT/sim"
  runner/run_gate_parallel.sh \
      --macro "$PLI_SIM/macros/nema_iq.mac" \
      --total-time-sec 2.2 --threads 16 --bounded-max-gb 4 \
      --batch-id $b --label c1_nema388 --output-dir "$OUT/sim" >> "$L/gate.log" 2>&1
  rc=$?; [ $rc -ne 0 ] && { say "ABORT batch $bd gate rc=$rc"; echo C1_FAILED >> "$L/run.log"; exit 1; }
  python "$PLI_ANALYSIS/results/2026-06-10_brain_tumor_1B_truetime/scripts/gate_root_to_singles_h5.py" \
      --root "$OUT/sim/sim-merged.root" --out "$OUT/singles_batch_$bd.h5" --batch-id $b \
      --keep-511 460 562 --keep-prompt 1080 1220 \
      --tight-511 481 541 --tight-prompt 1100 1200 >> "$L/convert.log" 2>&1
  rc=$?; [ $rc -ne 0 ] && { say "ABORT batch $bd convert rc=$rc"; echo C1_FAILED >> "$L/run.log"; exit 1; }
  python - "$OUT/singles_batch_$bd.h5" <<'GUARD' >> "$L/run.log" 2>&1
import sys, h5py
with h5py.File(sys.argv[1], "r") as h:
    g = h["singles"]; n = g["in_511"].shape[0]
    n5 = int(g["in_511"][:].sum()); np_ = int(g["in_prompt"][:].sum())
    print(f"  guard: {n:,} singles, in_511 {n5:,}, in_prompt {np_:,}")
    if n5 == 0 or np_ == 0:
        print("  GUARD FAILED: an energy flag is identically zero"); sys.exit(2)
GUARD
  rc=$?; [ $rc -ne 0 ] && { say "ABORT batch $bd guard rc=$rc"; echo C1_FAILED >> "$L/run.log"; exit 1; }
  rm -rf "$OUT/sim"
  say "  batch $bd/15 done -> $(du -h "$OUT/singles_batch_$bd.h5" | cut -f1) free=$(df -BG --output=avail "$PLI_OUT" | tail -1 | tr -dc '0-9')GB"
done
say "concat"
python "$PLI_ANALYSIS/results/2026-08-16_T6_nema/concat_singles_u16.py" \
    --in-dir "$OUT" --out "$OUT/singles_c1_nema388.h5" >> "$L/concat.log" 2>&1
rc=$?; say "concat rc=$rc"
if [ $rc -ne 0 ]; then say "concat FAILED; per-batch files KEPT"; echo C1_FAILED >> "$L/run.log"; exit 1; fi
say "per-batch files KEPT until a downstream pass verifies the merge"
echo C1_DONE >> "$L/run.log"
