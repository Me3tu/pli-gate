#!/bin/bash
# SPDX-License-Identifier: LGPL-3.0-only
# Copyright (C) 2026 Dhritman Dhar and German Cancer Research Center (DKFZ)
# C2: Sim-1 eight tubes under PLI_DIRECT_TAU_PS=388.
# Same macro, activity, digitizer and thread pattern as the frozen 0 ps arm;
# the two arms differ ONLY in the direct-annihilation delay.
#
# CONVERTER WINDOW (decision, 2026-08-24): keep the FULL digitizer band with no
# islands. The macro frames [435,1300] keV, so keep-511 [435,870] and
# keep-prompt [870,1300] tile it exactly -- no gap, no overlap, nothing
# discarded. Tight flags stay [481,541] / [1100,1200] so the comparison with
# the frozen arm is exact. This is the direct consequence of retraction 2:
# the ROOT dies at conversion and there is no second chance.
#
# ISOTOPE CHECK: Sim-1 is Sc-44, prompt 1157 keV -> tight prompt [1100,1200].
set -u
source "$(dirname "${BASH_SOURCE[0]}")/../runner/paths.sh"
L=$PLI_OUT/results/2026-08-24_C2_sim1_388
OUT=$PLI_OUT/c2_sim1_388
MACRO=$PLI_SIM/macros/eight_tube_sc44.mac
# PRIMARIES mode. The Sim-1 macro carries setTotalNumberOfPrimaries, so the
# runner must be driven with --total-primaries; --total-time-sec patches only
# setTimeStop/setTimeSlice and leaves the primary count intact, which is how
# six launches ended up asking each shard for 250,000,000 primaries.
NBATCH=${NBATCH:-10}
THREADS=${THREADS:-12}
PER_SHARD=${PER_SHARD:-1953125}       # = 250,000,000 / 128, the frozen arm's per-shard exposure
MAXGB=${MAXGB:-5}
PER_SHARD_SEC=${PER_SHARD_SEC:-0.137933}
TSEC=$(python3 -c "print($PER_SHARD_SEC*$THREADS)")
cd "$PLI_SIM"; source "$PLI_ENV_SCRIPT" >/dev/null 2>&1; [ -z "${PLI_VENV:-}" ] || source "$PLI_VENV/bin/activate"
export PLI_MATERIALS_JSON=${PLI_MATERIALS_JSON:-$PLI_PKG/materials/pli_materials.json}
export PLI_DIRECT_TAU_PS=388
say(){ echo "[$(date +%H:%M:%S)] $*" >> "$L/run.log"; }
mkdir -p "$OUT" "$L"
say "=== C2 START: Sim-1 8 tubes, Sc-44, direct at 388 ps, $NBATCH batches x $THREADS shards x ${PER_SHARD_SEC}s (TIME-driven, Poisson), ${MAXGB}G/proc ==="
say "    converter keeps FULL band [435,1300] keV (no islands)"
for b in $(seq 0 $((NBATCH-1))); do
  bd=$(printf %02d $b)
  [ -f "$OUT/singles_batch_$bd.h5" ] && [ -f "$OUT/singles_batch_$bd.json" ] && { say "  batch $bd present, skip"; continue; }
  rm -rf "$OUT/sim"; mkdir -p "$OUT/sim"
  runner/run_gate_parallel.sh \
      --macro "$MACRO" \
      --total-time-sec $TSEC --threads $THREADS --bounded-max-gb $MAXGB \
      --batch-id $b --label c2_sim1_388 --output-dir "$OUT/sim" >> "$L/gate.log" 2>&1
  rc=$?; [ $rc -ne 0 ] && { say "ABORT batch $bd gate rc=$rc"; echo C2_FAILED >> "$L/run.log"; exit 1; }
  python "$PLI_ANALYSIS/results/2026-06-10_brain_tumor_1B_truetime/scripts/gate_root_to_singles_h5.py" \
      --root "$OUT/sim/sim-merged.root" --out "$OUT/singles_batch_$bd.h5" --batch-id $b \
      --keep-511 435 870 --keep-prompt 870 1300 \
      --tight-511 481 541 --tight-prompt 1100 1200 >> "$L/convert.log" 2>&1
  rc=$?; [ $rc -ne 0 ] && { say "ABORT batch $bd convert rc=$rc"; echo C2_FAILED >> "$L/run.log"; exit 1; }
  python - "$OUT/singles_batch_$bd.h5" <<'GUARD' >> "$L/run.log" 2>&1
import sys, h5py, numpy as np
with h5py.File(sys.argv[1], "r") as h:
    g = h["singles"]; n = g["in_511"].shape[0]
    n5 = int(g["in_511"][:].sum()); np_ = int(g["in_prompt"][:].sum())
    e = g["energy_keV"][:]
    print(f"  guard: {n:,} singles, in_511 {n5:,}, in_prompt {np_:,}, "
          f"E range {e.min():.1f}-{e.max():.1f} keV, "
          f"in 562-1080 gap band {(int(((e>=562)&(e<=1080)).sum())):,}")
    if n5 == 0 or np_ == 0:
        print("  GUARD FAILED: an energy flag is identically zero"); sys.exit(2)
    # the whole point of the widened window: the old gap must NOT be empty
    if int(((e >= 562) & (e <= 1080)).sum()) == 0:
        print("  GUARD FAILED: 562-1080 keV is empty, the widened window did not take"); sys.exit(3)
GUARD
  rc=$?; [ $rc -ne 0 ] && { say "ABORT batch $bd guard rc=$rc"; echo C2_FAILED >> "$L/run.log"; exit 1; }
  rm -rf "$OUT/sim"
  say "  batch $bd/$((NBATCH-1)) done -> $(du -h "$OUT/singles_batch_$bd.h5" | cut -f1) free=$(df -BG --output=avail "$PLI_OUT" | tail -1 | tr -dc '0-9')GB"
done
say "concat"
python "$PLI_ANALYSIS/results/2026-08-16_T6_nema/concat_singles_u16.py" \
    --in-dir "$OUT" --out "$OUT/singles_c2_sim1_388.h5" >> "$L/concat.log" 2>&1
rc=$?; say "concat rc=$rc"
[ $rc -ne 0 ] && { say "concat FAILED; per-batch files KEPT"; echo C2_FAILED >> "$L/run.log"; exit 1; }
say "per-batch files KEPT until a downstream pass verifies the merge"
echo C2_DONE >> "$L/run.log"
