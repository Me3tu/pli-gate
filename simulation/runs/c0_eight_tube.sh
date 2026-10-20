#!/bin/bash
# SPDX-License-Identifier: LGPL-3.0-only
# Copyright (C) 2026 Dhritman Dhar and German Cancer Research Center (DKFZ)
# Chain: wait for the in-flight hybrid run, run its analysis, then run C0 (corrected zero-delay
# eight-tube arm, identical to C2 except PLI_DIRECT_TAU_PS unset) and its standard extraction.
# See REGISTERED_C0.md. Sequential on purpose: one GATE runner on the box at a time.
set -u
source "$(dirname "${BASH_SOURCE[0]}")/../runner/paths.sh"
F3=$PLI_OUT/results/2026-09-30_feedback3
I=$F3/inflight; L=$F3/c0_sim1
OUT=$PLI_OUT/feedback3_c0/c0_sim1_0ps
MACRO=$PLI_SIM/macros/eight_tube_sc44.mac
S=$PLI_ANALYSIS/results/2026-08-17_derivation_test
PER_SHARD_SEC=0.137933
say(){ echo "[$(date +%H:%M:%S)] $*" >> "$L/run.log"; }
mkdir -p "$L" "$OUT"
# v1.1.0: the wait for the in-flight arm and its analysis are disabled. C0 does not depend on them.
#[v1.1.0 disabled] say "chain started; waiting for HYB_DONE"
#[v1.1.0 disabled] while ! grep -q "HYB_DONE\|HYB_FAILED" "$I/run.log"; do sleep 60; done
#[v1.1.0 disabled] grep -q HYB_FAILED "$I/run.log" && { say "hybrid FAILED; C0 not started"; exit 1; }
#[v1.1.0 disabled] say "hybrid done; in-flight analysis"
#[v1.1.0 disabled] bash "$I/run_inflight_analysis.sh"; say "in-flight analysis rc=$?"
cd "$PLI_SIM"; source "$PLI_ENV_SCRIPT" >/dev/null 2>&1; [ -z "${PLI_VENV:-}" ] || source "$PLI_VENV/bin/activate"
export PLI_MATERIALS_JSON=${PLI_MATERIALS_JSON:-$PLI_PKG/materials/pli_materials.json}
unset PLI_DIRECT_TAU_PS PLI_INFLIGHT_STOCK GATE_BIN_DIR
say "=== C0 START: corrected binary, direct instantaneous, 10x12 + 1x8 shards x ${PER_SHARD_SEC}s ==="
for b in $(seq 0 10); do
  bd=$(printf %02d $b)
  THREADS=12; [ $b -eq 10 ] && THREADS=8
  TSEC=$(python3 -c "print($PER_SHARD_SEC*$THREADS)")
  [ -f "$OUT/singles_batch_$bd.h5" ] && [ -f "$OUT/singles_batch_$bd.json" ] && { say "  batch $bd present, skip"; continue; }
  rm -rf "$OUT/sim"; mkdir -p "$OUT/sim"
  runner/run_gate_parallel.sh --macro "$MACRO" --total-time-sec $TSEC --threads $THREADS --bounded-max-gb 5 \
      --batch-id $b --label c0_sim1_0ps --output-dir "$OUT/sim" >> "$L/gate.log" 2>&1
  rc=$?; [ $rc -ne 0 ] && { say "ABORT batch $bd gate rc=$rc"; echo C0_FAILED >> "$L/run.log"; exit 1; }
  grep -q "${GATE_INFLIGHT_BIN_DIR:-gate-inflight}" "$L/gate.log" && { say "ABORT: runner used the hybrid binary"; echo C0_FAILED >> "$L/run.log"; exit 1; }
  python "$PLI_ANALYSIS/results/2026-06-10_brain_tumor_1B_truetime/scripts/gate_root_to_singles_h5.py" \
      --root "$OUT/sim/sim-merged.root" --out "$OUT/singles_batch_$bd.h5" --batch-id $b \
      --keep-511 435 870 --keep-prompt 870 1300 --tight-511 481 541 --tight-prompt 1100 1200 >> "$L/convert.log" 2>&1
  rc=$?; [ $rc -ne 0 ] && { say "ABORT batch $bd convert rc=$rc"; echo C0_FAILED >> "$L/run.log"; exit 1; }
  rm -rf "$OUT/sim"
  say "  batch $bd done -> $(du -h "$OUT/singles_batch_$bd.h5" | cut -f1) free=$(df -BG --output=avail "$PLI_OUT" | tail -1 | tr -dc '0-9')GB"
done
say "concat"
python "$PLI_ANALYSIS/results/2026-08-16_T6_nema/concat_singles_u16.py" --in-dir "$OUT" --out "$OUT/singles_c0_sim1_0ps.h5" >> "$L/concat.log" 2>&1
rc=$?; say "concat rc=$rc"; [ $rc -ne 0 ] && { echo C0_FAILED >> "$L/run.log"; exit 1; }
echo C0_SIM_DONE >> "$L/run.log"
M=$OUT/singles_c0_sim1_0ps.h5
for RULE in closest all; do
  python "$S/sim1_extract.py" --singles "$M" --config "$S/sim1_voi_config.json" --phantom sim1_tubes --matching $RULE \
      --sigma-t-ps 90.9 --out "$L/out_$RULE" >> "$L/extract_$RULE.log" 2>&1; say "  extract/$RULE rc=$?"
  python "$S/measure_oversubtraction.py" --singles "$M" --config "$S/sim1_voi_config.json" --phantom sim1_tubes \
      --label sim1 --matching $RULE --out "$L/oversub_$RULE.json" >> "$L/oversub_$RULE.log" 2>&1; say "  oversub/$RULE rc=$?"
done
python "$PLI_ANALYSIS/results/2026-08-19_V4_stats/dump_shard_hists.py" --singles "$M" --config "$S/sim1_voi_config.json" \
    --phantom sim1_tubes --label c0_0ps --out "$L/hists_c0_0ps.npz" >> "$L/dump.log" 2>&1; say "  dump rc=$?"
echo C0_DONE >> "$L/run.log"
