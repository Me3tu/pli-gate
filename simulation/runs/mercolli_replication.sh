#!/bin/bash
# SPDX-License-Identifier: LGPL-3.0-only
# Copyright (C) 2026 Dhritman Dhar and German Cancer Research Center (DKFZ)
# G1 at TRUE exposure: Mercolli's 15 min separated acquisition, their per-tube
# activities, corrected densities, direct annihilation at 388 ps.
set -u
source "$(dirname "${BASH_SOURCE[0]}")/../runner/paths.sh"
L=$PLI_OUT/results/2026-08-22_G1_exposure
OUT=$PLI_OUT/g1_exposure   # 30 batches + concat need ~55 GB; / has only 47
cd "$PLI_SIM"; source "$PLI_ENV_SCRIPT" >/dev/null 2>&1; [ -z "${PLI_VENV:-}" ] || source "$PLI_VENV/bin/activate"
export PLI_MATERIALS_JSON=${PLI_MATERIALS_JSON:-$PLI_PKG/materials/pli_materials.json}
export PLI_DIRECT_TAU_PS=388          # their fitted tau_2
say(){ echo "[$(date +%H:%M:%S)] $*" >> "$L/run.log"; }
say "=== G1 true-exposure START (restart 2): 900 s, 4 tubes, densities corrected, direct 388 ps ==="
say "    I-124 energy windows passed explicitly; per-batch guard on in_511/in_prompt"
say "    total decays ~4.5e9; 30 batches x 30 s; output on /data"
mkdir -p "$OUT" "$L"
for b in $(seq 0 29); do
  bd=$(printf %02d $b)
  [ -f "$OUT/singles_batch_$bd.h5" ] && { say "  batch $bd present, skip"; continue; }
  rm -rf "$OUT/sim"; mkdir -p "$OUT/sim"
  runner/run_gate_parallel.sh \
      --macro "$PLI_SIM/macros/mercolli_xad4_i124.mac" \
      --total-time-sec 30 --threads 30 --bounded-max-gb 4 \
      --batch-id $b --label g1_exposure --output-dir "$OUT/sim" >> "$L/gate.log" 2>&1
  rc=$?; [ $rc -ne 0 ] && { say "ABORT batch $bd gate rc=$rc"; echo G1_FAILED >> "$L/run.log"; exit 1; }
  # ENERGY WINDOWS. The converter defaults are Sc-44 (keep [460,562] and
  # [1080,1220] keV, tight flags [481,541] and [1100,1200]). I-124's prompt is
  # the 603 keV line and the digitizer frames at [460,639], so those defaults
  # discard every prompt photon and leave in_prompt identically zero -- no triple
  # can be built. Batches 00-05 of the 2026-08-22 run were lost to exactly this.
  #   keep  : [460,556] + [556,639] = the full digitizer band, nothing dropped,
  #           so the 23 keV gap photons stay in the file and the gapped-vs-
  #           ungapped selection remains switchable in analysis.
  #   tight : Mercolli's own disjoint windows, LEW [460,545] and HEW [568,639].
  python "$PLI_ANALYSIS/results/2026-06-10_brain_tumor_1B_truetime/scripts/gate_root_to_singles_h5.py" \
      --root "$OUT/sim/sim-merged.root" --out "$OUT/singles_batch_$bd.h5" --batch-id $b \
      --keep-511 460 556 --keep-prompt 556 639 \
      --tight-511 460 545 --tight-prompt 568 639 \
      >> "$L/convert.log" 2>&1
  rc=$?; [ $rc -ne 0 ] && { say "ABORT batch $bd convert rc=$rc"; echo G1_FAILED >> "$L/run.log"; exit 1; }
  # GUARD: a batch with no prompts, or no 511s, cannot yield a triple. Fail on
  # the first batch rather than after 15 hours.
  python - "$OUT/singles_batch_$bd.h5" <<'GUARD' >> "$L/run.log" 2>&1
import sys, h5py
with h5py.File(sys.argv[1], "r") as h:
    g = h["singles"]
    n = g["in_511"].shape[0]
    n5 = int(g["in_511"][:].sum()); np_ = int(g["in_prompt"][:].sum())
    print(f"  guard: {n:,} singles, in_511 {n5:,}, in_prompt {np_:,}")
    if n5 == 0 or np_ == 0:
        print("  GUARD FAILED: an energy flag is identically zero"); sys.exit(2)
GUARD
  rc=$?; [ $rc -ne 0 ] && { say "ABORT batch $bd energy-window guard rc=$rc"; echo G1_FAILED >> "$L/run.log"; exit 1; }
  rm -rf "$OUT/sim"
  say "  batch $bd/29 done -> $(du -h "$OUT/singles_batch_$bd.h5" | cut -f1)  free=$(df -BG --output=avail "$PLI_OUT" | tail -1 | tr -dc '0-9')GB"
done
say "concat"
python "$PLI_ANALYSIS/results/2026-08-16_T6_nema/concat_singles_u16.py" \
    --in-dir "$OUT" --out "$OUT/singles_g1_exposure.h5" >> "$L/concat.log" 2>&1
rc=$?; say "concat rc=$rc"
# GUARD (added 2026-08-23, after the merge of that run was verified): this
# deletion used to run unconditionally. If the concat failed -- including on the
# block-count mismatch concat_singles_u16.py is written to REFUSE on -- 22 hours
# of simulation was deleted anyway. Delete only after a clean merge that also
# self-reports ok.
if [ $rc -ne 0 ]; then
  say "concat FAILED (rc=$rc); per-batch files KEPT"; echo G1_FAILED >> "$L/run.log"; exit 1
fi
if ! python - "$OUT/singles_g1_exposure.concat.json" <<'CHK'
import json, sys
d = json.load(open(sys.argv[1]))
sys.exit(0 if d.get("ok") and d["contiguous_blocks"] == d["global_blocks"] else 2)
CHK
then
  say "concat manifest not ok; per-batch files KEPT"; echo G1_FAILED >> "$L/run.log"; exit 1
fi
say "merge verified; removing per-batch files"
rm -f "$OUT"/singles_batch_*.h5
echo G1_EXPOSURE_DONE >> "$L/run.log"
