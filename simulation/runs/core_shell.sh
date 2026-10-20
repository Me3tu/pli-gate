#!/bin/bash
# SPDX-License-Identifier: LGPL-3.0-only
# Copyright (C) 2026 Dhritman Dhar and German Cancer Research Center (DKFZ)
# R2 transport-novelty experiment — 8 arms, launched after the R1 regression
# frees the CPU. Every arm uses the SAME macro and the SAME seed; arms differ
# only by PLI_LIFETIME_MODE (and, for R2a, by the materials JSON).
#
#   R2a (null control) : pair A geometry, equal-tau JSON -> transport, source
#   R2b pair A / pair B: production JSON              -> transport, source, fixed
#
# Because GetLifetime consumes no random numbers, transport and source arms
# run identical transport with an identical random stream and differ only in
# the tau applied to each annihilation: the comparison is paired event-by-event.
set -u
source "$(dirname "${BASH_SOURCE[0]}")/../runner/paths.sh"
HERE=$PLI_OUT/results/2026-08-13_referee_response
mkdir -p "$HERE"
LOG=$HERE/r2_driver.log
cd "$PLI_SIM"
source "$PLI_ENV_SCRIPT" >/dev/null 2>&1

# Sequenced by run_all_remaining.sh; no waiting needed.
echo "[r2] starting at $(date -Is)" >> $LOG

PRIMARIES=10000000
run_arm () {           # $1 label  $2 macro  $3 mode  $4 materials-json
    local label=$1 macro=$2 mode=$3 json=$4
    echo "[r2] START $label mode=$mode $(date -Is)" >> $LOG
    export PLI_MATERIALS_JSON=$json
    export PLI_LIFETIME_MODE=$mode
    export PLI_FIXED_TAU_NS=2.0
    export PLI_LOG_XMAT=1
    runner/run_gate_parallel.sh \
        --macro $macro \
        --total-primaries $PRIMARIES --threads 16 --bounded-max-gb 5 \
        --batch-id 0 --label r2_$label \
        --output-dir $HERE/r2_runs/$label >> $LOG 2>&1
    echo "[r2] END   $label rc=$? $(date -Is)" >> $LOG
    # keep only the PLI_R2 truth records; drop the bulky per-shard logs
    grep -h "^PLI_R2:" $HERE/r2_runs/$label/sim-*.log 2>/dev/null \
        | gzip -1 > $HERE/r2_runs/$label/pli_r2.txt.gz
    grep -h "^PLI_CHANNEL_TALLY:" $HERE/r2_runs/$label/sim-*.log 2>/dev/null \
        > $HERE/r2_runs/$label/channel_tally.txt
    rm -f $HERE/r2_runs/$label/sim-*.log $HERE/r2_runs/$label/sim-*.root
    rm -f $HERE/r2_runs/$label/sim-merged.root
}

PROD_JSON=${PLI_MATERIALS_JSON:-$PLI_PKG/materials/pli_materials.json}
EQ_JSON=$PLI_PKG/materials/variants/r2a_equal_tau_materials.json

# --- R2a null control (equal tau for both materials) ----------------------
run_arm r2a_transport $PLI_SIM/macros/core_shell_pairA.mac transport $EQ_JSON
run_arm r2a_source    $PLI_SIM/macros/core_shell_pairA.mac source    $EQ_JSON

# --- R2b pair A -----------------------------------------------------------
run_arm r2b_A_transport $PLI_SIM/macros/core_shell_pairA.mac transport $PROD_JSON
run_arm r2b_A_source    $PLI_SIM/macros/core_shell_pairA.mac source    $PROD_JSON
run_arm r2b_A_fixed     $PLI_SIM/macros/core_shell_pairA.mac fixed     $PROD_JSON

# --- R2b pair B (sensitivity to the material-pair choice) -----------------
run_arm r2b_B_transport $PLI_SIM/macros/core_shell_pairB.mac transport $PROD_JSON
run_arm r2b_B_source    $PLI_SIM/macros/core_shell_pairB.mac source    $PROD_JSON
run_arm r2b_B_fixed     $PLI_SIM/macros/core_shell_pairB.mac fixed     $PROD_JSON

echo "R2_ALL_ARMS_DONE $(date -Is)" >> $LOG
