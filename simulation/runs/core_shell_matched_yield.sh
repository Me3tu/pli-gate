#!/bin/bash
# SPDX-License-Identifier: LGPL-3.0-only
# Copyright (C) 2026 Dhritman Dhar and German Cancer Research Center (DKFZ)
# R2c matched-yield paired core-shell arms (REGISTERED_R2C.md). Waits until the oxygen v8 GATE stage has finished
# (one GATE runner on the box at a time), then runs four arms with the R2b recipe of
# 2026-08-13_referee_response/run_r2.sh: same macros, seeds, shards, primaries; database v8 with equal o-Ps yields;
# PLI_DISABLE_3GAMMA=1. Bulky per-shard logs and ROOT files of these arms are deleted after the PLI_R2 records are
# extracted, as in run_r2.sh (recovery: rerun this script; about 3 min per arm).
set -u
source "$(dirname "${BASH_SOURCE[0]}")/../runner/paths.sh"
HERE=$PLI_OUT/results/2026-10-08_r2c_matched_yield
R2=$PLI_SIM/macros
OXY=$PLI_OUT/results/2026-10-08_oxygen_v8/run.log
OUT=$PLI_OUT/r2c_runs
LOG=$HERE/run.log
say(){ echo "[$(date +%F' '%T)] $*" >> $LOG; }
mkdir -p $OUT "$HERE"
# v1.1.0: the wait for the oxygen run (one GATE runner per machine) is disabled.
#[v1.1.0 disabled] say "waiting for the oxygen GATE stage"
#[v1.1.0 disabled] while ! grep -q -E "GATE done|OXY_FAILED|OXY_DONE" $OXY; do sleep 60; done
#[v1.1.0 disabled] say "oxygen GATE stage finished; starting R2c"
cd "$PLI_SIM"
source "$PLI_ENV_SCRIPT" >/dev/null 2>&1
JSON=$PLI_PKG/materials/variants/r2c_equal_yield_materials.json
say "materials $(sha256sum $JSON | cut -c1-12); Gate $(command -v Gate) $(sha256sum $(command -v Gate) | cut -c1-12)"
run_arm () {           # $1 label  $2 macro  $3 mode
    local label=$1 macro=$2 mode=$3
    say "START $label mode=$mode"
    export PLI_MATERIALS_JSON=$JSON
    export PLI_LIFETIME_MODE=$mode
    export PLI_FIXED_TAU_NS=2.0
    export PLI_LOG_XMAT=1
    export PLI_DISABLE_3GAMMA=1
    unset PLI_DIRECT_TAU_PS PLI_INFLIGHT_STOCK GATE_BIN_DIR
    mkdir -p $OUT/$label $HERE/r2c_runs/$label
    runner/run_gate_parallel.sh --macro $macro --total-primaries 10000000 --threads 16 --bounded-max-gb 5 \
        --batch-id 0 --label r2c_$label --output-dir $OUT/$label >> $HERE/gate.log 2>&1
    local rc=$?
    say "END $label rc=$rc"
    grep -h "^PLI_R2:" $OUT/$label/sim-*.log 2>/dev/null | gzip -1 > $HERE/r2c_runs/$label/pli_r2.txt.gz
    grep -h "^PLI_CHANNEL_TALLY:" $OUT/$label/sim-*.log 2>/dev/null > $HERE/r2c_runs/$label/channel_tally.txt
    grep -h -i "3-gamma\|PLI_DISABLE_3GAMMA\|LIFETIME_MODE" $OUT/$label/sim-0001.log | head -5 > $HERE/r2c_runs/$label/switches.txt
    cp $OUT/$label/run.log $HERE/r2c_runs/$label/runner.log 2>/dev/null
    rm -f $OUT/$label/sim-*.log $OUT/$label/sim-*.root
}
run_arm A_transport $R2/core_shell_pairA.mac transport
run_arm A_source    $R2/core_shell_pairA.mac source
run_arm B_transport $R2/core_shell_pairB.mac transport
run_arm B_source    $R2/core_shell_pairB.mac source
say "R2C_DONE"
