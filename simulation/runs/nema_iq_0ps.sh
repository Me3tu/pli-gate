#!/bin/bash
# SPDX-License-Identifier: LGPL-3.0-only
# Copyright (C) 2026 Dhritman Dhar and German Cancer Research Center (DKFZ)
# NEMA IEC 2.5e8 benchmark: per-event Ps physics, standard Quadra.
# Pattern per scaling memory: sequential batches x 16 threads (VIPman prod).
# Total activity 7.1 MBq -> T_total = 2.5e8/7.1e6 = 35.2 s
# 16 batches x total-time-sec 2.2 (runner divides across 16 shards).
set -u
source "$(dirname "${BASH_SOURCE[0]}")/../runner/paths.sh"
cd "$PLI_SIM"
source "$PLI_ENV_SCRIPT" >/dev/null 2>&1 || true
export PLI_MATERIALS_JSON=${PLI_MATERIALS_JSON:-$PLI_PKG/materials/pli_materials.json}
OUT=$PLI_OUT/results/2026-08-11_nema_iec_250M
mkdir -p "$OUT"
for b in $(seq 0 15); do
  runner/run_gate_parallel.sh \
    --macro macros/nema_iq.mac \
    --total-time-sec 2.2 --threads 16 --bounded-max-gb 4 \
    --batch-id $b --label nema250M \
    --output-dir $OUT/batch_$(printf %02d $b) \
    >> $OUT/driver.log 2>&1
  echo "[driver] batch $b rc=$? $(date -Is)" >> $OUT/driver.log
done
echo NEMA_250M_ALL_BATCHES_DONE >> $OUT/driver.log
