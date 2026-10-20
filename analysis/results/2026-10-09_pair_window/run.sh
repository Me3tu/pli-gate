#!/bin/bash
# Round 9: pair-window arms on the corrected-yield zero-delay eight-tube stream (C0) and on the zero-delay NEMA stream,
# same inputs and timing smear as 2026-10-08_builder_sensitivity.
set -u
R=${PLI_RES:-$(cd "$(dirname "$0")/.." && pwd)}
H=$R/2026-10-09_pair_window
PY="env -i HOME=$HOME PATH=/usr/local/bin:/usr/bin:/bin PYTHONPATH=/usr/lib/python3.14/site-packages PLI_RAW=${PLI_RAW:-} ${PLI_REPO:-$(cd "$(dirname "$0")/../../.." && pwd)}/venv/bin/python"
cd $H
echo "[$(date +%F' '%T)] start C0" >> run.log
$PY pair_window.py --out $H/c0 > c0.log 2>&1 || { echo "C0 FAILED" >> run.log; exit 1; }
echo "[$(date +%F' '%T)] start NEMA" >> run.log
$PY pair_window.py --singles ${PLI_RAW:-<missing>}/2026-08-16_T6_nema/singles_nema_all.h5 --config $R/2026-08-16_T6_nema/nema_voi_config.json \
    --phantom nema_iq --out $H/nema > nema.log 2>&1 || { echo "NEMA FAILED" >> run.log; exit 1; }
echo "[$(date +%F' '%T)] PW_DONE" >> run.log
