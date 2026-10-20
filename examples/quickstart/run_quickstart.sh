#!/usr/bin/env bash
# SPDX-License-Identifier: LGPL-3.0-only
# Copyright (C) 2026 Dhritman Dhar and German Cancer Research Center (DKFZ)
# Run the quickstart and extract the lifetime. Two commands, one number.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG="$(cd "$HERE/../.." && pwd)"

: "${PLI_MATERIALS_JSON:=$PKG/materials/pli_materials.json}"
export PLI_MATERIALS_JSON
export PLI_LIFETIME_LOG=1

command -v Gate >/dev/null || { echo "ERROR: Gate not on PATH. Source your GATE environment first." >&2; exit 1; }

cd "$HERE"
echo "==> running quickstart.mac (a minute or two)"
# bash's SECONDS, not /usr/bin/time: the GNU time binary is not installed by
# default on every distribution and this script must not need it.
START=$SECONDS
Gate quickstart.mac > quickstart_run.log 2>&1 || {
    echo "GATE failed; last 20 lines:" >&2; tail -20 quickstart_run.log >&2; exit 1; }
echo "    wall $((SECONDS - START)) s"
echo
echo "==> extracting the lifetime"
python3 extract_tau.py quickstart_run.log
