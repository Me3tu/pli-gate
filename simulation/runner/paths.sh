#!/bin/bash
# SPDX-License-Identifier: LGPL-3.0-only
# Copyright (C) 2026 Dhritman Dhar and German Cancer Research Center (DKFZ)
#
# Path setup for the run scripts in simulation/runs/ (release v1.1.0). Sourced,
# not executed. Every location that the production scripts hard-coded for one
# workstation is one of the variables below; simulation/README.md lists them.
#
#   PLI_SIM            this simulation/ folder (set here)
#   PLI_PKG            the release root (set here)
#   PLI_OUT            output root, replaces the workstation data disk (required)
#   PLI_ANALYSIS       the analysis part of the release (default $PLI_PKG/analysis)
#   PLI_ENV_SCRIPT     script that puts Geant4, ROOT and Gate on the path
#                      (default runner/env_example.sh)
#   PLI_VENV           optional Python virtual environment to activate
#   PLI_SCANNER_MACRO  scanner geometry macro, NOT shipped (required)
#   PLI_GATE_MATERIALS_DB, PLI_VIPMAN_DIR
#                      further inputs that are NOT shipped, needed only by the
#                      NEMA, VIP-Man and stock-vs-custom macros; GATE stops with
#                      an alias error if a macro needs one that is not set
PLI_SIM="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PLI_PKG="$(cd "$PLI_SIM/.." && pwd)"
export PLI_SIM PLI_PKG
export PLI_ANALYSIS="${PLI_ANALYSIS:-$PLI_PKG/analysis}"
export PLI_ENV_SCRIPT="${PLI_ENV_SCRIPT:-$PLI_SIM/runner/env_example.sh}"
if [ -z "${PLI_OUT:-}" ]; then
    echo "ERROR: set PLI_OUT to a writable directory for the simulation output (simulation/README.md)" >&2
    exit 2
fi
mkdir -p "$PLI_OUT" && PLI_OUT="$(cd "$PLI_OUT" && pwd)" && export PLI_OUT
if [ -z "${PLI_SCANNER_MACRO:-}" ] || [ ! -f "$PLI_SCANNER_MACRO" ]; then
    echo "ERROR: PLI_SCANNER_MACRO must name the Quadra scanner geometry macro." >&2
    echo "       It is not part of this release; see 'Inputs you must supply' in simulation/README.md." >&2
    exit 2
fi
export PLI_SCANNER_MACRO
for v in PLI_GATE_MATERIALS_DB PLI_VIPMAN_DIR; do
    [ -n "${!v:-}" ] && export "$v"
done
