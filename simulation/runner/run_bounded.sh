#!/bin/bash
# SPDX-License-Identifier: LGPL-3.0-only
# Copyright (C) 2026 Dhritman Dhar and German Cancer Research Center (DKFZ)
# Wrap a command in a systemd user cgroup with hard memory limits.
#
# Origin: installed 2026-04-22 after session 1812 crash_diagnostic traced
# the Apr 21 19:07 OOM to reconstruction processes stacking in one tmux
# scope (peak 75 GB RSS + 32.6 GB swap across a 3d 6h session). Bounding
# each heavy job in its own cgroup keeps a runaway from taking down
# unrelated work.
#
# Usage:
#   scripts/run_bounded.sh <cmd> [args...]
#
# Env overrides:
#   BOUNDED_MAX=60G   MemoryMax        hard kernel-enforced ceiling (OOM here
#                                      kills inside the scope, not the host)
#   BOUNDED_HIGH=55G  MemoryHigh       soft throttle — kernel reclaims aggressively
#                                      above this level
#   BOUNDED_SWAP=16G  MemorySwapMax    prevents hours of swap-thrash death
#
# The scope inherits the invoking shell's stdin/stdout/stderr. Systemd logs
# the bounded unit to `systemctl --user status pli-bounded-<stamp>-<pid>`
# until the scope exits (~auto-cleanup on unit termination).
set -eu

MAX=${BOUNDED_MAX:-60G}
HIGH=${BOUNDED_HIGH:-55G}
SWAP=${BOUNDED_SWAP:-16G}

LOGDIR="${PLI_BOUNDED_LOGDIR:-${TMPDIR:-/tmp}/pli-bounded}"
mkdir -p "$LOGDIR"
STAMP=$(date +%Y%m%d_%H%M%S)
UNIT="pli-bounded-${STAMP}-$$"

echo "[run_bounded] unit=${UNIT}  max=${MAX}  high=${HIGH}  swap=${SWAP}" >&2
echo "[run_bounded] cmd: $*" >&2

exec systemd-run --user --scope \
    --unit="${UNIT}" \
    -p MemoryMax=${MAX} \
    -p MemoryHigh=${HIGH} \
    -p MemorySwapMax=${SWAP} \
    -- "$@"
