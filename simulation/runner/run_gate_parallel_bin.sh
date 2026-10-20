#!/bin/bash
# SPDX-License-Identifier: LGPL-3.0-only
# Copyright (C) 2026 Dhritman Dhar and German Cancer Research Center (DKFZ)
# ============================================================================
# run_gate_parallel_musire.sh: parallel time-fragmented GATE runner
#
# Pattern: J. Peter's Musiré (PetFragmentTimeSec + RunSimulation).
# Splits a time-based GATE macro into K time-fragmented per-proc macros and
# runs them concurrently under GNU parallel. Each fragment has its own seed,
# its own ROOT output, and its own systemd-run cgroup memory cap (via
# scripts/run_bounded.sh). The K shards are then hadd-merged.
#
# Usage:
#   scripts/run_gate_parallel_musire.sh \
#       --macro <path>                                \
#       (--total-time-sec <T> | --total-primaries <N>) \
#       --threads <K>                                 \
#       --output-dir <dir>                            \
#       [--batch-id <int>]                            \
#       [--bounded-max-gb <G>]                        \
#       [--label <tag>]                               \
#       [--keep-shards]
#
# Modes:
#   --total-time-sec  T : time-based macros (NEMA-style). Each thread runs
#                         T/K simulated seconds. Macro must have
#                         /gate/application/setTimeStop and setTimeSlice.
#   --total-primaries N : event-based macros (VIPman-style). Each thread
#                         runs floor(N/K) primaries. The (N mod K) extras
#                         are distributed one-per to the lowest-numbered
#                         shards. Macro must have setTotalNumberOfPrimaries.
#
# Behaviour (both modes):
#   1. setEngineSeed is forced to 1_000_000 + batch_id*1000 + thread_id
#      (no `auto`).  --batch-id defaults to 0. Pass distinct integers from
#      a multi-batch caller (e.g. run_nema_500M_batched.sh) so successive
#      batches do not reproduce identical events.
#   2. /gate/output/root/setFileName and any /gate/actor/.../save lines
#      are patched per shard so files don't collide.
#   3. After all complete, hadd merges into <output-dir>/sim-merged.root.
#   4. Triple-count consistency is verified across shards vs merged.
# ============================================================================
set -eEuo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"  # v1.1.0: this copy lives in runner/

# ---------------------------------------------------------------- arg parsing
MACRO=""
TOTAL_TIME_SEC=""
TOTAL_PRIMARIES=""
THREADS=""
OUTPUT_DIR=""
BOUNDED_MAX_GB=""
LABEL="run"
KEEP_SHARDS=0
BATCH_ID="0"

usage() {
    sed -n '4,34p' "${BASH_SOURCE[0]}"
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --macro)            MACRO="$2"; shift 2 ;;
        --total-time-sec)   TOTAL_TIME_SEC="$2"; shift 2 ;;
        --total-primaries)  TOTAL_PRIMARIES="$2"; shift 2 ;;
        --threads)          THREADS="$2"; shift 2 ;;
        --output-dir)       OUTPUT_DIR="$2"; shift 2 ;;
        --bounded-max-gb)   BOUNDED_MAX_GB="$2"; shift 2 ;;
        --batch-id)         BATCH_ID="$2"; shift 2 ;;
        --label)            LABEL="$2"; shift 2 ;;
        --keep-shards)      KEEP_SHARDS=1; shift ;;
        -h|--help)          usage ;;
        *) echo "ERROR: unknown arg: $1" >&2; usage ;;
    esac
done

[[ -z "$MACRO" ]]      && { echo "ERROR: --macro required" >&2; exit 2; }
[[ -z "$THREADS" ]]    && { echo "ERROR: --threads required" >&2; exit 2; }
[[ -z "$OUTPUT_DIR" ]] && { echo "ERROR: --output-dir required" >&2; exit 2; }
[[ ! -f "$MACRO" ]]    && { echo "ERROR: macro not found: $MACRO" >&2; exit 2; }

# --batch-id must be a non-negative integer, and thread_id ≤ 999 to keep the
# (batch_id*1000 + thread_id) seed offset collision-free across batches.
if ! [[ "$BATCH_ID" =~ ^[0-9]+$ ]]; then
    echo "ERROR: --batch-id must be a non-negative integer (got '$BATCH_ID')" >&2
    exit 2
fi
if [[ "$THREADS" -gt 999 ]]; then
    echo "ERROR: --threads must be ≤ 999 (got $THREADS)" >&2
    exit 2
fi

# Mutually-exclusive: exactly one of --total-time-sec / --total-primaries.
if [[ -n "$TOTAL_TIME_SEC" && -n "$TOTAL_PRIMARIES" ]]; then
    echo "ERROR: --total-time-sec and --total-primaries are mutually exclusive" >&2
    exit 2
fi
if [[ -z "$TOTAL_TIME_SEC" && -z "$TOTAL_PRIMARIES" ]]; then
    echo "ERROR: must give exactly one of --total-time-sec or --total-primaries" >&2
    exit 2
fi

if [[ -n "$TOTAL_PRIMARIES" ]]; then
    MODE="primaries"
else
    MODE="time"
fi

if [[ -z "$BOUNDED_MAX_GB" ]]; then
    # Default: floor(60 / K). 60 GB total cap is the project standard.
    BOUNDED_MAX_GB=$(( 60 / THREADS ))
    [[ $BOUNDED_MAX_GB -lt 2 ]] && BOUNDED_MAX_GB=2
fi

# ---------------------------------------------------------------- environment
if ! command -v parallel >/dev/null; then
    cat >&2 <<EOF
ERROR: GNU parallel not found on PATH.
  Install on Arch:  sudo pacman -S parallel
  Install on Debian/Ubuntu:  sudo apt-get install parallel
EOF
    exit 3
fi

if [[ -f "${PLI_ENV_SCRIPT:-$PROJECT_ROOT/runner/env_example.sh}" ]]; then
    # shellcheck disable=SC1091
    source "${PLI_ENV_SCRIPT:-$PROJECT_ROOT/runner/env_example.sh}"
fi
# feedback round 3: optional binary override (in-flight sensitivity arm), unset = unchanged
if [[ -n "${GATE_BIN_DIR:-}" ]]; then export PATH="${GATE_BIN_DIR}:${PATH}"; fi
if ! command -v Gate >/dev/null; then
    echo "ERROR: Gate not on PATH after sourcing GATE env" >&2; exit 4
fi
if ! command -v hadd >/dev/null; then
    echo "ERROR: hadd (ROOT) not on PATH" >&2; exit 4
fi

# PLI_MATERIALS_JSON: passthrough if exported, else auto-detect by macro path.
# VIPman macros need the extended materials JSON (17 atlas tissues + LSO/PTFE):
# without it, GATE falls back to scalar I3=0.35 for atlas tissues, which is OK
# for wallclock calibration but wrong for production physics.
if [[ -z "${PLI_MATERIALS_JSON:-}" ]]; then
    if [[ "$MACRO" == *vipman* || "$MACRO" == *VIPman* ]]; then
        VIPMAN_JSON="$PROJECT_ROOT/../materials/pli_materials_vipman.json"
        if [[ -f "$VIPMAN_JSON" ]]; then
            export PLI_MATERIALS_JSON="$VIPMAN_JSON"
            echo "[$(date -Iseconds)] auto-set PLI_MATERIALS_JSON=$VIPMAN_JSON" >&2
        else
            echo "WARNING: VIPman macro detected but $VIPMAN_JSON not found; using default" >&2
        fi
    fi
fi

mkdir -p "$OUTPUT_DIR"
RUN_LOG="$OUTPUT_DIR/run.log"
PARALLEL_LOG="$OUTPUT_DIR/parallel.log"
HADD_LOG="$OUTPUT_DIR/hadd.log"
MERGE_CHECK="$OUTPUT_DIR/merge_check.log"

# ---------------------------------------------------------------- per-shard split
# Time mode: per-shard time = TOTAL_TIME_SEC / THREADS (float division).
# Primary mode: per-shard primaries = floor(N / K), with (N mod K) extras
#   distributed one-per to the lowest-numbered shards (1..(N mod K)).
PER_SHARD_TIME=""
PER_SHARD_BASE=""
PER_SHARD_REMAINDER=""
if [[ "$MODE" == "time" ]]; then
    PER_SHARD_TIME=$(python3 -c "print(${TOTAL_TIME_SEC} / ${THREADS})")
else
    PER_SHARD_BASE=$(( TOTAL_PRIMARIES / THREADS ))
    PER_SHARD_REMAINDER=$(( TOTAL_PRIMARIES - PER_SHARD_BASE * THREADS ))
fi

{
    echo "[$(date -Iseconds)] run_gate_parallel_musire.sh"
    echo "  macro          = $MACRO"
    echo "  mode           = $MODE"
    if [[ "$MODE" == "time" ]]; then
        echo "  total_time_sec = $TOTAL_TIME_SEC"
        echo "  per_shard_time = $PER_SHARD_TIME s"
    else
        echo "  total_primaries  = $TOTAL_PRIMARIES"
        echo "  per_shard_base   = $PER_SHARD_BASE"
        echo "  per_shard_extras = $PER_SHARD_REMAINDER (first $PER_SHARD_REMAINDER shards get +1)"
    fi
    echo "  threads        = $THREADS"
    echo "  output_dir     = $OUTPUT_DIR"
    echo "  bounded_max_gb = $BOUNDED_MAX_GB (per proc)"
    echo "  total cap      = $((BOUNDED_MAX_GB * THREADS)) GB"
    echo "  batch_id       = $BATCH_ID  (seed = 1000000 + batch_id*1000 + thread_id)"
    echo "  label          = $LABEL"
    echo "  Gate           = $(command -v Gate)"
    echo "  parallel       = $(command -v parallel)"
    echo "  hadd           = $(command -v hadd)"
} | tee "$RUN_LOG"

# ---------------------------------------------------------------- build per-shard macros
make_shard_macro() {
    local i="$1"          # 1-indexed
    local tag             # zero-padded
    tag=$(printf "%04d" "$i")
    # Seed = 1_000_000 + batch_id*1000 + thread_id.  --batch-id varies
    # across calls in a multi-batch loop, and thread_id (i, 1..K) varies inside
    # one call.  Each batch therefore gets an independent random stream.
    # Without the batch_id term, repeated invocations re-use the same K
    # seeds and a time-based macro (NEMA-style with fixed setTimeStop)
    # produces 32 byte-identical copies. See findings.md of session
    # 2026-04-26_0956 for the original failure mode.
    local seed=$((1000000 + BATCH_ID * 1000 + i))
    local out="$OUTPUT_DIR/sim-${tag}"
    local mac="$OUTPUT_DIR/sim-${tag}.mac"

    cp "$MACRO" "$mac"

    # Patch ROOT output filename (GATE auto-appends .root)
    sed -i "s|^/gate/output/root/setFileName.*|/gate/output/root/setFileName ${out}|" "$mac"

    # Patch any /gate/actor/.../save <name> lines so per-actor outputs don't
    # collide between shards. Replace just the trailing filename.
    sed -i -E "s|^(/gate/actor/[^[:space:]]+/save)[[:space:]]+.*|\\1 ${out}_actor|" "$mac"

    if [[ "$MODE" == "time" ]]; then
        # Time-fragment: setTimeStop and setTimeSlice both = per-shard time.
        if grep -q "^/gate/application/setTimeStop" "$mac"; then
            sed -i "s|^/gate/application/setTimeStop.*|/gate/application/setTimeStop ${PER_SHARD_TIME} s|" "$mac"
        else
            echo "ERROR: macro missing /gate/application/setTimeStop line (--total-time-sec mode)" >&2; exit 5
        fi
        if grep -q "^/gate/application/setTimeSlice" "$mac"; then
            sed -i "s|^/gate/application/setTimeSlice.*|/gate/application/setTimeSlice ${PER_SHARD_TIME} s|" "$mac"
        else
            echo "/gate/application/setTimeSlice ${PER_SHARD_TIME} s" >> "$mac"
        fi
    else
        # Primary mode: setTotalNumberOfPrimaries = base + 1 if i <= remainder.
        local per_shard_n=$PER_SHARD_BASE
        if [[ $i -le $PER_SHARD_REMAINDER ]]; then
            per_shard_n=$((PER_SHARD_BASE + 1))
        fi
        if grep -q "^/gate/application/setTotalNumberOfPrimaries" "$mac"; then
            sed -i "s|^/gate/application/setTotalNumberOfPrimaries.*|/gate/application/setTotalNumberOfPrimaries ${per_shard_n}|" "$mac"
        else
            echo "ERROR: macro missing /gate/application/setTotalNumberOfPrimaries line (--total-primaries mode)" >&2; exit 5
        fi
        # Leave any setTimeStop/setTimeSlice lines as-is (inert when
        # setTotalNumberOfPrimaries is set).
    fi

    # Force per-shard seed (no `auto`).
    if grep -q "^/gate/random/setEngineSeed" "$mac"; then
        sed -i "s|^/gate/random/setEngineSeed.*|/gate/random/setEngineSeed ${seed}|" "$mac"
    else
        echo "ERROR: macro missing /gate/random/setEngineSeed line" >&2; exit 5
    fi
}

echo "[$(date -Iseconds)] building $THREADS shard macros..." | tee -a "$RUN_LOG"
for i in $(seq 1 "$THREADS"); do
    make_shard_macro "$i"
done

# ---------------------------------------------------------------- launch
# Worker function: invoke run_bounded.sh + Gate for one shard.
shard_run() {
    local i="$1"
    local tag
    tag=$(printf "%04d" "$i")
    local mac="$OUTPUT_DIR/sim-${tag}.mac"
    local log="$OUTPUT_DIR/sim-${tag}.log"
    # IMPORTANT: BOUNDED_HIGH = BOUNDED_MAX (no soft throttle).
    # The default 90% HIGH causes catastrophic slowdown when per-proc RSS
    # nudges above HIGH but below MAX: kernel enters aggressive
    # __mem_cgroup_handle_over_high reclaim and locks procs in state D
    # for hours. Setting HIGH=MAX disables soft throttle. OOM-kill at
    # MAX is the only memory cap. Lesson: session 2026-04-25_1057 lost
    # ~5h to this on Run B (K=5, M=10M, MAX=11G, HIGH=9G).
    # PLI_LIFETIME_LOG defaults to 1 (preserves prior production behaviour).
    # Caller may export PLI_LIFETIME_LOG=0 to skip oracle log emission,
    # important when the simulation has a high-activity small-volume source
    # (e.g. Steinberger 2026-06-22 disks at 2.5 MBq × 60s × dense glass
    # → 2.5M PLI_LIFETIME lines/shard ≈ 170 MB log + 20-min stall on log I/O).
    local pli_log_flag="${PLI_LIFETIME_LOG:-1}"
    BOUNDED_MAX="${BOUNDED_MAX_GB}G" \
    BOUNDED_HIGH="${BOUNDED_MAX_GB}G" \
    BOUNDED_SWAP="1G" \
        "$PROJECT_ROOT/runner/run_bounded.sh" \
        bash -c "PLI_LIFETIME_LOG=${pli_log_flag} Gate '${mac}' > '${log}' 2>&1"
}
export -f shard_run
export OUTPUT_DIR PROJECT_ROOT BOUNDED_MAX_GB PLI_LIFETIME_LOG

T_START=$(date +%s)
echo "[$(date -Iseconds)] launching $THREADS parallel shards..." | tee -a "$RUN_LOG"

# Run K parallel shards. --halt now,fail=1 makes the entire job abort if any
# shard exits non-zero.
seq 1 "$THREADS" | parallel -j "$THREADS" \
    --joblog "$PARALLEL_LOG" \
    --halt now,fail=1 \
    shard_run {} || {
    echo "[$(date -Iseconds)] ERROR: parallel reported failure (see $PARALLEL_LOG)" | tee -a "$RUN_LOG"
    exit 6
}

T_END=$(date +%s)
ELAPSED=$((T_END - T_START))
echo "[$(date -Iseconds)] all $THREADS shards finished. wall=${ELAPSED}s" | tee -a "$RUN_LOG"

# ---------------------------------------------------------------- merge
SHARDS=("$OUTPUT_DIR"/sim-[0-9][0-9][0-9][0-9].root)
if [[ ${#SHARDS[@]} -lt $THREADS ]]; then
    echo "ERROR: only ${#SHARDS[@]} shard ROOT files (expected $THREADS)" \
        | tee -a "$RUN_LOG"
    exit 7
fi

MERGED="$OUTPUT_DIR/sim-merged.root"
echo "[$(date -Iseconds)] merging ${#SHARDS[@]} shards via hadd ..." | tee -a "$RUN_LOG"
hadd -k -n 0 "$MERGED" "${SHARDS[@]}" > "$HADD_LOG" 2>&1
ls -la "$MERGED" | tee -a "$RUN_LOG"

# ---------------------------------------------------------------- verify
# Count Coincidences entries in each shard and the merged file (using python
# uproot, lightweight and present in the project venv).
verify_python() {
    python3 - "$@" <<'PY'
import sys
import uproot

paths = sys.argv[1:]
counts = {}
for p in paths:
    with uproot.open(p) as f:
        nc = int(f["Coincidences"].num_entries) if "Coincidences" in f else 0
        ns = int(f["Singles"].num_entries) if "Singles" in f else 0
    counts[p] = {"Coincidences": nc, "Singles": ns}

# print per-file
for p, c in counts.items():
    print(f"{p}\tSingles={c['Singles']}\tCoincidences={c['Coincidences']}")

# Sums vs merged. Last path is merged.
merged = paths[-1]
shards = paths[:-1]
sum_coinc = sum(counts[s]["Coincidences"] for s in shards)
sum_singles = sum(counts[s]["Singles"] for s in shards)
m_c = counts[merged]["Coincidences"]
m_s = counts[merged]["Singles"]

print(f"\nshard sum: Singles={sum_singles}  Coincidences={sum_coinc}")
print(f"merged   : Singles={m_s}  Coincidences={m_c}")

ok_s = (sum_singles == 0 and m_s == 0) or abs(m_s - sum_singles) / max(sum_singles, 1) < 0.001
ok_c = (sum_coinc == 0 and m_c == 0) or abs(m_c - sum_coinc) / max(sum_coinc, 1) < 0.001
print(f"singles match (±0.1%): {ok_s}")
print(f"coincs  match (±0.1%): {ok_c}")
sys.exit(0 if (ok_s and ok_c) else 8)
PY
}

# Use the project venv if available (uproot is installed there).
if [[ -f "${PLI_VENV:-$PROJECT_ROOT/venv}/bin/activate" ]]; then
    # shellcheck disable=SC1091
    source "${PLI_VENV:-$PROJECT_ROOT/venv}/bin/activate"
fi

set +e
verify_python "${SHARDS[@]}" "$MERGED" > "$MERGE_CHECK" 2>&1
RC=$?
set -e
cat "$MERGE_CHECK" | tee -a "$RUN_LOG"
if [[ $RC -ne 0 ]]; then
    echo "[$(date -Iseconds)] ERROR: merge_check failed (rc=$RC)" | tee -a "$RUN_LOG"
    exit $RC
fi

# ---------------------------------------------------------------- cleanup shards
if [[ $KEEP_SHARDS -eq 0 ]]; then
    echo "[$(date -Iseconds)] removing per-shard ROOT files (--keep-shards not set)" | tee -a "$RUN_LOG"
    rm -f "${SHARDS[@]}"
fi

echo "[$(date -Iseconds)] DONE. wall=${ELAPSED}s. merged ROOT: $MERGED" | tee -a "$RUN_LOG"
