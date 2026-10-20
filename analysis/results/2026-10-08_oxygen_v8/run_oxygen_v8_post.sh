#!/bin/bash
# Post-GATE chain of run_oxygen_v8.sh, rerun in the environment of the Phase 3s-2 orchestrator (env_with: venv python,
# PYTHONPATH=/usr/lib/python3.14/site-packages, no project ROOT on the path). The first attempt (01:57) sourced
# setup_environment_st.sh as well, which loads a second ROOT/cppyy and segfaulted in PyROOT before writing triples.
set -u
R=${PLI_REPO:-$(cd "$(dirname "$0")/../../.." && pwd)}
AN=$(cd "$(dirname "$0")/../.." && pwd)
L=${PLI_RES:-$AN/results}/2026-10-08_oxygen_v8
OUT=${PLI_RAW:-<missing>}/oxygen_v8
[ -d "$OUT/sim" ] || { echo "$OUT/sim: raw data not included in the release, see analysis/README.md" >&2; exit 1; }
A=$AN/pli-quadra-archive-2026-08-12/2026-05-06_2258_phase3s_2_paired_reruns
say(){ echo "[$(date +%F' '%T)] $*" >> "$L/run.log"; }
PY="env -i HOME=$HOME PATH=/usr/local/bin:/usr/bin:/bin PYTHONPATH=/usr/lib/python3.14/site-packages PLI_RAW=${PLI_RAW:-} PARALLELPROJ_C_LIB=/usr/local/lib/libparallelproj_c.so PARALLELPROJ_CUDA_LIB=/usr/local/lib/libparallelproj_cuda.so $R/venv/bin/python"
cd "$R"
say "post-GATE chain restarted in the orchestrator environment"
$PY "$AN/reconstruction/extract_triples_chunked.py" --input-format root --root-dir "$OUT/sim" --out-dir "$OUT/extracted" \
    --tag _hypoxia_te_v8 > "$L/extract2.log" 2>&1 || { say "ABORT extract (2)"; echo OXY_FAILED >> "$L/run.log"; exit 1; }
say "extracted"
$PY "$AN/reconstruction/filter_triples_tof_pann.py" --triples-detector "$OUT/extracted/triples_detector_hypoxia_te_v8.h5" \
    --triples-truth "$OUT/extracted/triples_truth_hypoxia_te_v8.h5" \
    --output-detector "$OUT/triples_detector_hypoxia_te_v8_rcut075mm.h5" \
    --output-truth "$OUT/triples_truth_hypoxia_te_v8_rcut075mm.h5" --r-cut 75.0 > "$L/filter.log" 2>&1 \
    || { say "ABORT filter"; echo OXY_FAILED >> "$L/run.log"; exit 1; }
say "filtered"
$PY "$A/analyzers/per_insert_hypoxia.py" --detector "$OUT/triples_detector_hypoxia_te_v8_rcut075mm.h5" \
    --truth "$OUT/triples_truth_hypoxia_te_v8_rcut075mm.h5" --tau-min 1.0 --materials-json "$R/materials/pli_materials.json" \
    --out-md "$L/hypoxia_te_v8_per_insert.md" --out-csv "$L/hypoxia_te_v8_per_insert.csv" \
    --label "review round 8 oxygen rerun (database v8, 250M)" > "$L/per_insert.log" 2>&1 \
    || { say "ABORT per-insert"; echo OXY_FAILED >> "$L/run.log"; exit 1; }
say "per-insert done"
$PY "$L/oxygen_boot_v8.py" > "$L/oxygen_boot_v8.log" 2>&1 || { say "ABORT boot"; echo OXY_FAILED >> "$L/run.log"; exit 1; }
say "done"; echo OXY_DONE >> "$L/run.log"
