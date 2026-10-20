#!/bin/bash
# B3: extract the four Mercolli tubes from the MERGED file, both matching rules.
set -u
R=${PLI_REPO:-$(cd "$(dirname "$0")/../../.." && pwd)}
L=${PLI_RES:-$(cd "$(dirname "$0")/.." && pwd)}/2026-08-23_B3
M=${PLI_RAW:-<missing>}/g1_exposure/singles_g1_exposure.h5
[ -f "$M" ] || { echo "$M: raw data not included in the release, see analysis/README.md" >&2; exit 1; }
cd "$R"; source venv/bin/activate
say(){ echo "[$(date +%H:%M:%S)] $*" >> "$L/extract.log"; }
say "=== B3 extraction START from the merged file ($(du -h $M | cut -f1)) ==="
for rule in closest all; do
  say "rule=$rule"
  python "$L/../2026-08-12_steinberger_allprompts/steinberger_extract_ap.py" \
      --singles "$M" --config "$L/voi_mercolli.json" --phantom mercolli_tubes \
      --out "$L/out_$rule" --matching "$rule" --sigma-t-ps 90.9 \
      >> "$L/extract_$rule.log" 2>&1
  say "  rule=$rule rc=$?"
done
echo B3_EXTRACT_DONE >> "$L/extract.log"
