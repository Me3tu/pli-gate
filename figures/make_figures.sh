#!/usr/bin/env bash
# Regenerate all 14 figures of the paper (6 main, 8 supplementary) from the shipped source_data/.
#
#   ./make_figures.sh [OUTDIR]        default OUTDIR: ./out
#   PYTHON=/path/to/python3 ./make_figures.sh
#
# Writes <stem>.pdf and <stem>.png for every figure, fig2_closure_benchmark.json, fig1_methods.pdf,
# and one log per generator in OUTDIR/logs/. Inputs are always taken from source_data/ next to this
# script (PLI_RES and PLI_FIGOUT are set here, PLI_INFLIGHT_JSON is unset). See README.md.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="${1:-$HERE/out}"
mkdir -p "$OUT/logs"
OUT="$(cd "$OUT" && pwd)"
PY="${PYTHON:-python3}"

export PLI_RES="$HERE/source_data"
export PLI_FIGOUT="$OUT"
unset PLI_INFLIGHT_JSON
export MPLBACKEND=Agg
export PYTHONDONTWRITEBYTECODE=1

echo "python: $("$PY" -c 'import sys, numpy, scipy, matplotlib; print(sys.version.split()[0], "numpy", numpy.__version__, "scipy", scipy.__version__, "matplotlib", matplotlib.__version__)')"
echo "inputs: $PLI_RES"
echo "output: $OUT"

# paper figure order. f08_pedestal.py writes Fig. 4 and Fig. S3
SCRIPTS=(
  f02b_closure_benchmark.py   # Fig. 2   fig2_closure_benchmark
  f09_mechanism.py            # Fig. 3   fig6_clinical_mechanism_v2
  f08_pedestal.py             # Fig. 4   fig9_pedestal_derivation, Fig. S3 figS2_pedestal_residuals (about 1 min)
  f11_realistic.py            # Fig. 5   fig11_realistic_forest
  f10_rate.py                 # Fig. 6   fig10_rate_law
  fS_benchspec.py             # Fig. S1  figS_benchspec
  f04_transport.py            # Fig. S2  fig_r2_transport_novelty
  f07_oxygen.py               # Fig. S4  fig5_sensitivity
  f06_brain.py                # Fig. S5  fig4_brain_tumor
  f05_vipman.py               # Fig. S6  fig3_vipman_axial
  fS1_refutation.py           # Fig. S7  figS3_eq2_refutation
  fS3_inflight.py             # Fig. S8  figS4_inflight
)

cd "$HERE"
for s in "${SCRIPTS[@]}"; do
  echo "== $s"
  if ! "$PY" "$s" > "$OUT/logs/${s%.py}.log" 2>&1; then
    echo "FAILED: $s (see $OUT/logs/${s%.py}.log)" >&2
    tail -n 20 "$OUT/logs/${s%.py}.log" >&2
    exit 1
  fi
done

# Fig. 1: standalone TikZ, compiled in a temporary directory
echo "== fig1_methods.tex"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
cp fig1_methods.tex "$TMP/"
if ! ( cd "$TMP" && pdflatex -interaction=nonstopmode -halt-on-error fig1_methods.tex > fig1_methods.pdflatex.log 2>&1 ); then
  cp "$TMP/fig1_methods.pdflatex.log" "$OUT/logs/" || true
  echo "FAILED: fig1_methods.tex (see $OUT/logs/fig1_methods.pdflatex.log)" >&2
  exit 1
fi
cp "$TMP/fig1_methods.pdf" "$OUT/"
cp "$TMP/fig1_methods.pdflatex.log" "$OUT/logs/"
if grep -q 'stix.sty' "$TMP/fig1_methods.log"; then
  echo "   note: fig1_methods compiled with stix.sty; the paper's Fig. 1 was compiled without it (mathptmx fallback), see README.md"
fi

echo "done: $(ls "$OUT"/*.pdf | wc -l) PDFs in $OUT"
