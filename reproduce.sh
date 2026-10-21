#!/usr/bin/env bash
# Regenerate the results. Default: the fast part (about 15 minutes on one M1 core) and the tests.
# FULL=1 also reruns the long exact searches (K = 12..14, hours) and the phase diagram (about 1 hour on 3 cores).
set -euo pipefail
cd "$(dirname "$0")"
PY="${PYTHON:-python}"
OUT=results/exact_search.jsonl
"$PY" scripts/find_fan_tunings.py                       # results/fan_tunings.json
for key in polyak100:3:8 ex32:3:6 convgap:3:11 prop13:3:11 fan7:3:11 fan9:3:11 p1:10:10; do
  IFS=: read -r t a b <<< "$key"
  "$PY" scripts/run_exact.py "$t" "$a" "$b" --out "$OUT"
done
if [[ "${FULL:-0}" == "1" ]]; then
  "$PY" scripts/run_exact.py convgap 12 12 --workers 4 --out "$OUT"
  "$PY" scripts/run_exact.py prop13 12 14 --workers 4 --out "$OUT"
  "$PY" scripts/run_exact.py fan7 12 14 --workers 4 --out "$OUT"
  "$PY" scripts/run_exact.py fan9 12 14 --workers 4 --out "$OUT"
  "$PY" scripts/run_exact.py p1 13 13 --workers 4 --out "$OUT"
  "$PY" scripts/phase_diagram.py --nrho 60 --normalized --workers 3   # results/phase_diagram_k100.json
fi
"$PY" scripts/plot_phase_diagram.py                     # figures/phase_diagram_1d_k100.png
"$PY" scripts/summarize.py                              # results/summary.{json,md}
"$PY" -m pytest -q
