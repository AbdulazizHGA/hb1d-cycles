#!/usr/bin/env bash
# Regenerate the results. Default part (about 22 minutes on an M1): exact searches with K <= 11, the rigorous
# export and its interval check, the Farkas-certificate check for K <= 12, the figure, the summary and the tests.
# FULL=1 also reruns the long exact searches (K = 12..14, a few hours on 4 cores), the phase diagram, and checks
# the K = 13, 14 certificates (1 to 2 hours more).
# Each run REPLACES the earlier record for the same (tuning, K) in results/exact_search.jsonl and rewrites
# certificates/<tuning>_K<K>.farkas.gz.
set -euo pipefail
cd "$(dirname "$0")"
# Interpreter: $PYTHON if set, else the repo's .venv (see README), else python3.
if [[ -n "${PYTHON:-}" ]]; then PY="$PYTHON"
elif [[ -x .venv/bin/python ]]; then PY=.venv/bin/python
else PY=python3; fi
echo "reproduce.sh: using $("$PY" -c 'import sys; print(sys.executable, sys.version.split()[0])')"
"$PY" scripts/find_fan_tunings.py                       # results/fan_tunings.json
for key in polyak100:3:8 ex32:3:6 convgap:3:11 prop13:3:11 fan7:3:11 fan9:3:11 p1:10:10; do
  IFS=: read -r t a b <<< "$key"
  "$PY" scripts/run_exact.py "$t" "$a" "$b"
done
if [[ "${FULL:-0}" == "1" ]]; then
  "$PY" scripts/run_exact.py convgap 12 12 --workers 4
  "$PY" scripts/run_exact.py prop13 12 14 --workers 4
  "$PY" scripts/run_exact.py fan7 12 14 --workers 4
  "$PY" scripts/run_exact.py fan9 12 14 --workers 4
  "$PY" scripts/run_exact.py p1 13 13 --workers 4
  "$PY" scripts/phase_diagram.py --nrho 60 --normalized --workers 3   # results/phase_diagram_k100.json
fi
"$PY" scripts/export_rigorous.py                        # certificates/mixed_and_omega_aug.json
"$PY" scripts/verify_rigorous.py                        # interval check of that file (independent of hb1d)
if [[ "${FULL:-0}" == "1" ]]; then                     # Farkas certificates (independent of hb1d)
  "$PY" scripts/verify_certificates.py certificates/*.farkas.gz certificates/grid/*.farkas.gz
else
  "$PY" scripts/verify_certificates.py --max-k 12 certificates/*.farkas.gz certificates/grid/*.farkas.gz
fi
"$PY" scripts/plot_phase_diagram.py                     # figures/phase_diagram_1d_k100.png
"$PY" scripts/summarize.py                              # results/summary.{json,md}
"$PY" -m pytest -q
