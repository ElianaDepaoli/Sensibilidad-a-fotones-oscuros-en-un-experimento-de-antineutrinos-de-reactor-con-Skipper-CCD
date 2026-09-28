#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.."
g++ -O2 -std=c++17 work/reactor_sensitivity/sensitivity.cpp $(root-config --cflags --libs) -o work/reactor_sensitivity/sensitivity
cd work/reactor_sensitivity
./sensitivity config.conf > output/run_receipt.txt 2>&1
cat output/run_receipt.txt
cd ../..
python3 work/reactor_sensitivity/build_reports.py
for language in en es; do
  TECTONIC_CACHE_DIR=code/darkphoton_v3/tex-cache XDG_CACHE_HOME=code/darkphoton_v3/tex-cache code/park_production_v2/tools/tectonic --only-cached --keep-logs "report_reactor_sensitivity_${language}.tex" > "work/reactor_sensitivity/output/latex_${language}.txt" 2>&1
done
PYTHONDONTWRITEBYTECODE=1 /home/eliana/.marimo/bin/marimo export html work/reactor_sensitivity/figures.py -o work/reactor_sensitivity/output/figures.html
