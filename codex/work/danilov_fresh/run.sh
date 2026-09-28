#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
python3 verify_algebra.py > output/algebra_receipt.txt 2>&1
g++ -O2 -std=c++17 calculate.cpp $(root-config --cflags --libs) -o calculate
./calculate experiment.conf > output/run_receipt.txt 2>&1
cat output/algebra_receipt.txt output/run_receipt.txt
cd ../..
python3 work/danilov_fresh/build_reports.py
for lang in en es; do
  pdflatex -interaction=nonstopmode -halt-on-error "report_danilov_fresh_${lang}.tex" > "work/danilov_fresh/output/latex_${lang}.txt" 2>&1
  pdflatex -interaction=nonstopmode -halt-on-error "report_danilov_fresh_${lang}.tex" >> "work/danilov_fresh/output/latex_${lang}.txt" 2>&1
done
PYTHONDONTWRITEBYTECODE=1 /home/eliana/.marimo/bin/marimo export html work/danilov_fresh/figures.py -o work/danilov_fresh/output/figures.html
