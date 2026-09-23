#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
mkdir -p output
pdftocairo -f 2 -l 2 -svg ../../papers_init/park2017.pdf output/park_page2.svg
python3 inspect_sources.py | tee output/digitization_receipt.txt
python3 derive.py | tee output/symbolic_receipt.txt
g++ -O2 -std=c++17 park_texono.cpp $(root-config --cflags --libs) -o park_texono
./park_texono "${1:-texono.conf}" 2>&1 | tee output/run_receipt.txt
