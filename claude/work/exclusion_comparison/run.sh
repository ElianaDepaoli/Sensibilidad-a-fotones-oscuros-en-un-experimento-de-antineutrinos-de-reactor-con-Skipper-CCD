#!/bin/bash
# Extracts the TEXONO epsilon_95 curves from the user's ROOT macros (mycodes/DP3_Danilov.C and
# Claude's mycodes/DP3_Danilov_new.C) for two flag settings each, copies the Codex exclusion CSVs,
# and draws the combined figure. The macros run on a temporary copy so mycodes/ is not touched.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
top=$(cd "$here/../.." && pwd)
gpt=$(cd "$top/../Proyecto_Final_DarkPhoton_gpt" && pwd)
tmp=$(mktemp -d)
cp "$top"/mycodes/DP3_Danilov.C "$top"/mycodes/DP3_Danilov_new.C \
   "$top"/mycodes/Uranium-photon_cross_section_Nist_gov.txt "$here"/dump_cSk.C "$tmp"/
cd "$tmp"
OFF='DANILOV_A=false;DANILOV_B=false;DANILOV_D=false;'
ON='DANILOV_A=true;DANILOV_B=true;DANILOV_D=true;'
root -l -b -q "dump_cSk.C(\"DP3_Danilov.C\",\"$OFF\",\"orig_none.txt\")"      > log_orig_none.txt 2>&1
root -l -b -q "dump_cSk.C(\"DP3_Danilov.C\",\"$ON\",\"orig_ABD.txt\")"        > log_orig_ABD.txt 2>&1
root -l -b -q "dump_cSk.C(\"DP3_Danilov_new.C\",\"${OFF}save_plots_new=false;\",\"new_none.txt\")" > log_new_none.txt 2>&1
root -l -b -q "dump_cSk.C(\"DP3_Danilov_new.C\",\"${ON}save_plots_new=false;\",\"new_ABD.txt\")"   > log_new_ABD.txt 2>&1
# Grid convergence: the same four runs with wpsize = 200 instead of the macro's 10 (~1-2 min each).
root -l -b -q "dump_cSk.C(\"DP3_Danilov.C\",\"${OFF}wpsize=200;\",\"orig_none_wp200.txt\")"  > log_orig_none_wp200.txt 2>&1 &
root -l -b -q "dump_cSk.C(\"DP3_Danilov.C\",\"${ON}wpsize=200;\",\"orig_ABD_wp200.txt\")"    > log_orig_ABD_wp200.txt 2>&1 &
root -l -b -q "dump_cSk.C(\"DP3_Danilov_new.C\",\"${OFF}save_plots_new=false;wpsize=200;\",\"new_none_wp200.txt\")" > log_new_none_wp200.txt 2>&1 &
root -l -b -q "dump_cSk.C(\"DP3_Danilov_new.C\",\"${ON}save_plots_new=false;wpsize=200;\",\"new_ABD_wp200.txt\")"   > log_new_ABD_wp200.txt 2>&1 &
wait
cp orig_*.txt new_*.txt "$here"/data/
cp "$gpt"/code/park_texono/output/figure2.csv             "$here"/data/codex_park_texono_figure2.csv
cp "$gpt"/code/darkphoton_v3/output/texono_limit.csv      "$here"/data/codex_v3_texono_limit.csv
cp "$gpt"/code/darkphoton_v3/output/danilov_digitized.csv "$here"/data/codex_v3_danilov_digitized.csv
cp "$gpt"/work/reactor_sensitivity/output/limits.csv      "$here"/data/codex_reactor_sensitivity_limits.csv
cp "$gpt"/work/danilov_fresh/output/exclusion.csv         "$here"/data/codex_danilov_fresh_exclusion.csv
(cd "$here"/data && sha256sum *.txt *.csv > sha256.txt)
rm -rf "$tmp"
cd "$here"
root -l -b -q exclusion_comparison.C
