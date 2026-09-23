# Fresh Park production derivation

New ROOT/C++ calculation under `code/park_production_v2/`; no prior implementation or results imported. Both exchange diagrams, exact invariant trace, analytic integral and independent numerical matrices agree. Reference 15 only validates the total. The full Figure 1 is not reproduced; finite-mass disagreement remains. No later corrections.

Reports: `report_park_production_v2_en.pdf`, `report_park_production_v2_es.pdf`, `report_park_production_v2_agent.html`. Task-local provenance: `code/park_production_v2/provenance/`. Links: [[park-production-v2-source]].

## Run receipt

- Date: 2026-09-17T23:54:05Z
- Environment: ROOT 6.36.000, Ubuntu g++ 13.3.0. Task-local Tectonic; no Python or base installs.
- Commands in task directory: `g++ -O2 -std=c++17 production.cpp $(root-config --cflags --libs) -lXMLIO -o production`, `./production > output/run_receipt.txt 2>&1`; similarly build and run `reports.cpp`.
- Actual physics output:

```text
ROOT 6.36.000; compiler 13.3.0
PASS: exact polynomial compact formula, analytic primitive derivative, photon limit, interference equality and crossing symmetry
PASS: Clifford algebra, on-shell kinematics, both Ward identities, physical polarizations and independent matrix traces
PASS: derived analytic totals, Reference 15 totals, Klein-Nishina totals, quadrature, gamma tail and printed photon count
reference15_max_relative_error=7.06101843662e-14; KN_max_relative_error=7.20534742982e-14; spectrum_nodes=1.80078174594e-13
PASS: independent source integration order; maximum relative difference=1.1990408666e-14
Extracted 907 horizontal PDF segments from 3 spectrum curves
Figure 1 M=0.1 available median residual=0.0543079495417; lower cutoff=0.0256612042407
Figure 1 M=0.1 central median residual=0.0414967524505; lower cutoff=0.0414967524505
Figure 1 M=0.5 available median residual=0.152838280092; lower cutoff=0.152838280092
Figure 1 M=0.5 central median residual=0.115551202919; lower cutoff=0.115551202919
Figure 1 M=1 available median residual=0.370145090123; lower cutoff=0.370145090123
Figure 1 M=1 central median residual=0.37854708634; lower cutoff=0.37854708634
Info in <TCanvas::Print>: pdf file output/figure1_comparison.pdf has been created
Info in <TCanvas::Print>: SVG file output/figure1_comparison.svg has been created
Info in <TCanvas::Print>: pdf file output/cross_section.pdf has been created
Info in <TCanvas::Print>: SVG file output/cross_section.svg has been created
Info in <TCanvas::Print>: pdf file output/compton_validity.pdf has been created
Info in <TCanvas::Print>: SVG file output/compton_validity.svg has been created
STATUS: source-equation computation completed; full Figure 1 reproduction not established.
```

Derived: trace, lab support, differential and total production. Supplied: QED rules, Park spectrum, constants, Reference 15 validation formula, NIST material data. Not checked: original source measurement, transport, atomic effects, loops, author's generator or event statistics.

Additional receipt: `output/validation_receipt.txt`

```text
PASS: report links, figure files and task-local provenance references
```

Additional receipt: `output/latex_en_receipt.txt`

```text
note: using only cached resource files
note: Running TeX ...
warning: code/park_production_v2/lecture.tex:348: Underfull \hbox (badness 2443) in paragraph at lines 347--348
warning: code/park_production_v2/lecture.tex:348: Underfull \hbox (badness 3492) in paragraph at lines 347--348
warning: code/park_production_v2/lecture.tex:349: Underfull \hbox (badness 3449) in paragraph at lines 348--349
warning: code/park_production_v2/lecture.tex:349: Underfull \hbox (badness 1043) in paragraph at lines 348--349
note: Rerunning TeX because "report_park_production_v2_en.aux" changed ...
warning: code/park_production_v2/lecture.tex:348: Underfull \hbox (badness 2443) in paragraph at lines 347--348
warning: code/park_production_v2/lecture.tex:348: Underfull \hbox (badness 3492) in paragraph at lines 347--348
warning: code/park_production_v2/lecture.tex:349: Underfull \hbox (badness 3449) in paragraph at lines 348--349
warning: code/park_production_v2/lecture.tex:349: Underfull \hbox (badness 1043) in paragraph at lines 348--349
warning: warnings were issued by the TeX engine; use --print and/or --keep-logs for details.
note: Running xdvipdfmx ...
note: Writing `report_park_production_v2_en.log` (13.5810546875 KiB)
note: Writing `report_park_production_v2_en.pdf` (156.134765625 KiB)
note: Skipped writing 2 intermediate files (use --keep-intermediates to keep them)
```

Additional receipt: `output/latex_es_receipt.txt`

```text
note: using only cached resource files
note: Running TeX ...
warning: code/park_production_v2/lecture.tex:353: Underfull \hbox (badness 2443) in paragraph at lines 352--353
warning: code/park_production_v2/lecture.tex:353: Underfull \hbox (badness 3895) in paragraph at lines 352--353
warning: code/park_production_v2/lecture.tex:354: Underfull \hbox (badness 3449) in paragraph at lines 353--354
warning: code/park_production_v2/lecture.tex:354: Underfull \hbox (badness 1163) in paragraph at lines 353--354
warning: code/park_production_v2/lecture.tex:355: Underfull \hbox (badness 2903) in paragraph at lines 354--355
note: Rerunning TeX because "report_park_production_v2_es.out" changed ...
warning: code/park_production_v2/lecture.tex:353: Underfull \hbox (badness 2443) in paragraph at lines 352--353
warning: code/park_production_v2/lecture.tex:353: Underfull \hbox (badness 3895) in paragraph at lines 352--353
warning: code/park_production_v2/lecture.tex:354: Underfull \hbox (badness 3449) in paragraph at lines 353--354
warning: code/park_production_v2/lecture.tex:354: Underfull \hbox (badness 1163) in paragraph at lines 353--354
warning: code/park_production_v2/lecture.tex:355: Underfull \hbox (badness 2903) in paragraph at lines 354--355
warning: warnings were issued by the TeX engine; use --print and/or --keep-logs for details.
note: Running xdvipdfmx ...
note: Writing `report_park_production_v2_es.log` (13.7958984375 KiB)
note: Writing `report_park_production_v2_es.pdf` (156.7626953125 KiB)
note: Skipped writing 2 intermediate files (use --keep-intermediates to keep them)
```
