# Park Figure 1: four-column plot data

Output: `work/park_figure1_digitization/park_figure1.txt`.
Columns: energy (MeV), then spectra for masses 0.1, 0.5, 1.0 MeV, in the plotted units of 10^21 MeV^-1 s^-1. Use a logarithmic ordinate and steps-mid. Comment lines describe the units; `nan` denotes absent or clipped strokes, not zero.

Read directly from vector paths on page 2 of `papers_init/park2017.pdf`. No physics calculation performed. New work: horizontal-segment extraction and shared bin-center sampling. Reused: axis calibration from `code/park_texono/inspect_sources.py`, checked visually against the page axes. Libraries: Poppler, NumPy, Python XML parser. Approximate plot-derived data, not the author's simulation output. No smoothing or extrapolation.

Related: [[park-2017-reactor-dark-photons]], [[park-figure1-normalization-diagnostic]]. Provenance: `provenance/numbers.json` entry `park-figure1-plot-digitization`; source hash and calibration in `work/park_figure1_digitization/receipt.json`.

## Run receipt

- Date: 2026-09-24
- Command: `python3 work/park_figure1_digitization/extract.py`
- Environment: Linux, Python 3.12.3, NumPy 2.5.2, system pdftocairo; no installation.
- Actual output:

```text
PASS: extracted three plotted vector paths; no physics calculation.
PASS: text readback has 430 rows and exactly 4 columns; energies increase; visible spectra are positive.
Visible samples by mass (MeV): {"1.0": 255, "0.1": 349, "0.5": 315}
UNCHECKED: underlying simulation and values clipped outside the published axes.
```

No new figure produced. The exported table was checked by readback; a rendered replot and digitization uncertainty were not independently validated.


## ROOT plot

Requested ROOT output: `work/park_figure1_digitization/park_figure1_root.png` and `.pdf`; macro: `plot_park_figure1.C` in the same folder. No physics derived; table values reused. ROOT supplies rendering. Steps preserve the sampled bins; nan values break the paths. Visually inspected the PNG. No independent pixel-level comparison or uncertainty validation.

Run date: 2026-09-24. Environment: existing ROOT 6.36.000, Linux; no installs.
Command: `root -l -b -q work/park_figure1_digitization/plot_park_figure1.C > work/park_figure1_digitization/root_run_receipt.txt 2>&1`.

Actual output:

```text
Processing work/park_figure1_digitization/plot_park_figure1.C...
PASS: spectrum column 2: 349 visible bins plotted
PASS: spectrum column 3: 315 visible bins plotted
PASS: spectrum column 4: 255 visible bins plotted
Info in <TCanvas::Print>: png file work/park_figure1_digitization/park_figure1_root.png has been created
Info in <TCanvas::Print>: pdf file work/park_figure1_digitization/park_figure1_root.pdf has been created
ROOT 6.36.000: plotted text data only; missing bins remain gaps.
```

Initial attempt failed because stream reads inside assertions were skipped with disabled assertions. Fixed by using unconditional reads and explicit runtime checks; the successful run above is after that correction.

## Macro path correction (2026-09-24)

Input and output paths now resolve relative to the macro via `std::filesystem::absolute(__FILE__).parent_path()`. Previously they depended on launching ROOT from the project root. Existing ROOT environment; no physical calculation changed. Verified interpreted execution from both directories; compiled ACLiC execution was not checked.

From the macro directory: `root -l -b -q plot_park_figure1.C`

```text
Processing plot_park_figure1.C...
PASS: spectrum column 2: 349 visible bins plotted
PASS: spectrum column 3: 315 visible bins plotted
PASS: spectrum column 4: 255 visible bins plotted
Info in <TCanvas::Print>: png file /home/eliana/Documentos/Materias/IA_OG/GW-AI-course/Proyecto_Final_DarkPhoton_gpt/work/park_figure1_digitization/park_figure1_root.png has been created
Info in <TCanvas::Print>: pdf file /home/eliana/Documentos/Materias/IA_OG/GW-AI-course/Proyecto_Final_DarkPhoton_gpt/work/park_figure1_digitization/park_figure1_root.pdf has been created
ROOT 6.36.000: plotted text data only; missing bins remain gaps.
```

From the project root: `root -l -b -q work/park_figure1_digitization/plot_park_figure1.C`

```text
Processing work/park_figure1_digitization/plot_park_figure1.C...
PASS: spectrum column 2: 349 visible bins plotted
PASS: spectrum column 3: 315 visible bins plotted
PASS: spectrum column 4: 255 visible bins plotted
Info in <TCanvas::Print>: png file /home/eliana/Documentos/Materias/IA_OG/GW-AI-course/Proyecto_Final_DarkPhoton_gpt/work/park_figure1_digitization/park_figure1_root.png has been created
Info in <TCanvas::Print>: pdf file /home/eliana/Documentos/Materias/IA_OG/GW-AI-course/Proyecto_Final_DarkPhoton_gpt/work/park_figure1_digitization/park_figure1_root.pdf has been created
ROOT 6.36.000: plotted text data only; missing bins remain gaps.
```
