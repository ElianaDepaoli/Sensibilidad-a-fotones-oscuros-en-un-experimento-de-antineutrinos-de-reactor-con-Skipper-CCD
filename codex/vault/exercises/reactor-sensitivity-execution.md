# Reactor sensitivity prompt execution

Executed `prompts/reproducesensitivitytodarkphotoninreactors.txt`. New ROOT/C++ survival and event-boundary calculation reuses the verified massive scattering kernels. English/Spanish LaTeX lectures and self-contained HTML: `report_reactor_sensitivity*`. Session skill: [[../skills/reactor-sensitivity/SKILL]]. Related: [[darkphoton-v3]], [[park-texono-reproduction]].

Park Figure 1 and the printed TEXONO limit remain unreproduced. The full Figure 2 compilation is not reproduced. The three-photon width is borrowed, and is not an exact near-threshold calculation. No new microscopic amplitude derived. Detector response and raw likelihood remain unchecked.

## Run receipt

Date: 2026-09-25T17:32:57.969925+00:00
Command: `bash work/reactor_sensitivity/run.sh` (compiles ROOT/C++, runs calculation, assembles and compiles both lectures, exports marimo). Environment: existing ROOT 6.36.000, g++ 13.3.0, Python 3.12, project Tectonic cache and existing marimo installation; no installation.

```text
PASS: decay scaling, zero-distance/mixing survival, analytic and missing roots, doubled vacuum quadrature, medium zero-decay comparison
VISIBLE medium survival: 12 points
VISIBLE medium survival: 23 points
Info in <TCanvas::Print>: pdf file output/conditional_limits.pdf has been created
Info in <TCanvas::Print>: SVG file output/conditional_limits.svg has been created
Info in <TCanvas::Print>: png file output/conditional_limits.png has been created
maximum_EH_optical_depth_at_vacuum_limit = 1.09746078459e-11
PASS: every reported boundary closes the configured event cap
UNCHECKED: exact three-photon loop near threshold; detector response; full Park Figure 2 compilation; published normalization discrepancy
ROOT 6.36.000; compiler 13.3.0
compute_seconds = 0.893734288
```

Independent earlier reruns in this execution: `bash code/park_texono/run.sh` and `make -C code/darkphoton_v3 reports`; actual outputs are embedded in the HTML. Provenance: `provenance/numbers.json`, `provenance/claims.yaml`, and the hashed output manifest.

## Export and visual review

The first run compiled both PDFs, then marimo failed because the sandbox blocked its local socket. With explicit escalation approval, reran:

`PYTHONDONTWRITEBYTECODE=1 /home/eliana/.marimo/bin/marimo export html work/reactor_sensitivity/figures.py -o work/reactor_sensitivity/output/figures.html`

Actual output: empty stdout/stderr; exit code 0. Both PDFs were rendered with `pdftoppm -scale-to 850 -png` and every page was visually inspected. No missing plots or clipped equations observed. TeX logs contain underfull boxes, but no overfull boxes or undefined references.
