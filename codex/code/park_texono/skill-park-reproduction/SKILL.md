---
name: park-texono-reproduction
description: Reproduce the TEXONO part of Park's reactor dark-photon paper with ROOT C++, symbolic cross-section checks, original-figure comparisons, and explicit accounting of failed paper matches.
---

# Park TEXONO reproduction

Use `../texono.conf` for experiment inputs and `../run.sh` for the numerical run.
Keep Park's production normalization, unpolarized detection average, point source,
count-limit convention, and ideal containment explicit. The requested scope excludes
NEOS and corrections from Danilov or later work.

Derive the two tree amplitudes before consulting reference formulas. Use Reference 15
only as an independent total-production check and Reference 20 as an amplitude and
phase-space check. Do not fit a normalization to recover a paper's reported limit.

Distinguish the low-energy Thomson approximation, the relativistic free-electron
Klein-Nishina cross section, and a material's total photon interaction cross section.
An NIST comparison diagnoses Park's approximation but must not change the baseline
rate unless the user extends the scope.

For figure comparisons, extract coordinates from the supplied PDF and compare on the
same axes. State when a figure or calculation fails reproduction. Preserve failures
alongside passing symbolic and numerical checks.

Use the existing `/home/eliana/.marimo/bin/python` environment for the presentation
notebook. ROOT C++ performs the numerical physics; the notebook displays its outputs.
Write concise English and Spanish HTML reports, plus provenance and a vault receipt.
Report runtime from timestamps. State token usage and reasoning effort as unavailable
when the runtime does not expose those fields; never estimate them as measurements.

Continuation diagnostics are in `../continuation.C` and `../continuation_report.py`.
Run the ROOT macro from its directory to reuse existing functions without rebuilding
the baseline executable. Conditional production normalization is an implementation
hypothesis, not Park's stated source equation. Compare horizontal PDF segments as
bins, preserve the baseline, and keep continuation timing and receipts separate.
