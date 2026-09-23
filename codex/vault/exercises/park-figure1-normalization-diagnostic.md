# Park Figure 1 normalization diagnostic

Continued [[park-texono-reproduction]] using its existing C++ implementation. Conditional production normalization plus a lower gamma cutoff closely follows the retained Figure 1 bins. This changes the published rate equation and remains a hypothesis about the author's implementation. It does not resolve the TEXONO limit. No later correction was applied.

Outputs: `report_park_texono_continuation.html`, `code/park_texono/output/continuation/summary.json`. Derived here: diagnostic normalization and bin averaging. Reused: existing amplitudes and kinematics. Supplied: [[park-2017-reactor-dark-photons]] and Reference 15 inputs.

## Run receipt

- Report date: 2026-09-17T22:43:42.210314+00:00
- Command from `code/park_texono/`: `root -l -b -q continuation.C > output/continuation/run_receipt.txt 2>&1`, then `python3 continuation_report.py`.
- Environment: existing ROOT/Cling, Python with PyYAML; no package installs. Baseline source and executable retained.
- Actual ROOT output:

```text
Processing continuation.C...
Continuation: ROOT 6.36.000; existing source included unchanged
PASS: finite-mass trace integrals match Reference 15; production KN and detection 2/3 KN limits; finite-mass detailed balance
PASS finite-mass rate M=0.001 nodes=2.109423747e-13 order=2.869926519e-13
PASS finite-mass rate M=0.1 nodes=2.117195308e-13 order=8.468337143e-12
PASS finite-mass rate M=0.5 nodes=2.119415754e-13 order=1.665845218e-07
PASS finite-mass rate M=1 nodes=2.114974862e-13 order=3.562469401e-06
PASS: conditional massless limit; small-mass normalization diagnostic leaves TEXONO mismatch
Bin residual M=0.1 baseline=0.0420647644 lower-cutoff=0.0420647644 conditional=0.04242800484 conditional+lower-cutoff=0.04242800484
Bin residual M=0.5 baseline=0.1161444674 lower-cutoff=0.1161444674 conditional=0.04660430614 conditional+lower-cutoff=0.04660430614
Bin residual M=1 baseline=0.3787300449 lower-cutoff=0.3787300449 conditional=0.0107603018 conditional+lower-cutoff=0.0107603018
Info in <TCanvas::Print>: SVG file output/continuation/figure1_diagnostic.svg has been created
Info in <TCanvas::Print>: pdf file output/continuation/figure1_diagnostic.pdf has been created
```

Not checked: author's generator, reactor transport, detector response, raw-event likelihood, prompt-fission spectrum, decay loop. No full-figure or printed-limit reproduction claim. Original run receipt remains in [[park-texono-reproduction]].

Post-run verification: `python3 code/park_texono/verify_outputs.py`; marimo HTML export completed successfully and contains the continuation panel. Baseline hashes matched the continuation snapshot.

```text
PASS: report links and embedded figures; provenance references; exported notebook exists
PASS: numerical-convergence checks; failed paper-limit comparison remains explicitly recorded
UNCOVERED: reactor transport, detector response, raw-event likelihood, decay loop, prompt-fission spectrum
```
