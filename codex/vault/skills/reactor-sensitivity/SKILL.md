---
name: reactor-sensitivity
description: Execute the reactor TEXONO sensitivity prompt using the checked ROOT scattering kernels, with decay and medium-model validity recorded separately from figure reproduction.
---

Read `prompts/reproducesensitivitytodarkphotoninreactors.txt` and the vault reproduction receipts. Use ROOT 6.36.000 and C++ for physical calculations; marimo presents their exported figures. Keep the historical Park convention separate from a one-sided selected-event cap.

Reuse `code/darkphoton_v3/physics_base.hxx` only after rerunning its matrix, Ward, detailed-balance and quadrature checks, plus the SymPy derivation in `code/park_texono/derive.py`. Never fit a normalization to declare reproduction.

Izaguirre Eq. (9) provides survival/decay geometry; its electron-pair channel is closed below the mass range endpoint in this prompt. Park Eq. (4) supplies a low-mass three-photon width, not an exact near-threshold loop result. Label that approximation. Do not multiply the Compton source by the alternative Danilov oscillation source. Reject the unmodeled material resonance.

Every figure in `papers_init` needs an explicit reproduction status, including third-party curves for which the underlying data are absent. A PDF extraction is a reference, not a physics reproduction. Write bilingual LaTeX/PDF lectures, a self-contained follow-up HTML, provenance and a vault receipt. Report unavailable token/effort telemetry as unavailable.
