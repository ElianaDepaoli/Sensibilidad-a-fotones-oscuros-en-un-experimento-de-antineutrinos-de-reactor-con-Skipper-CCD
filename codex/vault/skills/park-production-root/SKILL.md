---
name: park-production-root
description: Derive and check reactor Compton-like dark-photon production in ROOT C++, with bilingual LaTeX lectures and explicit figure-reproduction status.
---

Use `prompts/park_implementation_2.txt` as the task scope. All new computational,
extraction, validation and reporting utilities must be C++/ROOT. Do not add Python
or marimo code. The ROOT-only instruction overrides the general notebook convention.

For this fresh task use `code/park_production_v2/` and source documents in
`papers/park_production_v2/`; do not import prior implementations or outputs.
Derive both electron-exchange production diagrams, state spin averages and units,
and verify Ward identities, the photon limit and numerical convergence.
Reference 15 may independently check the integrated production result; identify
that borrowed formula separately from the derivation.

Keep Park's stated source equation distinct from a conditional energy distribution.
Never call a visually similar curve a reproduction without comparing actual PDF
coordinates and accounting for its input normalization. Do not apply later corrections.

Write English and Spanish lectures as LaTeX compiled to PDF, and an HTML handoff.
Separate Thomson validity from free-electron Klein–Nishina and from total material
attenuation. Record untested assumptions, actual command output, provenance and
known runtime metadata; mark unavailable token/effort fields as unavailable.
