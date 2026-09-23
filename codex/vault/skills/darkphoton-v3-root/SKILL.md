---
name: darkphoton-v3-root
description: ROOT/C++ production, medium photon-hidden-photon oscillations, TEXONO detection and exclusion limits.
---

Use this skill for the v3 task from `prompts/darkphoton_production_and_detection_3.txt`.
All calculations, plots, extraction, validation, report assembly and provenance
generation must be ROOT/C++. Do not add Python or marimo utilities.

Use Danilov, Demidov and Gorbunov (2019) for the medium oscillation formulas:
reactor and detector plasma masses, absorption width, production probability,
detection probability, and the TEXONO 95% event cap. Keep the tree-level
Compton cross section from the reviewed v2 implementation as a checked input;
Reference 15 is an independent total-production validation, not the derivation.

Use coefficient 2 in Danilov Eq. (5), as printed in the supplied PDF. These are
leading-order formulas, not exact all-orders probabilities. Do not use detector
Eq. (9) at resonance; it neglects absorption there. Bound the solver to small
mixing and mark absent roots explicitly. Include the assumed selection efficiency.

Keep ROOT graphs alive through canvas export. Assert finite in-frame points for
each curve and check pad ownership before saving. Render and visually inspect
every final plot; a nonzero PDF file size is not evidence of a nonempty plot.
Compare actual source curves with predictions and report unresolved discrepancies.
