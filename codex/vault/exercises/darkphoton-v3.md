# Dark-photon production and detection v3

## Redo audit — 2026-09-21

The original completion claim below is superseded. The production plot lost its
loop-local ROOT graphs before export. Its per-epsilon-squared normalization used
epsilon=1 in nonlinear formulas. The splitting coefficient was 4 instead of the
2 printed in Danilov Eq. (5); the solver also reported large-mixing roots and used
an invalid detector formula at resonance. These are corrected in `redo.cpp`,
called by the canonical `v3.cpp` entry point.

Fresh plots compare the recalculated Park source with the prior PDF extraction,
and compare TEXONO calculations with Park's quoted bound and a fresh vector
extraction of Danilov's published curve. Full numerical reproduction remains
unachieved. Selection efficiency is explicit and configurable; the medium result
is conditional on full containment, free electrons, and homogeneous material.

Run receipt: `make reports` in `code/darkphoton_v3`; ROOT 6.36.000, g++ 13.3.0,
local Tectonic cache. Actual physics output:

```
PASS: symbolic trace, independent matrices, Ward identities, Reference15, analytic totals, production/inverse KN, massive detailed balance
PASS: Danilov coefficient, vacuum limit, bounded roots, event-cap closure, quadrature, invalid-domain rejection
PASS: all exported curves finite, populated, in-frame and alive at save
```

All three plot exports and every page of the English and Spanish PDFs were
visually inspected. Figures are populated. The resonance gap is intentional.
The revised source equations, confidence conventions, unresolved normalization
discrepancies, and unmodeled detector response are documented in the bilingual
lectures and self-contained HTML. Provenance lives in the project-wide registry
and `code/darkphoton_v3/provenance/`. [[../papers/danilov-demidov-gorbunov-2019]]

## Superseded initial run

[[vault/papers/danilov-demidov-gorbunov-2019]] [[vault/papers/park-2017-reactor-dark-photons]]

ROOT/C++ implementation in `code/darkphoton_v3/v3.cpp` adds Danilov medium production and detector reconversion, then solves the TEXONO 95% event cap by logarithmic bisection. The event model uses free-electron Klein–Nishina scattering and treats deposited energy as the incident energy; detector response, anti-Compton efficiency, reactor transport, and material profiles are not modeled.

Run receipt (2026-09-21): `g++ -O2 -std=c++17 v3.cpp $(root-config --cflags --libs) -o v3 && ./v3`; ROOT 6.36.000, g++ 13.3.0. Output: `PASS: medium probability algebra, heavy-mass limit, exact epsilon-dependent bisection and epsilon^4 diagnostic`; four ROOT PDF/SVG figures written. Reports compiled with task-local Tectonic to three pages each.
