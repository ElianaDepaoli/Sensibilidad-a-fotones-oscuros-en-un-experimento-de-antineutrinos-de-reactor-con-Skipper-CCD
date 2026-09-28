# Fresh Danilov-guided TEXONO sensitivity

No prior project physical calculation or plot was reused. New code: `work/danilov_fresh/`; reports: `report_danilov_fresh*`. Session skill: [[../skills/danilov-fresh/SKILL]]. Related sources: [[../papers/danilov-fresh-corpus]].

New derivations: covariant total trace, analytic longitudinal absorption, transverse density matrix, damped source integral, event rate including optional detector decays. Borrowed inputs: photon spectrum, exposure, selection hypothesis and published exact-loop enhancement. The loop itself is not rederived. Published normalization comparison fails; no tuning. Material resonance, detector response and massive/longitudinal nuclear source remain unvalidated.

## Run receipt

Date: 2026-09-25T18:53:02.112743+00:00
Command: `bash work/danilov_fresh/run.sh`. Environment: ROOT 6.36.000, Ubuntu g++ 13.3.0, existing SymPy, LaTeX and marimo. No installation.

```text
PASS: fresh covariant trace equals compact F; Klein-Nishina differential identity; s/u symmetry
PASS: Hamiltonian discriminant including epsilon^4; damped-amplitude primitive; thick-source probability
PASS: absorption longitudinal trace derived; vanishes in the massless limit
F_L = 8*d*(a*b - 4*a*m**2 + 2*d*m**2 - 8*m**4)*(a**2*b + a**2*m**2 + a*b**2 - a*b*d + 2*a*b*m**2 + b**2*m**2)/(a**2*b**2*(a**2 - 2*a*d + d**2 - 4*d*m**2))

PASS: fresh polarized amplitudes, Ward identities, symbolic trace comparison, both KN limits, massive detailed balance, doubled angular quadrature
PASS: Danilov low-mass coefficient and rounded event cap; published exact-loop enhancement knots; decay mixing scaling
Info in <TCanvas::Print>: pdf file output/polarized_cross_sections.pdf has been created
Info in <TCanvas::Print>: SVG file output/polarized_cross_sections.svg has been created
Info in <TCanvas::Print>: png file output/polarized_cross_sections.png has been created
Info in <TCanvas::Print>: pdf file output/oscillation_source.pdf has been created
Info in <TCanvas::Print>: SVG file output/oscillation_source.svg has been created
Info in <TCanvas::Print>: png file output/oscillation_source.png has been created
Info in <TCanvas::Print>: pdf file output/texono_exclusion.pdf has been created
Info in <TCanvas::Print>: SVG file output/texono_exclusion.svg has been created
Info in <TCanvas::Print>: png file output/texono_exclusion.png has been created
Info in <TCanvas::Print>: pdf file output/decay_enhancement.pdf has been created
Info in <TCanvas::Print>: SVG file output/decay_enhancement.svg has been created
Info in <TCanvas::Print>: png file output/decay_enhancement.png has been created
selected_cap = 196.172274871
fresh_plateau_epsilon95 = 1.33627026745e-05
paper_rescaling_epsilon95 = 2.96412165583e-05
max_decay_optical_depth = 1.08388793272e-08
PUBLISHED NORMALIZATION COMPARISON: FAIL (no fitted factor)
PASS: event-cap closure, doubled energy/angle quadrature, decay-volume scaling, monotonic-domain guard, epsilon-fourth-power check
UNCHECKED: material profiles and resonance; full detector response; longitudinal nuclear emission; original likelihood; exact massive nuclear-source correction
ROOT 6.36.000; C++ 13.3.0
compute_seconds = 2.708228471
```
