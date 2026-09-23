# Cross-section audit of apuntes/ against papers_reactor_dark_photon/ (21 Sep 2026)
Prompt: `prompts/dark_photon_cross_section_1.txt`.

**Attempted.** Independent |M|² by explicit Dirac matrices + explicit polarization vectors at random on-shell points; comparison with v2 eq 36, deNiverville eq 18–20, Smirnov eq 3, and the Spanish notes' eq 40/41, 50/54, 52, 64; exact kinematics vs notes; polarization-resolved inverse cross sections; Park-channel vs oscillation-channel counts in 3–8 MeV.
**Came out.** See [[dark-compton-cross-section]], [[two-thirds-detection-factor]], [[oscillation-vs-compton-production]], [[longitudinal-fraction]], [[apuntes-spanish-notes]]. Reports: `reports/report_cross_section_{en,es}.pdf`, `reports/report_cross_section_agent.html`.
**Revision (same day).** Note 2 dropped from the review at the user's request; fixes added for note 1; figure 4 (dN/dE of both mechanisms) added; found that Park's Fig 1 tail = S_γ, not his eq 1 (see [[park-2017]]).
**Not checked.** Gondolo–Raffelt total, Park Fig 1, XCOM numbers, Du/Demidov MCs, resonance region, bound electrons.

## Run receipt
- date: 2026-09-21
- command: `work/.venv/bin/python work/check_cross_section.py work/check_cross_section_out.json` then `work/.venv/bin/marimo export html work/figures_cross_section.py -o work/figures_cross_section.html --no-include-code`
- environment: venv `work/.venv` (system site-packages: numpy 2.5.2, sympy 1.12, scipy 1.18.1, matplotlib 3.11.1; marimo 0.24.2 installed in the venv), Linux 6.8, pdflatex/latexmk for the reports
- output: `work/check_cross_section.log` (last line `ALL ASSERTIONS PASSED`), `work/check_cross_section_out.json`, `work/fig_*.pdf`, `work/figures_cross_section.html`
