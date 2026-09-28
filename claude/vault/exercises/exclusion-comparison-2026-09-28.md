# All TEXONO exclusion limits on one plot + final summary (28 Sep 2026)
**Attempted.** One figure with every eps_95(m_A') of the project: Park (quoted 2.1e-5), Danilov (Codex digitization),
the user's `mycodes/DP3_Danilov.C` and Claude's `DP3_Danilov_new.C` with DANILOV A/B/D off and on, and the four Codex
curves (`park_texono`, `reactor_sensitivity`, `darkphoton_v3`, `danilov_fresh`). Plus a 4-5 page report EN + ES separating
LLM-only work (Track 1) from work on the user's code and inputs (Track 2): `reports/final_summary_{en,es}.pdf`.
Pipeline: `work/exclusion_comparison/run.sh` runs both macros on a temp copy (flags set as globals by `dump_cSk.C`, no
file in `mycodes/` edited), copies the Codex CSVs to `data/`, draws `texono_exclusion_all.{pdf,png}`.

**Came out.**
- The macros' default `wpsize = 10` is not converged: flags-off plateau 1.26e-5 (10) -> 1.44 (40) -> 1.48 (100) -> 1.50 (200)
  -> 1.504 (400). The rise to 1 MeV in the default curves is grid noise.
- Converged, exact sigma (`_new.C`): 1.51e-5 plateau, 1.525e-5 at 1 MeV; Codex `park_texono` 1.565e-5 and 1.522e-5.
  Human chain + corrected sigma = independent Codex code, to 3% / 0.3%. See [[dp3-danilov-new-2026-09-21]].
- Exact vs eq. (50), converged: 1% on the plateau, 18% at 1 MeV (1.525 vs 1.29) = 1.99^(1/4). See [[cross-section-audit-2026-09-21]].
- Flags on, converged: 2.03e-5 (close to Park by coincidence of A x B). Codex oscillation 1.34e-5; Danilov digitized 3.47e-5;
  Park x 6^(1/4) = 3.29e-5 is within 5% of Danilov's plateau. See [[oscillation-vs-compton-production]], [[danilov-demidov-gorbunov-2019]].
- The dip of the human curve at 19.2 eV (1.8e-6) is the cap at 100 in `danilov_suppression`, not physics.
- `README_DP3_Danilov_new.md` sec. 4 says DANILOV_A is "off by default"; it is ON (and ON means factor 1). Not edited.

**Not checked.** Codex CSVs/digitization read as is, Codex codes not re-run. wpsize convergence of the A+B+D runs only at 200.

## Run receipt
- date: 2026-09-28; env: ROOT 6.36.000 (`~/root/bin/root`), Linux 6.8, latexmk/pdflatex (TeX Live)
- `work/exclusion_comparison/run.sh` -> real 2m01s; table printed by the plot macro, saved in `work/exclusion_comparison/table.log`:
  `DP3_Danilov.C off 1.259e-05 (1e-4 MeV) ... 1.465e-05 (1 MeV)`; `off wp200 1.496e-05 ... 1.291e-05`;
  `_new.C exact off wp200 1.513e-05 ... 1.525e-05`; `Codex park_texono 1.565e-05 ... 1.522e-05`; `Codex danilov_fresh 1.283e-05 ... 1.345e-05`
- by hand (scratch copy, same `dump_cSk.C`), flags off, plateau / 1 MeV: `DP3_Danilov.C` wpsize 40: 1.4371e-05 / 1.2444e-05;
  100: 1.4810e-05 / 1.3031e-05; 400: 1.5040e-05 / 1.2845e-05. `_new.C` wpsize 400: 1.5209e-05 / 1.5183e-05
- `cd reports && latexmk -pdf final_summary_en.tex` -> 4 pages; `final_summary_es.tex` -> 5 pages, no errors, no overfull boxes
- `echo '{}' | python3 .claude/hooks/provenance_gate.py` -> exit 0
