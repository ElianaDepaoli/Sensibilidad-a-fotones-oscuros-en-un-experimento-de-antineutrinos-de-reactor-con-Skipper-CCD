# Where the work stands — 23 Sep 2026

Read this first, then [[index]]. Task prompt: `prompts/dark_photon_cross_section_1.txt`.

## Done and verified
1. **Audit of `apuntes/`** against `papers_reactor_dark_photon/`. Result: the exact tree-level cross section is
   confirmed (explicit Dirac matrices vs v2 eq 36, deNiverville, Smirnov, to 2e-13); the Spanish note's eq (40), (50),
   kinematics and range are wrong; the v3 report uses a different, dominant mechanism. See
   [[cross-section-audit-2026-09-21]] and `reports/report_cross_section_{en,es}.pdf` (9 pp each).
2. **`mycodes/DP3_Danilov_new.C`** — the corrections applied as `*_new_f` functions, new curves on both canvases.
   See [[dp3-danilov-new-2026-09-21]] and `mycodes/README_DP3_Danilov_new.md`.
3. **Two normalisation factors in the original macro**: `h c` -> `hbar c` (applied by the user to
   `DP3_Danilov.C`, backup `DP3_Danilov.C.bak_20260922`) and the per-atom/per-electron `Z_U = 92`
   (applied by the user in `differential_number_DP`; gives eps_95 = 1.71205e-5, matching `epsilon95_oldxs_f`).
   They had been cancelling each other — that is why `h c` looked right against Park's Fig 1.

## In progress — resume here
**`reports/explicaciones_es.pdf`** (source `reports/explicaciones_es.tex`) is a *living* Spanish Q&A document: one
section per question, format "respuesta corta -> desarrollo -> números -> verificación", with a table of contents.
7 pp as of 23 Sep 2026.
- Section 1: why the mA' = 1 MeV curve of the Fig 1 reproduction breaks (the u-channel pole at w' = mA'^2/2me; the
  `|k'| ~ w'` approximation supplies the *range*, not the pole).
- Section 2 (23 Sep): the derivation of `x_pm` (eq (5) of the report), which the report stated without deriving. Two
  routes — CM two-body momentum boosted to the lab, and cos(theta) = +-1 in the *exact* relation (the note's own
  argument), which gives a quadratic with discriminant 16 E^2 lambda. With `|k'| ~ w'` that quadratic degenerates to
  two linear equations: that is where the second root and the threshold are lost. Also: x_+ is always forward
  emission, x_- is backward only when 2E(me-M) > M(2me-M) — at reactor energies with M >~ 0.5 MeV the A' is confined
  to a forward cone (67.4 deg at E=3, M=0.5). Verified in `work/check_xpm.py` (sympy identities + 296-point numeric
  scan; `work/check_xpm.log`).

- Section 3 (23 Sep): why their own `DP3_Danilov.C` spikes at w' -> mA' after taking only the `*_new_f` kinematics
  (see "Diagnosed, fix not applied" below). Includes the one-minute `delta_frac` diagnostic and the three fix options.

**The user will keep asking questions and each one becomes a new section of that PDF.** Candidates left, in the
order offered to them (a comment at the end of the .tex repeats this list):
- why dropping the `me^2` terms in eq (35)-(36) breaks Klein-Nishina;
- the jacobian / why `dsigma/dt` with `dt = 2 me dw'` avoids it;
- the per-electron normalisation.

Covered along the way: the threshold `E_th = mA' + mA'^2/2me` (section 2, from lambda >= 0;
`reports/ampliacion_del_reporte.pdf` still has the longer version) and the exact range `x_pm` (section 2).

## How to pick up
```
cd Proyecto_Final_DarkPhoton_claudia
work/.venv/bin/python work/check_cross_section.py work/check_cross_section_out.json   # ends in ALL ASSERTIONS PASSED
cd mycodes && root -l DP3_Danilov_new.C                                               # ~5 s, writes the two png
work/.venv/bin/python work/check_xpm.py work/check_xpm_out.json                       # ends in ALL ASSERTIONS PASSED
cd reports && latexmk -pdf explicaciones_es.tex                                       # after adding a section
```
Every number quoted anywhere is in `provenance/numbers.json`; every assertion in `provenance/claims.yaml`. The Stop
hook checks them (`echo '{}' | python3 .claude/hooks/provenance_gate.py` to dry-run; it passes as of 23 Sep 2026).

## Git
Nothing has been committed yet, by their instruction. The plan (stated 22 Sep 2026) is to commit this work
**together with** the Codex work in `~/Documentos/Materias/IA_OG/GW-AI-course/Proyecto_Final_DarkPhoton_gpt`, as one
consolidated final project. Do not commit this directory on its own; ask before touching the other tree.

## Conventions settled with the user
- Reports: LaTeX -> PDF. `.md` files in English (project rule), but `explicaciones_es` is Spanish because they asked.
- Never edit `apuntes/`; never edit the original functions in `mycodes/DP3_Danilov.C` — new work goes in
  `*_new_f` / `*_oldxs_f` functions inside `DP3_Danilov_new.C`, keeping the note's style and variable names.
- The anonymous second note `apuntes/Dark_Photon_cross_Section.pdf` was excluded from the review at their request.
- They want the *fix* for each error, not just the diagnosis.

## Diagnosed, fix not applied (23 Sep 2026)
**The spike at w' -> mA' in the user's `mycodes/DP3_Danilov.C`.** They took only the `*_new_f` kinematics into their
own macro (their `wp_production_f` is now the new grid, starting at w' = mA'), but kept eq (50) `ds_dwp_f` and the old
domain `Rw_f`. Eq (50) carries cos^2(theta) explicitly and the corrected `cos_theta_f` goes like 1/|k'|, so eq (50)
has a simple pole at w' = mA'; `Rw_f` never returns an empty photon range there (`invwpmin_f` = +inf for w' >= me/2).
The old grid started at mA'/2, below the mass, where `cos_theta_f` is NaN and the guard fired — that is why nothing
showed before. `work/check_wp_production_spike.py` (ends in ALL ASSERTIONS PASSED). Claim
`wp-production-new-grid-spike`; written up as section 3 of `reports/explicaciones_es.pdf`. The fix is to take `Rw_new_f` and `ds_dwp_new_f` as well: the three go together.
Offered to them, not applied — their macro is theirs to edit.

## Open questions not yet settled
- Park's Fig 1 tail above 2 MeV equals the photon spectrum, not his eq (1); the v2 report's Fig 2 follows Park.
  Neither generator is available. See the claim `park-figure1-tail-is-the-photon-spectrum`.
- `wpsize = 10` in the macro (line 170) still gives ~5 A' energies inside the 3-8 MeV window: the steps and the
  narrow dips near 0.1 MeV in the exclusion curves. On 23 Sep only `wpsize_fig1` was raised to 300 (Figure 1 block),
  and `cFig1_new.png` / `cSk_new.png` were regenerated; the global `wpsize` is untouched.
- The whole macro is production channel A (Compton-like); channel B (oscillation) is 7.7x larger in that window.
