---
name: formula-audit
description: Use when asked to compare theoretical derivations or published formulas (cross sections, rates, probabilities) across several documents and decide which differences are errors. Covers reading PDFs, re-deriving the target quantity from scratch with an assertion, comparing every candidate formula to it, and delivering the lecture-style PDF (EN + ES), the agent-facing HTML, the provenance record and the vault pages this project requires.
---

# Formula audit

The deliverable is a verdict per formula, backed by one independent computation, not a summary of who says what.

## 1. Read everything before computing
- `pdftotext -layout` every PDF in the target folders into the scratchpad; read the notes in full, grep the papers for the equations that matter (cross section, |M|², threshold, normalization, detection factor). Note equation numbers as you go — every verdict cites one.
- Build one table early: document → what it computes → its final formula → its assumptions. Different *physics* (e.g. oscillation vs scattering) is not a different derivation of the same thing; say so explicitly.

## 2. Derive the reference once, independently of all of them
- Prefer a method none of the documents used. For a 2→2 QED amplitude: explicit 4×4 Dirac matrices, explicit polarization vectors (all of them, longitudinal included), random on-shell kinematic points, spin sum by matrix multiplication. No trace identities, so nothing is shared with the trace-based derivations under review.
- Every block ends in an `assert`. Tolerances: 1e-9 relative for numeric identity, exact `sp.simplify(...) == 0` for symbolic.
- Transcribe each candidate formula into its own function with a docstring naming source + equation. Compare ratios, not differences.
- When a candidate fails, isolate *why* (e.g. subtract at a limit and simplify) so the report can name the missing term.
- Check kinematics separately from the matrix element: threshold, endpoints, angle relation. Errors compound.
- When integrating a candidate over its own claimed range, first look for poles inside it; if there is one, say the integral does not exist rather than reporting a number.

## 3. Quantify differences that are physics choices, not algebra
- Put a number on each modeling difference (e.g. events in the experiment's window under mechanism A vs B). "They differ" is not a finding; "×7.7 in 3–8 MeV" is.
- Say which mechanism is leading and cite the paper that settles it; do not adjudicate Monte Carlo disputes you did not rerun — report consistency of your simpler estimate with one side.

## 4. Deliverables for this project (all of them, every time)
- `work/check_*.py` (assertions, `main()` guarded so the notebook can import it), `work/check_*.log`, `work/check_*_out.json`.
- Figures from a **marimo** notebook `work/figures_*.py` exported with `marimo export html … --no-include-code`; PDF figures saved by the notebook cells. Per-exercise venv `work/.venv` (`--system-site-packages`, then `pip install marimo`).
- Human report: LaTeX → PDF, lecture structure (setup → kinematics → amplitude → cross sections → per-document verdict table → mechanism comparison → the formula to use → checked/not checked + run receipt). One `.tex` per language; the Spanish one is a translation, same equations, same numbers. Use `longtable` for the verdict table, `slashed` for Feynman slashes. Style: short declarative sentences, no hedging, no politeness; every number traceable to the JSON.
- Agent report: single self-contained `reports/report_*_agent.html` with: file map, the verdict table, the numbers, what is not checked, model + reasoning effort, wall time, tokens (from the `<total_tokens>` counter: 15 000 000 − remaining).
- `provenance/numbers.json` (six fields each) and `provenance/claims.yaml` (claims + every `.pdf/.png` newer than session start, including report PDFs and any handed-in file with a late mtime). Dry-run: `echo '{}' | python3 .claude/hooks/provenance_gate.py`. Delete scratch PNGs before finishing.
- Vault: one page per paper, per note under review, per concept; one exercise page with a run receipt; `vault/index.md`.

## 5. Traps met on 2026-09-21
- `sympy.Abs` inside `log` blocks `simplify`; use the sign known from the physical region (`log(-b/a)` for `b<0`).
- pdftotext garbles nested fractions in equations — reconstruct from the surrounding text and *verify the transcription* (e.g. `dG/db = F`) before judging the source.
- The provenance hook treats every `.pdf`/`.png` newer than session start as a figure, including reports and files the user dropped in after the session began.
- Units: cross sections computed in MeV⁻² need `(ħc)² = (1.973269804e-11 cm·MeV)²` before plotting in cm² or barn.
