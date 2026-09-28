# Sensitivity to dark photons in a reactor antineutrino experiment with Skipper CCDs

Final project of the course *Artificial Intelligence for Gravitational-Wave Astronomy* (IA_OG),
done as a physics study in its own right: dark-photon production at a nuclear reactor and its
detection with a Skipper CCD, following Park (2017) and its successors, and re-deriving the
production cross section from scratch to audit the calculations the project was handed.

The same problem was worked twice, in parallel, with two agent harnesses. Both trees are kept
side by side rather than merged, because each carries its own notes and provenance record and the
comparison between them is part of the result.

```
claude/   worked with Claude Code
codex/    worked with Codex
```

Each tree keeps the same internal layout:

| directory                                                       | contents                                                                                                                                                                                                                           |
| --------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `vault/`                                                        | one short page per paper, per exercise, per concept. **Start at `claude/vault/status.md`**                                                                                                                                         |
| `work/`, `code/`                                                | scripts and checks; every derivation is verified with `sympy` or against a published formula                                                                                                                                       |
| `reports/`, `report_*`                                          | LaTeX write-ups (EN + ES) and the agent-facing HTML report per task                                                                                                                                                                |
| `provenance/`                                                   | `numbers.json` and `claims.yaml`: what each number asserts, what produced it, what was written from scratch and what came from a library. A Stop hook refuses to end a session that produced a figure or a number without a record |
| `papers_reactor_dark_photon/`, `dark_photon_papers/`, `papers/` | the bibliography actually used, with its manifest                                                                                                                                                                                  |
| `prompts/`                                                      | the task prompts, verbatim                                                                                                                                                                                                         |
| `.claude/`, `.codex/`, `AGENTS.md`                              | the working agreement, hooks and skills the agents ran under                                                                                                                                                                       |

## Main results so far

- The exact tree-level cross section for `gamma e -> A' e` is confirmed against three independent
  published forms and against explicit Dirac matrices (agreement to 2e-13). See
  `claude/reports/report_cross_section_en.pdf`.
- The kinematics used in the material under review is wrong in a way that matters: the threshold is
  `E = m_A' + m_A'^2 / 2 m_e`, not `E = m_A'`, and the exact energy range `[x_-, x_+]` excludes a
  u-channel pole that the approximate range integrates across.
- `claude/reports/explicaciones_es.pdf` is a living Q&A document (in Spanish, by request): one
  section per question, each ending in an explicit verification and a statement of what was *not*
  checked.
- **Final summary** (in Spanish): `claude/reports/final_summary_es.pdf`. It says what was done by the
  agents alone and what was done on the author's own code, and puts every TEXONO 95% CL limit of the project
  on one plot (`claude/work/exclusion_comparison/texono_exclusion_all.pdf`). Once the cross section and the
  integration grid are corrected, the author's macro and Codex's independent code agree to 3%
  (eps_95 ~ 1.5e-5). Neither reproduces Park's 2.1e-5.

## What is deliberately not in this repository

- **The documents under review** (`apuntes/`): unpublished notes by third parties. The reports that
  audit them are here; the notes themselves are not ours to redistribute.
- **The ROOT macro `DP3_Danilov.C`** and its variants, pending publication. Its figures, run logs and
  input tables are included, so the work it produced is visible.
- Virtual environments, LaTeX caches and compiled binaries.

## Conventions

- Write-ups are LaTeX compiled to PDF; documentation is in English, except where a document was
  explicitly asked for in Spanish.
- A result with no assertion is a hypothesis, and is labelled as one. Numbers quoted anywhere in the
  reports are registered in the corresponding `provenance/numbers.json`.

Author: E. L. De Paoli. Agent-assisted work; each tree's `provenance/` says which parts were derived,
which were recalled, and which were verified.
