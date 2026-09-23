# Reactor dark-photon literature search

## Outcome

- Eight primary PDFs collected in `papers_reactor_dark_photon/`.
- Technical report: `report_reactor_dark_photon.html`.
- Human follow-up: `rfh_reactor_dark_photon_en.html` and `rfh_reactor_dark_photon_es.html`.
- Main correction chain: [[park-2017-reactor-dark-photons]] → [[danilov-demidov-gorbunov-2019]]; [[du-et-al-2024-scattering-enhancement]] → [[demidov-gorbunov-polonski-2025]].
- Direct result kept separate: [[neon-2025-light-dark-matter]].

## Not reproduced

No exclusion curve was digitized or recomputed. No collaboration likelihood was rerun. Numerical results are sourced claims; see `provenance/`.

## Run receipt

- Date: 2026-09-17
- Command: `python3 /home/eliana/.codex/skills/.system/skill-creator/scripts/quick_validate.py skills/reactor-dark-photon-literature && python3 work/verify_reactor_dark_photon_outputs.py`
- Environment: bash; Python 3 standard library; `pdftotext` 22.12 used during source inspection; local workspace only for verification.
- Actual output:

```text
Skill is valid!
PASS: 8 PDFs match the manifest
PASS: 3 HTML files parsed; all local links resolve
PASS: 42 quantitative provenance records are structurally complete
```

The stop-hook initially classified the eight preserved source PDFs as figure artifacts. Added one `figures:` provenance entry per PDF, recording acquisition, source status, selection choice, and supported claims. Final gate command: `python3 .claude/hooks/provenance_gate.py`; actual output was empty with exit status 0.
