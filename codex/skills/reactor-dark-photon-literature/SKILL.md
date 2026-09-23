---
name: reactor-dark-photon-literature
description: Research reactor-based dark-photon searches and produce a locally auditable paper set plus concise English and Spanish HTML reports. Use for literature updates on dark photons produced or detected at reactor-neutrino experiments; do not merge mediator-only neutrino-scattering bounds with reactor-produced dark-photon searches.
---

# Reactor dark-photon literature

Start from any seed paper in `papers_init/`, then trace its references, papers that cite or correct it, and later searches using reactor-on/off data.

Preserve every paper used in a task-specific local folder. Prefer the author or publisher PDF, record the canonical URL and SHA-256 digest, and distinguish these evidence classes:

- direct collaboration result using reactor data;
- phenomenological recast of published experimental data;
- projected sensitivity;
- theory or transport correction.

Keep physically different signal models separate: visible kinetically mixed dark photons, invisible dark-photon decay to light dark matter, and dark-axion-portal production do not constrain the same parameter plane.

For each search, extract the reactor production mechanism, detector signature, analysis energy window, mass reach, reported confidence level, and dominant limitation. If later work disputes a limit, place the original and criticism together and state what changed.

Produce a technical English HTML report and short human follow-up pages in English and Spanish. Link every claim to a local PDF and its canonical source. Record quantitative claims in `provenance/` and leave a run receipt in `vault/exercises/`.
