---
name: danilov-fresh
description: Calculate a new transverse hidden-photon reactor sensitivity following Danilov, with source provenance, polarization-resolved scattering and explicit validity limits.
---

For this session, build physical code only in `work/danilov_fresh`. Do not import any earlier project calculation, plot or extracted curve. Source PDFs and published supplementary numerical inputs are allowed; identify them as borrowed inputs.

Derive the damped two-state source from Danilov and Redondo. Define whether a damping coefficient attenuates amplitude or intensity. Oscillation production is not a Compton-convolved source and must not be counted twice. Use a transverse incident density matrix for detection. Derive the Compton production amplitude as a diagnostic, including its longitudinal contribution, without using it as a second source.

Use freshly implemented ROOT/C++ physics with comments and configurable experimental parameters. Independently check algebra with SymPy, Ward identities, massless Klein–Nishina and massive detailed balance. Read original TEXONO exposure and selection information. An assumed photon-selection efficiency is not a measured acceptance.

Use Izaguirre for propagation/decay geometry, but close the electron-pair channel below threshold. Check the three-photon width against McDermott's published numerical enhancement table. Do not claim to have derived that one-loop table. Exclude material resonance unless actual material profiles and geometry support its calculation.

Produce new bilingual LaTeX/PDF lectures, an HTML follow-up, a marimo presentation of the new ROOT figures, provenance and vault receipts. Record the measured wall time; unavailable model-effort/token telemetry must remain unavailable. Compare quoted source results without calibrating the prediction to them.
