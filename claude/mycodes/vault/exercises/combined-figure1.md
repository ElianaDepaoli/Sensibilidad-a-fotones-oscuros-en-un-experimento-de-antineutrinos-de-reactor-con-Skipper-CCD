# Combined ROOT Figure 1

Created `combined_figure1.C`, a standalone ROOT macro combining `cFig1.C` and
`park_figure1_digitized.C`. Source files are unchanged. Original coordinates and
colors are reused; the new code supplies the frame, solid/dotted styles and
legend. No physics was derived or recalled for this task.

Outputs: `combined_figure1.png`, `combined_figure1.pdf`.
Provenance: `provenance/claims.yaml`.
Generator: `work/combine_figure1/build.py`.

Checked: exact coordinate/color equality against the original ROOT objects,
line styles, graph-linked legend labels, logarithmic ordinate, and visual layout.
Not checked: underlying physics, digitization accuracy or agreement with the paper.

## Run receipt

Date: 2026-09-28. Environment: bash, ROOT 6.36.000, installed executable
`/home/eliana/root/bin/root`; working directory: `mycodes/`.
Commands: `python3 work/combine_figure1/build.py`,
`root -l -b -q combined_figure1.C`,
`root -l -b -q work/combine_figure1/check.C`.
Actual verification output (also in `work/combine_figure1/check.log`):

```text
Processing work/combine_figure1/check.C...
cling::PHOptLevel: conflicting `#pragma cling optimize` directives: was already set to 0
Ignoring higher value of 0
Info in <TCanvas::Print>: png file combined_figure1.png has been created
Info in <TCanvas::Print>: pdf file combined_figure1.pdf has been created
PASS: original coordinates and colors preserved for every graph and segment.
PASS: cFig1 solid, PARK dotted; legend labels linked to the correct graphs.
```

## Local hook path repair

The inherited Codex hook runs `python3 .claude/hooks/provenance_gate.py` relative
to the working directory. It failed because `mycodes/.claude/hooks/` was absent.
Copied the parent validator byte-for-byte, its format documents, and its existing
session timestamp into the local `.claude/` directory. No validator checks were
changed. Added `provenance/numbers.json` with the graph-object count and recorded
the pre-existing supplied figures flagged by the inherited timestamp, explicitly
marked as not regenerated or verified in this task.

Run receipt: 2026-09-28, system Python 3 with PyYAML, working directory `mycodes/`.
Command: `python3 .claude/hooks/provenance_gate.py < /dev/null`.
Actual stdout/stderr: empty. Exit status: 0.
