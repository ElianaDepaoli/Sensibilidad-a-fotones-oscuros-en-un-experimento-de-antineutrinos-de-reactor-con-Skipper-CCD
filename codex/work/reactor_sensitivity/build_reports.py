"""Assemble reports and provenance; no physics calculations are performed here."""
from pathlib import Path
import csv
import html
import json
import re
import hashlib
from datetime import datetime, timezone
import yaml

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "work/reactor_sensitivity"

def read(path):
    return (ROOT / path).read_text()

def assemble():
    rows = list(csv.DictReader((BASE / "output/limits.csv").open()))
    receipt = (BASE / "output/run_receipt.txt").read_text()
    peak = float(re.search(r"maximum_EH_optical_depth_at_vacuum_limit = (\S+)", receipt)[1])
    elapsed = float(re.search(r"compute_seconds = (\S+)", receipt)[1])
    now = datetime.now(timezone.utc).isoformat()
    summary = r"\begin{center}\begin{tabular}{lr}\hline " + "\n"
    summary += r"\both{Largest approximate optical depth}{Mayor profundidad óptica aproximada} & $" + f"{peak:.5g}".replace("e-",r"\times10^{-") + r"}$ \\" + "\n"
    summary += r"\both{Upper scan mass (MeV)}{Masa superior del barrido (MeV)} & " + rows[-1]["mass_MeV"] + r" \\ \hline\end{tabular}\end{center}" + "\n"
    (BASE / "output/summary.tex").write_text(summary)
    lecture = read("code/darkphoton_v3/lecture_redo.tex")
    lecture = lecture.replace("21 September 2026", "25 September 2026")
    lecture = lecture.replace(r"\begin{thebibliography}", read("work/reactor_sensitivity/addendum.tex") + "\n" + r"\begin{thebibliography}")
    (BASE / "lecture.tex").write_text(lecture)
    for lang, switch in [("en", "false"), ("es", "true")]:
        wrapper = read(f"report_darkphoton_v3_{lang}.tex").replace("code/darkphoton_v3/lecture_redo.tex", "work/reactor_sensitivity/lecture.tex")
        (ROOT / f"report_reactor_sensitivity_{lang}.tex").write_text(wrapper)
    figures = [
        ("Scattering", "code/park_texono/output/cross_sections.svg"),
        ("Park Figure 1 comparison: not fully reproduced", "code/darkphoton_v3/output/figure1_comparison.svg"),
        ("Medium source", "code/darkphoton_v3/output/figure1_medium.svg"),
        ("TEXONO comparison: not the full Park Figure 2", "code/darkphoton_v3/output/figure2_texono_medium.svg"),
        ("Decay-aware conditional sensitivity", "work/reactor_sensitivity/output/conditional_limits.svg")]
    page = ['<!doctype html><html lang="en"><meta charset="utf-8"><title>Reactor sensitivity execution</title><style>body{max-width:1100px;margin:3em auto;padding:1em;font:18px/1.5 system-ui;color:#18232b}svg{width:100%;height:auto}pre{white-space:pre-wrap;font-size:13px;background:#eef2f5;padding:1em}table{border-collapse:collapse}td,th{border:1px solid #bbb;padding:.6em}</style><h1>Reactor dark-photon sensitivity</h1>',
        '<p>The ROOT calculation and bilingual lectures were regenerated. Full numerical reproduction of Park is not achieved. These are conditional sensitivity calculations, not a validated detector-level exclusion.</p>',
        '<p><a href="report_reactor_sensitivity_en.pdf">English lecture</a> · <a href="report_reactor_sensitivity_es.pdf">Clase en español</a> · <a href="work/reactor_sensitivity/output/figures.html">Marimo figures</a></p>',
        '<p>New work: decay propagation, cached event integration and bounded first-crossing solver. Reused and rerun: QED traces, exact massive production/inverse scattering and Danilov medium formulas. The three-photon loop width is borrowed from Park, not derived here. Izaguirre Eq. (9) supplies the survival/decay geometry; its pair decay is closed in this mass range.</p>',
        '<table><tr><th>Requested reference</th><th>Status</th></tr><tr><td>Park Figure 1</td><td>Calculated and compared with extracted vector data. Heavier-mass curves disagree: not reproduced.</td></tr><tr><td>Park Figure 2</td><td>TEXONO component calculated; printed bound disagrees. NEOS/private analysis and external compilation are not reproduced.</td></tr></table>',
        '<p>The approximate three-photon width has not been validated near the electron-pair threshold. Material resonance, detector response, Compton-beam polarization transport, reactor transport and the original event likelihood remain unchecked. No fitted normalization is used.</p>',
        '<p>Run <code>bash work/reactor_sensitivity/run.sh</code>. The configuration exposes power, baseline, exposure, mass, composition, ROI, efficiency, medium inputs and an external event cap or Gaussian excess/error prescription. Baseline checks: <code>bash code/park_texono/run.sh</code> and <code>make -C code/darkphoton_v3 reports</code>.</p>']
    for title, path in figures:
        svg = read(path)
        page += [f'<h2>{html.escape(title)}</h2>', svg[svg.index('<svg'):]]
    for title, path in [("Fresh symbolic checks", "code/park_texono/output/symbolic_receipt.txt"), ("Fresh Park run", "code/park_texono/output/run_receipt.txt"), ("Fresh medium run", "code/darkphoton_v3/output/run_receipt.txt"), ("Decay-aware run", "work/reactor_sensitivity/output/run_receipt.txt")]:
        page += [f'<h2>{title}</h2><pre>{html.escape(read(path))}</pre>']
    page += ['<p>Sources: local Park, Danilov, Gondolo–Raffelt and Izaguirre PDFs listed in the lecture. Implementation reference: <a href="https://root.cern.ch/doc/master/classTGraph.html">ROOT TGraph documentation</a>. Bibliographic cross-check: <a href="https://arxiv.org/abs/1507.02681">Izaguirre et al.</a></p>',
             f'<footer><h2>Execution metadata</h2><p>Generated {now}. Model: GPT-6 (Codex). Reasoning effort: unavailable. Tokens used: unavailable. End-to-end resolution time: not instrumented. Measured C++ computation: {elapsed} seconds; this excludes reading, derivation, checks, compilation and report preparation.</p></footer></html>']
    (ROOT / "report_reactor_sensitivity.html").write_text('\n'.join(page))
    numbers_path = ROOT / "provenance/numbers.json"
    numbers = json.loads(numbers_path.read_text())
    numbers["reactor-sensitivity-execution"] = dict(
        value={"limits": rows, "maximum_approximate_optical_depth": peak, "compute_seconds": elapsed,
               "configuration": read("work/reactor_sensitivity/config.conf"),
               "park_run": read("code/park_texono/output/run_receipt.txt")},
        statement="Conditional decay-aware boundaries, runtime, configured inputs and fresh Park comparison; not measured detector limits.",
        produced_by="work/reactor_sensitivity/sensitivity.cpp::Calculate,Run; code/park_texono/park_texono.cpp",
        from_scratch="Survival, cached convolution weights and first upward crossing added in this execution; existing checked scattering kernels reused.",
        from_library="ROOT Gauss-Legendre quadrature, TGraph and TCanvas; source and medium parameters from Park and Danilov.",
        choices=["Same selected cap and efficiency for comparing source hypotheses", "Park unpolarized inverse convention retained for vacuum", "Low-mass three-photon approximation is a diagnostic near the upper endpoint", "No material resonance claim", "No raw-likelihood reconstruction or normalization fit"])
    numbers_path.write_text(json.dumps(numbers, indent=2) + '\n')
    claims_path = ROOT / "provenance/claims.yaml"
    claims = yaml.safe_load(claims_path.read_text())
    identifier = "reactor-sensitivity-execution"
    claims["claims"] = [c for c in claims["claims"] if c["id"] != identifier]
    claims["claims"].append(dict(id=identifier, statement="Decay-aware conditional boundaries satisfy their configured event caps. Park figure reproduction remains incomplete; the width near threshold and detector response are not validated.", evidence=["work/reactor_sensitivity/sensitivity.cpp::Checks,Calculate", "work/reactor_sensitivity/output/run_receipt.txt", "code/park_texono/output/run_receipt.txt", "papers/park_reproduction/reference20.pdf::Eq.9", "papers/darkphoton_v3/park2017.pdf::Eq.4"], numbers=[identifier]))
    files = [f"work/reactor_sensitivity/output/conditional_limits.{ext}" for ext in ("svg", "pdf", "png")]
    files += [f"report_reactor_sensitivity_{lang}.pdf" for lang in ("en", "es")]
    files += [str(p.relative_to(ROOT)) for p in (BASE / "output").glob("lecture_*.png")]
    claims["figures"] = [f for f in claims["figures"] if f["file"] not in files]
    for path in files:
        claims["figures"].append(dict(file=path, produced_by="work/reactor_sensitivity/sensitivity.cpp::Calculate; work/reactor_sensitivity/build_reports.py::assemble", shows="Mixing boundary versus dark-photon mass, with separate vacuum and medium hypotheses; PDFs also explain the reused source-comparison plots.", from_scratch="Decay propagation and conditional event-boundary calculation", from_library="ROOT rendering; Tectonic PDF layout; checked project kernels", choices=["Separate medium graphs across excluded resonance", "Label approximate decay and unvalidated response"], supports=[identifier]))
    claims_path.write_text(yaml.safe_dump(claims, sort_keys=False, allow_unicode=True))
    hashes = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in [BASE / "sensitivity.cpp", BASE / "config.conf", BASE / "output/limits.csv"]}
    (BASE / "output/manifest.json").write_text(json.dumps({"generated_utc": now, "sha256": hashes}, indent=2)+'\n')
    vault = '# Reactor sensitivity prompt execution\n\nExecuted `prompts/reproducesensitivitytodarkphotoninreactors.txt`. New ROOT/C++ survival and event-boundary calculation reuses the verified massive scattering kernels. English/Spanish LaTeX lectures and self-contained HTML: `report_reactor_sensitivity*`. Session skill: [[../skills/reactor-sensitivity/SKILL]]. Related: [[darkphoton-v3]], [[park-texono-reproduction]].\n\nPark Figure 1 and the printed TEXONO limit remain unreproduced. The full Figure 2 compilation is not reproduced. The three-photon width is borrowed, and is not an exact near-threshold calculation. No new microscopic amplitude derived. Detector response and raw likelihood remain unchecked.\n\n## Run receipt\n\nDate: '+now+'\nCommand: `bash work/reactor_sensitivity/run.sh` (compiles ROOT/C++, runs calculation, assembles and compiles both lectures, exports marimo). Environment: existing ROOT 6.36.000, g++ 13.3.0, Python 3.12, project Tectonic cache and existing marimo installation; no installation.\n\n```text\n'+receipt+'```\n\nIndependent earlier reruns in this execution: `bash code/park_texono/run.sh` and `make -C code/darkphoton_v3 reports`; actual outputs are embedded in the HTML. Provenance: `provenance/numbers.json`, `provenance/claims.yaml`, and the hashed output manifest.\n'
    if (BASE / "output/export_receipt.txt").exists():
        vault += '\n## Export and visual review\n\n' + (BASE / "output/export_receipt.txt").read_text()
    (ROOT / "vault/exercises/reactor-sensitivity-execution.md").write_text(vault)
    print("Reports assembled; provenance and vault receipt updated.")

if __name__ == '__main__':
    assemble()
