"""Report the continuation diagnostics without regenerating baseline physics."""
import csv
from datetime import datetime, timezone
import hashlib
from html import escape
import json
from pathlib import Path
import statistics
import yaml
from make_reports import CSS, fig, table

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[1]
OUT = BASE / 'output/continuation'
START = '2026-09-17T22:37:52+00:00'  # First recorded diagnostic timestamp, not session start.
VARIANTS = ('baseline', 'low_cutoff', 'conditional', 'conditional_low_cutoff')


def read_csv(name):
    with (OUT / name).open() as stream:
        return [{k: float(v) for k, v in row.items()} for row in csv.DictReader(stream)]


def summarize():
    """Use each horizontal PDF segment once, comparing its bin-average prediction."""
    rows = read_csv('bin_comparison.csv')
    result = {}
    for mass in (.1, .5, 1.):
        result[str(mass)] = {}
        for name, lo, hi in [('central', 1.8, 3.2), ('below1', 0, 1), ('available', 0, 4.1)]:
            selected = [r for r in rows if r['mass_MeV'] == mass and r['low_MeV'] >= lo and r['high_MeV'] <= hi]
            if not selected:
                continue
            metrics = {'bin_count': len(selected), 'energy_low_MeV': min(r['low_MeV'] for r in selected),
                       'energy_high_MeV': max(r['high_MeV'] for r in selected)}
            for variant in VARIANTS:
                residuals = [abs(r[variant] / r['paper_flux'] - 1) for r in selected]
                metrics[variant] = {'median_relative_error': statistics.median(residuals),
                                    'median_percent': 100 * statistics.median(residuals),
                                    'maximum_relative_error': max(residuals),
                                    'fraction_above_10percent': sum(e > .1 for e in residuals) / len(residuals)}
            result[str(mass)][name] = metrics
        # Actual assertions about the diagnostic, not assertions of reproduction.
        assert result[str(mass)]['central']['conditional']['median_relative_error'] < .1
        assert result[str(mass)]['central']['baseline'] == result[str(mass)]['central']['low_cutoff']
    for mass in ('0.5', '1.0'):
        assert result[mass]['central']['conditional']['median_relative_error'] < result[mass]['central']['baseline']['median_relative_error']
    limit = read_csv('conditional_limit.csv')[0]
    baseline = json.loads((BASE / 'output/results.json').read_text())
    assert abs(limit['epsilon_limit'] / baseline['epsilon_limit'] - 1) < 1e-5
    assert abs(limit['epsilon_limit'] / baseline['paper_epsilon_limit'] - 1) > .025
    checks = read_csv('finite_mass_checks.csv')
    assert all(r['node_relative_difference'] < 2e-5 and r['order_relative_difference'] < 2e-4 for r in checks)
    summary = {'regions': result, 'conditional_limit': limit, 'finite_mass_checks': checks,
               'choices': {'gamma_low_diagnostic_MeV': .2, 'central_low_MeV': 1.8, 'central_high_MeV': 3.2,
                           'diagnostic_tolerance': .1, 'massless_tolerance': 1e-8, 'small_mass_rate_tolerance': 1e-5,
                           'bin_quadrature_nodes': 8}}
    (OUT / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    return summary


def notices():
    """Add current findings to existing reports and preserve their historical bodies."""
    en = '''<aside id="continuation-update" class="note"><strong>Continuation:</strong>
    A conditional production normalization and a lower source cutoff closely track the Figure 1 curves.
    This is a hypothesis about an undocumented implementation choice; it changes the stated equation.
    The TEXONO limit remains unreproduced.
    <a href="report_park_texono_continuation.html">New diagnostics, plot and run receipt</a>.</aside>'''
    es = '''<aside id="continuation-update" class="note"><strong>Continuación:</strong>
    Normalizar el espectro de producción por su propia sección eficaz total y bajar el corte de la fuente
    aproxima las curvas de la figura original. Es una hipótesis sobre una elección de implementación no documentada;
    modifica la ecuación publicada. El límite de TEXONO sigue sin reproducirse.
    <a href="report_park_texono_continuation.html">Diagnóstico, figura y registro de ejecución, en inglés</a>.</aside>'''
    for lang, content in [('en', en), ('es', es)]:
        (OUT / f'notice_{lang}.html').write_text(content)
    for filename in ('report_park_texono_human.html', 'report_park_texono_agent.html',
                     'rfh_park_texono_en.html', 'rfh_park_texono_es.html'):
        path = ROOT / filename
        content = path.read_text()
        if '<aside id="continuation-update"' in content:
            start = content.index('<aside id="continuation-update"')
            end = content.index('</aside>', start) + len('</aside>')
            content = content[:start] + content[end:]
        notice = es if filename.endswith('_es.html') else en
        path.write_text(content.replace('</h1>', '</h1>' + notice, 1))


def report(summary, metadata):
    rows = []
    for mass, regions in summary['regions'].items():
        for region in ('central', 'below1', 'available'):
            if region not in regions:
                continue
            values = regions[region]
            rows.append([mass, region, values['bin_count']] +
                        [f"{values[v]['median_percent']:.2f}%" for v in VARIANTS])
    body = '''<p>The existing implementation was retained. This continuation tests why the source spectra differ,
    using the existing cross sections, digitized PDF paths and ROOT integration functions.</p>
    <div class="note"><strong>New result:</strong> a conditional normalization closely follows the finite-mass Figure 1 curves.
    Including lower-energy source photons also improves their low-energy portion. This is an implementation hypothesis,
    not evidence that Park's stated rate equation has been reproduced. The TEXONO limit still disagrees.</div>
    <h2>Normalization hypothesis</h2>
    <p>The baseline uses the ratio (dσ<sub>prod</sub>/dE′)/σ<sub>KN</sub>, with mixing factored out.
    The diagnostic replaces the denominator with σ<sub>prod,total</sub> for the same mass and incident energy.
    Its integral over outgoing energy is unity: it samples the conditional energy distribution of a produced A′,
    while dropping the finite-mass production branching suppression σ<sub>prod,total</sub>/σ<sub>KN</sub>.</p>
    <p>No amplitude or fitted scale was introduced. Reference 15 supplies the independent total-production expression
    already tested against the existing trace integral. In the massless limit the two denominators agree, as checked numerically.</p>
    <p>The additional source-cutoff diagnostic uses the lower boundary of the empirical spectrum's stated domain.
    These are alternatives to the baseline, not Danilov or later corrections. Neither is applied to the retained event calculation.</p>'''
    body += fig('continuation/figure1_diagnostic.svg',
                'Figure 1 diagnostic: PDF staircase-bin values, baseline bin averages, and the conditional-normalization/lower-cutoff hypothesis. Each panel uses the same energy and rate units as the paper.')
    body += table(['Mass [MeV]', 'Interval', 'Bins', 'Baseline', 'Lower cutoff only',
                   'Conditional only', 'Both'], rows)
    body += '''<p>Entries are median absolute relative residuals. “Central” uses the original comparison interval;
    “below1” uses available bins below one MeV; “available” uses all retained horizontal segments of the earlier PDF extraction.
    See the <a href="code/park_texono/output/continuation/summary.json">summary JSON</a> for interval boundaries,
    maximum residuals and the fraction of bins outside the diagnostic tolerance. The original extraction clips at the plot floor;
    these are not complete published event data. Bin widths come from the PDF, not an assumed Monte Carlo sample.</p>
    <p>The baseline discrepancy survives bin averaging: duplicate vertical staircase endpoints were not its cause.
    The cutoff change has no effect in the central interval. The agreement of the conditional hypothesis is a clue about
    possible event generation, but does not establish the author's algorithm. No statistical goodness-of-fit claim is possible
    without the original sample and its uncertainties; full figure reproduction is still not established.</p>'''
    eps = summary['conditional_limit']['epsilon_limit']
    body += f'''<h2>TEXONO and verification</h2><p>The small-mass conditional diagnostic gives ε₉₅={eps:.10g}.
    It agrees with the baseline within the checked tolerance and still fails the printed limit.
    A source cutoff below the deposited-energy ROI cannot explain that discrepancy in this calculation.</p>
    <p>New checks compare doubled quadrature order and reversed source/detector integration order across the finite-mass
    benchmarks. Existing finite-mass production/reference and detailed-balance checks were rerun through the included source.</p>
    <pre>{escape((OUT / 'run_receipt.txt').read_text())}</pre>
    <p>Derived here: the conditional-normalization diagnostic and bin-level residuals. Reused: the previously derived
    production/detection amplitudes and kinematics. Supplied: Park's spectrum, curves and reported limit; Reference 15's total cross section.
    ROOT supplies quadrature and graphics. Python only assembles metrics and reports.</p>
    <p>Not checked: the author's original event generator, detector response, reactor transport, raw-event likelihood,
    prompt-fission spectrum, and the three-photon loop amplitude. The next decisive input is the original source/rate implementation.</p>
    <h2>Files and execution</h2>
    <p><a href="report_park_texono_human.html">Original derivation and baseline report</a> ·
    <a href="code/park_texono/continuation.C">ROOT diagnostic macro</a> ·
    <a href="code/park_texono/output/continuation/bin_comparison.csv">Bin comparisons</a> ·
    <a href="code/park_texono/output/figures.html">Marimo export</a> ·
    <a href="papers/park_reproduction/park2017.pdf">Park, Figure 1 and Eqs. (1)–(3)</a> ·
    <a href="papers/park_reproduction/reference15.pdf">Reference 15, Appendix</a>.</p>
    <p>The macro includes the existing C++ file and runs through Cling; the baseline executable is not rebuilt.
    This uses <a href="https://root.cern.ch/manual/root_macros_and_shared_libraries/">ROOT's documented macro workflow</a>.</p>
    <pre>cd code/park_texono
root -l -b -q continuation.C &gt; output/continuation/run_receipt.txt 2&gt;&amp;1
python3 continuation_report.py</pre>
    <footer><h2>Execution record</h2><p>Model: GPT-6, as identified by the session. Reasoning effort and exact token usage:
    not exposed. Recorded diagnostic start: {metadata['diagnostic_start_utc']}; report build: {metadata['report_build_utc']}.
    Diagnostic wall time: {metadata['diagnostic_elapsed_seconds']:.1f} s. This excludes the earlier context-recovery work;
    total time across sessions is unavailable. The original report retains its original runtime record.</p>
    <p>ROOT 6.36.000 with Cling; no installation and no baseline executable rebuild. Baseline compilation was previously
    recorded with g++ 13.3.0. The new diagnostic is interpreted, not a fresh g++ compilation.</p></footer>'''
    document = f'<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Park / TEXONO continuation</title><style>{CSS}</style></head><body><main><h1>Park / TEXONO: Figure 1 normalization diagnostic</h1>{body}</main></body></html>'
    (ROOT / 'report_park_texono_continuation.html').write_text(document)


def provenance(summary, metadata):
    numbers_path = ROOT / 'provenance/numbers.json'
    numbers = json.loads(numbers_path.read_text())
    ids = []

    def record(path, value):
        if isinstance(value, dict):
            for key, item in value.items():
                record(path + '_' + key, item)
        elif isinstance(value, list):
            for key, item in enumerate(value):
                record(path + '_' + str(key), item)
        elif isinstance(value, (int, float)):
            slug = 'park_continuation_' + path
            ids.append(slug)
            numbers[slug] = {'value': value, 'statement': path.replace('_', ' '),
                'produced_by': 'code/park_texono/continuation.C; continuation_report.py::summarize',
                'from_scratch': 'Conditional source normalization and PDF-bin comparison using the existing trace implementation.',
                'from_library': 'ROOT Gauss-Legendre nodes; Python statistics.median and CSV arithmetic.',
                'choices': ['No normalization fit; preserve the baseline; compare horizontal PDF segments once; report all retained-bin and central-interval metrics.',
                            'Cutoff follows the empirical spectrum domain; conditional denominator is a hypothesis, not the stated equation.',
                            'The local residual tolerance is diagnostic, not an experimental uncertainty.']}
    record('summary', summary)
    record('diagnostic_elapsed_seconds', metadata['diagnostic_elapsed_seconds'])
    numbers_path.write_text(json.dumps(numbers, ensure_ascii=False, indent=2) + '\n')
    claims_path = ROOT / 'provenance/claims.yaml'
    claims = yaml.safe_load(claims_path.read_text())
    claim = {'id': 'park-continuation-normalization',
             'statement': 'A conditional-production normalization and lower gamma cutoff reduce Figure 1 residuals but change the stated equation; they do not establish the author implementation or resolve the TEXONO limit. Finite-mass rate convergence checks pass.',
             'evidence': ['code/park_texono/continuation.C', 'code/park_texono/continuation_report.py::summarize',
                          'code/park_texono/output/continuation/run_receipt.txt', 'papers/park_reproduction/park2017.pdf::Figure 1 and Eqs. 1-3'],
             'numbers': ids}
    claims['claims'] = [c for c in claims['claims'] if c['id'] != claim['id']] + [claim]
    paths = [str(p.relative_to(ROOT)) for p in OUT.glob('figure1_diagnostic*') if p.suffix in ('.svg', '.pdf', '.png')]
    claims['figures'] = [f for f in claims['figures'] if f['file'] not in paths]
    for path in paths:
        producer = 'code/park_texono/continuation.C::PlotComparison'
        if path.endswith('.png'):producer += '; pdftoppm::PDF preview'
        claims['figures'].append({'file': path, 'produced_by': producer,
            'shows': 'Dark-photon energy versus source rate: original PDF bin heights, unchanged baseline, and conditional-normalization/lower-cutoff hypothesis.',
            'from_scratch': 'Bin averages from reused tree-level amplitudes and a diagnostic source denominator.',
            'from_library': 'ROOT quadrature and graphics; Poppler PDF-path extraction and preview rendering.',
            'choices': ['Separate mass panels; same axes and units as the original; no fitted scale; clipped PDF bins are not complete source data.'],
            'supports': [claim['id']]})
    claims_path.write_text(yaml.safe_dump(claims, sort_keys=False, allow_unicode=True))


def vault(metadata):
    receipt = (OUT / 'run_receipt.txt').read_text().strip()
    text = f'''# Park Figure 1 normalization diagnostic

Continued [[park-texono-reproduction]] using its existing C++ implementation. Conditional production normalization plus a lower gamma cutoff closely follows the retained Figure 1 bins. This changes the published rate equation and remains a hypothesis about the author's implementation. It does not resolve the TEXONO limit. No later correction was applied.

Outputs: `report_park_texono_continuation.html`, `code/park_texono/output/continuation/summary.json`. Derived here: diagnostic normalization and bin averaging. Reused: existing amplitudes and kinematics. Supplied: [[park-2017-reactor-dark-photons]] and Reference 15 inputs.

## Run receipt

- Report date: {metadata['report_build_utc']}
- Command from `code/park_texono/`: `root -l -b -q continuation.C > output/continuation/run_receipt.txt 2>&1`, then `python3 continuation_report.py`.
- Environment: existing ROOT/Cling, Python with PyYAML; no package installs. Baseline source and executable retained.
- Actual ROOT output:

```text
{receipt}
```

Not checked: author's generator, reactor transport, detector response, raw-event likelihood, prompt-fission spectrum, decay loop. No full-figure or printed-limit reproduction claim. Original run receipt remains in [[park-texono-reproduction]].
'''
    (ROOT / 'vault/exercises/park-figure1-normalization-diagnostic.md').write_text(text)
    original = ROOT / 'vault/exercises/park-texono-reproduction.md'
    note = '\nContinuation: [[park-figure1-normalization-diagnostic]] tests a conditional production normalization without changing this baseline.\n'
    if note.strip() not in original.read_text():
        original.write_text(original.read_text() + note)


def main():
    summary = summarize()
    now = datetime.now(timezone.utc)
    metadata = {'model': 'GPT-6', 'reasoning_effort': None, 'tokens_used': None,
                'diagnostic_start_utc': START, 'report_build_utc': now.isoformat(),
                'diagnostic_elapsed_seconds': (now - datetime.fromisoformat(START)).total_seconds(),
                'timing_scope': 'From first recorded diagnostic timestamp, excluding context recovery and prior sessions.'}
    (OUT / 'runtime_metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
    preserved = ['park_texono.cpp', 'park_texono', 'matrix_element_generated.h', 'texono.conf', 'output/results.json']
    hashes = {name: hashlib.sha256((BASE / name).read_bytes()).hexdigest() for name in preserved}
    manifest = OUT / 'baseline_sha256.json'
    if manifest.exists():
        assert json.loads(manifest.read_text()) == hashes, 'Baseline changed after continuation snapshot'
    else:
        manifest.write_text(json.dumps(hashes, indent=2) + '\n')
    report(summary, metadata)
    notices()
    provenance(summary, metadata)
    vault(metadata)
    print('PASS: bin-average diagnostics; conditional massless/rate checks; finite-mass convergence receipts')
    print('Wrote continuation report, bilingual notices, provenance, and separate vault receipt; baseline retained.')
    print('UNRESOLVED: original Figure 1 implementation and printed TEXONO limit; no full reproduction claim.')


if __name__ == '__main__':
    main()
