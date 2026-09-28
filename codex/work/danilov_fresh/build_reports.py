"""Assemble fresh results, source manifest, bilingual lectures and provenance.

This file performs formatting, not physical calculations. Physical quantities
come from calculate.cpp, experiment.conf, or explicitly identified source data.
"""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import html
import json
import re
import yaml

ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'work/danilov_fresh'
START=datetime.fromisoformat('2026-09-25T18:04:04+00:00')
PAPERS=ROOT/'papers/danilov_fresh'

def scientific(x):
    mantissa, exponent=f'{float(x):.5e}'.split('e')
    return mantissa.rstrip('0').rstrip('.')+r'\times10^{'+str(int(exponent))+'}'

def assemble():
    now=datetime.now(timezone.utc)
    rows=list(csv.DictReader((BASE/'output/exclusion.csv').open()))
    raw=(BASE/'output/benchmarks.txt').read_text()
    scalars={line.split('=')[0]:float(line.split('=')[1]) for line in raw.splitlines() if line.count('=')==1}
    numerical=(BASE/'output/run_receipt.txt').read_text()
    symbolic=(BASE/'output/algebra_receipt.txt').read_text()
    assert 'Traceback' not in symbolic and 'FAIL:' not in numerical
    macros={'Plateau':scientific(scalars['fresh_plateau']), 'SelectedCap':f"{scalars['selected_cap']:.5f}",
            'PaperRescale':scientific(scalars['paper_algebraic_rescaling']), 'MaxTau':scientific(scalars['max_decay_optical_depth'])}
    (BASE/'output/results_macros.tex').write_text('\n'.join('\\newcommand{\\'+key+'}{'+value+'}' for key,value in macros.items())+'\n')
    selected=[]
    for m in [1,10,1000,100000,999000]:
        selected.append(next(row for row in rows if float(row['mass_eV'])==m))
    table=r'\begin{center}\begin{tabular}{rrr}\toprule $M\,[\mathrm{eV}]$ & $\epsilon_{95}\ (\eta_D=0)$ & $\epsilon_{95}\ (\eta_D=\eta_A)$ \\ \midrule'+'\n'
    for r in selected:
        table+=f"{float(r['mass_eV']):g} & ${scientific(r['epsilon95'])}$ & ${scientific(r['including_selected_decays'])}$ \\\\\n"
    table+=r'\bottomrule\end{tabular}\end{center}'+'\n'
    (BASE/'output/benchmark_table.tex').write_text(table)
    for language,switch in [('en','false'),('es','true')]:
        (ROOT/f'report_danilov_fresh_{language}.tex').write_text('\\newif\\ifspanish\\spanish'+switch+'\n\\input{work/danilov_fresh/lecture.tex}\n')
    sources=[
        ('danilov2019.pdf','https://arxiv.org/abs/1804.10777','Existing source PDF; guiding Hamiltonian, probabilities and TEXONO assumptions.'),
        ('redondo2015.pdf','https://arxiv.org/abs/1501.07292','New download; Danilov Ref. 6, damping convention and homogeneous-source amplitude.'),
        ('texono2010.pdf','https://arxiv.org/abs/0911.1597','New download; Danilov Ref. 16, original exposure and event selection.'),
        ('an2013.pdf','https://arxiv.org/abs/1302.3884','Existing corpus PDF; Danilov Ref. 14, polarization and longitudinal source distinction.'),
        ('redondo_raffelt2013.pdf','https://arxiv.org/abs/1305.2920','Existing corpus PDF; Danilov Ref. 8, medium polarization formalism.'),
        ('izaguirre2015.pdf','https://arxiv.org/abs/1507.02681','Existing source PDF; survival and decay-in-detector geometry, absorption amplitude.'),
        ('mcdermott2018.pdf','https://arxiv.org/abs/1705.00619','New download; full-loop decay correction below pair threshold.'),
        ('decay_enhancement.txt','https://arxiv.org/src/1705.00619','New ancillary table; published numerical input, not calculated here.'),
        ('park2017.pdf','https://arxiv.org/abs/1705.02470','Existing source PDF; printed-bound comparison only.'),
        ('nist_csi.html','https://physics.nist.gov/cgi-bin/Star/compos.pl?matno=141','New material-data download; CsI density for decay volume.')]
    manifest=[]
    for filename,url,purpose in sources:
        p=PAPERS/filename
        assert p.exists()
        if filename.endswith('.pdf'): assert p.read_bytes().startswith(b'%PDF'),filename
        manifest.append(dict(file=str(p.relative_to(ROOT)),url=url,purpose=purpose,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
    (PAPERS/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    metadata=dict(model='GPT-6 (Codex)',reasoning_effort='not exposed',tokens_used='not exposed',
                  first_recorded_command_utc=START.isoformat(),assembled_utc=now.isoformat(),
                  elapsed_seconds_since_first_command=(now-START).total_seconds(),
                  compute_seconds=float(re.search(r'compute_seconds = (\S+)',numerical)[1]))
    (BASE/'output/metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    figs=[('oscillation_source','Transverse oscillation production'),('polarized_cross_sections','Polarization-resolved cross sections'),('texono_exclusion','Conditional TEXONO exclusion'),('decay_enhancement','Published decay correction, freshly interpolated')]
    body=['<!doctype html><html lang="en"><meta charset="utf-8"><title>Fresh Danilov-guided TEXONO calculation</title><style>body{max-width:1100px;margin:3em auto;padding:1em;font:18px/1.55 system-ui;color:#14212b}svg{width:100%;height:auto}pre{white-space:pre-wrap;font-size:13px;padding:1em;background:#eef2f6}td,th{border:1px solid #b8c0c8;padding:.5em}table{border-collapse:collapse}a{color:#13568d}</style><h1>Fresh Danilov-guided TEXONO calculation</h1>',
          '<p>New code and new plots. No earlier project calculation or curve is imported. The numerical one-loop decay enhancement is published source data; the loop amplitude itself is not rederived.</p>',
          '<p><a href="report_danilov_fresh_en.pdf">English lecture</a> · <a href="report_danilov_fresh_es.pdf">Clase en español</a> · <a href="work/danilov_fresh/output/figures.html">Marimo figure notebook</a></p>',
          f'<p>The absorption-only plateau boundary is ε₉₅ = {scalars["fresh_plateau"]:.8g} at 1 keV. This uses a one-sided Gaussian cap of {scalars["selected_cap"]:.8g} selected events, transverse detection and assumed acceptance 1/6. It is a conditional reinterpretation, not a detector-calibrated exclusion.</p>',
          '<p>Danilov’s source is coherent transverse conversion in an absorbing medium. It preserves the original photon energy. No Compton-convolved source is multiplied into it. Both Compton amplitudes are derived independently, including a new analytic longitudinal absorption trace. The transverse result has the ordinary Klein–Nishina massless limit.</p>',
          '<p>Baseline decay losses are small at the boundary, but accepted detector decays are treated separately. The main result sets their acceptance to zero; a second curve assumes the same acceptance as absorption. Resonances and roots outside the declared small-mixing domain are omitted.</p>',
          f'<p>Published-normalization comparison remains unresolved: rescaling Park’s printed bound gives {scalars["paper_algebraic_rescaling"]:.8g}; the direct source integral gives {scalars["fresh_plateau"]:.8g}. No normalization is fitted and the published figure is not claimed as reproduced.</p>',
          '<p>Not established: actual energy-dependent signal response, material profiles, resonance transport, exact nuclear mass corrections and longitudinal nuclear emission, secondary reactor cascade, or the original TEXONO likelihood. At the upper mass endpoint the source’s relativistic expansion parameter can reach about 0.11 in the ROI.</p>',
          '<table><tr><th>Mass [eV]</th><th>Absorption only</th><th>Equal decay/absorption acceptance</th></tr>']
    for r in selected: body.append(f'<tr><td>{float(r["mass_eV"]):g}</td><td>{float(r["epsilon95"]):.8g}</td><td>{float(r["including_selected_decays"]):.8g}</td></tr>')
    body.append('</table>')
    for name,title in figs:
        svg=(BASE/'output'/f'{name}.svg').read_text()
        body.extend([f'<h2>{title}</h2>',svg[svg.index('<svg'):]])
    body+=['<h2>New corpus additions and retained primary sources</h2><ul>']
    for s in manifest:body.append(f'<li><a href="{html.escape(s["url"])}">{html.escape(Path(s["file"]).name)}</a>: {html.escape(s["purpose"])}</li>')
    body+=['</ul><p>The original FRJ-1 annual report was located at the Jülich archive, but download attempts returned an HTML challenge, not a PDF. Those responses are retained as HTML, not mislabelled as papers. The source spectrum is cited through Danilov. The APS supplement endpoint returned HTTP 403; the decay data were successfully obtained from the arXiv source archive instead.</p>',
           '<h2>Reproduction commands and actual output</h2><pre>bash work/danilov_fresh/run.sh</pre>',
           '<p>ROOT 6.36.000; g++ (Ubuntu 13.3.0-6ubuntu2~24.04) 13.3.0. Configure experiment.conf for exposure, composition, density, ROI, background excess/error or a supplied event cap, efficiencies and medium inputs. No installation was needed. ROOT provides quadrature/graphics; SymPy supplies independent algebra; marimo presents new figures.</p>',
           '<pre>'+html.escape(symbolic+'\n'+numerical)+'</pre>',
           '<p>Implementation reference: <a href="https://root.cern.ch/doc/v636/classTF1.html">ROOT 6.36 TF1 documentation</a>. Source hashes, derivations and numeric provenance are stored with this run.</p>',
           '<h2>Follow-up needed for an experimental exclusion</h2><p>Obtain photon/three-photon selection response, energy migration and containment, reactor material profiles and nuclear source information. Resolve the normalization difference with the author implementation; do not absorb it into an efficiency fit. Replace the homogeneous resonance omission with transport only when the necessary geometry is available.</p>',
           '<footer><h2>Execution metadata</h2><pre>'+html.escape(json.dumps(metadata,indent=2))+'</pre><p>The elapsed wall time starts at the first recorded command and ends at report assembly; compute time is measured separately. Token and reasoning-effort telemetry were not exposed.</p></footer></html>']
    (ROOT/'report_danilov_fresh.html').write_text('\n'.join(body))
    values={**scalars,'selected_rows':selected,'full_scan':rows,'configuration':(BASE/'experiment.conf').read_text(),
            'metadata':metadata,'reported_source_constants':{'photon_normalization_per_MW':.58e18,'photon_scale_MeV':.91,'alpha_inverse':137.035999084,'electron_MeV':.510998950,'hbarc_MeV_cm':1.973269804e-11,'source_min_MeV':.2,'paper_cap_selected':195.7,'paper_cap_preselection':1174.1,'park_bound':2.1e-5,'park_cap':197.2,'source_expansion_parameters':[.11,.0011]},
            'sources':manifest}
    registry=ROOT/'provenance/numbers.json';numbers=json.loads(registry.read_text())
    key='danilov-fresh-2026-09-25'
    numbers[key]=dict(value=values,statement='Fresh transverse-source conditional TEXONO boundaries, source inputs, decay hypotheses and actual execution metadata.',produced_by='work/danilov_fresh/calculate.cpp::Limits,PhysicsChecks; physics.hpp::EventKernel,Events,Boundary; experiment.conf',from_scratch='New covariant traces, polarization-resolved C++ matrix amplitudes, homogeneous conversion derivation, event integration and boundary solver.',from_library='ROOT quadrature/graphics; SymPy identities; published reactor spectrum and TEXONO inputs; McDermott ancillary decay table.',choices=['No reuse of previous project calculations','Danilov energy-preserving source instead of Compton convolution','Two transverse incident states','Original fiducial exposure','Exact Gaussian quantile instead of rounded 1.64','Main absorption acceptance 1/6, zero decay acceptance; separate equal-acceptance curve','Benchmark plasma masses instead of reconstructed profiles','Resonance and invalid weak-mixing roots omitted','Published decay enhancement, no claim of a new loop evaluation','No fitted normalization to published exclusion'])
    registry.write_text(json.dumps(numbers,indent=2)+'\n')
    cp=ROOT/'provenance/claims.yaml';claims=yaml.safe_load(cp.read_text());claims['claims']=[x for x in claims['claims'] if x['id']!=key]
    claims['claims'].append(dict(id=key,statement='A fresh polarization-resolved calculation passes symbolic and numerical consistency checks and produces conditional TEXONO boundaries. Its direct source normalization does not reproduce a rescaling of the published bound; detector response and resonance remain unvalidated.',evidence=['work/danilov_fresh/verify_algebra.py::derive_trace,derive_absorption_polarization,derive_mixing','work/danilov_fresh/calculate.cpp::PhysicsChecks,Limits','work/danilov_fresh/output/algebra_receipt.txt','work/danilov_fresh/output/run_receipt.txt','papers/danilov_fresh/danilov2019.pdf::Eqs.6,9','papers/danilov_fresh/texono2010.pdf::Tables II,III','papers/danilov_fresh/mcdermott2018.pdf::Eqs.9,10 and ancillary table'],numbers=[key]))
    figure_files=[BASE/'output'/f'{name}.{ext}' for name,_ in figs for ext in ('pdf','svg','png')]
    figure_files += [ROOT/f'report_danilov_fresh_{lang}.pdf' for lang in ('en','es')]
    figure_files += list((BASE/'output').glob('review_*.png'))
    relpaths=[str(p.relative_to(ROOT)) for p in figure_files]
    claims['figures']=[x for x in claims['figures'] if x['file'] not in relpaths]
    for path in relpaths:
        claims['figures'].append(dict(file=path,produced_by='work/danilov_fresh/calculate.cpp::DrawSource,DrawCrossSections,Limits,DecayPlot; lecture.tex',shows='Fresh source, polarized cross sections, conditional exclusion or published decay enhancement; lecture PDFs assemble these with derivations.',from_scratch='New amplitudes, conversion and event calculations; rendering of downloaded decay input is distinguished in captions.',from_library='ROOT graphics and quadrature; published decay table; LaTeX layout',choices=['No inherited curves','Gaps for omitted resonance/invalid domain','Absorption-only and decay-acceptance hypotheses shown separately'],supports=[key]))
    cp.write_text(yaml.safe_dump(claims,sort_keys=False,allow_unicode=True))
    receipt='# Fresh Danilov-guided TEXONO sensitivity\n\nNo prior project physical calculation or plot was reused. New code: `work/danilov_fresh/`; reports: `report_danilov_fresh*`. Session skill: [[../skills/danilov-fresh/SKILL]]. Related sources: [[../papers/danilov-fresh-corpus]].\n\nNew derivations: covariant total trace, analytic longitudinal absorption, transverse density matrix, damped source integral, event rate including optional detector decays. Borrowed inputs: photon spectrum, exposure, selection hypothesis and published exact-loop enhancement. The loop itself is not rederived. Published normalization comparison fails; no tuning. Material resonance, detector response and massive/longitudinal nuclear source remain unvalidated.\n\n## Run receipt\n\nDate: '+now.isoformat()+'\nCommand: `bash work/danilov_fresh/run.sh`. Environment: ROOT 6.36.000, Ubuntu g++ 13.3.0, existing SymPy, LaTeX and marimo. No installation.\n\n```text\n'+symbolic+'\n'+numerical+'```\n'
    if (BASE/'output/final_review.txt').exists():receipt+='\n## Artifact verification\n\n'+(BASE/'output/final_review.txt').read_text()
    (ROOT/'vault/exercises/danilov-fresh-texono.md').write_text(receipt)
    corpus='# Sources for the fresh Danilov calculation\n\nGuiding paper: Danilov arXiv:1804.10777v2. New essential references downloaded: Redondo arXiv:1501.07292 (Ref. 6), Deniz/TEXONO arXiv:0911.1597 (Ref. 16). Added McDermott arXiv:1705.00619 and its ancillary enhancement table to avoid a low-mass-only decay approximation. Existing primary PDFs for An, Redondo–Raffelt, Izaguirre and Park were copied as sources, not as calculation outputs.\n\nOriginal FRJ-1 report retrieval returned an HTML challenge; it is not recorded as a downloaded PDF. The spectrum remains a Danilov-quoted empirical input. The APS decay supplement was blocked by HTTP 403; arXiv ancillary data succeeded. NIST CsI material data downloaded for decay-volume normalization. URLs, purposes and hashes: `papers/danilov_fresh/manifest.json`. Results: [[../exercises/danilov-fresh-texono]].\n'
    (ROOT/'vault/papers/danilov-fresh-corpus.md').write_text(corpus)
    print('PASS: fresh reports assembled, source PDF signatures verified, provenance and vault written')

if __name__=='__main__':assemble()
