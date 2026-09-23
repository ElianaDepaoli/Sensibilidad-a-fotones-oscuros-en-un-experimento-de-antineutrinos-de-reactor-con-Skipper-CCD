"""Assemble HTML reports and provenance from actual ROOT outputs.

Presentation and bookkeeping only. All numerical physics is in ROOT C++;
symbolic algebra is in derive.py. Existing project provenance is preserved.
"""
import base64
import csv
from datetime import datetime, timezone
import hashlib
from html import escape
import json
from pathlib import Path
import platform
import subprocess
import yaml

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]
OUT=BASE/'output'
START=datetime.fromisoformat('2026-09-17T15:45:41+00:00')
CSS='''body{margin:0;background:#f7f6f2;color:#172331;font:17px/1.55 system-ui,sans-serif}main{max-width:1100px;margin:auto;padding:40px 25px 80px}h1{font-size:clamp(2rem,4vw,3.2rem);line-height:1.1}h2{margin-top:2em;border-bottom:1px solid #bbb}h3{margin-top:1.5em}a{color:#145a80}table{border-collapse:collapse;width:100%;font-size:.91rem}td,th{border:1px solid #ccc;padding:9px;text-align:left;vertical-align:top}th{background:#e6edf1}pre{overflow:auto;background:#eaf0f3;padding:15px;font-size:.8rem}code{font-size:.9em}.note{background:#fff1d9;border-left:5px solid #b77715;padding:16px}.pass{color:#126849}.fail{color:#a12d2d}.figure img{width:100%;max-height:650px}figcaption{font-size:.88rem;color:#475669}.meta{font-size:.88rem;color:#536473}.equation{padding:12px;background:#fff;overflow:auto}.scroll{overflow-x:auto}li{margin:.3em 0}'''

def load():
    """Read stored physics results and their independent-comparison summaries."""
    values=json.loads((OUT/'results.json').read_text())
    comparison=json.loads((OUT/'figure1_status.json').read_text())
    with (OUT/'compton_validity.csv').open() as stream:
        material=[{k:float(v) for k,v in row.items()} for row in csv.DictReader(stream)]
    return values,comparison,material

def fig(name,caption):
    """Embed the ROOT vector output for a standalone HTML document."""
    src=base64.b64encode((OUT/name).read_bytes()).decode()
    return f'<figure class="figure"><img alt="{escape(caption)}" src="data:image/svg+xml;base64,{src}"><figcaption>{caption}</figcaption></figure>'

def table(head,rows):
    """One source of table formatting across the reports."""
    return '<div class="scroll"><table><tr>'+''.join('<th>'+str(h)+'</th>' for h in head)+'</tr>'+''.join('<tr>'+''.join('<td>'+str(v)+'</td>' for v in row)+'</tr>' for row in rows)+'</table></div>'

def footer(now,lang='en'):
    """Expose known runtime metadata and explicitly identify unavailable fields."""
    elapsed=(now-START).total_seconds()/60
    if lang=='es':
        return f'''<footer><h2>Registro de ejecución</h2><p>Modelo: GPT-6, según la sesión.
        Nivel de razonamiento: no disponible. Uso exacto de tokens: no disponible; no se lo reemplaza por una estimación.
        Tiempo transcurrido hasta generar este informe: {elapsed:.2f} minutos, incluida la interrupción y las esperas de autorización.
        Inicio: {START.isoformat()}; generación: {now.isoformat()}.</p><p>ROOT 6.36.000;
        g++ (Ubuntu 13.3.0-6ubuntu2~24.04) 13.3.0. Entorno existente:
        /home/eliana/.marimo/bin/python, marimo 0.24.0. No se instalaron paquetes en el entorno base.</p></footer>'''
    return f'''<footer><h2>Execution record</h2><p>Model: GPT-6 (reported by the session).
    Reasoning effort: not exposed. Exact token usage: not exposed; no estimate is substituted.
    Elapsed wall time to this report build: {elapsed:.2f} minutes, including the interruption and approval waits.
    Start: {START.isoformat()}; report build: {now.isoformat()}.</p>
    <p>ROOT 6.36.000; g++ (Ubuntu 13.3.0-6ubuntu2~24.04) 13.3.0.
    Existing notebook environment: /home/eliana/.marimo/bin/python (marimo 0.24.0).
    No packages were installed in the base environment.</p></footer>'''

def page(name,title,body,now,lang='en'):
    # Preserve the continuation notice when the baseline reports are regenerated.
    notice=OUT/'continuation'/f'notice_{lang}.html'
    if notice.exists():body=notice.read_text()+body
    document=f'<!doctype html><html lang="{lang}"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>{title}</title><style>{CSS}</style></head><body><main><h1>{title}</h1>{body}{footer(now,lang)}</main></body></html>'
    (ROOT/name).write_text(document,encoding='utf-8')

def build_human(v,comparison,material,now):
    """Derivation, reproduction status, and uncertainty choices for a physicist."""
    rows=[
        ('Reactor thermal power','2900 MW','Park TEXONO input'),
        ('Baseline / target / live time','28 m / 187 kg CsI(Tl) / 160 days','Park TEXONO input'),
        ('Deposited-energy ROI','3–8 MeV','Full photon + electron containment'),
        ('Target Z / molar mass A','108 / 259.8099219 g mol⁻¹','CsI formula unit; Tl trace dopant neglected'),
        ('Electron count',f"{v['electron_count']:.9g}",'Calculated from mass, Z/A, Avogadro constant'),
        ('Event-count cap',f"1.96 × 100.6 = {v['N95']:.3f}",'Park convention; not a new likelihood fit'),
        ('Acceptance','1','Park idealized baseline; no later anti-Compton correction'),
        ('Incident γ integration','1–40 MeV','Extended to 60 MeV as a convergence check')]
    comparisons=[(f'{m} MeV',f"{100*r['median_absolute_relative_error']:.2f}%",'passes the local 10% criterion' if r['comparison_pass_10percent'] else 'fails the local 10% criterion') for m,r in comparison.items()]
    matrows=[(f"{r['E_MeV']:g}",f"{r['KN_cm2_g']:.5f}",f"{r['total_cm2_g']:.5f}",f"{r['KN_fraction_total']:.3f}",f"{r['Thomson_over_KN']:.2f}") for r in material if r['E_MeV'] in (.2,1,2,3,5,8)]
    body=f'''<p class="meta">Park 2017 · TEXONO only · ROOT/C++ calculation with symbolic checks</p>
    <div class="note"><strong>The cross sections are checked; the published TEXONO limit is not reproduced.</strong>
    Direct integration gives ε₉₅ = {v['epsilon_limit']:.8g}; Park reports {v['paper_epsilon_limit']:.3g}.
    Figure 1 also disagrees at finite mass. No normalization was adjusted to force agreement.</div>
    <p><a href="rfh_park_texono_en.html">Short English follow-up</a> · <a href="rfh_park_texono_es.html">Resumen en español</a> ·
    <a href="report_park_texono_agent.html">Technical handoff</a> · <a href="code/park_texono/output/figures.html">Exported marimo notebook</a></p>
    <h2>Inputs and scope</h2>{table(['Input','Value','Meaning'],rows)}
    <p>The run uses Park's single-interaction source normalization, point-source dilution, stationary free electrons,
    the unpolarized detection factor 2/3, and his event-count prescription. No Danilov or later correction is implemented.
    The source flux is an ε² coefficient: evaluating its coefficient at ε=1 follows Park's plotting convention;
    it is not a physical prediction at strong mixing. General theory constraints in Park's introduction and the non-TEXONO
    exclusion regions are outside this TEXONO calculation.</p>
    <h2>Production cross section derived from the diagrams</h2>
    <p>Write p for the incoming electron, k for the incoming photon, q for the outgoing A′, and p′=p+k−q.
    Let m=mₑ, M=m<sub>A′</sub>, s=(p+k)², u=(p−q)². From the interaction −eA·J−εeA′·J:</p>
    <div class="equation">𝓜 = εe² ū(p′) [γ<sup>μ</sup>(p̸+k̸+m)γ<sup>ν</sup>/(s−m²)
    + γ<sup>ν</sup>(p̸−q̸+m)γ<sup>μ</sup>/(u−m²)] u(p) ε′*<sub>μ</sub> ε<sub>ν</sub>.</div>
    <p>The two terms are the electron propagators in the s and u channels. Average over the initial electron spin
    and photon polarization, and sum over final states. The current is conserved on shell, so the q<sub>μ</sub>q<sub>ν</sub>/M²
    term in the massive polarization sum vanishes. The script expands the Dirac trace using anticommutators; it does not call
    a matrix-element library.</p>
    <p>Define a=s−m², b=u−m² and Q=M². The resulting dimensionless factor is</p>
    <div class="equation">⟨|𝓜|²⟩ = ε²e⁴ F,<br>
    F = −2 [a/b + b/a + 2Q²/(ab) − 2Qm²(a⁻²+b⁻²)<br>
    − (2Q+4m²)(a⁻¹+b⁻¹) − 4m⁴(a⁻¹+b⁻¹)²].</div>
    <p>In the electron-rest laboratory frame, with incident energy E and outgoing dark-photon energy E′:</p>
    <div class="equation">s=m²+2mE, u=m²+M²−2mE′,<br>
    dσ<sub>prod</sub>/dE′ = ε²e⁴ F / (32πmE²).</div>
    <p>Natural units are used in this equation; multiplication by (ℏc)² converts MeV⁻³ to cm²/MeV.
    Production requires E ≥ M+M²/(2m). With λ=[s−(m+M)²][s−(m−M)²], the exact lab endpoints are</p>
    <div class="equation">E′<sub>±</sub> = [(E+m)(s+M²−m²) ± E√λ]/(2s).</div>
    <p><strong>Reference 15 is not needed to obtain the amplitude.</strong> It is used as an independent test:
    the integral of the derived differential cross section agrees with its vector total cross section, Appendix A.1–A.3,
    at the checked incident energies and finite masses. At M→0, SymPy verifies the exact Klein–Nishina expression.
    Explicit Dirac matrices at an independent massive noncollinear point also agree, and both Ward contractions vanish.</p>
    <h2>Detection derived by reversing the reaction</h2>
    <p>For A′e→γe, let E′ be the incoming A′ energy and w the outgoing photon energy.
    The squared amplitude is the same invariant function, with s=m²+M²+2mE′ and u=m²−2mw.
    Park averages over three incoming A′ polarizations instead of two incoming photon polarizations, hence</p>
    <div class="equation">dσ<sub>det</sub>/dw = (2/3) ε²e⁴ F / [32πm(E′²−M²)],<br>
    w<sub>±</sub> = (M²+2mE′)/[2(m+E′∓√(E′²−M²))].</div>
    <p>The flux factor and lab Jacobian follow from dσ/dt=⟨|𝓜|²⟩/(16πλ) and |dt/dw|=2m.
    <strong>Reference 20 is not needed as a numerical dependency.</strong> Its two amplitudes and invariant phase-space
    expression provide a check on the reversed reaction. In the light-mass limit, the integrated detection cross section
    equals (2/3)ε²σ<sub>KN</sub>. The 2/3 factor is deliberately retained as requested.</p>
    <div class="note">A printed-bound ambiguity matters. Park's text below Eq. (7) gives the lower photon energy
    as 2mE′/(m+2E′). Exact Compton kinematics gives mE′/(m+2E′). Both interpretations are evaluated below.
    This is a check of Park's own equations, not a later in-medium correction.</div>
    {fig('cross_sections.svg','Derived production and unpolarized detection total cross sections per electron. ε² is factored out.')}
    <h2>From source to TEXONO events</h2>
    <div class="equation">dṄ<sub>γ</sub>/dE = 0.58×10¹⁸ (P/MW) exp[−E/(0.91 MeV)] s⁻¹ MeV⁻¹,<br>
    dṄ<sub>A′</sub>/dE′ = ∫ (dσ<sub>prod</sub>/dE′)/σ<sub>KN</sub>(E) × (dṄ<sub>γ</sub>/dE) dE,<br>
    N = NₑT η/(4πR²) ∫<sub>ROI</sub> dE′ (dṄ<sub>A′</sub>/dE′) ∫<sub>w−</sub><sup>w+</sup> dw dσ<sub>det</sub>/dw.</div>
    <p>The ROI is on the total deposited energy E′=w+Tₑ when both secondaries are contained.
    A and Z describe the detector's electron count. They cancel from the production ratio when both numerator and denominator
    are multiplied by the same number of free electrons. The gamma spectrum is an empirical input, not derived here.</p>
    {table(['Calculation','This run','Status'],[
      ('∫γ spectrum above 1 MeV at 1 GW',f"{v['photon_integral_above_1']:.8g} s⁻¹",'Matches 1.76×10²⁰ within its printed rounding'),
      ('Prompt-photon number above 1 MeV','6.82×10¹⁹ s⁻¹','Quoted input; underlying fission spectrum is not supplied'),
      ('Flux ratio / fourth-root limit factor',f"{v['flux_ratio']:.6f} / {v['epsilon_flux_factor']:.6f}",'Recovers the reported approximate 30% weakening'),
      ('Eq. (5) prefactor from Eq. (4)',f"{v['decay_prefactor_from_eq4_m']:.6f} m",'Matches the printed 505 m'),
      ('N / ε⁴ at M=0.001 MeV',f"{v['coefficient_eps4']:.9g}",'Independent integral; gamma tail and order checks pass'),
      ('N at Park’s printed ε',f"{v['events_at_paper_epsilon']:.6f}",f"Does not equal the cap {v['N95']:.3f}"),
      ('ε limit at M=0.001 MeV',f"{v['epsilon_limit']:.9g}",'Fails the printed-limit comparison'),
      ('ε limit at M=0.1, 0.5, 1 MeV',f"{v['epsilon_mass_01']:.7g}; {v['epsilon_mass_05']:.7g}; {v['epsilon_mass_1']:.7g}",'Finite-mass result, not a reproduced flat limit')])}
    <p>Eq. (4), Γ<sub>3γ</sub>=2.16×10⁻¹⁶ e⁴ε²M⁹/m⁸, is adopted from Park; the electron-loop decay calculation
    is not rederived. From it, L=(ℏc)p/(MΓ), and for E′≫M the coefficient in Eq. (5) follows.
    For E′=3 MeV, M=0.1 MeV and the calculated small-mass ε limit, L={v['decay_length_E3_M01_epslimit_m']:.6g} m.
    This verifies negligible decay in that benchmark; it does not validate the low-mass width expansion near the pair threshold.</p>
    <h2>Why the published limit still differs</h2>
    {table(['Photon normalization','Detection lower bound','ε₉₅'],[
       ('Eq. (3)','Exact kinematics',f"{v['epsilon_limit']:.9g}"),
       ('Eq. (3)','Literal printed Eq. (7)',f"{v['epsilon_printed_bound']:.9g}"),
       ('Quoted prompt flux / same assumed shape','Exact kinematics',f"{v['epsilon_prompt_flux']:.9g}"),
       ('Quoted prompt flux / same assumed shape','Literal printed Eq. (7)',f"{v['epsilon_prompt_flux_printed_bound']:.9g}")])}
    <p>The last case is close to 2.1×10⁻⁵. It changes the normalization that Park says was used, so it does not reproduce
    the stated Eq. (3) calculation. <strong>The hypothesis is that an undocumented normalization or integration choice
    explains the difference.</strong> The author's code and likelihood are unavailable here. A fitted acceptance would be
    {v['implicit_acceptance_to_match']:.6f}; this is only the algebraic factor needed to force agreement and is not applied.</p>
    <h2>Figure reproduction</h2>
    {fig('figure1_comparison.svg','Figure 1: ROOT calculation from Park’s gamma input versus vector coordinates extracted from the supplied PDF. Dotted curves are the published data, not new calculations.')}
    {table(['A′ mass','Median absolute fractional residual','Comparison'],comparisons)}
    <p>The residual metric covers 1.8–3.2 MeV; its 10% tolerance is an explicitly chosen diagnostic, not a paper uncertainty.
    It avoids the threshold edges and noisy high-energy bins. The PDF staircase coordinates are preserved in CSV.
    No published Monte Carlo seed or event sample was supplied. The full figure is not reproduced within this criterion;
    passing the low-mass central interval does not certify the entire curve.</p>
    {fig('figure2_texono.svg','TEXONO component of Figure 2 on logarithmic mass/mixing axes. The red reference is Park’s reported line; the blue curve is the calculated result. Other experiments and NEOS are omitted by scope.')}
    <p>The flat original line can be drawn from the printed number, but that is a reference overlay. The independent calculation
    does not reproduce it. The ultra-low-mass extension retains Park's assumptions and is not advertised as a modern physical exclusion.</p>
    <h2>Where “Compton total” is not an acceptable substitution</h2>
    <p><strong>Classical Thomson:</strong> σ<sub>T</sub>=8πrₑ²/3 requires E≪mₑ. It cannot replace the relativistic
    cross section in the MeV range. <strong>Klein–Nishina:</strong> it describes scattering by a free electron and is appropriate
    for that channel, but it is not the total photon interaction cross section of uranium or CsI.
    At low energy, atomic binding, photoelectric absorption and coherent scattering matter. Pair production is allowed above
    2mₑ and becomes increasingly important at higher energy; there is no universal material-independent boundary.</p>
    {table(['Eγ [MeV]','Free-electron KN × NₐZ/A [cm²/g]','NIST uranium total [cm²/g]','KN proxy / total','Thomson / KN'],matrows)}
    <p>For uranium, the free-electron KN proxy is only about 0.60 of the total at 3 MeV and 0.29 at 8 MeV.
    Thus σ<sub>tot</sub>=σ<sub>KN</sub> is not a precision approximation across the TEXONO ROI.
    No sampled point from 0.1 to 10 MeV satisfies a 10% total-attenuation criterion. This statement is limited to the tabulated
    uranium diagnostic; the actual reactor composition requires its own transport input.
    The baseline code retains Park's substitution because this task asks for his calculation.</p>
    {fig('compton_validity.svg','Diagnostic only: free-electron Klein–Nishina attenuation versus NIST total uranium attenuation. These data do not modify the event calculation.')}
    <h2>What has and has not been established</h2>
    <p>Derived here: both tree amplitudes, trace contractions, differential cross sections, lab endpoints, electron count,
    flux integrals, event convolution, and inversion for ε. Recalled or supplied: QED Feynman rules, physical constants,
    reactor and detector inputs, the gamma parameterization, the three-photon width coefficient, and the reported limits.
    ROOT provides quadrature nodes and plotting; SymPy provides symbolic simplification.</p>
    <p>Checked: Ward identities and explicit Dirac matrices at a noncollinear point; symbolic photon limit and crossing symmetry;
    finite-mass total production against Reference 15; detection’s 2/3 limit; convergence, swapped integration order,
    power scaling, photon count and decay prefactor. Not checked: reactor-specific transport, detector response and containment,
    the three-photon loop amplitude, the quoted prompt-fission spectrum, or a statistical likelihood from event data.
    The complete Figure 1 and Park's TEXONO limit failed reproduction. No later correction is included.</p>
    <h2>Sources and runnable files</h2>
    <ul><li><a href="papers/park_reproduction/park2017.pdf">Park, supplied paper</a> · <a href="https://arxiv.org/abs/1705.02470">arXiv</a>.</li>
    <li><a href="papers/park_reproduction/reference15.pdf">Gondolo–Raffelt, Reference 15</a> · <a href="https://arxiv.org/abs/0807.2926">arXiv</a>.</li>
    <li><a href="papers/park_reproduction/reference20.pdf">Izaguirre–Krnjaic–Pospelov, Reference 20</a> · <a href="https://arxiv.org/abs/1507.02681">arXiv</a>.</li>
    <li><a href="https://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z92.html">NIST uranium table</a> · <a href="https://physics.nist.gov/PhysRefData/XrayMassCoef/tab1.html">NIST Z/A</a> · <a href="code/park_texono/nist_uranium.csv">stored numeric input</a>.</li>
    <li><a href="https://root.cern.ch/doc/v636/classROOT_1_1Math_1_1GaussLegendreIntegrator.html">ROOT numerical integration documentation</a>.</li>
    <li><a href="code/park_texono/park_texono.cpp">ROOT C++</a> · <a href="code/park_texono/derive.py">symbolic derivation</a> · <a href="code/park_texono/texono.conf">experiment configuration</a>.</li></ul>
    <pre>bash code/park_texono/run.sh
PYTHONDONTWRITEBYTECODE=1 /home/eliana/.marimo/bin/marimo export html code/park_texono/figures.py -o code/park_texono/output/figures.html
python3 code/park_texono/make_reports.py</pre>'''
    page('report_park_texono_human.html','Park’s TEXONO calculation: derivation and reproduction audit',body,now)

def build_followups(v,now):
    en=f'''<p>The requested ROOT/C++ calculation is implemented. The cross-section checks pass.
    Park’s published TEXONO limit is not reproduced.</p>
    <p>The direct result is ε₉₅={v['epsilon_limit']:.7g}; the paper prints {v['paper_epsilon_limit']:.3g}.
    The original Figure 1 also differs from the finite-mass calculation.</p>
    <ul><li>Production and detection were derived from the QED diagrams. Reference 15 checks the production integral;
    Reference 20 checks the detection amplitudes and phase space.</li>
    <li>Park’s 2/3 detection factor is retained. No Danilov or later correction is used.</li>
    <li>Park’s printed lower outgoing-photon bound differs from Compton kinematics. Testing both bounds exposes part of the discrepancy.</li>
    <li>Using the quoted prompt-photon normalization and the literal bound gives {v['epsilon_prompt_flux_printed_bound']:.7g}.
    This is a diagnostic coincidence; the author's actual normalization remains unverified.</li>
    <li>Thomson is unsuitable for MeV photons. Klein–Nishina is valid for free-electron scattering, but does not equal total attenuation in reactor material.</li></ul>
    <p>The next useful input is Park’s original rate code, including the source normalization and detector integration limits.
    Until then, treat the unexplained normalization as a hypothesis.</p>
    <p><a href="report_park_texono_human.html">Derivation, figures and checks</a> · <a href="report_park_texono_agent.html">Technical handoff</a> ·
    <a href="code/park_texono/output/figures.html">Notebook export</a> · <a href="rfh_park_texono_es.html">Español</a>.</p>'''
    es=f'''<p>El cálculo pedido está implementado en ROOT/C++. Los controles de las secciones eficaces pasan.
    El límite de TEXONO publicado por Park no se reproduce.</p>
    <p>El resultado directo es ε₉₅={v['epsilon_limit']:.7g}; el artículo indica {v['paper_epsilon_limit']:.3g}.
    La figura original también difiere del cálculo a masa finita.</p>
    <ul><li>Las secciones eficaces de producción y detección se derivaron de los diagramas de QED.
    La referencia 15 verifica la integral de producción; la referencia 20 permite verificar las amplitudes y el espacio de fases de detección.</li>
    <li>Se conserva el factor 2/3 de Park para detección. No se aplican correcciones de Danilov ni posteriores.</li>
    <li>El límite inferior impreso para la energía del fotón saliente difiere del permitido por la cinemática de Compton.
    Se evaluaron ambas opciones.</li>
    <li>Usar la normalización citada para fotones prontos junto con el límite inferior literal da {v['epsilon_prompt_flux_printed_bound']:.7g}.
    La coincidencia sirve como diagnóstico; no demuestra qué normalización usó el autor.</li>
    <li>Thomson no sirve para fotones de MeV. Klein–Nishina describe la dispersión en electrones libres,
    pero no equivale a la atenuación total del material del reactor.</li></ul>
    <p>El siguiente dato útil es el código original de Park, con la normalización de la fuente y los límites de integración del detector.
    Hasta entonces, la explicación de la diferencia de normalización es una hipótesis.</p>
    <p><a href="report_park_texono_human.html">Derivación, figuras y controles, en inglés</a> ·
    <a href="code/park_texono/output/figures.html">Notebook exportado</a> · <a href="rfh_park_texono_en.html">English</a>.</p>'''
    page('rfh_park_texono_en.html','Park / TEXONO: follow-up',en,now)
    page('rfh_park_texono_es.html','Park / TEXONO: seguimiento',es,now,'es')

def build_agent(v,now):
    rows=[('derive.py','Dirac trace recursion, symbolic photon limit, Ward identities, header generation'),
          ('matrix_element_generated.h','Generated compact invariant squared amplitude'),
          ('park_texono.cpp','All numerical physics and ROOT plots; comments at each function'),
          ('texono.conf','Power, baseline, exposure, target mass, A/Z, ROI, counts, uncertainty and acceptance'),
          ('inspect_sources.py','PDF vector path extraction; fixed, documented axis calibration'),
          ('figures.py','Marimo display of ROOT outputs; no replacement physics calculation'),
          ('make_reports.py','HTML generation, provenance merge, and vault receipt')]
    body=f'''<p>The remaining issue is the paper comparison, not compilation or quadrature.
    The result ε₉₅={v['epsilon_limit']:.10g} does not match the printed {v['paper_epsilon_limit']:.3g}.</p>
    {table(['File in code/park_texono/','Responsibility'],rows)}
    <p>Change experiment inputs in texono.conf. For a compound, A is grams per mole of formula units and Z is the electron count
    per formula unit. For a mixture, use the weighted electron-per-mass ratio. baseline_m is metres; the code converts it to centimetres.
    sigma_events is the total event-count uncertainty. background=observed preserves Park's uncertainty-only convention;
    a different background enables the documented positive-residual count cap. It is not a Poisson likelihood.</p>
    <p>efficiency is a constant full-containment acceptance. A response matrix or energy-dependent efficiency is not implemented.
    The numerical implementation is intended for the paper's sub-MeV masses and MeV energies; porting it to lower energies
    requires atomic physics and a validated source spectrum. NIST data are diagnostic only.</p>
    <h2>Numerical consistency</h2>
    <p>Doubling quadrature nodes changes the event coefficient by {v['quadrature_relative_difference']:.6g} fractionally.
    Swapping the order of source and detector integration changes it by {v['integration_order_relative_difference']:.6g}.
    The incident gamma-tail extension also passes its stored tolerance.</p>
    <pre>{escape((OUT/'symbolic_receipt.txt').read_text())}</pre><pre>{escape((OUT/'run_receipt.txt').read_text())}</pre>
    <h2>Continuation constraints</h2><ul>
    <li>Do not overwrite the original corpus, reports, or previous provenance records.</li>
    <li>Do not apply a fitted acceptance to claim reproduction. The prompt-flux/literal-bound diagnostic is not the stated Eq. (3) computation.</li>
    <li>Do not claim the complete exclusion compilation was reproduced. The requested scope is TEXONO.</li>
    <li>The photon source and the three-photon width are literature inputs. A tree-level derivation does not verify them.</li>
    <li>Existing marimo environment: /home/eliana/.marimo/bin/python. No PyMuPDF dependency is needed; pdftocairo supplies SVG paths.</li></ul>
    <p><a href="report_park_texono_human.html">Complete human report</a> · <a href="code/park_texono/output/results.json">Machine-readable results</a> ·
    <a href="provenance/claims.yaml">Claims and figure provenance</a>.</p>'''
    page('report_park_texono_agent.html','Park / TEXONO: implementation handoff',body,now)

def record_provenance(v,comparison,material,now):
    """Merge this task's assertions without replacing earlier project records."""
    numberfile=ROOT/'provenance/numbers.json';claimfile=ROOT/'provenance/claims.yaml'
    numbers=json.loads(numberfile.read_text())
    claims=yaml.safe_load(claimfile.read_text())
    ids=[]
    def record(key,value,statement,producer='code/park_texono/park_texono.cpp::Results',choices=None):
        slug='park_impl_'+key;ids.append(slug)
        numbers[slug]={'value':value,'statement':statement,'produced_by':producer,
                      'from_scratch':'See the named function: derived cross sections and rates; cited physical constants and experimental inputs are not new measurements.',
                      'from_library':'ROOT quadrature nodes and plotting; SymPy simplification for symbolic identities; Python standard library for report arithmetic.',
                      'choices':choices or ['Retain Park assumptions, report failed matches, and do not fit an overall normalization.']}
    for k,value in v.items():
        if isinstance(value,(int,float)) and not isinstance(value,bool):
            record(k,value,k.replace('_',' ')+'; actual run value in output/results.json')
    for line in (BASE/'texono.conf').read_text().splitlines():
        if '=' not in line or line.startswith('#'):continue
        k,val=line.split('=',1);record('input_'+k,float(val),'Configuration input '+k,'code/park_texono/texono.conf::'+k)
    for m,r in comparison.items():
        record('fig1_residual_'+m,r['median_absolute_relative_error'],'Median absolute relative residual to original Figure 1 at mass '+m+' MeV','code/park_texono/park_texono.cpp::CompareFigure1')
        record('fig1_percent_'+m,100*r['median_absolute_relative_error'],'Same residual expressed as percent','code/park_texono/make_reports.py::build_human')
    for row in material:
        for k,value in row.items():
            record('material_'+str(row['E_MeV'])+'_'+k,value,'Uranium diagnostic '+k,'code/park_texono/park_texono.cpp::MaterialDiagnostic',
                   ['Total attenuation taken from NIST elemental uranium table; KN proxy uses free electrons and NIST Z/A; not used in the baseline yield.'])
    constants={'electron_mass_MeV':.510998950,'alpha_inverse':137.035999084,'hbarc_MeV_cm':1.973269804e-11,
      'Avogadro':6.02214076e23,'gamma_normalization_MW':.58e18,'gamma_scale_MeV':.91,'uranium_Z_over_A':.38651,
      'three_gamma_width_coefficient':2.16e-16,'paper_flux_above1':1.76e20,'paper_N95_rounded':197.2,
      'figure1_tolerance':.1,'figure1_comparison_low_MeV':1.8,'figure1_comparison_high_MeV':3.2,
      'comparison_mass_MeV':.001,'gamma_tail_test_MeV':60,'figure2_mass_low_eV':1e-18,'figure2_mass_high_eV':1e6,
      'quadrature_relative_tolerance':2e-5,'order_relative_tolerance':2e-4,'paper_limit_tolerance':.025,
      'decay_benchmark_E_MeV':3,'decay_benchmark_mass_MeV':.1,'elapsed_wall_minutes':(now-START).total_seconds()/60}
    for k,value in constants.items():record(k,value,k.replace('_',' '),'code/park_texono/park_texono.cpp; derive.py; make_reports.py::explicit convention/input')
    taskclaims=[
      {'id':'park-impl-cross-sections','statement':'The two tree-level cross sections were derived from Dirac traces; symbolic photon limits, explicit gamma matrices, Ward identities, and finite-mass reference-production integrals pass.',
       'evidence':['code/park_texono/derive.py::derive','code/park_texono/derive.py::check_explicit_matrices','code/park_texono/park_texono.cpp::Tests','papers/park_reproduction/reference15.pdf::Appendix A.1-A.3','papers/park_reproduction/reference20.pdf::Eqs. 11-13'], 'numbers':[]},
      {'id':'park-impl-limit-failure','statement':'Independent integration using Park source normalization and exact free-electron kinematics does not reproduce the printed TEXONO limit. The diagnostic prompt-flux/literal-bound coincidence does not establish the original implementation.',
       'evidence':['code/park_texono/park_texono.cpp::Results','code/park_texono/output/results.json','papers/park_reproduction/park2017.pdf::Eq. 7 and TEXONO paragraph'], 'numbers':ids},
      {'id':'park-impl-figure-failure','statement':'The original Figure 1 vector paths were extracted and compared without normalization fitting. The full finite-mass figure fails the selected local residual criterion; only the low-mass central interval passes.',
       'evidence':['code/park_texono/inspect_sources.py::extract_figure1','code/park_texono/park_texono.cpp::CompareFigure1','code/park_texono/output/figure1_status.json'], 'numbers':[]},
      {'id':'park-impl-material','statement':'Free-electron Klein-Nishina cannot stand in for precision total uranium attenuation in the TEXONO energy range; the diagnostic is not incorporated into Park event rates.',
       'evidence':['code/park_texono/park_texono.cpp::MaterialDiagnostic','code/park_texono/nist_uranium.csv','https://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z92.html'], 'numbers':[]},
      {'id':'park-impl-environment','statement':'The calculation ran with the requested ROOT and compiler; the existing marimo environment was found under /home/eliana/.marimo and used for the figure notebook.',
       'evidence':['code/park_texono/output/run_receipt.txt','code/park_texono/output/figures.html'], 'numbers':[]}]
    claims['claims']=[c for c in claims.get('claims',[]) if not c['id'].startswith('park-impl-')]+taskclaims
    figrecords=[]
    for path in sorted(list(OUT.glob('*.svg'))+list(OUT.glob('*.pdf'))+list(OUT.glob('*.png'))+list((ROOT/'papers/park_reproduction').glob('*.pdf'))):
        rel=str(path.relative_to(ROOT));source='papers/' in rel or path.name.startswith('park_page2')
        if source:
            producer='pdftocairo / pdftoppm::render supplied Park PDF' if path.name.startswith('park_page2') else 'local source preservation::copy pre-existing PDF without modification'
            shows='Preserved source document or rendered original paper page. This is not a new calculation or a reproduced figure.'
            support=['park-impl-cross-sections','park-impl-figure-failure']
        else:
            stem=path.stem.removesuffix('_preview')
            fname={'figure1':'SaveFigure1','figure1_comparison':'CompareFigure1','figure2_texono':'SaveFigure2','compton_validity':'MaterialDiagnostic','cross_sections':'SaveCrossSections'}.get(stem,'presentation')
            producer='code/park_texono/park_texono.cpp::'+fname
            shows={'figure1':'Dark-photon energy in MeV versus source rate per MeV-second at the paper reference power and stripped mixing.',
              'figure1_comparison':'Calculated source spectra and original PDF curve coordinates on the same axes.',
              'figure2_texono':'Log dark-photon mass versus log epsilon: calculated TEXONO result and reported reference line; disagreement is explicit.',
              'compton_validity':'Photon energy versus mass attenuation in uranium; free-electron KN compared to NIST total.',
              'cross_sections':'Incident energy versus production and detection cross sections per electron with mixing factored out.'}.get(stem,'Diagnostic report rendering.')
            if path.suffix=='.png':producer+='; pdftoppm::render PDF preview'
            support=['park-impl-material'] if path.stem=='compton_validity' else ['park-impl-cross-sections','park-impl-limit-failure','park-impl-figure-failure']
        figrecords.append({'file':rel,'produced_by':producer,'shows':shows,
                          'from_scratch':'None for source PDFs; our QED calculation for ROOT curves. See producing function.',
                          'from_library':'ROOT graphics for calculated plots; Poppler for source rendering.',
                          'choices':['Keep original source data distinct from calculated curves; disclose all failed comparisons.'], 'supports':support})
    oldfigs=[f for f in claims.get('figures',[]) if not f['file'].startswith(('code/park_texono/','papers/park_reproduction/'))]
    claims['figures']=oldfigs+figrecords
    numberfile.write_text(json.dumps(numbers,ensure_ascii=False,indent=2)+'\n')
    claimfile.write_text(yaml.safe_dump(claims,sort_keys=False,allow_unicode=True))
    manifest={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((ROOT/'papers/park_reproduction').glob('*.pdf'))}
    (OUT/'source_sha256.json').write_text(json.dumps(manifest,indent=2)+'\n')

def vault_receipt(now):
    """Leave the actual checks and unresolved scientific issue for the next run."""
    text=f'''# Park / TEXONO reproduction

Derived production and detection from Dirac traces; retained Park's 2/3 detection average and source normalization. No later correction is applied. The printed TEXONO limit and the full Figure 1 are not reproduced. The prompt-flux/literal-bound alternative nearly matches the printed limit, but the author's implementation is unknown.

Use `report_park_texono_human.html`, `report_park_texono_agent.html`, and `code/park_texono/`. Links: [[park-2017-reactor-dark-photons]]. All numerical outputs and caveats are recorded in `provenance/`.

## Run receipt

- Date: {now.isoformat()}
- Command: `bash code/park_texono/run.sh`
- Environment: ROOT 6.36.000; Ubuntu g++ 13.3.0; Python {platform.python_version()} with SymPy; notebook environment `/home/eliana/.marimo/bin/python`.
- Actual symbolic output:

```text
{(OUT/'symbolic_receipt.txt').read_text().strip()}
```

- Actual numerical output:

```text
{(OUT/'run_receipt.txt').read_text().strip()}
```

No new base-environment packages installed. Notebook export uses existing marimo. Not checked: source reactor transport, detector response, raw-event likelihood, the three-photon loop amplitude, and the prompt-fission spectrum. Complete figure and limit reproduction failed; numerical and symbolic consistency checks passed.
'''
    if (ROOT/'vault/exercises/park-figure1-normalization-diagnostic.md').exists():
        text+='\nContinuation: [[park-figure1-normalization-diagnostic]] tests a conditional production normalization without changing this baseline.\n'
    (ROOT/'vault/exercises/park-texono-reproduction.md').write_text(text)

def main():
    v,comparison,material=load();now=datetime.now(timezone.utc)
    build_human(v,comparison,material,now);build_followups(v,now);build_agent(v,now)
    record_provenance(v,comparison,material,now);vault_receipt(now)
    (OUT/'runtime_metadata.json').write_text(json.dumps({'model':'GPT-6','reasoning_effort':None,'tokens_used':None,
       'start_utc':START.isoformat(),'report_build_utc':now.isoformat(),'elapsed_seconds':(now-START).total_seconds()},indent=2)+'\n')
    print('Wrote human and technical reports, English/Spanish follow-ups, provenance, and vault receipt.')

if __name__=='__main__': main()
