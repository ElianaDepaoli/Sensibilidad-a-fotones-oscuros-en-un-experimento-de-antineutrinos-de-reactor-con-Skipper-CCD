# Park / TEXONO reproduction

Derived production and detection from Dirac traces; retained Park's 2/3 detection average and source normalization. No later correction is applied. The printed TEXONO limit and the full Figure 1 are not reproduced. The prompt-flux/literal-bound alternative nearly matches the printed limit, but the author's implementation is unknown.

Use `report_park_texono_human.html`, `report_park_texono_agent.html`, and `code/park_texono/`. Links: [[park-2017-reactor-dark-photons]]. All numerical outputs and caveats are recorded in `provenance/`.

## Run receipt

- Date: 2026-09-17T16:06:50.125460+00:00
- Command: `bash code/park_texono/run.sh`
- Environment: ROOT 6.36.000; Ubuntu g++ 13.3.0; Python 3.12.3 with SymPy; notebook environment `/home/eliana/.marimo/bin/python`.
- Actual symbolic output:

```text
PASS: exact explicit Dirac matrices agree; both Ward identities vanish on shell
PASS: interference traces agree; photon limit equals Klein-Nishina exactly
PASS: s/u symmetry and external electron mass shell
F(s,u)= -2*(2*M2**2*m**4 - 2*M2**2*m**2*s - 2*M2**2*m**2*u + 2*M2**2*s*u - 2*M2*m**4*s - 2*M2*m**4*u + 8*M2*m**2*s*u - 2*M2*s**2*u - 2*M2*s*u**2 - 6*m**8 + 3*m**4*s**2 + 14*m**4*s*u + 3*m**4*u**2 - m**2*s**3 - 7*m**2*s**2*u - 7*m**2*s*u**2 - m**2*u**3 + s**3*u + s*u**3)/((-m**2 + s)**2*(-m**2 + u)**2)
```

- Actual numerical output:

```text
ROOT 6.36.000; C++ 13.3.0
PASS: finite-mass trace integrals match Reference 15; production KN and detection 2/3 KN limits
N95 = 197.176
N_e = 4.681235955e+28
N(epsilon=1) = 3.282783995e+21
epsilon95 = 1.565499511e-05
Park printed epsilon95 = 2.1e-05
Paper limit comparison = FAIL (reported, no normalization fitted)
PASS: quadrature, gamma-tail, photon integral, count cap, and power scaling checks
Info in <TCanvas::Print>: SVG file output/figure1.svg has been created
Info in <TCanvas::Print>: pdf file output/figure1.pdf has been created
Info in <TCanvas::Print>: SVG file output/figure2_texono.svg has been created
Info in <TCanvas::Print>: pdf file output/figure2_texono.pdf has been created
Info in <TCanvas::Print>: SVG file output/cross_sections.svg has been created
Info in <TCanvas::Print>: pdf file output/cross_sections.pdf has been created
Info in <TCanvas::Print>: SVG file output/compton_validity.svg has been created
Info in <TCanvas::Print>: pdf file output/compton_validity.pdf has been created
Info in <TCanvas::Print>: SVG file output/figure1_comparison.svg has been created
Info in <TCanvas::Print>: pdf file output/figure1_comparison.pdf has been created
Figure 1 m=0.1 MeV: median relative residual=0.04223486988; PASS (10% diagnostic tolerance)
Figure 1 m=0.5 MeV: median relative residual=0.1141071908; FAIL (10% diagnostic tolerance)
Figure 1 m=1 MeV: median relative residual=0.3812906341; FAIL (10% diagnostic tolerance)
```

No new base-environment packages installed. Notebook export uses existing marimo. Not checked: source reactor transport, detector response, raw-event likelihood, the three-photon loop amplitude, and the prompt-fission spectrum. Complete figure and limit reproduction failed; numerical and symbolic consistency checks passed.

Continuation: [[park-figure1-normalization-diagnostic]] tests a conditional production normalization without changing this baseline.
