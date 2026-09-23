"""Marimo presentation notebook: all physical curves are computed in ROOT C++."""
import marimo

app=marimo.App(width="full")

@app.cell
def _():
    import marimo as mo
    from pathlib import Path
    import json
    import base64
    base=Path(__file__).resolve().parent
    values=json.loads((base/'output/results.json').read_text())
    return mo,base,json,base64,values

@app.cell
def _(mo,values):
    mo.md(f"""# Park / TEXONO: reproduction audit

ROOT C++ computes the cross sections and rates. The stored result is
ε₉₅ = **{values['epsilon_limit']:.8g}**, while Park prints
**{values['paper_epsilon_limit']:.8g}**. The printed limit is not reproduced.

The code keeps Park's unpolarized detection average and omits later corrections.
The original PDF curves are used only for comparison, not normalization.
""")
    return

@app.cell
def _(mo,base,base64):
    def show_svg(filename,title):
        encoded=base64.b64encode((base/'output'/filename).read_bytes()).decode()
        return mo.vstack([mo.md('## '+title),mo.image('data:image/svg+xml;base64,'+encoded)])
    mo.vstack([
        show_svg('cross_sections.svg','Derived cross sections'),
        show_svg('figure1_comparison.svg','Figure 1: calculated curves and original PDF'),
        show_svg('figure2_texono.svg','Figure 2: TEXONO component only'),
        show_svg('compton_validity.svg','Klein–Nishina versus uranium total attenuation')
    ])
    return

@app.cell
def _(mo,base):
    mo.md('## Actual numerical run\n```text\n'+(base/'output/run_receipt.txt').read_text()+'\n```')
    return

@app.cell
def _(mo,base,base64):
    diagnostic=base/'output/continuation/figure1_diagnostic.svg'
    if diagnostic.exists():
        continuation_panel=mo.vstack([
            mo.md('''## Continuation: normalization hypothesis

The conditional production distribution and lower source cutoff closely follow
the retained Figure 1 bins. This changes the stated source equation and does not
resolve the TEXONO limit. It is a hypothesis about the original implementation,
not a reproduced physical rate. The baseline above is retained.
'''),
            mo.image('data:image/svg+xml;base64,'+base64.b64encode(diagnostic.read_bytes()).decode()),
            mo.md('### Actual continuation run\n```text\n'+(base/'output/continuation/run_receipt.txt').read_text()+'\n```')
        ])
    else:
        continuation_panel=mo.md('Continuation diagnostics have not been run.')
    continuation_panel
    return

if __name__=='__main__': app.run()
