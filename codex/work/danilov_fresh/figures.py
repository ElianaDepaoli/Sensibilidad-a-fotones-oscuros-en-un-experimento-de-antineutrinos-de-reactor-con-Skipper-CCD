"""Fresh ROOT/C++ outputs, presented in marimo; no physics is calculated here."""
import marimo
app = marimo.App(width="full")

@app.cell
def _():
    import marimo as mo
    from pathlib import Path
    import base64
    folder = Path(__file__).resolve().parent
    return mo, folder, base64

@app.cell
def _(mo, folder, base64):
    mo.vstack([
        mo.md("# Fresh Danilov-guided TEXONO calculation\nConditional response; resonance omitted. All plots are newly computed in ROOT/C++. The decay enhancement is published input, freshly interpolated."),
        *[mo.vstack([mo.md("## " + title),mo.image("data:image/svg+xml;base64,"+base64.b64encode((folder/"output"/(name+".svg")).read_bytes()).decode())]) for name,title in [
            ("oscillation_source","Transverse oscillation source"),
            ("polarized_cross_sections","Polarized cross sections"),
            ("texono_exclusion","Conditional exclusion boundaries"),
            ("decay_enhancement","Published decay enhancement")]],
        mo.md("## Actual execution\n```text\n"+(folder/"output/run_receipt.txt").read_text()+"\n```")])
    return

if __name__ == "__main__":
    app.run()
