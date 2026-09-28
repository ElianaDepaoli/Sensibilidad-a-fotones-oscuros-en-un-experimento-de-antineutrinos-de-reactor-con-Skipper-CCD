"""Presentation only; ROOT/C++ produces all physical curves."""
import marimo
app = marimo.App(width="full")

@app.cell
def _():
    import marimo as mo
    from pathlib import Path
    import base64
    root = Path(__file__).resolve().parents[2]
    return mo, root, base64

@app.cell
def _(mo, root, base64):
    paths = [
        ("Scattering cross sections", "code/park_texono/output/cross_sections.svg"),
        ("Park Figure 1: comparison, not full reproduction", "code/darkphoton_v3/output/figure1_comparison.svg"),
        ("Medium source", "code/darkphoton_v3/output/figure1_medium.svg"),
        ("Published TEXONO comparisons", "code/darkphoton_v3/output/figure2_texono_medium.svg"),
        ("Conditional limits with approximate three-photon survival", "work/reactor_sensitivity/output/conditional_limits.svg")
    ]
    mo.vstack([mo.vstack([mo.md("## " + title), mo.image("data:image/svg+xml;base64," + base64.b64encode((root / path).read_bytes()).decode())]) for title, path in paths])
    return

if __name__ == "__main__":
    app.run()
