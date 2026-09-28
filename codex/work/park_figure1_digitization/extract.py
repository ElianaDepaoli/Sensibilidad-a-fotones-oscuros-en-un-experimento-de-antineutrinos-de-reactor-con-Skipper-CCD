"""Read Figure 1 vector strokes; no production-rate calculation or fitted model."""
from pathlib import Path
import hashlib
import json
import re
import subprocess
import xml.etree.ElementTree as ET
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent


def extract():
    source = ROOT / 'papers_init/park2017.pdf'
    subprocess.run(['pdftocairo', '-f', '2', '-l', '2', '-svg',
                    str(source), str(HERE / 'page2.svg')], check=True)
    paths = [p for p in ET.parse(HERE / 'page2.svg').iter()
             if p.tag.endswith('path') and p.get('stroke')
             and len(p.get('d', '')) > 500]
    assert len(paths) == 3
    # Axis calibration reused from code/park_texono/inspect_sources.py,
    # checked against the labeled axes in the page image.
    x0, x4, y_one, decade = 380.160, 522.250, 64.363, 31.282
    energy = np.arange(430) * .01 + .005
    table = np.full((len(energy), 4), np.nan)
    table[:, 0] = energy
    counts = {}
    for p in paths:
        mass = .1 if '100%' in p.get('stroke') else (1. if p.get('stroke-dasharray') else .5)
        col = {.1: 1, .5: 2, 1.: 3}[mass]
        a, b, c, d, e, f = map(float, re.findall(r'[-+]?\d*\.?\d+', p.get('transform')))
        tokens = re.findall(r'[MLCZ]|[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', p.get('d'))
        i, previous = 0, None
        for_command = {'M': 2, 'L': 2, 'C': 6, 'Z': 0}
        while i < len(tokens):
            command = tokens[i]
            n = for_command[command]
            values = list(map(float, tokens[i+1:i+1+n]))
            i += n + 1
            if not n:
                continue
            x, y = values[-2:]
            point = (a*x+c*y+e, b*x+d*y+f)
            if command == 'L' and previous is not None:
                px, py = previous
                # Only horizontal histogram strokes inside the visible frame.
                # The bottom-frame paths encode clipping, not measured zeros.
                if abs(point[1]-py) < 1e-4 and point[0] > px and 49.438 < py < 152.60:
                    left = (px-x0)*4/(x4-x0)
                    right = (point[0]-x0)*4/(x4-x0)
                    selected = (energy >= left) & (energy < right)
                    assert np.all(np.isnan(table[selected, col]))
                    table[selected, col] = 10**((y_one-py)/decade)
            previous = point
        counts[str(mass)] = int(np.isfinite(table[:, col]).sum())
        assert counts[str(mass)] > 200
    header = '\n'.join([
        'Source: papers_init/park2017.pdf, page 2, Figure 1; digitized PDF vector strokes.',
        'Columns: E_Aprime_MeV  spectrum_m0.1MeV  spectrum_m0.5MeV  spectrum_m1.0MeV',
        'All spectra in units of 10^21 MeV^-1 s^-1, as on the plotted ordinate.',
        'Caption normalization: epsilon=1; reactor thermal power=1 GW.',
        'Energy samples are bin centers; use steps-mid and a logarithmic y axis.',
        'nan means no visible curve at that energy (including clipping); not zero.',
        'No physics calculation, smoothing, extrapolation, or normalization fit.',
        'Digitization precision is limited by the published plot and axis calibration.'
    ])
    out = HERE / 'park_figure1.txt'
    np.savetxt(out, table, fmt=['%.3f', '%.6e', '%.6e', '%.6e'], header=header)
    readback = np.loadtxt(out)
    assert readback.shape == (430, 4)
    assert np.allclose(readback, table, equal_nan=True)
    assert np.all(np.diff(readback[:, 0]) > 0)
    assert np.all(readback[:, 1:][np.isfinite(readback[:, 1:])] > 0)
    metadata = dict(source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                    rows=430, columns=4, visible_samples_by_mass_MeV=counts,
                    axis_calibration=dict(x0=x0, x4=x4, y_one=y_one, decade=decade))
    (HERE / 'receipt.json').write_text(json.dumps(metadata, indent=2)+'\n')
    print('PASS: extracted three plotted vector paths; no physics calculation.')
    print('PASS: text readback has 430 rows and exactly 4 columns; energies increase; visible spectra are positive.')
    print('Visible samples by mass (MeV): ' + json.dumps(counts))
    print('UNCHECKED: underlying simulation and values clipped outside the published axes.')


if __name__ == '__main__':
    extract()
