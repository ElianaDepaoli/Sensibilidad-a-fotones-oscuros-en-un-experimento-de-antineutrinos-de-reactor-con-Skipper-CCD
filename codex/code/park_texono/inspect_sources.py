"""Extract original PDF vector-curve coordinates; do not use them to fit rates."""
import csv
import json
import re
import xml.etree.ElementTree as ET
from pathlib import Path

BASE=Path(__file__).resolve().parent

def extract_figure1():
    """PDF-to-SVG preserves actual staircase coordinates, without raster picking."""
    tree=ET.parse(BASE/'output/park_page2.svg')
    paths=[e for e in tree.iter() if e.tag.endswith('path') and e.get('stroke') and len(e.get('d',''))>500]
    assert len(paths)==3
    out=[]
    for path in paths:
        mass=.1 if '100%' in path.get('stroke') else (1. if path.get('stroke-dasharray') else .5)
        a,b,c,d,e,f=map(float,re.findall(r'[-+]?\d*\.?\d+',path.get('transform')))
        tokens=re.findall(r'[MLCZ]|[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?',path.get('d'))
        i=0; points=[]
        while i<len(tokens):
            command=tokens[i];i+=1
            n={'M':2,'L':2,'C':6,'Z':0}[command]
            values=list(map(float,tokens[i:i+n]));i+=n
            if n:
                x,y=values[-2:];points.append((a*x+c*y+e,b*x+d*y+f))
        # Tick calibration from this PDF: x=0 at 380.160, x=4 at 522.250.
        # y=1 at 64.363; logarithmic decades are 31.282 PDF points apart.
        for x,y in points:
            energy=(x-380.160)/(142.090/4)
            flux=10**((64.363-y)/31.282)*1e21
            if .15<energy<4.1 and y<152.60:
                out.append((energy,mass,flux))
    with (BASE/'output/park_figure1_digitized.csv').open('w') as stream:
        writer=csv.writer(stream);writer.writerow(['EA_MeV','mass_MeV','flux_per_MeV_s']);writer.writerows(out)
    for mass in (.1,.5,1.):
        row=min((r for r in out if r[1]==mass),key=lambda r:abs(r[0]-2))
        print('PDF vector point near 2 MeV:',row)

if __name__=='__main__': extract_figure1()
