"""Check artifact links and provenance, including explicitly failed comparisons."""
import json
from html.parser import HTMLParser
from pathlib import Path
import yaml

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]

class Links(HTMLParser):
    def __init__(self):
        super().__init__();self.hrefs=[];self.images=[]
    def handle_starttag(self,tag,attrs):
        attributes=dict(attrs)
        if tag=='a' and attributes.get('href'):self.hrefs.append(attributes['href'])
        if tag=='img':self.images.append(attributes.get('src',''))

def main():
    reports=['report_park_texono_human.html','report_park_texono_agent.html','rfh_park_texono_en.html','rfh_park_texono_es.html']
    if (ROOT/'report_park_texono_continuation.html').exists():reports.append('report_park_texono_continuation.html')
    for filename in reports:
        parser=Links();parser.feed((ROOT/filename).read_text())
        for link in parser.hrefs:
            if not link.startswith(('http:','https:','#')):assert (ROOT/link).exists(),link
        assert all(src.startswith('data:image/svg+xml;base64,') for src in parser.images)
    numbers=json.loads((ROOT/'provenance/numbers.json').read_text())
    claims=yaml.safe_load((ROOT/'provenance/claims.yaml').read_text())
    ids={c['id'] for c in claims['claims']}
    for claim in claims['claims']:
        assert claim['evidence']
        for slug in claim.get('numbers',[]):assert slug in numbers,slug
    for fig in claims['figures']:
        assert (ROOT/fig['file']).exists(),fig['file']
        assert set(fig['supports'])<=ids
    vals=json.loads((BASE/'output/results.json').read_text())
    # A failed paper match must be visible; it must never be hidden as success.
    assert vals['paper_limit_reproduced_within_2p5_percent'] is False
    assert abs(vals['epsilon_limit']/vals['paper_epsilon_limit']-1)>.025
    assert vals['quadrature_relative_difference']<2e-5
    assert vals['integration_order_relative_difference']<2e-4
    assert (BASE/'output/figures.html').stat().st_size>1000
    print('PASS: report links and embedded figures; provenance references; exported notebook exists')
    print('PASS: numerical-convergence checks; failed paper-limit comparison remains explicitly recorded')
    print('UNCOVERED: reactor transport, detector response, raw-event likelihood, decay loop, prompt-fission spectrum')

if __name__=='__main__':main()
