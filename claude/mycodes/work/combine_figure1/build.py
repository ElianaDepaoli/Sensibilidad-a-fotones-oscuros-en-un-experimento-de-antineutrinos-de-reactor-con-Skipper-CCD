from pathlib import Path
import re

base = Path(__file__).resolve().parents[2]
header = '''// Combined ROOT plot. Data copied verbatim from cFig1.C and park_figure1_digitized.C.
// Run: root -l 'combined_figure1.C()'
// No rescaling, smoothing, interpolation or joining of digitized segments.
#include <TCanvas.h>
#include <TColor.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>
#include <vector>

void combined_figure1()
{
   if (auto *old = gROOT->FindObject("combined_figure1")) delete old;
   auto *canvas = new TCanvas("combined_figure1", "cFig1 + PARK", 1100, 750);
   canvas->SetLogy();
   canvas->SetTicks(1, 1);
   canvas->SetLeftMargin(0.15);
   canvas->SetRightMargin(0.04);
   canvas->SetBottomMargin(0.13);
   canvas->SetTopMargin(0.07);
   auto *frame = canvas->DrawFrame(0., 0.0001, 4.3, 5.);
   frame->SetTitle("");
   frame->GetXaxis()->SetTitle("E_{A'} [MeV]");
   frame->GetYaxis()->SetTitle("dN_{A'}/dE_{A'} [10^{21} MeV^{-1} s^{-1}]");
   frame->GetXaxis()->SetTitleSize(0.045);
   frame->GetYaxis()->SetTitleSize(0.045);
   frame->GetYaxis()->SetTitleOffset(1.5);
   auto *legend = new TLegend(0.56, 0.62, 0.94, 0.91);
   legend->SetBorderSize(0);
   legend->SetFillStyle(0);
   legend->SetTextFont(42);
   legend->SetTextSize(0.029);
'''
parts = [header]
for filename, tag, style, representatives in [
    ('cFig1.C', 'cFig1', 1, {0: '0.1', 1: '0.5', 2: '1.0'}),
    ('park_figure1_digitized.C', 'PARK', 3, {0: '0.1', 6: '0.5', 16: '1.0'}),
]:
    source = (base / filename).read_text()
    vectors = re.findall(r'   std::vector<Double_t> graph_[xy]_vect\d+\{.*?\};', source, re.S)
    graphs = re.findall(r'(?:TGraph \*)?graph = new TGraph\((.*?)\);(.*?)(?=graph->Draw)', source, re.S)
    assert len(vectors) == 2 * len(graphs)
    parts.append(f'\n   // {filename}: original coordinates and colors.\n   {{\n')
    for i, (constructor, attrs) in enumerate(graphs):
        parts.extend([vectors[2*i] + '\n', vectors[2*i+1] + '\n'])
        color = re.search(r'graph->SetLineColor\((.*?)\);', attrs)
        color = color.group(1) if color else 'kBlack'
        name = f'{tag}_{i}'
        parts.append(f'''   auto *{name} = new TGraph({constructor});
   {name}->SetName("{name}");
   {name}->SetLineColor({color});
   {name}->SetLineStyle({style});
   {name}->SetLineWidth(2);
   {name}->Draw("L SAME");
''')
        if i in representatives:
            label = f"m_{{A'}} = {representatives[i]} MeV ({tag})"
            parts.append(f'   legend->AddEntry({name}, "{label}", "l");\n')
    parts.append('   }\n')
parts.append('''
   legend->Draw();
   canvas->RedrawAxis();
   canvas->Modified();
   canvas->Update();
   canvas->SaveAs("combined_figure1.png");
   canvas->SaveAs("combined_figure1.pdf");
}
''')
(base / 'combined_figure1.C').write_text(''.join(parts))
print('Created standalone combined_figure1.C; original coordinate vectors copied verbatim.')
