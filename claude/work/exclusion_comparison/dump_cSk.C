// Runs one of the user's macros with given DANILOV flags and dumps every TGraph in the cSk
// canvas (label from the legend) to a text file. The macro itself is not modified.
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TList.h>
#include <TROOT.h>
#include <fstream>
void dump_cSk(const char *macro, const char *flags, const char *out)
{
   gROOT->ProcessLine(Form(".L %s", macro));
   gROOT->ProcessLine(flags);
   TString fn(macro); fn.ReplaceAll(".C", "");
   gROOT->ProcessLine(Form("%s();", fn.Data()));
   auto *c = (TCanvas *)gROOT->FindObject("cSk");
   TLegend *leg = nullptr;
   for (auto *o : *c->GetListOfPrimitives()) if (o->InheritsFrom("TLegend")) leg = (TLegend *)o;
   std::ofstream f(out);
   for (auto *e : *leg->GetListOfPrimitives()) {
      auto *le = (TLegendEntry *)e;
      auto *g = dynamic_cast<TGraph *>(le->GetObject());
      if (!g) continue;
      f << "# " << le->GetLabel() << "  flags: " << flags << "\n";
      for (int i = 0; i < g->GetN(); ++i) f << g->GetX()[i] << " " << g->GetY()[i] << "\n";
      f << "\n\n";
   }
}
