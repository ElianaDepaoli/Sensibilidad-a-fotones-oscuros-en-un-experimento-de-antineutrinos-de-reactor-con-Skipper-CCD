#include "../../combined_figure1.C"
#include "../../cFig1.C"
#include "../../park_figure1_digitized.C"
#include <TLegendEntry.h>
#include <cassert>
#include <iostream>
#include <string>

std::vector<TGraph *> graphs_on(TCanvas *canvas)
{
   std::vector<TGraph *> result;
   TIter next(canvas->GetListOfPrimitives());
   while (auto *obj = next())
      if (auto *graph = dynamic_cast<TGraph *>(obj)) result.push_back(graph);
   return result;
}

void check()
{
   cFig1();
   auto original = graphs_on(static_cast<TCanvas *>(gROOT->FindObject("cFig1")));
   park_figure1_digitized();
   auto park = graphs_on(static_cast<TCanvas *>(gROOT->FindObject("park_figure1")));
   const auto solid_count = original.size();
   original.insert(original.end(), park.begin(), park.end());
   combined_figure1();
   auto *canvas = static_cast<TCanvas *>(gROOT->FindObject("combined_figure1"));
   auto combined = graphs_on(canvas);
   assert(combined.size() == original.size());
   for (size_t i = 0; i < combined.size(); ++i) {
      assert(combined[i]->GetN() == original[i]->GetN());
      assert(combined[i]->GetLineColor() == original[i]->GetLineColor());
      assert(combined[i]->GetLineStyle() == (i < solid_count ? 1 : 3));
      for (int j = 0; j < combined[i]->GetN(); ++j) {
         assert(combined[i]->GetPointX(j) == original[i]->GetPointX(j));
         assert(combined[i]->GetPointY(j) == original[i]->GetPointY(j));
      }
   }
   auto *legend = static_cast<TLegend *>(canvas->GetListOfPrimitives()->FindObject("TPave"));
   assert(legend);
   assert(legend->GetListOfPrimitives()->GetSize() == 6);
   TIter entries(legend->GetListOfPrimitives());
   int entry_index = 0;
   while (auto *entry = static_cast<TLegendEntry *>(entries())) {
      auto *graph = dynamic_cast<TGraph *>(entry->GetObject());
      assert(graph);
      const bool is_park = entry_index++ >= 3;
      assert(graph->GetLineStyle() == (is_park ? 3 : 1));
      assert((std::string(entry->GetLabel()).find("PARK") != std::string::npos) == is_park);
   }
   assert(canvas->GetLogy());
   std::cout << "PASS: original coordinates and colors preserved for every graph and segment.\n"
             << "PASS: cFig1 solid, PARK dotted; legend labels linked to the correct graphs.\n";
}
