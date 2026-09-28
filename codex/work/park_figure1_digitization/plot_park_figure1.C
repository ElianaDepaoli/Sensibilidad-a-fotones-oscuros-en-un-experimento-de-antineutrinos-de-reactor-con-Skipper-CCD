#include <TCanvas.h>
#include <TGraph.h>
#include <TH1.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>
#include <array>
#include <cassert>
#include <cmath>
#include <fstream>
#include <filesystem>
#include <iostream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <vector>

void plot_park_figure1() {
  const std::string base = std::filesystem::absolute(__FILE__).parent_path().string() + "/";
  std::ifstream input(base + "park_figure1.txt");
  if (!input.good()) throw std::runtime_error("Cannot open input table: " + base + "park_figure1.txt");
  std::vector<std::array<double, 4>> rows;
  std::string line;
  while (std::getline(input, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream stream(line);
    std::array<double, 4> row;
    std::string token;
    for (auto &value : row) {
      if (!(stream >> token)) throw std::runtime_error("Missing column");
      value = std::stod(token);
    }
    if (stream >> token) throw std::runtime_error("Extra column");
    if (!rows.empty() && row[0] <= rows.back()[0]) throw std::runtime_error("Unsorted energy");
    rows.push_back(row);
  }
  if (rows.size() != 430) throw std::runtime_error("Unexpected row count");
  gStyle->SetOptStat(0);
  auto canvas = new TCanvas("park_figure1", "Park Figure 1: digitized curves", 1000, 700);
  canvas->SetLogy();
  canvas->SetLeftMargin(.15);
  canvas->SetBottomMargin(.14);
  canvas->SetTopMargin(.06);
  canvas->SetRightMargin(.04);
  canvas->SetTicks(1, 1);
  auto frame = canvas->DrawFrame(0., 1.5e-3, 4.3, 3.);
  frame->GetXaxis()->SetTitle("E_{A'} (MeV)");
  frame->GetYaxis()->SetTitle("dN_{A'}/dE_{A'} [10^{21} MeV^{-1} s^{-1}]");
  frame->GetXaxis()->SetTitleSize(.045);
  frame->GetYaxis()->SetTitleSize(.045);
  frame->GetYaxis()->SetTitleOffset(1.5);
  auto legend = new TLegend(.64, .71, .94, .90);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  const char *labels[] = {"m_{A'} = 0.1 MeV", "m_{A'} = 0.5 MeV", "m_{A'} = 1.0 MeV"};
  for (int column = 1; column <= 3; ++column) {
    TGraph *segment = nullptr;
    bool labeled = false;
    int count = 0;
    for (const auto &row : rows) {
      if (!std::isfinite(row[column])) {
        segment = nullptr; // Never connect across missing or clipped bins.
        continue;
      }
      if (row[column] <= 0) throw std::runtime_error("Nonpositive spectrum");
      if (!segment) {
        segment = new TGraph();
        segment->SetLineColor(column == 1 ? kBlue : kBlack);
        segment->SetLineStyle(column == 3 ? 3 : 1);
        segment->SetLineWidth(1);
        segment->Draw("L SAME");
        if (!labeled) {
          legend->AddEntry(segment, labels[column-1], "l");
          labeled = true;
        }
      }
      // Each input row is a bin center; draw the original histogram steps.
      segment->SetPoint(segment->GetN(), row[0]-.005, row[column]);
      segment->SetPoint(segment->GetN(), row[0]+.005, row[column]);
      ++count;
    }
    if (count == 0) throw std::runtime_error("Empty spectrum");
    std::cout << "PASS: spectrum column " << column+1 << ": " << count << " visible bins plotted\n";
  }
  legend->Draw();
  canvas->RedrawAxis();
  canvas->SaveAs((base + "park_figure1_root.png").c_str());
  canvas->SaveAs((base + "park_figure1_root.pdf").c_str());
  std::cout << "ROOT " << gROOT->GetVersion() << ": plotted text data only; missing bins remain gaps.\n";
}
