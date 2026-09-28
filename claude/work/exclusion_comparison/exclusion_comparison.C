// TEXONO 95% CL exclusion in (m_A', epsilon): every approach of the project on one canvas.
// Run: root -l -b -q exclusion_comparison.C   (after ./run.sh has filled data/)
// Inputs (data/): curves dumped from the cSk canvas of mycodes/DP3_Danilov.C and DP3_Danilov_new.C
// (dump_cSk.C), and the exclusion CSVs written by the Codex projects, copied verbatim.
// No curve is rescaled, smoothed or interpolated; rows that are nan or non-positive are skipped.
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TColor.h>
#include <TH1F.h>
#include <TPad.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TStyle.h>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

// Reads one block of a dump_cSk.C file: the block whose '#' header contains `label`.
TGraph *read_dump(const std::string &file, const std::string &label)
{
   std::ifstream f("data/" + file);
   std::string line;
   bool on = false;
   auto *g = new TGraph();
   while (std::getline(f, line)) {
      if (line.rfind("#", 0) == 0) { on = line.find(label) != std::string::npos; continue; }
      if (!on || line.empty()) continue;
      double x, y;
      std::istringstream(line) >> x >> y;
      if (std::isfinite(y) && y > 0) g->AddPoint(x, y);
   }
   return g;
}

// Reads column `ycol` against column `xcol` of a CSV with a header; x is multiplied by `xscale`
// (1e-6 turns eV into MeV). Rows where y is nan, empty or non-positive are skipped.
TGraph *read_csv(const std::string &file, const std::string &xcol, const std::string &ycol, double xscale)
{
   std::ifstream f("data/" + file);
   std::string line;
   std::getline(f, line);
   std::vector<std::string> head;
   {
      std::stringstream ss(line);
      std::string c;
      while (std::getline(ss, c, ',')) head.push_back(c);
   }
   int ix = -1, iy = -1;
   for (int i = 0; i < (int)head.size(); ++i) {
      if (head[i] == xcol) ix = i;
      if (head[i] == ycol) iy = i;
   }
   auto *g = new TGraph();
   if (ix < 0 || iy < 0) { printf("missing column in %s\n", file.c_str()); return g; }
   while (std::getline(f, line)) {
      std::vector<std::string> v;
      std::stringstream ss(line);
      std::string c;
      while (std::getline(ss, c, ',')) v.push_back(c);
      if ((int)v.size() <= std::max(ix, iy)) continue;
      double x = std::atof(v[ix].c_str()) * xscale, y = std::atof(v[iy].c_str());
      if (v[iy] == "nan" || v[iy].empty() || !std::isfinite(y) || y <= 0) continue;
      g->AddPoint(x, y);
   }
   return g;
}

void style(TGraph *g, const char *hex, int ls, int lw)
{
   g->SetLineColor(TColor::GetColor(hex));
   g->SetLineStyle(ls);
   g->SetLineWidth(lw);
}

// Prints epsilon at the grid point closest (in log m) to each benchmark mass, for the report table.
void table_row(const char *name, TGraph *g)
{
   printf("%-36s", name);
   for (double m : {1e-6, 1e-4, 1e-2, 0.1, 0.5, 1.0}) {
      double best = 1e9, y = NAN;
      for (int i = 0; i < g->GetN(); ++i) {
         double d = std::fabs(std::log(g->GetX()[i] / m));
         if (d < best) { best = d; y = g->GetY()[i]; }
      }
      if (best < 0.25) printf("  %9.3e", y); else printf("  %9s", "-");
   }
   printf("\n");
}

void exclusion_comparison()
{
   gStyle->SetOptStat(0);
   gStyle->SetPaperSize(20., 20.9);   // same aspect as the canvas, so the pdf is not letterboxed
   auto *c = new TCanvas("cExcl", "TEXONO exclusion, all approaches", 1100, 1150);
   // Top pad: the full plane. Bottom pad: zoom on the plateau, where the approaches differ by < 3x.
   auto *top = new TPad("top", "", 0, 0.34, 1, 1);
   auto *bot = new TPad("bot", "", 0, 0, 1, 0.34);
   for (auto *p : {top, bot}) {
      p->SetLogx();
      p->SetLogy();
      p->SetTicks(1, 1);
      p->SetLeftMargin(0.13);
      p->SetRightMargin(0.03);
      p->Draw();
   }
   top->SetTopMargin(0.03);
   top->SetBottomMargin(0.09);
   bot->SetTopMargin(0.04);
   bot->SetBottomMargin(0.2);
   top->cd();
   const double xmin = 1e-7, xmax = 1.2, ymin = 1e-6, ymax = 1e-1;
   auto *frame = top->DrawFrame(xmin, ymin, xmax, ymax);
   frame->GetXaxis()->SetTitle("m_{A'} [MeV]");
   frame->GetYaxis()->SetTitle("#varepsilon_{95}  (TEXONO, 95% CL)");
   frame->GetXaxis()->SetTitleSize(0.045);
   frame->GetYaxis()->SetTitleSize(0.045);
   frame->GetYaxis()->SetTitleOffset(1.2);

   // Published. Park (2017) states eps < 2.1e-5 for m_A' < 1 MeV (Fig. 2 is flat there).
   auto *park = new TGraph();
   park->AddPoint(xmin, 2.1e-5);
   park->AddPoint(1.0, 2.1e-5);
   style(park, "#000000", 1, 4);
   auto *danilov = read_csv("codex_v3_danilov_digitized.csv", "mass_eV", "epsilon", 1e-6);
   style(danilov, "#000000", 2, 4);

   // Human macro (DP3_Danilov.C) and Claude's corrected cross section (DP3_Danilov_new.C).
   auto *h_off = read_dump("orig_none.txt", "TEXONO");
   auto *h_on = read_dump("orig_ABD.txt", "TEXONO");
   auto *n_off = read_dump("new_none.txt", "exact");
   auto *n_on = read_dump("new_ABD.txt", "exact");
   style(h_off, "#1f5fbf", 1, 3);
   style(h_on, "#1f5fbf", 2, 3);
   style(n_off, "#e8710a", 1, 3);
   style(n_on, "#e8710a", 2, 3);
   // Same macros with wpsize = 200 (converged E_A' grid); drawn only in the zoom.
   auto *h_off2 = read_dump("orig_none_wp200.txt", "TEXONO");
   auto *h_on2 = read_dump("orig_ABD_wp200.txt", "TEXONO");
   auto *n_off2 = read_dump("new_none_wp200.txt", "exact");
   auto *n_on2 = read_dump("new_ABD_wp200.txt", "exact");
   style(h_off2, "#8fb8ef", 1, 4);
   style(h_on2, "#8fb8ef", 2, 4);
   style(n_off2, "#f7b97a", 1, 4);
   style(n_on2, "#f7b97a", 2, 4);

   // Codex, from scratch.
   auto *x_park = read_csv("codex_park_texono_figure2.csv", "mass_MeV", "epsilon_limit", 1.0);
   auto *x_rs_vac = read_csv("codex_reactor_sensitivity_limits.csv", "mass_MeV", "vacuum_no_decay", 1.0);
   auto *x_v3 = read_csv("codex_v3_texono_limit.csv", "mass_eV", "epsilon95", 1e-6);
   auto *x_fresh = read_csv("codex_danilov_fresh_exclusion.csv", "mass_eV", "epsilon95", 1e-6);
   style(x_park, "#1b9e77", 1, 2);
   style(x_rs_vac, "#7b3294", 1, 2);
   style(x_v3, "#1b9e77", 2, 2);
   style(x_fresh, "#d01c8b", 2, 2);

   for (auto *g : {x_park, x_rs_vac, x_v3, x_fresh, h_off, h_on, n_off, n_on, danilov, park}) g->Draw("L SAME");

   auto *leg = new TLegend(0.36, 0.50, 0.97, 0.965);
   leg->SetBorderSize(0);
   leg->SetFillStyle(1001);
   leg->SetFillColorAlpha(kWhite, 0.85);
   leg->SetTextFont(42);
   leg->SetTextSize(0.027);
   leg->SetHeader("solid: vacuum / Park-like      dashed: with Danilov-type corrections");
   leg->AddEntry(park, "Park 2017, published (2.1#times10^{-5})", "l");
   leg->AddEntry(danilov, "Danilov et al. 2019, published (digitized by Codex)", "l");
   leg->AddEntry(h_off, "DP3_Danilov.C (human), note eq. (50), flags off", "l");
   leg->AddEntry(h_on, "DP3_Danilov.C (human), DANILOV A+B+D on", "l");
   leg->AddEntry(n_off, "DP3_Danilov_new.C (human+Claude), exact #sigma, flags off", "l");
   leg->AddEntry(n_on, "DP3_Danilov_new.C (human+Claude), exact #sigma, A+B+D on", "l");
   leg->AddEntry(x_park, "Codex park_texono (vacuum Compton-like)", "l");
   leg->AddEntry(x_rs_vac, "Codex reactor_sensitivity (vacuum, #eta = 1/6)", "l");
   leg->AddEntry(x_v3, "Codex darkphoton_v3 (in-medium #gamma#leftrightarrowA', #eta = 1/6)", "l");
   leg->AddEntry(x_fresh, "Codex danilov_fresh (oscillation, pol., decay, #eta = 1/6)", "l");
   leg->Draw();

   top->RedrawAxis();

   bot->cd();
   auto *zoom = bot->DrawFrame(5e-5, 1.1e-5, xmax, 3.9e-5);
   zoom->GetYaxis()->SetNdivisions(505);
   zoom->GetXaxis()->SetTitle("m_{A'} [MeV]");
   zoom->GetYaxis()->SetTitle("#varepsilon_{95}  (zoom)");
   for (auto *ax : {zoom->GetXaxis(), zoom->GetYaxis()}) {
      ax->SetTitleSize(0.085);
      ax->SetLabelSize(0.07);
   }
   zoom->GetYaxis()->SetTitleOffset(0.72);
   zoom->GetYaxis()->SetNoExponent(false);
   zoom->GetXaxis()->SetTitleOffset(1.0);
   zoom->GetYaxis()->SetMoreLogLabels();
   for (auto *g : {h_off2, h_on2, n_off2, n_on2}) g->Draw("L SAME");
   for (auto *g : {x_park, x_rs_vac, x_v3, x_fresh, h_off, h_on, n_off, n_on, danilov, park}) g->Draw("L SAME");
   auto *leg2 = new TLegend(0.30, 0.66, 0.97, 0.92);
   leg2->SetBorderSize(0);
   leg2->SetFillStyle(0);
   leg2->SetTextFont(42);
   leg2->SetTextSize(0.055);
   leg2->AddEntry(h_off2, "light: same macros, wpsize = 200 (converged grid)", "l");
   leg2->Draw();
   bot->RedrawAxis();
   c->SaveAs("texono_exclusion_all.eps");   // eps -> pdf with epstopdf: ROOT's own pdf letterboxes the two pads
   gSystem->Exec("epstopdf texono_exclusion_all.eps && rm texono_exclusion_all.eps");
   c->SaveAs("texono_exclusion_all.png");

   printf("\n%-34s  %9s  %9s  %9s  %9s  %9s  %9s\n", "eps95 at m_A' [MeV] =", "1e-6", "1e-4", "1e-2", "0.1", "0.5", "1.0");
   table_row("Park 2017 published", park);
   table_row("Danilov 2019 digitized", danilov);
   table_row("DP3_Danilov.C off", h_off);
   table_row("DP3_Danilov.C A+B+D", h_on);
   table_row("DP3_Danilov_new.C exact off", n_off);
   table_row("DP3_Danilov_new.C exact A+B+D", n_on);
   table_row("DP3_Danilov.C off wp200", h_off2);
   table_row("DP3_Danilov.C A+B+D wp200", h_on2);
   table_row("DP3_Danilov_new.C exact off wp200", n_off2);
   table_row("DP3_Danilov_new.C exact ABD wp200", n_on2);
   table_row("Codex park_texono", x_park);
   table_row("Codex reactor_sensitivity vac", x_rs_vac);
   table_row("Codex darkphoton_v3", x_v3);
   table_row("Codex danilov_fresh", x_fresh);
}
