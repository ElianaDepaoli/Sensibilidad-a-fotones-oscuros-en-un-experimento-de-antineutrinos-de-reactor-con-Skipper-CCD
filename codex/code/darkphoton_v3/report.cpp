// ROOT/C++ self-contained report assembly. Run from code/darkphoton_v3.
#include <TROOT.h>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <ctime>
#include <iomanip>
std::string Read(const std::string& p){std::ifstream f(p);if(!f)throw std::runtime_error("Cannot read "+p);return {std::istreambuf_iterator<char>(f),{}};}
std::string Escape(std::string s){std::string o;for(char c:s){if(c=='&')o+="&amp;";else if(c=='<')o+="&lt;";else if(c=='>')o+="&gt;";else o+=c;}return o;}
std::string Svg(const std::string& p){auto s=Read(p);auto a=s.find("<svg");if(a==std::string::npos)throw std::runtime_error("missing SVG");return s.substr(a);}
void WriteHtml(){
 std::ofstream h("../../report_darkphoton_v3.html");
 h<<R"(<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>TEXONO: corrected ROOT calculation</title><style>body{font:17px/1.55 system-ui;max-width:1120px;margin:35px auto;padding:0 24px;color:#182633;background:#faf9f5}svg{width:100%;height:auto}pre{white-space:pre-wrap;background:#edf1f3;padding:16px;font-size:13px}.note{background:#fff0d5;padding:18px}table{border-collapse:collapse}td,th{padding:9px;border-bottom:1px solid #ddd}</style><h1>TEXONO: production, detection and comparison</h1>
 <p><a href="report_darkphoton_v3_en.pdf">English lecture</a> · <a href="report_darkphoton_v3_es.pdf">Clase en español</a></p>
 <p class="note">The old production plot was empty because ROOT graphs died before canvas export. The corrected run keeps them alive and checks every curve for finite points inside the axes. All three figures have been visually inspected. The independent calculation does not reproduce either published result completely; the differences are shown.</p>
 <h2>Corrections and source audit</h2><ul><li>Danilov Eq. (5) uses coefficient 2, not 4.</li><li>The medium source is normalized at small mixing, not at epsilon = 1.</li><li>The detector conversion formula is not valid at resonance. A ten-percent mass window around each plasma mass is omitted.</li><li>Roots are bounded to epsilon at most 0.01. No root means no reported limit.</li><li>The default selection efficiency is 1/6, following Danilov's estimate. The unit-efficiency result is also plotted.</li><li>Park and Danilov reference results now appear explicitly.</li></ul>
 <h2>Recalculation of Park Figure 1</h2><p>The production trace kernel is inherited from v2 and checked again against independent Dirac matrices, Ward identities, the analytic integral, and Reference 15. Open points are the v2 vector extraction from the supplied Park PDF; lines are freshly recalculated. No normalization is fitted. The full published spectrum is not reproduced.</p>)";
 h<<Svg("output/figure1_comparison.svg");
 h<<R"(<h2>Medium source</h2><p>This source follows Danilov's direct oscillation prescription and preserves photon energy. It is different from the energy-redistributing Compton convolution above. The probabilities are leading-order medium approximations, not exact all-orders results.</p>)"<<Svg("output/figure1_medium.svg");
 h<<R"(<h2>TEXONO comparison</h2><p>The red line is Park's quoted bound. The black line is digitized from the actual Danilov vector PDF, calibrated using its log-axis endpoints. The magenta curve recalculates the Park convention with massive inverse scattering and a 1.96-sigma cap. Blue/cyan curves use direct medium conversion with estimated/unit efficiency and the 195.7-event one-sided cap. No medium claim is made above 10 keV or in the resonance gap. The normalization disagreement with Danilov remains unresolved; no tuning was applied.</p>)"<<Svg("output/figure2_texono_medium.svg");
 h<<R"(<h2>Regeneration and experiment settings</h2><pre>cd code/darkphoton_v3
g++ -O2 -std=c++17 v3.cpp $(root-config --cflags --libs) -o v3
./v3 config.conf &gt; output/run_receipt.txt 2&gt;&amp;1
g++ -O2 -std=c++17 report.cpp $(root-config --cflags --libs) -o report
./report</pre><p>Edit config.conf for power, baseline, mass, live time, A, Z, ROI, plasma masses, absorption, event cap and efficiency. A different background can be supplied as observed_minus_background and excess_sigma with use_gaussian_cap=1; this uses a simple one-sided Gaussian approximation. For Poisson or nuisance-parameter models supply an externally justified N95_events instead. The third command-line argument can override the constant signal efficiency.</p>
<h2>Actual execution receipt</h2><pre>)"<<Escape(Read("output/run_receipt.txt"))<<"</pre>";
 h<<R"(<h2>Unchecked physics</h2><p>No detector Monte Carlo, finite detector propagation, material-dependent attenuation, energy-dependent efficiency, photon escape, or reactor rescattering is supplied. Full energy containment is assumed. These are conditional model limits, not a validated detector-level exclusion. The Park statistical convention is not the same one-sided confidence construction as Danilov's.</p>
<h2>Sources and implementation references</h2><p>Local physics sources: papers/darkphoton_v3/{park2017,danilov2019,reference15}.pdf. The English and Spanish lectures state the equations and validity domains. ROOT graph lifetime guidance: <a href="https://root.cern.ch/doc/master/classTGraph.html">TGraph documentation</a> and <a href="https://root-forum.cern.ch/t/no-display-when-drawing-graph/63217">ROOT forum scope discussion</a>.</p><footer><h2>Execution metadata</h2><pre>)"<<Escape(Read("output/runtime.txt"))<<"</pre>";
 std::tm tm{};std::istringstream start("2026-09-21 16:04:12");start>>std::get_time(&tm,"%Y-%m-%d %H:%M:%S");auto stamp=timegm(&tm);
 h<<"<p>Measured audit interval since the first recorded UTC clock sample: "<<std::difftime(std::time(nullptr),stamp)<<" seconds. Work began before that sample, so this is a lower bound on task time. Full task duration, exact token usage, and reasoning-effort metadata are not exposed.</p></footer></html>";
}
int main(){WriteHtml();}
