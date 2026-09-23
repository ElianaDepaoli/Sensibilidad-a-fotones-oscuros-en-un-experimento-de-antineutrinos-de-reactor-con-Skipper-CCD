// ROOT/C++ report and provenance assembly. No Python utilities are used.
#include "physics.hxx"
#include <TROOT.h>
#include <nlohmann/json.hpp>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <regex>

namespace Reports {
using namespace Production;
using Json=nlohmann::json;
using Row=std::map<std::string,std::string>;
Json numbers=Json::object();
std::string Read(const std::string&p){std::ifstream f(p);Require(bool(f),"Missing "+p);return std::string(std::istreambuf_iterator<char>(f),{});}
std::vector<std::string> Split(const std::string&s){std::istringstream in(s);std::string item;std::vector<std::string> fields;while(std::getline(in,item,','))fields.push_back(item);return fields;}
std::vector<Row> CSV(const std::string&file){std::ifstream f(file);Require(bool(f),"Missing CSV "+file);std::string line;std::getline(f,line);auto header=Split(line);std::vector<Row> rows;
    while(std::getline(f,line)){auto fields=Split(line);Require(fields.size()==header.size(),"CSV column mismatch");Row r;for(size_t i=0;i<header.size();++i)r[header[i]]=fields[i];rows.push_back(r);}return rows;}
std::string Fixed(double v,int n=3){std::ostringstream s;s<<std::fixed<<std::setprecision(n)<<v;return s.str();}
std::string Scientific(double v){std::ostringstream s;s<<std::scientific<<std::setprecision(5)<<v;auto value=s.str();auto pos=value.find('e');return value.substr(0,pos)+"\\times10^{"+std::to_string(std::stoi(value.substr(pos+1)))+"}";}
std::string Escape(std::string s){for(auto pair:std::vector<std::pair<std::string,std::string>>{{"&","&amp;"},{"<","&lt;"},{">","&gt;"}}){size_t pos=0;while((pos=s.find(pair.first,pos))!=std::string::npos){s.replace(pos,pair.first.size(),pair.second);pos+=pair.second.size();}}return s;}
void Record(const std::string&id,double value,const std::string&statement,const std::string&producer){
    numbers[id]={{"value",value},{"statement",statement},{"produced_by",producer},
    {"from_scratch","New C++ invariant-trace algebra, phase space, source integration and comparison; empirical inputs and constants are supplied, not new measurements."},
    {"from_library","ROOT quadrature, XML reader and graphics; C++ standard library; nlohmann JSON serialization."},
    {"choices",{"Stationary free electrons and tree level; Park spectrum and Compton-only source denominator retained; no later correction or fitted normalization.","Figure comparisons use fresh horizontal PDF segments and bin averages; residual tolerance is a diagnostic choice, not a statistical uncertainty."}}};
}
void RecordCSV(const std::string&name,const std::vector<Row>&rows){
    for(size_t i=0;i<rows.size();++i)for(auto&[key,text]:rows[i]){try{size_t used;double value=std::stod(text,&used);if(used==text.size())Record(name+"_"+std::to_string(i)+"_"+key,value,name+" row "+std::to_string(i)+" "+key,"output/"+name+".csv");}catch(const std::exception&) {}}
}

double ThomsonBoundary(){
    // Boundary for a specified error of Thomson relative to the exact KN total.
    double lo=1e-5,hi=.1;
    for(int i=0;i<70;++i){double mid=(lo+hi)/2;if(Thomson()/KleinNishina(mid)-1>.05)hi=mid;else lo=mid;}
    double energy=(lo+hi)/2;Close(Thomson()/KleinNishina(energy)-1,.05,1e-8,"Thomson boundary failed");return energy;
}

void Tables(){
    auto checks=CSV("output/checks.csv"),summary=CSV("output/comparison_summary.csv"),material=CSV("output/material.csv");
    RecordCSV("checks",checks);RecordCSV("comparison_summary",summary);RecordCSV("material",material);
    std::map<std::string,double> values;for(auto&r:checks)values[r.at("quantity")]=std::stod(r.at("value"));
    std::ofstream macros("output/numbers.tex");
    macros<<"\\newcommand{\\ReferenceError}{"<<Scientific(values["reference15_max_relative_error"])<<"}\n"
      <<"\\newcommand{\\KNError}{"<<Scientific(values["KN_max_relative_error"])<<"}\n"
      <<"\\newcommand{\\NodeError}{"<<Scientific(values["spectrum_node_max_relative_error"])<<"}\n"
      <<"\\newcommand{\\PhotonIntegral}{"<<Scientific(values["photon_rate_above1_per_s"])<<"}\n";
    double boundary=ThomsonBoundary()*1000;
    macros<<"\\newcommand{\\ThomsonBoundary}{"<<Fixed(boundary,3)<<"}\n";Record("Thomson_5percent_boundary_keV",boundary,"Thomson/KN minus one equals 0.05", "reports.cpp::ThomsonBoundary");
    std::ofstream thresholds("output/threshold_table.tex");thresholds<<"\\begin{tabular}{rr}\\toprule $M$ [MeV] & $E_{\\rm th}$ [MeV]\\\\\\midrule\n";
    for(double mass:{.1,.5,1.}){double threshold=Threshold(mass);thresholds<<mass<<" & "<<Fixed(threshold,6)<<"\\\\\n";Record("threshold_M_"+Fixed(mass,1),threshold,"Production threshold in MeV","physics.hxx::Threshold");}
    thresholds<<"\\bottomrule\\end{tabular}\n";
    std::ofstream residual("output/residual_table.tex");residual<<"\\begin{tabular}{rrrr}\\toprule $M$ [MeV] & \\both{Bins}{Bins} & \\both{Median residual}{Residuo mediano} & \\both{Maximum}{Máximo}\\\\\\midrule\n";
    for(auto&r:summary)if(r.at("region")=="central"){
        double med=100*std::stod(r.at("baseline_median_relative_error")),maximum=100*std::stod(r.at("baseline_max_relative_error"));
        residual<<r.at("mass_MeV")<<" & "<<r.at("bin_count")<<" & "<<Fixed(med,2)<<"\\% & "<<Fixed(maximum,2)<<"\\%\\\\\n";
        Record("central_median_percent_M_"+r.at("mass_MeV"),med,"Median bin residual percent","production.cpp::PlotSpectra; reports.cpp::Tables");
        Record("central_max_percent_M_"+r.at("mass_MeV"),maximum,"Maximum bin residual percent","production.cpp::PlotSpectra; reports.cpp::Tables");
    }residual<<"\\bottomrule\\end{tabular}\n";
    std::ofstream mat("output/material_table.tex");mat<<"\\begin{tabular}{rrrr}\\toprule $E$ [MeV] & $(\\mu/\\rho)_{\\rm NIST}$ [cm$^2$/g] & KN/NIST & $\\sigma_T/\\sigma_{\\rm KN}$\\\\\\midrule\n";
    double maxFraction=0;
    for(auto&r:material){double E=std::stod(r.at("energy_MeV")),fraction=std::stod(r.at("KN_fraction"));maxFraction=std::max(maxFraction,fraction);
        if(E==.2||E==1||E==2||E==3||E==5||E==8)mat<<E<<" & "<<Fixed(std::stod(r.at("total_cm2_g")),5)<<" & "<<Fixed(fraction)<<" & "<<Fixed(std::stod(r.at("Thomson_over_KN")),2)<<"\\\\\n";
    }mat<<"\\bottomrule\\end{tabular}\n";
    Require(maxFraction<.9,"Material diagnostic no longer supports reported criterion");
    Record("material_max_KN_fraction",maxFraction,"Largest KN proxy fraction among sampled uranium energies","reports.cpp::Tables");
}

Json Runtime(){
    std::time_t now=std::time(nullptr);std::tm start{};std::istringstream input("2026-09-17 23:08:41");input>>std::get_time(&start,"%Y-%m-%d %H:%M:%S");
    std::time_t epoch=timegm(&start);std::ostringstream end;end<<std::put_time(std::gmtime(&now),"%Y-%m-%dT%H:%M:%SZ");
    Json result={{"model","GPT-6"},{"reasoning_effort",nullptr},{"tokens_used",nullptr},{"first_recorded_environment_check_utc","2026-09-17T23:08:41Z"},
        {"report_build_utc",end.str()},{"elapsed_seconds",std::difftime(now,epoch)},{"scope","From first recorded environment check; includes dependency downloads and permission waits. Exact earlier prompt-reading time is not measured."}};
    std::ofstream("output/runtime.json")<<result.dump(2)<<'\n';return result;
}

void Provenance(const Json&runtime){
    std::filesystem::create_directories("provenance");
    std::map<std::string,double> constants={{"electron_mass_MeV",electron},{"alpha_inverse",1/alpha},{"hbarc_MeV_cm",hbarc},{"reference_power_MW",1000},
      {"gamma_normalization",.58e18},{"gamma_scale_MeV",.91},{"gamma_cutoff_MeV",1},{"alternative_cutoff_MeV",.2},{"gamma_max_MeV",40},{"gamma_tail_check_MeV",60},
      {"quadrature_nodes",96},{"fine_nodes",192},{"bin_nodes",12},{"relative_comparison_tolerance",.1},{"central_low_MeV",1.8},{"central_high_MeV",3.2},
      {"uranium_Z",92},{"uranium_A_g_mol",238.02891},{"Avogadro_per_mol",6.02214076e23},{"pair_threshold_MeV",2*electron},{"Thomson_error_criterion",.05},
      {"calibration_x0_pt",380.16015625},{"calibration_x_per_MeV",35.5224609375},{"calibration_y_rate1_pt",64.36328125},{"calibration_y_per_decade",31.28125}};
    for(auto&[key,value]:constants)Record(key,value,"Explicit constant or diagnostic convention: "+key,"physics.hxx; production.cpp; reports.cpp::Provenance");
    Record("elapsed_seconds",runtime.at("elapsed_seconds"),"Wall time since first recorded environment check","reports.cpp::Runtime");
    std::ofstream("provenance/numbers.json")<<numbers.dump(2)<<'\n';
    std::vector<std::string> keys;for(auto i=numbers.begin();i!=numbers.end();++i)keys.push_back(i.key());
    Json claims=Json::array({
      {{"id","v2-production-derivation"},{"statement","The new invariant trace gives the compact tree-level production cross section; exact photon-limit and crossing identities and independent matrix, Ward and polarization checks pass."},
       {"evidence",{"algebra.hxx::Derive","physics.hxx::CheckMatrices","output/run_receipt.txt","output/trace_terms.tex"}},{"numbers",Json::array()}},
      {{"id","v2-total-check"},{"statement","The newly derived analytic total and quadrature agree with Reference 15 on the tested grid; Reference 15 is a validation formula, not the production kernel."},
       {"evidence",{"physics.hxx::AnalyticTotal","physics.hxx::Reference15","production.cpp::Checks","../../papers/park_production_v2/reference15.pdf::Appendix A.1-A.3"}},{"numbers",keys}},
      {{"id","v2-figure-failure"},{"statement","The stated-equation source calculation does not reproduce full Figure 1. Fresh PDF-bin comparisons expose finite-mass disagreement; no normalization has been fitted."},
       {"evidence",{"production.cpp::Digitize","production.cpp::PlotSpectra","output/comparison_summary.csv","../../papers/park_production_v2/park2017.pdf::Eqs. 1-3 and Figure 1"}},{"numbers",keys}},
      {{"id","v2-compton-validity"},{"statement","Thomson is unsuitable in the plotted MeV region; free-electron KN differs from sampled total uranium attenuation beyond the chosen precision criterion. Material data do not modify the source."},
       {"evidence",{"reports.cpp::ThomsonBoundary","production.cpp::Material","output/material.csv","../../papers/park_production_v2/nist_uranium.html"}},{"numbers",keys}}
    });
    Json figures=Json::array();
    for(auto stem:{"figure1_comparison","cross_section","compton_validity"})for(auto ext:{"pdf","svg"}){
        std::string name=stem,id=name=="figure1_comparison"?"v2-figure-failure":name=="cross_section"?"v2-production-derivation":"v2-compton-validity";
        figures.push_back({{"file","output/"+name+"."+ext},{"produced_by","production.cpp::"+(name=="figure1_comparison"?std::string("PlotSpectra"):name=="cross_section"?std::string("PlotCrossSection"):std::string("Material"))},
          {"shows",name=="figure1_comparison"?"Dark-photon energy versus reactor source rate; calculated curves versus fresh PDF bins.":name=="cross_section"?"Incident photon energy versus total production cross section per electron, with mixing squared stripped.":"Incident energy versus uranium mass attenuation; free-electron KN proxy and NIST total."},
          {"from_scratch","New QED trace and source implementation; empirical spectrum, paper curves and material table are supplied inputs."},{"from_library","ROOT quadrature, XML reader and graphics; Poppler source PDF to SVG."},
          {"choices",{"Fixed source convention; show failed comparison; no later correction or fitted normalization."}},{"supports",{id}}});
    }
    figures.push_back({{"file","output/paper_page2.svg"},{"produced_by","pdftocairo -f 2 -l 2 -svg ../../papers_init/park2017.pdf output/paper_page2.svg"},
        {"shows","Original source page with Figure 1; retained for vector-coordinate extraction, not a computed figure."},
        {"from_scratch","No scientific curve: original source page."},{"from_library","Poppler PDF-to-SVG conversion."},
        {"choices",{"Read tick positions and horizontal steps from the supplied journal PDF; no raster fitting."}},{"supports",{"v2-figure-failure"}}});
    for(auto&entry:std::filesystem::directory_iterator("output"))if(entry.path().extension()==".png")
        figures.push_back({{"file",entry.path().generic_string()},{"produced_by","pdftoppm::render task PDF for visual verification"},
          {"shows","Raster preview of a task source, plot or lecture page; not a new numerical calculation."},
          {"from_scratch","No additional scientific calculation."},{"from_library","Poppler PDF rendering."},
          {"choices",{"Preview used to check readable equations, axes and layout."}},{"supports",{"v2-production-derivation","v2-figure-failure"}}});
    // JSON is also valid YAML; this avoids introducing a separate YAML dependency.
    std::ofstream("provenance/claims.yaml")<<Json{{"claims",claims},{"figures",figures}}.dump(2)<<'\n';
}

void Handoff(const Json&runtime){
    std::ofstream html("../../report_park_production_v2_agent.html");
    html<<R"HTML(<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Fresh Park production derivation</title><style>body{font:17px/1.55 system-ui;background:#faf9f5;color:#182633;max-width:1050px;margin:40px auto;padding:0 24px}h1,h2{line-height:1.2}pre{background:#eef1f3;padding:16px;overflow:auto;font-size:12px}svg{width:100%;height:auto}.note{border-left:4px solid #ad711e;background:#fff0d5;padding:16px}a{color:#145b87}footer{border-top:1px solid #bbb;margin-top:30px}</style><h1>Dark-photon production: independent ROOT derivation</h1>
<p class="note">The new cross-section derivation passes the independent checks. The source equation does not reproduce the full published Figure 1. No normalization has been fitted and no later correction is applied.</p>
<p><a href="report_park_production_v2_en.pdf">English physics lecture</a> · <a href="report_park_production_v2_es.pdf">Apuntes de física en español</a></p>
<h2>What was built</h2><p>This task uses a new C++ implementation. It imports no previous code, symbolic expression, digitized data, or result. The supplied Park PDF is the requested source. Reference 15 and the NIST table were downloaded afresh. The derivation reduces Dirac traces to exact dyadic polynomials and checks the compact expression with independent complex matrices and physical polarizations.</p>
<p>Reference 15 is used only to test the integrated result. The total cross section is also integrated analytically from the new differential formula. The English and Spanish LaTeX lectures contain the background formulas, both amplitudes, Ward identities, trace terms, compact result, phase space, integration, figure comparison and Compton-validity discussion.</p>
<h2>Actual run</h2><pre>)HTML"<<Escape(Read("output/run_receipt.txt"))<<"</pre>";
    for(auto stem:{"figure1_comparison","cross_section","compton_validity"}){std::string svg=Read(std::string("output/")+stem+".svg");auto begin=svg.find("<svg");html<<"<section>"<<svg.substr(begin)<<"</section>";}
    html<<R"HTML(<h2>Reproduction status and uncovered work</h2><p>The central-interval median residual criterion passes only for the smallest mass. It does not certify the full curve. The lower incident-photon cutoff improves low-energy agreement but cannot resolve the central finite-mass mismatch. Raw simulated events, original generator settings and bin uncertainties are unavailable. An undocumented implementation choice is a hypothesis, not an established explanation.</p>
<p>Not checked: original reactor-spectrum measurement, reactor transport, bound-electron physics, loop corrections, the author's generator and event statistics. This task calculates production only. Background, detector mass, exposure, ROI, A and Z enter a subsequent detection model; adding them to this source cross section would create a false dependence.</p>
<h2>Files and continuation</h2><ul><li><a href="code/park_production_v2/production.cpp">ROOT calculation, extraction and figures</a></li><li><a href="code/park_production_v2/algebra.hxx">Exact polynomial trace derivation</a> and <a href="code/park_production_v2/physics.hxx">cross sections and independent matrices</a></li><li><a href="code/park_production_v2/reports.cpp">C++ report and provenance generator</a></li><li><a href="code/park_production_v2/source.conf">Source configuration</a></li><li><a href="code/park_production_v2/provenance/claims.yaml">Task-local claims and figure provenance</a> and <a href="code/park_production_v2/provenance/numbers.json">numbers</a></li><li><a href="vault/skills/park-production-root/SKILL.md">ROOT-only task skill</a></li></ul>
<p>Run from code/park_production_v2:</p><pre>g++ -O2 -std=c++17 production.cpp $(root-config --cflags --libs) -lXMLIO -o production
./production &gt; output/run_receipt.txt 2&gt;&amp;1
g++ -O2 -std=c++17 reports.cpp $(root-config --cflags --libs) -o reports
./reports</pre><p>The PDF lectures are compiled from the root LaTeX wrappers with the task-local Tectonic binary and cache. Source and artifact hashes and command receipts remain with the task.</p>
<footer><h2>Execution metadata</h2><p>Model: GPT-6, as identified by this session. Reasoning effort: not exposed. Exact tokens used: not exposed; no estimate substituted.</p><p>ROOT 6.36.000; g++ (Ubuntu 13.3.0-6ubuntu2~24.04) 13.3.0. No base-environment packages installed. A standalone LaTeX compiler and its cache are confined to this task directory.</p><p>)HTML";
    html<<"First recorded environment check: "<<runtime["first_recorded_environment_check_utc"].get<std::string>()<<". Report build: "<<runtime["report_build_utc"].get<std::string>()
        <<". Elapsed wall time: "<<runtime["elapsed_seconds"]<<" seconds, including dependency downloads and permission waits. Exact earlier prompt-reading time is not measured.</p></footer></html>";
}

void Validate(){
    auto values=Json::parse(Read("provenance/numbers.json")),claims=Json::parse(Read("provenance/claims.yaml"));
    for(auto&claim:claims["claims"]){Require(!claim["evidence"].empty(),"Missing evidence");for(auto&id:claim["numbers"])Require(values.contains(id.get<std::string>()),"Missing number record");}
    for(auto&figure:claims["figures"])Require(std::filesystem::exists(figure["file"].get<std::string>()),"Missing figure");
    auto html=Read("../../report_park_production_v2_agent.html");std::regex link("href=\"([^\"]+)\"");
    for(std::sregex_iterator i(html.begin(),html.end(),link),end;i!=end;++i){std::string path=(*i)[1];if(path.find(":")==std::string::npos)Require(std::filesystem::exists("../../"+path),"Missing report link "+path);}
    Require(html.find("full published Figure 1")!=std::string::npos,"Missing failure disclosure");
    std::cout<<"PASS: report links, figure files and task-local provenance references\n";
}

void Vault(const Json&runtime){
    std::ofstream vault("../../vault/exercises/park-production-v2.md");
    vault<<"# Fresh Park production derivation\n\nNew ROOT/C++ calculation under `code/park_production_v2/`; no prior implementation or results imported. Both exchange diagrams, exact invariant trace, analytic integral and independent numerical matrices agree. Reference 15 only validates the total. The full Figure 1 is not reproduced; finite-mass disagreement remains. No later corrections.\n\nReports: `report_park_production_v2_en.pdf`, `report_park_production_v2_es.pdf`, `report_park_production_v2_agent.html`. Task-local provenance: `code/park_production_v2/provenance/`. Links: [[park-production-v2-source]].\n\n## Run receipt\n\n- Date: "<<runtime["report_build_utc"].get<std::string>()<<"\n- Environment: ROOT 6.36.000, Ubuntu g++ 13.3.0. Task-local Tectonic; no Python or base installs.\n- Commands in task directory: `g++ -O2 -std=c++17 production.cpp $(root-config --cflags --libs) -lXMLIO -o production`, `./production > output/run_receipt.txt 2>&1`; similarly build and run `reports.cpp`.\n- Actual physics output:\n\n```text\n"<<Read("output/run_receipt.txt")<<"```\n\nDerived: trace, lab support, differential and total production. Supplied: QED rules, Park spectrum, constants, Reference 15 validation formula, NIST material data. Not checked: original source measurement, transport, atomic effects, loops, author's generator or event statistics.\n";
    std::ofstream paper("../../vault/papers/park-production-v2-source.md");paper<<"# Park Figure 1 — new source check\n\nSource: `papers/park_production_v2/park2017.pdf`, copied from the requested `papers_init/` source. Reference 15 and NIST were downloaded into the same new source directory. [[park-production-v2]] derives production from the two QED diagrams and checks its total independently. Park's stated source convolution does not reproduce the full finite-mass Figure 1; exact figure reproduction remains unresolved.\n";
    for(auto receipt:{"validation_receipt.txt","latex_en_receipt.txt","latex_es_receipt.txt"})if(std::filesystem::exists(std::string("output/")+receipt))
        vault<<"\nAdditional receipt: `output/"<<receipt<<"`\n\n```text\n"<<Read(std::string("output/")+receipt)<<"```\n";
}
}

int main(int argc,char**argv){try{
    if(argc>1&&std::string(argv[1])=="--validate"){Reports::Validate();return 0;}
    Reports::Tables();auto runtime=Reports::Runtime();Reports::Provenance(runtime);Reports::Handoff(runtime);Reports::Vault(runtime);
    std::cout<<"PASS: generated LaTeX tables, ROOT/C++ HTML handoff, isolated provenance and vault receipt\n";
}catch(const std::exception&e){std::cerr<<"FAIL: "<<e.what()<<'\n';return 1;}return 0;}
