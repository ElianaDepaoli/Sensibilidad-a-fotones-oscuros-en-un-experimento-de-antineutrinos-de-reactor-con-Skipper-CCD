// ROOT 6.36.000 / C++17. Cross sections and rates use MeV, cm, s, and kg.
// Build/run instructions and configurable inputs are in run.sh and texono.conf.
#include "matrix_element_generated.h"
#include <TF1.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TH1.h>
#include <TROOT.h>
#include <TStyle.h>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace Park {
constexpr double me=0.510998950, alpha=1/137.035999084;
constexpr double hbarc=1.973269804e-11; // MeV cm
constexpr double hbarc_m=1.973269804e-13; // MeV m
constexpr double NA=6.02214076e23, pi=3.14159265358979323846;
constexpr double e4=16*pi*pi*alpha*alpha;

struct Config {
    double power_MW=2900, baseline_m=28, mass_kg=187, days=160;
    double Z=108, A=259.8099219; // CsI formula-unit electron count and g/mol
    double roi_low=3, roi_high=8, observed=414, background=414;
    double sigma_events=100.6, cl_multiplier=1.96, efficiency=1;
    double gamma_low=1, gamma_high=40, reference_limit=2.1e-5;
    int quadrature=96;
};

Config ReadConfig(const std::string& path) {
    // Keep inputs in plain key=value form so another target needs no code edit.
    Config c;
    std::map<std::string,double*> fields={{"power_MW",&c.power_MW},{"baseline_m",&c.baseline_m},
      {"mass_kg",&c.mass_kg},{"days",&c.days},{"Z",&c.Z},{"A",&c.A},
      {"roi_low",&c.roi_low},{"roi_high",&c.roi_high},{"observed",&c.observed},
      {"background",&c.background},{"sigma_events",&c.sigma_events},
      {"cl_multiplier",&c.cl_multiplier},{"efficiency",&c.efficiency},
      {"gamma_low",&c.gamma_low},{"gamma_high",&c.gamma_high},
      {"reference_limit",&c.reference_limit}};
    std::ifstream in(path); if(!in) throw std::runtime_error("Missing config: "+path);
    std::string line;
    while(std::getline(in,line)) {
        line=line.substr(0,line.find('#')); auto at=line.find('=');
        if(at==std::string::npos) continue;
        std::string key=line.substr(0,at);
        key.erase(std::remove_if(key.begin(),key.end(),::isspace),key.end());
        if(key=="quadrature") c.quadrature=std::stoi(line.substr(at+1));
        else if(fields.count(key)) *fields.at(key)=std::stod(line.substr(at+1));
        else throw std::runtime_error("Unknown input: "+key);
    }
    if(c.baseline_m<=0||c.mass_kg<=0||c.days<=0||c.Z<=0||c.A<=0||c.power_MW<=0||
       c.roi_low<=0||c.roi_high<=c.roi_low||c.gamma_high<=c.gamma_low||
       c.efficiency<=0||c.efficiency>1||c.sigma_events<=0||c.quadrature<16)
        throw std::runtime_error("Unphysical configuration");
    return c;
}

struct Gauss {
    // ROOT provides the nodes/weights; this small wrapper integrates lambdas.
    std::vector<double> x,w;
    explicit Gauss(int n):x(n),w(n) { TF1::CalcGaussLegendreSamplingPoints(n,x.data(),w.data(),1e-14); }
    template<class F> double Integral(F f,double lo,double hi) const {
        if(hi<=lo) return 0;
        double sum=0, mid=(lo+hi)/2, half=(hi-lo)/2;
        for(size_t i=0;i<x.size();++i) sum+=w[i]*f(mid+half*x[i]);
        return half*sum;
    }
};

double KN(double E) {
    // Relativistic total Klein-Nishina cross section per free electron, cm^2.
    const double x=E/me, re2=alpha*alpha*hbarc*hbarc/(me*me);
    if(x<1e-4) return 8*pi*re2/3*(1-2*x+26*x*x/5);
    const double l=std::log1p(2*x);
    return 2*pi*re2*((1+x)/(x*x)*(2*(1+x)/(1+2*x)-l/x)+l/(2*x)-(1+3*x)/std::pow(1+2*x,2));
}

double GammaSpectrum(double E,double MW) {
    // Park Eq. (3), photons/(s MeV); empirical input, not derived here.
    return 0.58e18*MW*std::exp(-E/0.91);
}

std::pair<double,double> ProductionBounds(double E,double M) {
    // Lorentz-boost the two-body CM endpoints to the stationary-electron lab.
    const double s=me*me+2*me*E;
    if(E<=M+M*M/(2*me)) return {0,0};
    const double K=std::sqrt((s-std::pow(me+M,2))*(s-std::pow(me-M,2)));
    const double center=(E+me)*(s+M*M-me*me)/(2*s), half=E*K/(2*s);
    return {center-half,center+half};
}

double ProductionDS(double E,double EA,double M) {
    // d sigma(gamma e -> A' e)/d EA, cm^2/MeV, with epsilon stripped off.
    const auto [lo,hi]=ProductionBounds(E,M);
    if(EA<lo||EA>hi||hi<=lo) return 0;
    const double s=me*me+2*me*E, u=me*me+M*M-2*me*EA;
    return e4*TraceFactor(s,u,me,M*M)*hbarc*hbarc/(32*pi*me*E*E);
}

double Reference15Total(double E,double M) {
    // Independent test only: Ref. 15 Appendix (A.1-A.3), vector branch.
    double s=me*me+2*me*E, root=std::sqrt(s), d=s-me*me;
    if(E<=M+M*M/(2*me)) return 0;
    double p0=(s-me*me+M*M)/(2*root), p=std::sqrt(p0*p0-M*M);
    double k0=(s+me*me)/(2*root), k=root-k0;
    double A=2+2*(me*me-M*M)/s+16*(M*M+2*me*me)*s/(d*d);
    double B=2-4*(M*M+2*me*me)/d-4*(4*std::pow(me,4)-std::pow(M,4))/(d*d);
    return pi*alpha*alpha/(2*s)*(p/k)*(A+B*root/p*std::log((2*p0*k0+2*p*k-M*M)/(2*p0*k0-2*p*k-M*M)))*hbarc*hbarc;
}

std::pair<double,double> DetectionBounds(double EA,double M) {
    // Outgoing real-photon energy endpoints; they include the vector's rest mass.
    if(EA<=M) return {0,0};
    double p=std::sqrt(EA*EA-M*M), a=M*M+2*me*EA;
    return {a/(2*(me+EA+p)),a/(2*(me+EA-p))};
}

double DetectionDS(double EA,double w,double M) {
    // Crossed trace: spin average 1/6 instead of 1/4, exactly Park's 2/3.
    auto [lo,hi]=DetectionBounds(EA,M);
    if(w<lo||w>hi||hi<=lo) return 0;
    double s=me*me+M*M+2*me*EA,u=me*me-2*me*w;
    return (2./3)*e4*TraceFactor(s,u,me,M*M)*hbarc*hbarc/(32*pi*me*(EA*EA-M*M));
}

double DetectionTotal(double EA,double M,const Gauss& q,bool printed_bound=false) {
    auto [lo,hi]=DetectionBounds(EA,M);
    // Optional diagnostic: the printed Eq. (7) lower bound has an extra factor 2.
    if(printed_bound) lo=std::max(lo,2*me*EA/(me+2*EA));
    return q.Integral([&](double w){return DetectionDS(EA,w,M);},lo,hi);
}

double Flux(double EA,double M,const Config& c,const Gauss& q) {
    // Park Eq. (1). Integrate exactly over allowed incoming gamma energies.
    if(EA<=M) return 0;
    double p=std::sqrt(EA*EA-M*M), n=me*EA-M*M/2;
    double d=me-EA+p; if(d<=0) return 0;
    double low=std::max({c.gamma_low,n/d,M+M*M/(2*me)});
    double high=c.gamma_high;
    if(me-EA-p>0) high=std::min(high,n/(me-EA-p));
    return q.Integral([&](double E){return GammaSpectrum(E,c.power_MW)*ProductionDS(E,EA,M)/KN(E);},low,high);
}

double ElectronCount(const Config& c) { return c.mass_kg*1000*NA*c.Z/c.A; }
double CountLimit(const Config& c) {
    // Park prescription. It uses measured uncertainty, not sqrt(background).
    // Optional positive observed-minus-background residual allows a simple recast.
    return std::max(0.,c.observed-c.background)+c.cl_multiplier*c.sigma_events;
}
double Events(double M,const Config& c,const Gauss& q,bool printed=false) {
    // Full containment: ROI applies to EA = photon energy + recoil kinetic energy.
    double rate=q.Integral([&](double EA){return Flux(EA,M,c,q)*DetectionTotal(EA,M,q,printed);},
                          std::max(c.roi_low,M),c.roi_high);
    return ElectronCount(c)*(c.days*86400)*c.efficiency*rate/(4*pi*std::pow(c.baseline_m*100,2));
}
double EpsilonLimit(double coefficient,const Config& c) { return std::pow(CountLimit(c)/coefficient,0.25); }

double EventsReordered(double M,const Config& c,const Gauss& q) {
    // Independent integration order for checking the source-to-event convolution.
    auto kernel=[&](double E) {
        auto [lo,hi]=ProductionBounds(E,M);lo=std::max(lo,c.roi_low);hi=std::min(hi,c.roi_high);
        return GammaSpectrum(E,c.power_MW)/KN(E)*q.Integral([&](double EA) {
            return ProductionDS(E,EA,M)*DetectionTotal(EA,M,q);
        },lo,hi);
    };
    double rate=0;std::vector<double> edges={c.gamma_low,c.roi_low,c.roi_high,c.gamma_high};
    std::sort(edges.begin(),edges.end());
    for(size_t i=1;i<edges.size();++i)rate+=q.Integral(kernel,std::max(c.gamma_low,edges[i-1]),std::min(c.gamma_high,edges[i]));
    return ElectronCount(c)*c.days*86400*c.efficiency*rate/(4*pi*std::pow(c.baseline_m*100,2));
}

double Width3Gamma(double M,double eps) { return 2.16e-16*e4*eps*eps*std::pow(M,9)/std::pow(me,8); }
double DecayLength(double E,double M,double eps) { return hbarc_m*std::sqrt(E*E-M*M)/(M*Width3Gamma(M,eps)); }

void Near(double a,double b,double tol,const std::string& label) {
    if(std::abs(a/b-1)>tol) throw std::runtime_error("Failed "+label+": "+std::to_string(a/b));
}

void Tests(const Config& c,const Gauss& q) {
    // Test all physical normalization factors before attempting a paper match.
    for(double E:{0.1,1.,3.,8.,15.}) {
        auto [lo,hi]=ProductionBounds(E,0);
        Near(q.Integral([&](double w){return ProductionDS(E,w,0);},lo,hi),KN(E),1e-8,"KN production");
        Near(DetectionTotal(E,0,q),(2./3)*KN(E),1e-8,"Park detection 2/3");
        for(double M:{0.1,0.5,1.}) {
            auto [a,b]=ProductionBounds(E,M); if(b<=a) continue;
            double production=q.Integral([&](double w){return ProductionDS(E,w,M);},a,b);
            Near(production,Reference15Total(E,M),1e-8,"Reference 15");
            double s=me*me+2*me*E,lambda=(s-std::pow(me+M,2))*(s-std::pow(me-M,2));
            double inverse_energy=(s-me*me-M*M)/(2*me);
            Near(DetectionTotal(inverse_energy,M,q)/production,
                 (2./3)*std::pow(s-me*me,2)/lambda,1e-8,"finite-mass detailed balance");
        }
    }
    Near(CountLimit(c),c.cl_multiplier*c.sigma_events+std::max(0.,c.observed-c.background),1e-12,"count limit");
    std::cout<<"PASS: finite-mass trace integrals match Reference 15; production KN and detection 2/3 KN limits; finite-mass detailed balance\n";
}

void SaveFigure1(const Config& c,const Gauss& q) {
    Config one=c; one.power_MW=1000;
    TCanvas canvas("fig1","Park Figure 1: independent calculation",820,510);
    canvas.SetLogy(); canvas.SetLeftMargin(.13);
    auto frame=canvas.DrawFrame(0,.0015,4.3,3.);
    frame->SetTitle(";E_{A'} [MeV];dN_{A'}/dE_{A'} [10^{21} MeV^{-1} s^{-1}]");
    TLegend leg(.64,.70,.90,.89);leg.SetBorderSize(0);
    std::vector<TGraph*> curves;
    std::ofstream csv("output/figure1.csv");csv<<"EA_MeV,mass_MeV,flux_per_MeV_s\n"<<std::setprecision(15);
    int j=0;
    for(double M:{0.1,0.5,1.}) {
        auto g=new TGraph();curves.push_back(g);
        for(int i=0;i<=210;++i) {
            double EA=i*.02,v=Flux(EA,M,one,q);
            csv<<EA<<','<<M<<','<<v<<'\n'; if(v>0) g->SetPoint(g->GetN(),EA,v/1e21);
        }
        g->SetLineColor(j==0?kBlue:kBlack);g->SetLineStyle(j==2?3:1);g->SetLineWidth(2);g->Draw("L SAME");
        leg.AddEntry(g,Form("m_{A'} = %.1f MeV",M),"l");++j;
    }
    leg.Draw();canvas.SaveAs("output/figure1.svg");canvas.SaveAs("output/figure1.pdf");
    for(auto g:curves) delete g;
}

void SaveFigure2(const Config& c,const Gauss& q) {
    // Only the TEXONO component of Park Fig. 2 is in scope. No imported exclusions.
    TCanvas canvas("fig2","TEXONO component of Park Figure 2",900,540);
    auto frame=canvas.DrawFrame(-18,-18,6,0);
    frame->SetTitle(";log_{10}(m_{A'} / eV);log_{10}(#epsilon)");
    TGraph original,computed;
    std::ofstream csv("output/figure2.csv");csv<<"mass_MeV,epsilon_limit,park_limit\n"<<std::setprecision(15);
    for(int i=0;i<=60;++i) {
        double logm=-18+24.*i/60, M=std::pow(10.,logm-6);
        double limit=EpsilonLimit(Events(M,c,q),c);
        original.SetPoint(i,logm,std::log10(c.reference_limit));computed.SetPoint(i,logm,std::log10(limit));
        csv<<M<<','<<limit<<','<<c.reference_limit<<'\n';
    }
    original.SetLineColor(kRed);original.SetLineStyle(3);original.SetLineWidth(2);original.Draw("L SAME");
    computed.SetLineColor(kBlue+1);computed.SetLineWidth(2);computed.Draw("L SAME");
    TLegend leg(.13,.70,.68,.87);leg.SetBorderSize(0);
    leg.AddEntry(&original,"Park reported TEXONO line (reference)","l");
    leg.AddEntry(&computed,"Independent integration; no later corrections","l");leg.Draw();
    canvas.SaveAs("output/figure2_texono.svg");canvas.SaveAs("output/figure2_texono.pdf");
}

void SaveCrossSections(const Config& c,const Gauss& q) {
    std::ofstream csv("output/cross_sections.csv");
    csv<<"E_MeV,mass_MeV,production_cm2,detection_cm2,KN_cm2\n"<<std::setprecision(15);
    TCanvas cv("cross","Production and detection",1000,470);cv.Divide(2,1);
    std::vector<TGraph> prod(4),det(4);int j=0;
    for(double M:{0.,0.1,0.5,1.}) {for(int i=1;i<=120;++i) {
        double E=.1*i;auto [lo,hi]=ProductionBounds(E,M);
        double p=q.Integral([&](double w){return ProductionDS(E,w,M);},lo,hi),d=DetectionTotal(E,M,q);
        csv<<E<<','<<M<<','<<p<<','<<d<<','<<KN(E)<<'\n';
        if(p>0)prod[j].SetPoint(prod[j].GetN(),E,p/1e-24);
        if(d>0)det[j].SetPoint(det[j].GetN(),E,d/1e-24);
    }++j;}
    TLegend lp(.46,.61,.87,.88),ld(.46,.61,.87,.88);lp.SetBorderSize(0);ld.SetBorderSize(0);
    for(int panel=1;panel<=2;++panel) {
        cv.cd(panel)->SetLogy();auto frame=gPad->DrawFrame(.1,.001,12,1.);
        frame->SetTitle(panel==1?"Production;E_{#gamma} [MeV];#sigma/#epsilon^{2} [barn/electron]":"Detection;E_{A'} [MeV];#sigma/#epsilon^{2} [barn/electron]");
        auto& curves=panel==1?prod:det;auto& leg=panel==1?lp:ld;
        for(int k=0;k<4;++k) {
            curves[k].SetLineColor(k==0?kGray+2:(k==1?kBlue:(k==2?kBlack:kRed)));
            curves[k].SetLineWidth(2);curves[k].Draw("L SAME");
            leg.AddEntry(&curves[k],Form("m_{A'} = %.1f MeV",std::vector<double>{0,.1,.5,1}[k]),"l");
        }leg.Draw();
    }
    cv.SaveAs("output/cross_sections.svg");cv.SaveAs("output/cross_sections.pdf");
}

void MaterialDiagnostic() {
    // Diagnose sigma_tot ~= sigma_C. Never feed these data into Park's rate.
    std::ifstream source("nist_uranium.csv"); std::string line;
    std::ofstream csv("output/compton_validity.csv");
    csv<<"E_MeV,KN_cm2_g,total_cm2_g,KN_fraction_total,Thomson_over_KN\n"<<std::setprecision(15);
    const double thomson=8*pi*alpha*alpha*hbarc*hbarc/(3*me*me);
    TGraph kn,tot;
    while(std::getline(source,line)) {
        if(line.empty()||line[0]=='#'||line[0]=='E')continue;
        auto comma=line.find(',');double E=std::stod(line.substr(0,comma)),mu=std::stod(line.substr(comma+1));
        double comp=NA*.38651*KN(E);
        csv<<E<<','<<comp<<','<<mu<<','<<comp/mu<<','<<thomson/KN(E)<<'\n';
        kn.SetPoint(kn.GetN(),E,comp);tot.SetPoint(tot.GetN(),E,mu);
    }
    TCanvas cv("material","Uranium total attenuation diagnostic",850,520);cv.SetLogx();cv.SetLogy();
    auto frame=cv.DrawFrame(.1,.008,10,3);frame->SetTitle(";E_{#gamma} [MeV];#mu/#rho [cm^{2}/g], uranium");
    kn.SetLineColor(kBlue);kn.SetLineWidth(2);kn.Draw("L SAME");
    tot.SetLineColor(kRed);tot.SetLineWidth(2);tot.Draw("LP SAME");
    TLegend leg(.48,.65,.87,.87);leg.SetBorderSize(0);
    leg.AddEntry(&kn,"Free-electron KN times N_{A}Z/A","l");leg.AddEntry(&tot,"NIST total attenuation","l");leg.Draw();
    cv.SaveAs("output/compton_validity.svg");cv.SaveAs("output/compton_validity.pdf");
}

void CompareFigure1(const Config& c,const Gauss& q) {
    // Compare actual PDF staircase coordinates without a fitted rescaling.
    std::ifstream in("output/park_figure1_digitized.csv");if(!in)return;
    Config one=c;one.power_MW=1000;
    std::ofstream out("output/figure1_comparison.csv");
    out<<"EA_MeV,mass_MeV,paper_flux,calculated_flux,ratio\n"<<std::setprecision(15);
    std::string line;std::getline(in,line);
    std::map<double,std::vector<double>> errors;
    std::map<double,TGraph> original;
    while(std::getline(in,line)) {
        std::replace(line.begin(),line.end(),',',' ');std::istringstream row(line);
        double E,M,f;row>>E>>M>>f;
        if(E<1.8||E>3.2) continue; // Avoid thresholds and noisy last bins.
        double calc=Flux(E,M,one,q);out<<E<<','<<M<<','<<f<<','<<calc<<','<<calc/f<<'\n';
        errors[M].push_back(std::abs(calc/f-1));original[M].SetPoint(original[M].GetN(),E,f/1e21);
    }
    std::ofstream summary("output/figure1_status.json");summary<<"{\n";int n=0;
    for(auto& [mass,errs]:errors) {
        std::sort(errs.begin(),errs.end());double median=errs[errs.size()/2];
        if(n++)summary<<",\n";
        summary<<'"'<<mass<<"\": {\"median_absolute_relative_error\":"<<median
          <<",\"comparison_pass_10percent\":"<<(median<.1?"true":"false")<<'}';
        std::cout<<"Figure 1 m="<<mass<<" MeV: median relative residual="<<median<<"; "<<(median<.1?"PASS":"FAIL")<<" (10% diagnostic tolerance)\n";
    }
    summary<<"\n}\n";
    TCanvas cv("comparison","Park Figure 1 comparison",850,520);cv.SetLogy();
    auto frame=cv.DrawFrame(0,.0015,4.3,3);frame->SetTitle(";E_{A'} [MeV];dN_{A'}/dE_{A'} [10^{21} MeV^{-1} s^{-1}]");
    TLegend leg(.52,.61,.89,.88);leg.SetBorderSize(0);std::vector<TGraph> calculated(3);int j=0;
    for(double M:{.1,.5,1.}) {
        int color=j==0?kBlue:(j==1?kBlack:kRed);
        // All digitized curve coordinates are loaded again for the full overlay.
        std::ifstream src("output/park_figure1_digitized.csv");std::getline(src,line);TGraph& ref=original[M];ref.Set(0);
        while(std::getline(src,line)) {
            std::replace(line.begin(),line.end(),',',' ');std::istringstream row(line);double E,mi,f;row>>E>>mi>>f;
            if(std::abs(mi-M)<1e-9)ref.SetPoint(ref.GetN(),E,f/1e21);
        }
        ref.SetLineColor(color);ref.SetLineStyle(3);ref.Draw("L SAME");
        auto& g=calculated[j];for(int i=0;i<=215;++i) {
            double E=.02*i,v=Flux(E,M,one,q);if(v>0)g.SetPoint(g.GetN(),E,v/1e21);
        }
        g.SetLineColor(color);g.SetLineWidth(2);g.Draw("L SAME");
        leg.AddEntry(&g,Form("Derived %.1f MeV",M),"l");leg.AddEntry(&ref,Form("Park PDF %.1f MeV",M),"l");++j;
    }
    leg.Draw();cv.SaveAs("output/figure1_comparison.svg");cv.SaveAs("output/figure1_comparison.pdf");
}

void Results(const Config& c,const Gauss& q) {
    const double coeff=Events(.001,c,q), eps=EpsilonLimit(coeff,c);
    Gauss fine(2*c.quadrature);double finer=Events(.001,c,fine);
    Near(coeff,finer,2e-5,"quadrature convergence");
    double swapped=EventsReordered(.001,c,fine);
    Near(coeff,swapped,2e-4,"independent integration order");
    Config longer=c;longer.gamma_high=60;
    Near(coeff,Events(.001,longer,q),2e-5,"gamma tail convergence");
    Config scaled=c;scaled.power_MW*=16;
    Near(EpsilonLimit(Events(.001,scaled,q),scaled),eps/2,1e-10,"fourth-root power scaling");
    const double fluxint=0.58e21*.91*std::exp(-1/.91);
    Near(fluxint,1.76e20,.004,"printed photon integral");
    const double n95=1.96*100.6;Near(n95,197.2,.0002,"printed TEXONO count cap");
    Near(hbarc_m*std::pow(me,8)/(2.16e-16*e4),505,.0002,"decay-length prefactor");
    const bool match=std::abs(eps/c.reference_limit-1)<.025;
    std::ofstream json("output/results.json");json<<std::setprecision(16)<<"{\n";
    auto number=[&](const char* key,double value){json<<'"'<<key<<"\": "<<value<<",\n";};
    number("electron_count",ElectronCount(c));number("N95",CountLimit(c));number("coefficient_eps4",coeff);
    number("epsilon_limit",eps);number("paper_epsilon_limit",c.reference_limit);
    number("epsilon_printed_bound",EpsilonLimit(Events(.001,c,q,true),c));
    number("epsilon_prompt_flux",eps*std::pow(fluxint/6.82e19,.25));
    number("epsilon_prompt_flux_printed_bound",EpsilonLimit(Events(.001,c,q,true),c)*std::pow(fluxint/6.82e19,.25));
    number("events_at_paper_epsilon",coeff*std::pow(c.reference_limit,4));
    number("epsilon_mass_01",EpsilonLimit(Events(.1,c,q),c));
    number("epsilon_mass_05",EpsilonLimit(Events(.5,c,q),c));
    number("epsilon_mass_1",EpsilonLimit(Events(1,c,q),c));
    number("photon_integral_above_1",fluxint);number("prompt_flux_quoted",6.82e19);
    number("flux_ratio",fluxint/6.82e19);number("epsilon_flux_factor",std::pow(fluxint/6.82e19,.25));
    number("decay_prefactor_from_eq4_m",hbarc_m*std::pow(me,8)/(2.16e-16*e4));
    number("decay_prefactor_printed_m",505);
    number("decay_length_E3_M01_epslimit_m",DecayLength(3,.1,eps));
    number("quadrature_relative_difference",std::abs(finer/coeff-1));
    number("integration_order_relative_difference",std::abs(swapped/coeff-1));
    number("implicit_acceptance_to_match",CountLimit(c)/(coeff*std::pow(c.reference_limit,4)));
    json<<"\"paper_limit_reproduced_within_2p5_percent\": "<<(match?"true":"false")<<"\n}\n";
    std::cout<<std::setprecision(10)<<"N95 = "<<CountLimit(c)<<"\nN_e = "<<ElectronCount(c)
      <<"\nN(epsilon=1) = "<<coeff<<"\nepsilon95 = "<<eps<<"\nPark printed epsilon95 = "<<c.reference_limit
      <<"\nPaper limit comparison = "<<(match?"PASS":"FAIL (reported, no normalization fitted)")
      <<"\nPASS: quadrature, gamma-tail, photon integral, count cap, and power scaling checks\n";
}
}

int main(int argc,char** argv) {
    try {
        gROOT->SetBatch(true);gStyle->SetOptStat(0);
        auto c=Park::ReadConfig(argc>1?argv[1]:"texono.conf");Park::Gauss q(c.quadrature);
        std::cout<<"ROOT "<<gROOT->GetVersion()<<"; C++ "<<__VERSION__<<'\n';
        Park::Tests(c,q);Park::Results(c,q);Park::SaveFigure1(c,q);
        Park::SaveFigure2(c,q);Park::SaveCrossSections(c,q);Park::MaterialDiagnostic();Park::CompareFigure1(c,q);
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
