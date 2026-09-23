// Compile with g++ -O2 -std=c++17 production.cpp $(root-config --cflags --libs) -lXMLIO -o production
// Run from this directory. Every computation is a named, documented function.
#include "physics.hxx"
#include <TCanvas.h>
#include <TGraph.h>
#include <TH1.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TXMLEngine.h>
#include <fstream>
#include <iomanip>
#include <regex>
#include <filesystem>
#include <chrono>

namespace Task {
using namespace Production;
inline std::string Read(const std::string&path){std::ifstream f(path);Require(bool(f),"Cannot read "+path);return std::string(std::istreambuf_iterator<char>(f),{});}
void Checks(){
    const auto&r=Symbolic();
    std::ofstream algebra("output/trace_terms.tex");
    algebra<<"\\begin{align}\nT_{ss}&="<<r.ss.Tex()<<",\\\\\nT_{uu}&="<<r.uu.Tex()<<",\\\\\nT_{su}=T_{us}&="<<r.su.Tex()<<".\n\\end{align}\n";
    std::cout<<"PASS: exact polynomial compact formula, analytic primitive derivative, photon limit, interference equality and crossing symmetry\n";
    CheckMatrices();std::cout<<"PASS: Clifford algebra, on-shell kinematics, both Ward identities, physical polarizations and independent matrix traces\n";
    Quadrature q(96),fine(192);double worstReference=0,worstKN=0,worstConvergence=0,worstOrder=0;
    for(double energy:{.02,.2,1.,2.2,4.,8.,20.}){
        double kn=Total(energy,0,q);Close(kn,KleinNishina(energy),1e-9,"Integrated Klein-Nishina failed");
        worstKN=std::max(worstKN,std::abs(kn/KleinNishina(energy)-1));
        for(double mass:{.1,.5,1.})if(energy>Threshold(mass)){
            double total=Total(energy,mass,q),reference=Reference15(energy,mass);
            Close(total,reference,2e-8,"Reference 15 failed");
            Close(total,AnalyticTotal(energy,mass),2e-8,"Derived analytic integral failed");
            Close(total,Total(energy,mass,fine),2e-8,"Total convergence failed");
            worstReference=std::max(worstReference,std::abs(total/reference-1));
        }
    }
    Source source;Source extended=source;extended.gammaMax=60;
    for(double mass:{.1,.5,1.})for(double out:{1.5,2.,3.,4.}){
        double a=Spectrum(out,mass,source,q),b=Spectrum(out,mass,source,fine);
        Close(a,b,2e-7,"Spectrum quadrature failed");
        Close(a,Spectrum(out,mass,extended,q),2e-7,"Spectrum tail failed");
        worstConvergence=std::max(worstConvergence,std::abs(a/b-1));
    }
    for(double mass:{.1,.5,1.})for(double low:{1.8,2.5,3.1}){
        double inverse=fine.Integrate([&](double x){return Spectrum(x,mass,source,fine);},low,low+.1);
        double forward=ForwardBin(low,low+.1,mass,source,fine);
        Close(inverse,forward,1e-6,"Independent source integration order failed");
        worstOrder=std::max(worstOrder,std::abs(inverse/forward-1));
    }
    double photonIntegral=.58e21*.91*std::exp(-1/.91);
    Close(photonIntegral,1.76e20,.001,"Printed photon normalization failed");
    std::ofstream checks("output/checks.csv");checks<<std::setprecision(16)<<"quantity,value\nreference15_max_relative_error,"<<worstReference<<"\nKN_max_relative_error,"<<worstKN
        <<"\nspectrum_node_max_relative_error,"<<worstConvergence<<"\nsource_order_max_relative_error,"<<worstOrder<<"\nphoton_rate_above1_per_s,"<<photonIntegral<<'\n';
    std::cout<<"PASS: derived analytic totals, Reference 15 totals, Klein-Nishina totals, quadrature, gamma tail and printed photon count\n";
    std::cout<<"reference15_max_relative_error="<<worstReference<<"; KN_max_relative_error="<<worstKN<<"; spectrum_nodes="<<worstConvergence<<'\n';
    std::cout<<"PASS: independent source integration order; maximum relative difference="<<worstOrder<<'\n';
}

struct Bin{double lo,hi,mass,value;};
std::vector<Bin> Digitize(){
    // Fresh extraction from the supplied PDF's SVG; retain horizontal steps once.
    TXMLEngine xml;auto doc=xml.ParseFile("output/paper_page2.svg");Require(doc,"Cannot parse source SVG");
    std::vector<Bin> bins;int curves=0;
    std::function<void(XMLNodePointer_t)> walk=[&](XMLNodePointer_t node){
        for(auto n=node;n;n=xml.GetNext(n)){
            if(std::string(xml.GetNodeName(n))=="path"){
                const char* data=xml.GetAttr(n,"d"),*stroke=xml.GetAttr(n,"stroke"),*transform=xml.GetAttr(n,"transform");
                if(data&&stroke&&transform&&std::strlen(data)>1500){
                    ++curves;double mass=std::string(stroke)=="rgb(0%, 0%, 100%)"?.1:(xml.GetAttr(n,"stroke-dasharray")?1:.5);
                    double a,b,c,d,e,f;Require(std::sscanf(transform,"matrix(%lf, %lf, %lf, %lf, %lf, %lf)",&a,&b,&c,&d,&e,&f)==6,"Bad SVG transform");
                    std::istringstream path(data);char command;double previousX=0,previousY=0;bool have=false;
                    while(path>>command){
                        if(command=='Z'){have=false;continue;}
                        int count=command=='C'?6:2;Require(command=='M'||command=='L'||command=='C',"Unsupported path command");
                        double coords[6];for(int i=0;i<count;++i)path>>coords[i];
                        double x=a*coords[count-2]+c*coords[count-1]+e,y=b*coords[count-2]+d*coords[count-1]+f;
                        // Axis ticks read directly from the newly generated SVG.
                        if(have&&command!='M'&&x>previousX&&std::abs(y-previousY)<1e-5&&y<152.6){
                            double lo=(previousX-380.16015625)/35.5224609375,hi=(x-380.16015625)/35.5224609375;
                            if(lo>=0&&hi<=4.3)bins.push_back({lo,hi,mass,1e21*std::pow(10.,(64.36328125-y)/31.28125)});
                        }
                        previousX=x;previousY=y;have=true;
                    }
                }
            }walk(xml.GetChild(n));
        }
    };walk(xml.DocGetRootElement(doc));xml.FreeDoc(doc);
    Require(curves==3,"Expected exactly three PDF spectrum paths");
    std::ofstream csv("output/paper_bins.csv");csv<<"low_MeV,high_MeV,mass_MeV,rate_per_s_MeV\n"<<std::setprecision(16);
    for(auto r:bins)csv<<r.lo<<','<<r.hi<<','<<r.mass<<','<<r.value<<'\n';
    std::cout<<"Extracted "<<bins.size()<<" horizontal PDF segments from "<<curves<<" spectrum curves\n";return bins;
}

double Median(std::vector<double> v){Require(!v.empty(),"Empty residual set");std::sort(v.begin(),v.end());return v.size()%2?v[v.size()/2]:(v[v.size()/2-1]+v[v.size()/2])/2;}
void PlotSpectra(const std::vector<Bin>&bins,const Source&source){
    Quadrature q(source.nodes),binq(12);Source lower=source;lower.gammaMin=.2;
    TCanvas cv("production_figure","Figure 1 from first principles",1100,650);cv.SetLogy();cv.SetLeftMargin(.14);
    auto frame=cv.DrawFrame(0,.0015,4.3,3.);frame->SetTitle(";E_{A'} [MeV];dN_{A'}/dE_{A'} [10^{21} MeV^{-1} s^{-1}]");
    TLegend legend(.54,.60,.89,.89);legend.SetTextSize(.027);legend.SetBorderSize(0);
    std::array<TGraph,3> calculated,paper;
    std::ofstream csv("output/spectra.csv");csv<<"energy_MeV,mass_MeV,stated_cutoff_rate,lower_cutoff_rate\n"<<std::setprecision(16);
    std::ofstream comparisons("output/comparison.csv");comparisons<<"low_MeV,high_MeV,mass_MeV,paper_rate,stated_cutoff_rate,lower_cutoff_rate\n"<<std::setprecision(16);
    std::ofstream summary("output/comparison_summary.csv");summary<<"mass_MeV,region,bin_count,baseline_median_relative_error,lower_cutoff_median_relative_error,baseline_max_relative_error\n"<<std::setprecision(16);
    int index=0;for(double mass:{.1,.5,1.}){
        int color=index==0?kBlue:index==1?kBlack:kRed+1;
        for(int step=1;step<=430;++step){double energy=step*.01,v=Spectrum(energy,mass,source,q),w=Spectrum(energy,mass,lower,q);
            csv<<energy<<','<<mass<<','<<v<<','<<w<<'\n';if(v>0)calculated[index].SetPoint(calculated[index].GetN(),energy,v/1e21);}
        std::map<std::string,std::vector<double>> errors,lowerErrors;
        for(auto r:bins)if(r.mass==mass){
            double a=binq.Integrate([&](double E){return Spectrum(E,mass,source,q);},r.lo,r.hi)/(r.hi-r.lo);
            double b=binq.Integrate([&](double E){return Spectrum(E,mass,lower,q);},r.lo,r.hi)/(r.hi-r.lo);
            comparisons<<r.lo<<','<<r.hi<<','<<mass<<','<<r.value<<','<<a<<','<<b<<'\n';
            paper[index].SetPoint(paper[index].GetN(),r.lo,r.value/1e21);paper[index].SetPoint(paper[index].GetN(),r.hi,r.value/1e21);
            for(std::string region:{"available","central"})if(region=="available"||(r.lo>=1.8&&r.hi<=3.2)){
                errors[region].push_back(std::abs(a/r.value-1));lowerErrors[region].push_back(std::abs(b/r.value-1));}
        }
        for(auto&[region,residual]:errors){
            summary<<mass<<','<<region<<','<<residual.size()<<','<<Median(residual)<<','<<Median(lowerErrors[region])<<','<<*std::max_element(residual.begin(),residual.end())<<'\n';
            std::cout<<"Figure 1 M="<<mass<<" "<<region<<" median residual="<<Median(residual)<<"; lower cutoff="<<Median(lowerErrors[region])<<'\n';
        }
        calculated[index].SetLineColor(color);calculated[index].SetLineWidth(2);calculated[index].Draw("L SAME");
        paper[index].SetLineColor(color);paper[index].SetLineStyle(3);paper[index].Draw("L SAME");
        legend.AddEntry(&calculated[index],Form("QED M=%.1f MeV",mass),"l");legend.AddEntry(&paper[index],Form("Park PDF M=%.1f MeV",mass),"l");++index;
    }
    legend.Draw();cv.SaveAs("output/figure1_comparison.pdf");cv.SaveAs("output/figure1_comparison.svg");
}

void PlotCrossSection(const Source&source){
    Quadrature q(source.nodes);TCanvas cv("cross_section","Production cross section",900,570);cv.SetLogy();cv.SetLeftMargin(.13);
    auto frame=cv.DrawFrame(.01,1e-3,10.,1.);frame->SetTitle(";E_{#gamma} [MeV];#sigma_{prod}/#epsilon^{2} [barn / electron]");
    TLegend legend(.61,.62,.89,.88);legend.SetBorderSize(0);std::array<TGraph,4> lines;
    std::ofstream csv("output/cross_section.csv");csv<<"energy_MeV,mass_MeV,cross_section_cm2,branching_coefficient\n"<<std::setprecision(16);
    int index=0;for(double mass:{0.,.1,.5,1.}){for(int step=1;step<=250;++step){double energy=.04*step,v=Total(energy,mass,q);csv<<energy<<','<<mass<<','<<v<<','<<v/KleinNishina(energy)<<'\n';if(v>0)lines[index].SetPoint(lines[index].GetN(),energy,v/1e-24);}
        lines[index].SetLineColor(index==0?kGray+2:index==1?kBlue:index==2?kBlack:kRed+1);lines[index].SetLineWidth(2);lines[index].Draw("L SAME");legend.AddEntry(&lines[index],Form("M=%.1f MeV",mass),"l");++index;}
    legend.Draw();cv.SaveAs("output/cross_section.pdf");cv.SaveAs("output/cross_section.svg");
}

void Material(){
    // Fresh NIST table, diagnostic only. No material correction enters the rate.
    std::string page=Read("../../papers/park_production_v2/nist_uranium.html");
    auto pos=page.find("<PRE>");if(pos==std::string::npos)pos=page.find("<pre>");Require(pos!=std::string::npos,"Missing NIST ASCII table");
    page=page.substr(pos);std::regex row(R"((\d\.\d+E[+-]\d+)\s+(\d\.\d+E[+-]\d+)\s+(\d\.\d+E[+-]\d+))");
    TGraph kn,total;std::ofstream csv("output/material.csv");csv<<"energy_MeV,KN_cm2_g,total_cm2_g,KN_fraction,Thomson_over_KN\n"<<std::setprecision(16);
    for(std::sregex_iterator i(page.begin(),page.end(),row),end;i!=end;++i){double E=std::stod((*i)[1]),mu=std::stod((*i)[2]);if(E<.1||E>10)continue;
        double proxy=6.02214076e23*92/238.02891*KleinNishina(E);
        csv<<E<<','<<proxy<<','<<mu<<','<<proxy/mu<<','<<Thomson()/KleinNishina(E)<<'\n';
        kn.SetPoint(kn.GetN(),E,proxy);total.SetPoint(total.GetN(),E,mu);
    }
    Require(kn.GetN()>10,"NIST extraction incomplete");
    TCanvas cv("compton_material","Compton validity",900,560);cv.SetLogx();cv.SetLogy();cv.SetLeftMargin(.13);
    auto frame=cv.DrawFrame(.1,.008,10,3);frame->SetTitle(";E_{#gamma} [MeV];#mu/#rho [cm^{2}/g], uranium");
    kn.SetLineColor(kBlue);kn.SetLineWidth(2);kn.Draw("L SAME");total.SetLineColor(kRed+1);total.SetLineWidth(2);total.Draw("L SAME");
    TLegend legend(.49,.69,.89,.88);legend.SetBorderSize(0);legend.AddEntry(&kn,"N_{A}(Z/A)#sigma_{KN}","l");legend.AddEntry(&total,"NIST: #mu/#rho","l");legend.Draw();
    cv.SaveAs("output/compton_validity.pdf");cv.SaveAs("output/compton_validity.svg");
}

Source Configuration(const std::string&path){
    // Production depends on reactor/source inputs; detector settings do not enter.
    Source result;std::ifstream f(path);Require(bool(f),"Missing configuration");std::string line;
    while(std::getline(f,line)){if(line.empty()||line[0]=='#')continue;auto at=line.find('=');Require(at!=std::string::npos,"Bad config line");auto key=line.substr(0,at);double value=std::stod(line.substr(at+1));
        if(key=="power_MW")result.powerMW=value;else if(key=="gamma_min_MeV")result.gammaMin=value;else if(key=="gamma_max_MeV")result.gammaMax=value;else if(key=="quadrature_nodes")result.nodes=int(value);else if(key=="gamma_normalization")result.normalization=value;else if(key=="gamma_scale_MeV")result.scale=value;else throw std::runtime_error("Unknown configuration key "+key);}
    Require(result.powerMW>0&&result.gammaMin>0&&result.gammaMax>result.gammaMin&&result.nodes>=32&&result.normalization>0&&result.scale>0,"Invalid source configuration");return result;
}
}

int main(int argc,char**argv){
    try{gROOT->SetBatch(true);gStyle->SetOptStat(0);std::cout<<std::setprecision(12)<<"ROOT "<<gROOT->GetVersion()<<"; compiler "<<__VERSION__<<'\n';
        Task::Checks();if(argc>1&&std::string(argv[1])=="--checks")return 0;
        auto source=Task::Configuration(argc>1?argv[1]:"source.conf");auto bins=Task::Digitize();Task::PlotSpectra(bins,source);Task::PlotCrossSection(source);Task::Material();
        std::cout<<"STATUS: source-equation computation completed; full Figure 1 reproduction not established.\n";
    }catch(const std::exception&e){std::cerr<<"FAIL: "<<e.what()<<'\n';return 1;}return 0;
}
