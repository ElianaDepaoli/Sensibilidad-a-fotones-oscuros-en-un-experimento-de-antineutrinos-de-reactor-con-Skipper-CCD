// Run from work/reactor_sensitivity. Reuse the audited ROOT kernels without
// executing the historical entry point. All new physics is in named functions.
#define DARKPHOTON_KERNEL_ONLY
#include "../../code/darkphoton_v3/redo.cpp"
#undef DARKPHOTON_KERNEL_ONLY

namespace Sensitivity {
using namespace V3;
// Park Eq. (4): low-mass electron-loop approximation, energies in MeV.
double Width3Gamma(double mass,double eps){
    Require(mass>=0&&mass<2*electron&&eps>=0,"width outside sub-pair domain");
    return 2.16e-16*e4*eps*eps*std::pow(mass,9)/std::pow(electron,8);
}
// Exact beta*gamma in the propagation exponent, hbar*c in MeV m.
double OpticalDepth(double E,double mass,double eps,double distance){
    Require(E>mass&&distance>=0,"invalid propagation kinematics");
    return distance*mass*Width3Gamma(mass,eps)/(1.973269804e-13*std::sqrt(E*E-mass*mass));
}
double Survival(double E,double mass,double eps,double distance){return std::exp(-OpticalDepth(E,mass,eps,distance));}
// Cache the epsilon-independent quadrature weights for repeated limit solving.
struct VacuumKernel {std::vector<double> energies,weights;};
VacuumKernel BuildKernel(double mass,const Config& c,const Quadrature& q){
    Source source;source.powerMW=c.powerMW;source.gammaMin=1;source.gammaMax=c.gammaHi;
    VacuumKernel k;double half=(c.roiHi-c.roiLo)/2,mid=(c.roiHi+c.roiLo)/2;
    for(size_t i=0;i<q.nodes.size();++i){double E=mid+half*q.nodes[i];k.energies.push_back(E);
        k.weights.push_back(half*q.weights[i]*Redo::Geometry(c)*Spectrum(E,mass,source,q)*Redo::DetectionTotal(E,mass,q));}
    return k;
}
double VacuumEvents(const VacuumKernel& k,double mass,double eps,const Config& c,bool decay){
    double n=0;for(size_t i=0;i<k.weights.size();++i)n+=k.weights[i]*(decay?Survival(k.energies[i],mass,eps,c.baseline_m):1);
    return c.eta*std::pow(eps,4)*n;
}
double MediumEvents(double massEV,double eps,const Config& c,const Quadrature& q){
    return Redo::Geometry(c)*c.eta*q.Integrate([&](double E){return V3::PhotonSpectrum(E,c)*Redo::Prod(E,massEV,eps,c)*Redo::Det(massEV,eps,c)*KleinNishina(E)*Survival(E,massEV*1e-6,eps,c.baseline_m);},c.roiLo,c.roiHi);
}
// Find the first upward crossing. Decay can make the global signal nonmonotonic;
// no statement is made beyond the explicitly bounded small-mixing interval.
template<class F> double LowerBoundary(F events,double cap){
    double lo=1e-14,prev=events(lo);
    for(int j=1;j<=144;++j){double hi=std::pow(10.,-14+j/12.);double next=events(hi);
        if(prev<cap&&next>=cap){for(int z=0;z<70;++z){double mid=std::sqrt(lo*hi);if(events(mid)>=cap)hi=mid;else lo=mid;}return std::sqrt(lo*hi);}
        lo=hi;prev=next;
    }return NAN;
}
void Checks(const Config& c,const Quadrature& q){
    Close(Survival(3,.1,0,28),1,1e-15,"zero mixing survival");
    Close(Survival(3,.1,.001,0),1,1e-15,"zero distance survival");
    Close(Width3Gamma(.1,2e-5)/Width3Gamma(.1,1e-5),4,1e-14,"width mixing scaling");
    Close(Width3Gamma(.2,1e-5)/Width3Gamma(.1,1e-5),512,1e-14,"width mass scaling");
    Close(LowerBoundary([](double e){return 1e22*std::pow(e,4);},100),1e-5,1e-12,"analytic root");
    Require(std::isnan(LowerBoundary([](double){return 0.;},100)),"missing root not marked");
    auto k=BuildKernel(.001,c,q);auto fine=BuildKernel(.001,c,Quadrature(c.nodes*2));
    Close(VacuumEvents(k,.001,1e-5,c,true),VacuumEvents(fine,.001,1e-5,c,true),1e-7,"vacuum quadrature");
    Close(MediumEvents(1000,1e-5,c,q),Redo::Events(1000,1e-5,c,q,c.eta),1e-12,"negligible low-mass decay");
    std::cout<<"PASS: decay scaling, zero-distance/mixing survival, analytic and missing roots, doubled vacuum quadrature, medium zero-decay comparison\n";
}
void Calculate(const Config& c,const Quadrature& q){
    std::ofstream csv("output/limits.csv");csv<<"mass_MeV,vacuum_no_decay,vacuum_with_decay,medium_with_decay,max_optical_depth_at_vacuum_limit\n"<<std::setprecision(16);
    TGraph vacuum,mediumLow,mediumHigh;double maximum=0;
    for(int i=0;i<=70;++i){double m=std::pow(10.,-8+8.*i/70);if(i==70)m=.999;
        auto k=BuildKernel(m,c,q);double e0=LowerBoundary([&](double e){return VacuumEvents(k,m,e,c,false);},c.N95);
        double ed=LowerBoundary([&](double e){return VacuumEvents(k,m,e,c,true);},c.N95);
        double em=NAN;if(m<=.01&&!Redo::Resonance(m*1e6,c))em=LowerBoundary([&](double e){return MediumEvents(m*1e6,e,c,q);},c.N95);
        Require(std::isfinite(ed),"missing vacuum boundary");Close(VacuumEvents(k,m,ed,c,true),c.N95,1e-11,"vacuum root closure");
        double tau=OpticalDepth(c.roiLo,m,ed,c.baseline_m);maximum=std::max(maximum,tau);
        vacuum.SetPoint(vacuum.GetN(),m*1e6,ed);
        if(std::isfinite(em)){Close(MediumEvents(m*1e6,em,c,q),c.N95,1e-11,"medium root closure");auto& g=m*1e6<c.reactorM?mediumLow:mediumHigh;g.SetPoint(g.GetN(),m*1e6,em);}
        csv<<m<<','<<e0<<','<<ed<<','<<em<<','<<tau<<'\n';
    }
    TCanvas canvas("decay_limit","Conditional sensitivity",1000,650);canvas.SetLogx();canvas.SetLogy();canvas.SetLeftMargin(.13);
    canvas.DrawFrame(.01,1e-7,1e6,.02)->SetTitle(";m_{A'} [eV];#epsilon, conditional 95% lower boundary");
    vacuum.SetLineColor(kBlue+1);vacuum.SetLineWidth(2);vacuum.Draw("L SAME");
    for(auto* g:{&mediumLow,&mediumHigh}){g->SetLineColor(kGreen+2);g->SetLineWidth(2);g->Draw("L SAME");Redo::CheckGraph(*g,.01,1e6,1e-7,.02,"medium survival");}
    TLegend legend(.19,.16,.87,.35);legend.SetBorderSize(0);legend.AddEntry(&vacuum,"Vacuum Compton hypothesis + survival","l");legend.AddEntry(&mediumHigh,"Danilov medium + survival (validity restricted)","l");legend.SetHeader("Same selected-event cap and efficiency");legend.Draw();Redo::Save(canvas,"conditional_limits",3);
    std::cout<<std::setprecision(12)<<"maximum_EH_optical_depth_at_vacuum_limit = "<<maximum<<'\n';
    std::cout<<"PASS: every reported boundary closes the configured event cap\nUNCHECKED: exact three-photon loop near threshold; detector response; full Park Figure 2 compilation; published normalization discrepancy\n";
}
int Run(const std::string& path){auto start=std::chrono::steady_clock::now();gROOT->SetBatch(true);gStyle->SetOptStat(0);auto c=ReadConfig(path);Require(c.roiLo>1&&c.roiHi>c.roiLo,"this sub-MeV scan requires ROI lower edge above 1 MeV");Quadrature q(c.nodes);Checks(c,q);Calculate(c,q);std::cout<<"ROOT "<<gROOT->GetVersion()<<"; compiler "<<__VERSION__<<"\ncompute_seconds = "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<'\n';return 0;}
}
int main(int argc,char**argv){try{return Sensitivity::Run(argc>1?argv[1]:"config.conf");}catch(const std::exception& e){std::cerr<<"FAIL: "<<e.what()<<'\n';return 1;}}
