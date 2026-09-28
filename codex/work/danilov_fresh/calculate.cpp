#include "physics.hpp"
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TH1.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TLine.h>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <memory>

namespace Fresh {
// The formula is independently derived in verify_algebra.py, not used to compute the polarized curves.
double InvariantF(const Event&v,double mass){double a=2*Dot(v.p,v.in),b=mass*mass-2*Dot(v.p,v.out),c=me*me,d=mass*mass;return -2*(a/b+b/a)+4*(2*c+d)*(1/a+1/b+c*std::pow(1/a+1/b,2)-d/(a*b));}
// Independently derived closed longitudinal absorption trace, checked against matrices.
double AbsorptionLongitudinalF(const Event&v,double mass){double d=mass*mass,c=me*me,a=d+2*Dot(v.p,v.in),b=-2*Dot(v.p,v.out);
    return 8*d*(a*b-4*a*c+2*d*c-8*c*c)*(a*a*b+a*a*c+a*b*b-a*b*d+2*a*b*c+b*b*c)/(a*a*b*b*((a-d)*(a-d)-4*d*c));}
void WardCheck(const Event&v){auto ai=Transverse(v.in),ao=Transverse(v.out);for(auto p:ai)Check(std::abs(SpinTrace(v,p,v.out))<1e-8,"outgoing Ward identity");for(auto p:ao)Check(std::abs(SpinTrace(v,v.in,p))<1e-8,"incoming Ward identity");Near(Dot(v.recoil,v.recoil),me*me,1e-10,"electron mass shell");}
void PhysicsChecks(const Config&c,const DecayTable&t){Rule q(64),q2(128);
    for(double E:{1.,3.,8.}){Near(CrossSection(E,0,false,false,q),KleinNishina(E),1e-10,"production KN");Near(CrossSection(E,0,true,false,q),KleinNishina(E),1e-10,"transverse absorption KN");}
    for(double m:{.01,.3,.9})for(double E:{2.,4.}){
        for(double z:{-.7,.1,.8}){auto v=Kinematics(E,0,m,z);auto inverse=Kinematics(E,m,0,z);WardCheck(v);WardCheck(inverse);Near(Squared(v,m,false,false)+Squared(v,m,false,true),InvariantF(v,m),1e-9,"physical polarizations / invariant trace");Near(Squared(inverse,m,true,true),AbsorptionLongitudinalF(inverse,m),1e-8,"closed longitudinal absorption trace");}
        double s=me*me+m*m+2*me*E,Eg=(s-me*me)/(2*me),p2=E*E-m*m;
        double unpol=(2*CrossSection(E,m,true,false,q)+CrossSection(E,m,true,true,q))/3;
        double prod=CrossSection(Eg,m,false,false,q)+CrossSection(Eg,m,false,true,q);
        Near(unpol,2./3*Eg*Eg/p2*prod,1e-9,"massive detailed balance");
        Near(CrossSection(E,m,true,false,q),CrossSection(E,m,true,false,q2),1e-10,"angle quadrature");
    }
    Near(ConversionCoefficient(1,1,20,10),1./std::pow(20.,4),.006,"Danilov low-mass production estimate");
    Near(WidthEH(.1,2e-5)/WidthEH(.1,1e-5),4,1e-12,"width mixing scaling");
    for(size_t j=0;j<t.m.size();j+=7)if(t.m[j]<1)Near(t.Enhancement(t.m[j]),t.r[j],1e-13,"published decay-table knot");
    double printedCap=195.7;Near(30.7+1.64*100.6,printedCap,.0001,"Danilov rounded event cap");
    std::cout<<"PASS: fresh polarized amplitudes, Ward identities, symbolic trace comparison, both KN limits, massive detailed balance, doubled angular quadrature\nPASS: Danilov low-mass coefficient and rounded event cap; published exact-loop enhancement knots; decay mixing scaling\n";
}
void Save(TCanvas&canvas,const std::string&name){canvas.Modified();canvas.Update();for(auto ext:{".pdf",".svg",".png"})canvas.SaveAs(("output/"+name+ext).c_str());}
void GraphCheck(const TGraph&g){Check(g.GetN()>1,"empty plotted graph");for(int i=0;i<g.GetN();++i)Check(std::isfinite(g.GetPointX(i))&&std::isfinite(g.GetPointY(i))&&g.GetPointY(i)>0,"invalid plotted data");}
void DrawCrossSections(const Config&c){Rule q(c.na);std::array<TGraph,4> gs;std::ofstream csv("output/cross_sections.csv");csv<<"energy_MeV,production_T,production_L,absorption_T,absorption_L\n"<<std::setprecision(16);
    for(int i=0;i<=65;++i){double E=.7+9.3*i/65,m=.3;std::array<double,4> y={CrossSection(E,m,false,false,q),CrossSection(E,m,false,true,q),CrossSection(E,m,true,false,q),CrossSection(E,m,true,true,q)};csv<<E;for(int j=0;j<4;++j){csv<<','<<y[j];gs[j].SetPoint(i,E,y[j]/1e-24);}csv<<'\n';}
    TCanvas cv("cross","Fresh polarized scattering",1000,680);cv.SetLogy();cv.SetLeftMargin(.13);cv.DrawFrame(.7,1e-5,10,1)->SetTitle("m_{A'}=0.3 MeV, epsilon^{2} stripped;Incident energy [MeV];Cross section / #epsilon^{2} [barn/electron]");TLegend leg(.17,.13,.54,.34);leg.SetBorderSize(0);const char* names[]={"Production: sum T","Production: L","Absorption: average T","Absorption: L"};for(int j=0;j<4;++j){GraphCheck(gs[j]);gs[j].SetLineColor(j+1);gs[j].SetLineWidth(2);gs[j].SetLineStyle(j%2+1);gs[j].Draw("L SAME");leg.AddEntry(&gs[j],names[j],"l");}leg.Draw();Save(cv,"polarized_cross_sections");
}
void DrawSource(const Config&c){std::array<TGraph,3> gs;std::array<double,3> mass={1,10,1000};std::ofstream csv("output/source.csv");csv<<"energy_MeV,mass_eV,rate_per_epsilon2\n"<<std::setprecision(16);TCanvas cv("source","Oscillation production",1000,680);cv.SetLogy();cv.SetLeftMargin(.15);cv.SetRightMargin(.09);cv.SetBottomMargin(.14);cv.DrawFrame(.2,1e10,10,1e22)->SetTitle("Danilov transverse oscillation source;E_{A'}=E_{#gamma} [MeV];dN_{A'}/dE / #epsilon^{2} [s^{-1} MeV^{-1}]");TLegend leg(.57,.68,.87,.87);leg.SetBorderSize(0);for(int j=0;j<3;++j){for(int i=0;i<=98;++i){double E=.2+.1*i,y=PhotonSource(E,c)*ConversionCoefficient(E,mass[j],c.mr,c.lr);gs[j].SetPoint(i,E,y);csv<<E<<','<<mass[j]<<','<<y<<'\n';}GraphCheck(gs[j]);gs[j].SetLineColor(j+1);gs[j].SetLineWidth(2);gs[j].Draw("L SAME");leg.AddEntry(&gs[j],Form("m_{A'} = %.0f eV",mass[j]),"l");}leg.Draw();Save(cv,"oscillation_source");}

// Plot separate line segments across every invalid mass, never connect a resonance gap.
struct Segments {std::vector<std::unique_ptr<TGraph>> graphs;bool gap=true;void Add(double x,double y){if(!std::isfinite(y)){gap=true;return;}if(gap){graphs.push_back(std::make_unique<TGraph>());gap=false;}auto&g=*graphs.back();g.SetPoint(g.GetN(),x,y);}void Draw(int color,int style){for(auto&g:graphs)if(g->GetN()>1){GraphCheck(*g);g->SetLineColor(color);g->SetLineStyle(style);g->SetLineWidth(2);g->Draw("L SAME");}}TGraph* First(){return graphs.empty()?nullptr:graphs.front().get();}};
void Limits(const Config&c,const DecayTable&t){Rule qe(c.ne),qa(c.na);Segments main,alt,unit,withDecays;std::ofstream csv("output/exclusion.csv");csv<<"mass_eV,epsilon95,unit_efficiency,plasma60_epsilon95,no_decay_epsilon95,including_selected_decays,status\n"<<std::setprecision(16);std::ofstream benchmark("output/benchmarks.txt");benchmark<<std::setprecision(15);double maxLoss=0,maxDecayChange=0;std::vector<double> masses;for(int i=0;i<=115;++i)masses.push_back(std::pow(10.,-2+8.*i/115));for(double m:{1.,10.,18.,19.,20.,21.,22.,50.,58.,60.,62.,100.,1000.,10000.,100000.,500000.,999000.})masses.push_back(m);std::sort(masses.begin(),masses.end());
    for(double ev:masses){if(ev>=1e6)continue;double e=NAN,e1=NAN,e60=NAN,e0=NAN,eBoth=NAN;
        if(AwayFromResonance(ev,c)){auto k=EventKernel(ev*1e-6,c,qe,qa);e=Boundary(k,c,t);e0=Boundary(k,c,t,false);Config cd=c;cd.etaDecay=c.eta;eBoth=Boundary(k,cd,t);if(std::isfinite(eBoth))Near(Events(k,eBoth,cd,t),cd.cap,1e-10,"decay boundary closure");Config u=c;u.eta=1;auto ku=k;for(auto&w:ku.w)w/=c.eta;e1=Boundary(ku,u,t);if(std::isfinite(e1))Near(Events(ku,e1,u,t),u.cap,1e-10,"unit efficiency closure");
            if(std::isfinite(e)){Near(Events(k,e,c,t),c.cap,1e-10,"event cap closure");maxLoss=std::max(maxLoss,DecayDepth(c.lo,ev*1e-6,e,100*c.R,t));maxDecayChange=std::max(maxDecayChange,std::abs(e/e0-1));}
            Config d=c;d.mr=d.md=60;if(AwayFromResonance(ev,d)){auto kd=k;for(size_t j=0;j<k.w.size();++j){kd.sourceWeight[j]*=ConversionCoefficient(k.E[j],ev,d.mr,d.lr)/ConversionCoefficient(k.E[j],ev,c.mr,c.lr);kd.w[j]*=ConversionCoefficient(k.E[j],ev,d.mr,d.lr)/ConversionCoefficient(k.E[j],ev,c.mr,c.lr)*DetectorCoefficient(ev,d.md)/DetectorCoefficient(ev,c.md);}e60=Boundary(kd,d,t);if(std::isfinite(e60))Near(Events(kd,e60,d,t),d.cap,1e-10,"alternate medium closure");}
        }else {Config d=c;d.mr=d.md=60;if(AwayFromResonance(ev,d)){auto kd=EventKernel(ev*1e-6,d,qe,qa);e60=Boundary(kd,d,t);if(std::isfinite(e60))Near(Events(kd,e60,d,t),d.cap,1e-10,"alternate medium closure");}}
        std::string status=std::isfinite(e)?"conditional_limit":AwayFromResonance(ev,c)?"outside_small_mixing":"material_resonance";
        withDecays.Add(ev,eBoth);main.Add(ev,e);alt.Add(ev,e60);unit.Add(ev,e1);csv<<ev<<','<<e<<','<<e1<<','<<e60<<','<<e0<<','<<eBoth<<','<<status<<'\n';
        if(ev==1||ev==10||ev==1000||ev==100000||ev==999000)benchmark<<"mass_eV="<<ev<<" epsilon95="<<e<<" including_decays="<<eBoth<<"\n";
    }
    // Independent energy-refinement check at plateau and upper mass endpoint.
    for(double M:{.001,.999}){auto k=EventKernel(M,c,qe,qa),fine=EventKernel(M,c,Rule(c.ne*2),Rule(c.na*2));Near(Boundary(k,c,t),Boundary(fine,c,t),1e-8,"energy/angle convergence");Config both=c;both.etaDecay=c.eta;Near(Boundary(k,both,t),Boundary(fine,both,t),1e-8,"decay energy convergence");Config dense=both;dense.density*=2;double base=Events(k,1e-5,c,t),extra=Events(k,1e-5,both,t)-base;if(M>.5)Near(Events(k,1e-5,dense,t)-base,extra/2,1e-11,"decay volume/density scaling");double e=Boundary(k,c,t);Near(Events(k,e/2,c,t)*16,Events(k,e,c,t),1e-8,"epsilon fourth-power regime");}
    auto plateau=EventKernel(.001,c,qe,qa);double result=Boundary(plateau,c,t),reference=.000021*std::pow((1174.1/197.2)/(1.5),.25);
    benchmark<<"selected_cap="<<c.cap<<"\npreselection_cap="<<c.cap/c.eta<<"\nfiducial_kg_days="<<c.kg*c.days<<"\nmax_decay_optical_depth="<<maxLoss<<"\nmax_relative_boundary_shift_from_decay="<<maxDecayChange<<"\npaper_algebraic_rescaling="<<reference<<"\nfresh_plateau="<<result<<"\n";
    std::cout<<std::setprecision(12)<<"selected_cap = "<<c.cap<<"\nfresh_plateau_epsilon95 = "<<result<<"\npaper_rescaling_epsilon95 = "<<reference<<"\nmax_decay_optical_depth = "<<maxLoss<<"\n";
    std::cout<<"PUBLISHED NORMALIZATION COMPARISON: "<<(std::abs(result/reference-1)<.1?"PASS":"FAIL (no fitted factor)")<<"\nPASS: event-cap closure, doubled energy/angle quadrature, decay-volume scaling, monotonic-domain guard, epsilon-fourth-power check\n";
    TCanvas cv("limit","Fresh conditional TEXONO bound",1100,740);cv.SetLogx();cv.SetLogy();cv.SetLeftMargin(.13);cv.DrawFrame(.01,1e-7,1e6,.02)->SetTitle("TEXONO: one-sided Gaussian 95% CL, conditional response; m_{A'} [eV];#epsilon_{95}");main.Draw(kBlue+1,1);alt.Draw(kRed+1,2);unit.Draw(kGreen+2,3);withDecays.Draw(kMagenta+1,4);TLegend leg(.24,.14,.88,.38);leg.SetBorderSize(0);leg.AddEntry(main.First(),"20 eV media; efficiency 1/6","l");leg.AddEntry(alt.First(),"60 eV media; efficiency 1/6","l");leg.AddEntry(unit.First(),"20 eV media; unit efficiency","l");leg.AddEntry(withDecays.First(),"20 eV media; equal decay/absorption acceptance","l");leg.Draw();Save(cv,"texono_exclusion");
}
void DecayPlot(const DecayTable&t){TGraph g;std::ofstream f("output/decay.csv");f<<"mass_MeV,enhancement\n"<<std::setprecision(16);for(int i=0;i<=120;++i){double m=.01+.989*i/120,v=t.Enhancement(m);g.SetPoint(i,m,v);f<<m<<','<<v<<'\n';}TCanvas cv("decay","Published loop enhancement",950,650);cv.SetLogy();cv.DrawFrame(.01,1,.999,200)->SetTitle("Published exact-loop input, freshly interpolated; m_{A'} [MeV];#Gamma_{3#gamma}/#Gamma_{EH}");g.SetLineWidth(2);g.Draw("L SAME");GraphCheck(g);Save(cv,"decay_enhancement");}
int Run(const std::string&path){auto start=std::chrono::steady_clock::now();gROOT->SetBatch(true);gStyle->SetOptStat(0);auto c=ReadConfig(path);DecayTable t("../../papers/danilov_fresh/decay_enhancement.txt");PhysicsChecks(c,t);DrawCrossSections(c);DrawSource(c);Limits(c,t);DecayPlot(t);std::cout<<"UNCHECKED: material profiles and resonance; full detector response; longitudinal nuclear emission; original likelihood; exact massive nuclear-source correction\nROOT "<<gROOT->GetVersion()<<"; C++ "<<__VERSION__<<"\ncompute_seconds = "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<'\n';return 0;}
}
int main(int argc,char**argv){try{return Fresh::Run(argc>1?argv[1]:"experiment.conf");}catch(const std::exception&e){std::cerr<<"FAIL: "<<e.what()<<'\n';return 1;}}
