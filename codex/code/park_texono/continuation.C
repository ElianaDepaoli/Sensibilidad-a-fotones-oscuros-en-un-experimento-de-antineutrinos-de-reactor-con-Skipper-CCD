// Targeted diagnostics reusing the existing implementation, without rebuilding it.
// Run from this directory: root -l -b -q continuation.C
#include "park_texono.cpp"
#include <TSystem.h>

namespace Continuation {
double ConditionalFlux(double EA,double M,const Park::Config& c,const Park::Gauss& q) {
    // Hypothesis only: normalize each production spectrum to its own total.
    // This drops the massive production branching suppression of Park Eq. (1).
    using namespace Park;
    if(EA<=M)return 0;
    double p=std::sqrt(EA*EA-M*M),n=me*EA-M*M/2,d=me-EA+p;
    if(d<=0)return 0;
    double lo=std::max({c.gamma_low,n/d,M+M*M/(2*me)}),hi=c.gamma_high;
    if(me-EA-p>0)hi=std::min(hi,n/(me-EA-p));
    return q.Integral([&](double E){return GammaSpectrum(E,c.power_MW)*ProductionDS(E,EA,M)/Reference15Total(E,M);},lo,hi);
}

double Median(std::vector<double> v) {
    std::sort(v.begin(),v.end());size_t n=v.size();
    if(!n)throw std::runtime_error("Empty comparison interval");
    return n%2?v[n/2]:(v[n/2-1]+v[n/2])/2;
}

void CompareBins(const Park::Config& c,const Park::Gauss& q) {
    // Horizontal PDF segments are histogram bins; vertical edges are not data.
    using namespace Park;
    std::ifstream in("output/park_figure1_digitized.csv");
    if(!in)throw std::runtime_error("Missing stored digitization");
    std::ofstream out("output/continuation/bin_comparison.csv");
    out<<"low_MeV,high_MeV,mass_MeV,paper_flux,baseline,low_cutoff,conditional,conditional_low_cutoff\n"<<std::setprecision(16);
    std::string line;std::getline(in,line);double previousE=0,previousM=-1,previousF=0;
    Config low=c;low.gamma_low=.2;Gauss binq(8);
    std::map<double,std::vector<std::vector<double>>> residual;
    while(std::getline(in,line)) {
        std::replace(line.begin(),line.end(),',',' ');std::istringstream row(line);double E,M,f;row>>E>>M>>f;
        if(M==previousM&&E>previousE&&std::abs(f/previousF-1)<1e-10) {
            double a=previousE,b=E;
            double values[]={binq.Integral([&](double x){return Flux(x,M,c,q);},a,b)/(b-a),
                binq.Integral([&](double x){return Flux(x,M,low,q);},a,b)/(b-a),
                binq.Integral([&](double x){return ConditionalFlux(x,M,c,q);},a,b)/(b-a),
                binq.Integral([&](double x){return ConditionalFlux(x,M,low,q);},a,b)/(b-a)};
            out<<a<<','<<b<<','<<M<<','<<f;for(double v:values)out<<','<<v;out<<'\n';
            if(a>=1.8&&b<=3.2) {
                if(!residual.count(M))residual[M].resize(4);
                for(int j=0;j<4;++j)residual[M][j].push_back(std::abs(values[j]/f-1));
            }
        }
        previousE=E;previousM=M;previousF=f;
    }
    for(auto& [M,values]:residual) {
        std::cout<<"Bin residual M="<<M<<" baseline="<<Median(values[0])
          <<" lower-cutoff="<<Median(values[1])<<" conditional="<<Median(values[2])
          <<" conditional+lower-cutoff="<<Median(values[3])<<'\n';
    }
}

void CheckFiniteMass(const Park::Config& c,const Park::Gauss& q) {
    using namespace Park;
    Gauss fine(2*c.quadrature);
    std::ofstream out("output/continuation/finite_mass_checks.csv");
    out<<"mass_MeV,events_coefficient,node_relative_difference,order_relative_difference\n"<<std::setprecision(16);
    for(double M:{.001,.1,.5,1.}) {
        double a=Events(M,c,q),b=Events(M,c,fine),d=EventsReordered(M,c,fine);
        Near(a,b,2e-5,"finite mass nodes");Near(a,d,2e-4,"finite mass integration order");
        out<<M<<','<<a<<','<<std::abs(a/b-1)<<','<<std::abs(a/d-1)<<'\n';
        std::cout<<"PASS finite-mass rate M="<<M<<" nodes="<<std::abs(a/b-1)<<" order="<<std::abs(a/d-1)<<'\n';
    }
    Config one=c;one.power_MW=1000;
    for(double E:{.3,1.5,3.})
        Near(ConditionalFlux(E,0,one,q),Flux(E,0,one,q),1e-8,"conditional massless limit");
    double conditionalRate=q.Integral([&](double E){return ConditionalFlux(E,.001,c,q)*DetectionTotal(E,.001,q);},c.roi_low,c.roi_high)
        *ElectronCount(c)*c.days*86400*c.efficiency/(4*pi*std::pow(c.baseline_m*100,2));
    std::ofstream limit("output/continuation/conditional_limit.csv");
    limit<<std::setprecision(16)<<"events_coefficient,epsilon_limit\n"<<conditionalRate<<','<<EpsilonLimit(conditionalRate,c)<<'\n';
    Near(EpsilonLimit(conditionalRate,c),EpsilonLimit(Events(.001,c,q),c),1e-5,"small mass conditional rate");
    std::cout<<"PASS: conditional massless limit; small-mass normalization diagnostic leaves TEXONO mismatch\n";
}

void PlotComparison() {
    // Plot stored bin averages, with one panel per mass and no fitted scale.
    std::ifstream in("output/continuation/bin_comparison.csv");std::string line;std::getline(in,line);
    std::map<double,std::vector<TGraph>> graphs;
    while(std::getline(in,line)) {
        std::replace(line.begin(),line.end(),',',' ');std::istringstream row(line);
        double a,b,m,p,base,low,conditional,both;row>>a>>b>>m>>p>>base>>low>>conditional>>both;
        if(!graphs.count(m))graphs[m].resize(3);
        double values[]={p,base,both};
        for(int j=0;j<3;++j)if(values[j]>0)graphs[m][j].SetPoint(graphs[m][j].GetN(),(a+b)/2,values[j]/1e21);
    }
    TCanvas cv("continuation_comparison","Figure 1 normalization diagnostic",1380,500);cv.Divide(3,1);
    std::vector<TLegend*> legends;int panel=0;
    for(auto& [m,g]:graphs) {
        cv.cd(++panel)->SetLogy();gPad->SetLeftMargin(.15);
        auto frame=gPad->DrawFrame(0,.0015,4.3,3.);
        frame->SetTitle(Form("m_{A'} = %.1f MeV;E_{A'} [MeV];dN_{A'}/dE_{A'} [10^{21} MeV^{-1} s^{-1}]",m));
        auto leg=new TLegend(.39,.67,.89,.88);legends.push_back(leg);leg->SetBorderSize(0);leg->SetTextSize(.027);
        const char* names[]={"Park PDF bins","Stated-equation baseline","Conditional + 0.2 MeV cutoff"};
        int colors[]={kBlack,kBlue,kRed};
        for(int j=0;j<3;++j){g[j].SetLineColor(colors[j]);g[j].SetLineWidth(2);g[j].SetLineStyle(j==0?3:1);g[j].Draw("L SAME");leg->AddEntry(&g[j],names[j],"l");}
        leg->Draw();
    }
    cv.SaveAs("output/continuation/figure1_diagnostic.svg");cv.SaveAs("output/continuation/figure1_diagnostic.pdf");
    for(auto leg:legends)delete leg;
}
}

void continuation() {
    try {
        auto c=Park::ReadConfig("texono.conf");Park::Gauss q(c.quadrature);
        std::cout<<std::setprecision(10)<<"Continuation: ROOT "<<gROOT->GetVersion()<<"; existing source included unchanged\n";
        Park::Tests(c,q);Continuation::CheckFiniteMass(c,q);
        c.power_MW=1000;Continuation::CompareBins(c,q);Continuation::PlotComparison();
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';gSystem->Exit(1);}
}
