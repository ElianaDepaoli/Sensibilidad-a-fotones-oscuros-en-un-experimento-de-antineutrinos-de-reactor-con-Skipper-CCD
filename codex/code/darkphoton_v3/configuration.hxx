#pragma once
#include "physics_base.hxx"
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TH1.h>
#include <TROOT.h>
#include <TStyle.h>
#include <fstream>
#include <iomanip>
#include <filesystem>
namespace V3 {
using namespace Production;
constexpr double hbarc_eV_cm=1.973269804e-5, NA=6.02214076e23, pi=3.14159265358979323846;
struct Config {double powerMW=2900,baseline_m=28,masskg=187,days=160,Z=108,A=259.8099219,roiLo=3,roiHi=8,gammaLo=.2,gammaHi=40,reactorM=20,detectorM=20,absLength=10,N95=195.7,eta=1./6,excess=30.7,error=100.6;bool gaussian=false;int nodes=96;};
Config ReadConfig(const std::string&path){Config c;std::ifstream f(path);if(!f)throw std::runtime_error("missing config");std::string line;while(std::getline(f,line)){if(line.empty()||line[0]=='#')continue;auto at=line.find('=');if(at==std::string::npos)throw std::runtime_error("bad config");std::string k=line.substr(0,at);double v=std::stod(line.substr(at+1));
    if(k=="power_MW")c.powerMW=v;else if(k=="baseline_m")c.baseline_m=v;else if(k=="detector_mass_kg")c.masskg=v;else if(k=="live_days")c.days=v;else if(k=="Z")c.Z=v;else if(k=="A_g_mol")c.A=v;else if(k=="roi_low_MeV")c.roiLo=v;else if(k=="roi_high_MeV")c.roiHi=v;else if(k=="gamma_low_MeV")c.gammaLo=v;else if(k=="gamma_high_MeV")c.gammaHi=v;else if(k=="quadrature_nodes")c.nodes=int(v);else if(k=="reactor_mgamma_eV")c.reactorM=v;else if(k=="detector_mgamma_eV")c.detectorM=v;else if(k=="reactor_absorption_length_cm")c.absLength=v;else if(k=="N95_events")c.N95=v;else if(k=="signal_efficiency")c.eta=v;else if(k=="observed_minus_background")c.excess=v;else if(k=="excess_sigma")c.error=v;else if(k=="use_gaussian_cap")c.gaussian=v!=0;else throw std::runtime_error("unknown key "+k);}
    if(c.gaussian)c.N95=std::max(0.,c.excess+1.644853626951*c.error);
    if(c.Z<=0||c.A<=0||c.N95<=0||c.eta<=0||c.eta>1||c.error<=0||c.absLength<=0||c.reactorM<=0||c.detectorM<=0||c.roiLo<c.gammaLo||c.roiHi>c.gammaHi)throw std::runtime_error("unphysical response settings");
    if(c.powerMW<=0||c.baseline_m<=0||c.masskg<=0||c.days<=0||c.roiHi<=c.roiLo||c.gammaHi<=c.gammaLo||c.nodes<32)return throw std::runtime_error("unphysical config"),c;return c;}
inline double GammaWidthEV(double length){return hbarc_eV_cm/length;}
inline double PhotonSpectrum(double E,const Config&c){return .58e18*c.powerMW*std::exp(-E/.91);}
inline double ElectronCount(const Config&c){return c.masskg*1000*NA*c.Z/c.A;}
}
