// Fresh implementation: no earlier project headers or numerical results imported.
#pragma once
#include <TF1.h>
#include <Math/QuantFuncMathCore.h>
#include <array>
#include <complex>
#include <cmath>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>

namespace Fresh {
constexpr double pi=3.14159265358979323846, alpha=1/137.035999084;
constexpr double me=.510998950, hc=1.973269804e-11, hcEVcm=1.973269804e-5;
constexpr double NA=6.02214076e23, e4=16*pi*pi*alpha*alpha;
using Four=std::array<double,4>;
using Complex=std::complex<double>;
using Matrix=std::array<Complex,16>;
inline void Check(bool ok,const std::string& text){if(!ok)throw std::runtime_error(text);}
inline void Near(double a,double b,double r,const std::string& text){Check(std::isfinite(a)&&std::isfinite(b)&&std::abs(a-b)<=r*std::max(std::abs(b),1e-100),text);}
inline double Dot(const Four&a,const Four&b){return a[0]*b[0]-a[1]*b[1]-a[2]*b[2]-a[3]*b[3];}
inline Four Add(const Four&a,const Four&b,double sign=1){Four r;for(int i=0;i<4;++i)r[i]=a[i]+sign*b[i];return r;}
inline Matrix Plus(Matrix a,const Matrix&b){for(int i=0;i<16;++i)a[i]+=b[i];return a;}
inline Matrix Scale(Matrix a,Complex z){for(auto&v:a)v*=z;return a;}
inline Matrix Multiply(const Matrix&a,const Matrix&b){Matrix r{};for(int i=0;i<4;++i)for(int j=0;j<4;++j)for(int k=0;k<4;++k)r[4*i+j]+=a[4*i+k]*b[4*k+j];return r;}
inline Matrix Unit(){Matrix r{};for(int i=0;i<4;++i)r[5*i]=1;return r;}
inline const std::array<Matrix,4>& Gamma(){static auto g=[](){std::array<Matrix,4> v{};v[0][0]=v[0][5]=1;v[0][10]=v[0][15]=-1;std::array<std::array<Complex,4>,3> s={{{0,1,1,0},{0,Complex(0,-1),Complex(0,1),0},{1,0,0,-1}}};for(int k=0;k<3;++k)for(int i=0;i<2;++i)for(int j=0;j<2;++j){v[k+1][4*i+j+2]=s[k][2*i+j];v[k+1][4*(i+2)+j]=-s[k][2*i+j];}return v;}();return g;}
inline Matrix Slash(const Four&p){Matrix r{};for(int i=0;i<4;++i)r=Plus(r,Scale(Gamma()[i],(i?-1:1)*p[i]));return r;}
inline Matrix Adjoint(const Matrix&a){Matrix d{};for(int i=0;i<4;++i)for(int j=0;j<4;++j)d[4*i+j]=std::conj(a[4*j+i]);return Multiply(Multiply(Gamma()[0],d),Gamma()[0]);}
inline double TraceReal(const Matrix&a){Complex t=0;for(int i=0;i<4;++i)t+=a[5*i];Check(std::abs(t.imag())<1e-7,"complex spin trace");return t.real();}
inline Matrix SpinSum(const Four&p){return Plus(Slash(p),Scale(Unit(),me));}

struct Rule {std::vector<double>x,w;explicit Rule(int n):x(n),w(n){TF1::CalcGaussLegendreSamplingPoints(n,x.data(),w.data(),1e-14);}template<class F>double Integrate(F f,double a,double b)const{if(b<=a)return 0;double r=0;for(size_t i=0;i<x.size();++i)r+=w[i]*f((a+b)/2+(b-a)*x[i]/2);return r*(b-a)/2;}};
struct Event {Four p,in,out,recoil;double s,pinCM,poutCM;};
// Construct exact two-body lab momenta by boosting a CM scattering angle.
inline Event Kinematics(double energy,double incomingMass,double outgoingMass,double z){
    Event v;v.p={me,0,0,0};double p=std::sqrt(energy*energy-incomingMass*incomingMass);v.in={energy,0,0,p};v.s=me*me+incomingMass*incomingMass+2*me*energy;double root=std::sqrt(v.s);
    auto triangle=[](double x,double y,double z){return (x-y-z)*(x-y-z)-4*y*z;};
    v.pinCM=me*p/root;v.poutCM=std::sqrt(std::max(0.,triangle(v.s,me*me,outgoingMass*outgoingMass)))/(2*root);
    double ec=(v.s+outgoingMass*outgoingMass-me*me)/(2*root),beta=p/(energy+me),boost=(energy+me)/root;
    v.out={boost*(ec+beta*v.poutCM*z),v.poutCM*std::sqrt(1-z*z),0,boost*(v.poutCM*z+beta*ec)};v.recoil=Add(Add(v.p,v.in),v.out,-1);return v;
}
// Two transverse polarizations in the lab. The scattering plane is x-z.
inline std::array<Four,2> Transverse(const Four&q){double p=std::sqrt(q[1]*q[1]+q[3]*q[3]);return {Four{0,q[3]/p,0,-q[1]/p},Four{0,0,1,0}};}
// Subtract q/M from the longitudinal vector using the Ward identity, avoiding E/M cancellation.
inline Four Longitudinal(const Four&q,double mass){double p=std::sqrt(q[1]*q[1]+q[3]*q[3]),a=mass/(q[0]+p);return {-a,a*q[1]/p,0,a*q[3]/p};}
// s- and u-channel electron exchange, with the coupling epsilon*e^2 stripped.
inline Matrix Vertex(const Event&v,const Four&incomingPol,const Four&outgoingPol){
    Four a=Add(v.p,v.in),b=Add(v.p,v.out,-1);auto x=Slash(incomingPol),y=Slash(outgoingPol);
    auto s=Multiply(Multiply(y,SpinSum(a)),x),u=Multiply(Multiply(x,SpinSum(b)),y);
    return Plus(Scale(s,1/(Dot(a,a)-me*me)),Scale(u,1/(Dot(b,b)-me*me)));
}
inline double SpinTrace(const Event&v,const Four&a,const Four&b){auto vertex=Vertex(v,a,b);return TraceReal(Multiply(Multiply(Multiply(SpinSum(v.recoil),vertex),SpinSum(v.p)),Adjoint(vertex)));}
// Production sums outgoing T or L; absorption averages incoming T or L.
// Electron spins are averaged by 1/2. No artificial 1/3 average is applied to a T beam.
inline double Squared(const Event&v,double mass,bool absorption,bool longitudinal=false){
    auto ti=Transverse(v.in),to=Transverse(v.out);double sum=0;
    if(absorption){if(longitudinal){for(auto b:to)sum+=SpinTrace(v,Longitudinal(v.in,mass),b);return sum/2;}
        for(auto a:ti)for(auto b:to)sum+=SpinTrace(v,a,b);return sum/4;}
    if(longitudinal){for(auto a:ti)sum+=SpinTrace(v,a,Longitudinal(v.out,mass));return sum/4;}
    for(auto a:ti)for(auto b:to)sum+=SpinTrace(v,a,b);return sum/4;
}
inline double DifferentialCM(double energy,double mass,double z,bool absorption,bool longitudinal=false){
    if((!absorption&&energy<=mass+mass*mass/(2*me))||(absorption&&energy<=mass))return 0;
    auto v=Kinematics(energy,absorption?mass:0,absorption?0:mass,z);
    return e4*hc*hc*v.poutCM/v.pinCM*Squared(v,mass,absorption,longitudinal)/(32*pi*v.s);
}
inline double CrossSection(double E,double m,bool absorption,bool longitudinal,const Rule&q){return q.Integrate([&](double z){return DifferentialCM(E,m,z,absorption,longitudinal);},-1,1);}
inline double KleinNishina(double E){double x=E/me,l=std::log1p(2*x);return 2*pi*alpha*alpha*hc*hc/(me*me)*((1+x)/(x*x)*(2*(1+x)/(1+2*x)-l/x)+l/(2*x)-(1+3*x)/((1+2*x)*(1+2*x)));}

struct Config {double power=2900,R=28,kg=187,days=160,Z=108,A=259.8099219,lo=3,hi=8,excess=30.7,error=100.6,eta=1./6,etaDecay=0,density=4.51,cap=-1,mr=20,md=20,lr=10,ld=10,epsmax=.01;int ne=48,na=48;};
inline Config ReadConfig(const std::string&name){Config c;std::map<std::string,double*> fields={{"power_MW",&c.power},{"distance_m",&c.R},{"detector_mass_kg",&c.kg},{"live_days",&c.days},{"Z",&c.Z},{"A_g_mol",&c.A},{"roi_low_MeV",&c.lo},{"roi_high_MeV",&c.hi},{"excess_events",&c.excess},{"error_events",&c.error},{"efficiency",&c.eta},{"decay_efficiency",&c.etaDecay},{"detector_density_g_cm3",&c.density},{"cap_events",&c.cap},{"reactor_plasma_eV",&c.mr},{"detector_plasma_eV",&c.md},{"reactor_absorption_cm",&c.lr},{"detector_absorption_cm",&c.ld},{"epsilon_max",&c.epsmax}};std::ifstream f(name);Check(bool(f),"missing configuration");std::string line;while(std::getline(f,line)){if(line.empty()||line[0]=='#')continue;auto pos=line.find('=');Check(pos!=line.npos,"bad configuration line");auto key=line.substr(0,pos);double val=std::stod(line.substr(pos+1));Check(std::isfinite(val),"nonfinite config");if(key=="energy_nodes")c.ne=int(val);else if(key=="angle_nodes")c.na=int(val);else {Check(fields.count(key),"unknown config key: "+key);*fields.at(key)=val;}}
    Check(c.power>0&&c.R>0&&c.kg>0&&c.days>0&&c.A>0&&c.Z>0&&c.lo>1&&c.hi>c.lo&&c.error>0&&c.eta>0&&c.eta<=1&&c.etaDecay>=0&&c.etaDecay<=1&&c.density>0&&c.mr>0&&c.md>0&&c.lr>0&&c.ld>0&&c.epsmax>0&&c.epsmax<=.01&&c.ne>=24&&c.na>=24,"unsupported configuration");if(c.cap<0)c.cap=std::max(0.,c.excess+ROOT::Math::normal_quantile(.95,1)*c.error);Check(c.cap>0,"nonpositive cap needs a physical-boundary likelihood");return c;}
inline double PhotonSource(double E,const Config&c){return .58e18*c.power*std::exp(-E/.91);}
inline double Exposure(const Config&c){return c.kg*1000*NA*c.Z/c.A*c.days*86400/(4*pi*std::pow(100*c.R,2));}
// Probability / epsilon^2 for a transverse photon in a thick absorbing homogeneous source.
inline double ConversionCoefficient(double E,double massEV,double plasmaEV,double absCM){double delta=massEV*massEV-plasmaEV*plasmaEV;return std::pow(massEV,4)/(delta*delta+std::pow(E*1e6*hcEVcm/absCM,2));}
// Eq.(9): detector absorption neglected only away from material resonance.
inline double DetectorCoefficient(double massEV,double plasmaEV){return std::pow(massEV,4)/std::pow(massEV*massEV-plasmaEV*plasmaEV,2);}
inline bool AwayFromResonance(double mEV,const Config&c){return std::abs(mEV*mEV-c.mr*c.mr)>10*c.hi*1e6*hcEVcm/c.lr&&std::abs(mEV*mEV-c.md*c.md)>10*c.hi*1e6*hcEVcm/c.ld;}
inline bool WeakMixing(double mEV,double eps,const Config&c){return eps*mEV*mEV/std::min(std::abs(mEV*mEV-c.mr*c.mr),std::abs(mEV*mEV-c.md*c.md))<.1;}

struct DecayTable {std::vector<double>m,r;explicit DecayTable(const std::string&path){std::ifstream f(path);Check(bool(f),"missing published decay table");std::string s;while(std::getline(f,s)){if(s.empty()||s[0]=='#')continue;std::replace(s.begin(),s.end(),',',' ');std::istringstream row(s);double x,y;Check(bool(row>>x>>y)&&x>0&&y>=1,"bad decay data");m.push_back(x);r.push_back(y);}Check(m.size()>10&&std::is_sorted(m.begin(),m.end()),"invalid table order");}
    // Low-mass series from Table I, not an extrapolation of coarse tabulated points.
    double Enhancement(double M)const{if(M<m.front()){double x=M*M/(me*me);return 1+x*(335./714+x*(128941./839664+x*(44787./1026256+x*(1249649333./108064756800+x*(36494147./12382420050+x*867635449./1614300688000)))));}
        auto it=std::lower_bound(m.begin(),m.end(),M);Check(it!=m.end(),"decay-table extrapolation forbidden");size_t j=it-m.begin();if(j==0)return r[0];double w=(M-m[j-1])/(m[j]-m[j-1]);return std::exp((1-w)*std::log(r[j-1])+w*std::log(r[j]));}
};
inline double WidthEH(double M,double eps){return 17*std::pow(alpha,4)*eps*eps*std::pow(M,9)/(11664000*std::pow(pi,3)*std::pow(me,8));}
inline double DecayDepth(double E,double M,double eps,double lengthCM,const DecayTable&table){return lengthCM*M*WidthEH(M,eps)*table.Enhancement(M)/(hc*std::sqrt(E*E-M*M));}
struct Kernel {std::vector<double>E,w,sourceWeight;double mass;};
inline Kernel EventKernel(double M,const Config&c,const Rule&energy,const Rule&angle){Kernel k;k.mass=M;for(size_t j=0;j<energy.x.size();++j){double E=(c.lo+c.hi)/2+(c.hi-c.lo)*energy.x[j]/2;
    double sigma=CrossSection(E,M,true,false,angle);double weight=Exposure(c)*c.eta*(c.hi-c.lo)/2*energy.w[j]*PhotonSource(E,c)*ConversionCoefficient(E,M*1e6,c.mr,c.lr)*DetectorCoefficient(M*1e6,c.md)*sigma;
    k.E.push_back(E);k.w.push_back(weight);k.sourceWeight.push_back(Exposure(c)*(c.hi-c.lo)/2*energy.w[j]*PhotonSource(E,c)*ConversionCoefficient(E,M*1e6,c.mr,c.lr));}return k;}
inline double Events(const Kernel&k,double eps,const Config&c,const DecayTable&t,bool decay=true){double n=0;for(size_t j=0;j<k.E.size();++j)n+=(k.w[j]+c.etaDecay*k.sourceWeight[j]*DecayDepth(k.E[j],k.mass,1,1,t)/(c.density*NA*c.Z/c.A))*(decay?std::exp(-DecayDepth(k.E[j],k.mass,eps,100*c.R,t)):1);return std::pow(eps,4)*n;}
inline double Boundary(const Kernel&k,const Config&c,const DecayTable&t,bool decay=true){if(!AwayFromResonance(k.mass*1e6,c))return NAN;double lo=1e-14,hi=c.epsmax;Check(!decay||DecayDepth(c.lo,k.mass,hi,100*c.R,t)<2,"nonmonotonic decay regime: requires a two-boundary solver");if(Events(k,hi,c,t,decay)<c.cap)return NAN;for(int i=0;i<70;++i){double mid=std::sqrt(lo*hi);if(Events(k,mid,c,t,decay)>c.cap)hi=mid;else lo=mid;}double eps=std::sqrt(lo*hi);return WeakMixing(k.mass*1e6,eps,c)?eps:NAN;}
}
