// New production implementation: exact trace polynomials plus independent matrices.
#pragma once
#include "algebra.hxx"
#include <TF1.h>
#include <complex>
#include <iostream>

namespace Production {
constexpr double electron=0.510998950,alpha=1.0/137.035999084,hbarc=1.973269804e-11;
constexpr double pi=3.14159265358979323846,e4=16*pi*pi*alpha*alpha;
using Vector=std::array<double,4>;
using Complex=std::complex<double>;
struct Matrix {std::array<Complex,16>x{};Complex& at(int i,int j){return x[4*i+j];}Complex at(int i,int j)const{return x[4*i+j];}};
inline Matrix operator+(Matrix a,const Matrix&b){for(int i=0;i<16;++i)a.x[i]+=b.x[i];return a;}
inline Matrix operator*(Complex k,Matrix a){for(auto&v:a.x)v*=k;return a;}
inline Matrix operator*(const Matrix&a,const Matrix&b){Matrix r;for(int i=0;i<4;++i)for(int j=0;j<4;++j)for(int k=0;k<4;++k)r.at(i,j)+=a.at(i,k)*b.at(k,j);return r;}
inline Matrix Identity(){Matrix a;for(int i=0;i<4;++i)a.at(i,i)=1;return a;}
inline const std::array<Matrix,4>& Gamma(){
    static const auto gamma=[](){
        std::array<Matrix,4> g;g[0].at(0,0)=g[0].at(1,1)=1;g[0].at(2,2)=g[0].at(3,3)=-1;
        Complex pauli[3][2][2]={{{0,1},{1,0}},{{0,Complex(0,-1)},{Complex(0,1),0}},{{1,0},{0,-1}}};
        for(int k=0;k<3;++k)for(int i=0;i<2;++i)for(int j=0;j<2;++j){g[k+1].at(i,j+2)=pauli[k][i][j];g[k+1].at(i+2,j)=-pauli[k][i][j];}
        return g;
    }();return gamma;
}
inline Matrix Slash(const Vector&p){Matrix a;for(int mu=0;mu<4;++mu)a=a+Complex((mu? -1:1)*p[mu])*Gamma()[mu];return a;}
inline Matrix Bar(const Matrix&a){Matrix dagger;for(int i=0;i<4;++i)for(int j=0;j<4;++j)dagger.at(i,j)=std::conj(a.at(j,i));return Gamma()[0]*dagger*Gamma()[0];}
inline Complex Trace(const Matrix&a){Complex r=0;for(int i=0;i<4;++i)r+=a.at(i,i);return r;}
inline double Norm(const Matrix&a){double r=0;for(auto z:a.x)r=std::max(r,std::abs(z));return r;}
inline Vector Sum(const Vector&a,const Vector&b,double sign=1){Vector r;for(int i=0;i<4;++i)r[i]=a[i]+sign*b[i];return r;}
inline double Dot(const Vector&a,const Vector&b){double r=a[0]*b[0];for(int i=1;i<4;++i)r-=a[i]*b[i];return r;}
inline void Require(bool ok,const std::string&message){if(!ok)throw std::runtime_error(message);}
inline void Close(double x,double y,double tolerance,const std::string&message){Require(std::isfinite(x)&&std::isfinite(y)&&std::abs(x-y)<=tolerance*std::max(std::abs(y),1e-100),message);}

struct Quadrature {
    std::vector<double> nodes,weights;
    explicit Quadrature(int n):nodes(n),weights(n){TF1::CalcGaussLegendreSamplingPoints(n,nodes.data(),weights.data(),1e-14);}
    template<class F>double Integrate(F f,double lo,double hi)const{
        if(hi<=lo)return 0;double r=0,half=(hi-lo)/2,mid=(hi+lo)/2;
        for(size_t i=0;i<nodes.size();++i)r+=weights[i]*f(mid+half*nodes[i]);return half*r;
    }
};

inline double Threshold(double mass){return mass+mass*mass/(2*electron);}
inline double Kallen(double x,double y,double z){return (x-y-z)*(x-y-z)-4*y*z;}
inline std::pair<double,double> Endpoints(double energy,double mass){
    if(energy<=Threshold(mass))return {0,0};
    double s=electron*electron+2*electron*energy;
    double center=(energy+electron)*(s+mass*mass-electron*electron)/(2*s);
    double half=energy*std::sqrt(Kallen(s,electron*electron,mass*mass))/(2*s);
    return {center-half,center+half};
}
inline const Algebra::Result& Symbolic(){static const auto result=Algebra::Derive();return result;}
inline double TracePolynomial(double energy,double outgoing,double mass){
    double a=2*electron*energy,b=mass*mass-2*electron*outgoing;
    return Symbolic().numerator.Value(a,b,electron*electron,mass*mass)/(a*a*b*b);
}
inline double Differential(double energy,double outgoing,double mass){
    // d sigma / d outgoing energy, cm^2 / MeV, per electron, epsilon^2 stripped.
    auto bounds=Endpoints(energy,mass);
    if(bounds.second<=bounds.first||outgoing<bounds.first||outgoing>bounds.second)return 0;
    return e4*hbarc*hbarc*TracePolynomial(energy,outgoing,mass)/(32*pi*electron*energy*energy);
}
inline double Total(double energy,double mass,const Quadrature&q){auto [lo,hi]=Endpoints(energy,mass);return q.Integrate([&](double out){return Differential(energy,out,mass);},lo,hi);}
inline double Primitive(double a,double b,double mass){
    // Antiderivative of F with respect to b; the logarithm has a dimensionless argument.
    double c=electron*electron,d=mass*mass,H=2*c+d;
    return -b*b/a+4*H*(1/a+c/(a*a))*b+(-2*a+4*H*(1+(2*c-d)/a))*std::log(std::abs(b)/a)-4*H*c/b;
}
inline double AnalyticTotal(double energy,double mass){
    auto [lo,hi]=Endpoints(energy,mass);if(hi<=lo)return 0;
    double a=2*electron*energy,blo=mass*mass-2*electron*lo,bhi=mass*mass-2*electron*hi;
    return e4*hbarc*hbarc/(64*pi*electron*electron*energy*energy)*(Primitive(a,blo,mass)-Primitive(a,bhi,mass));
}
inline double Thomson(){return 8*pi*alpha*alpha*hbarc*hbarc/(3*electron*electron);}
inline double KleinNishina(double energy){
    double x=energy/electron;
    if(x<1e-4)return Thomson()*(1-2*x+5.2*x*x);
    double l=std::log1p(2*x);
    return 2*pi*alpha*alpha*hbarc*hbarc/(electron*electron)*
      ((1+x)/(x*x)*(2*(1+x)/(1+2*x)-l/x)+l/(2*x)-(1+3*x)/((1+2*x)*(1+2*x)));
}
inline double Reference15(double energy,double mass){
    // Borrowed only for validation: Gondolo & Raffelt Appendix, vector branch.
    if(energy<=Threshold(mass))return 0;
    double m2=electron*electron,M2=mass*mass,s=m2+2*electron*energy,root=std::sqrt(s),a=s-m2;
    double p0=(s-m2+M2)/(2*root),p=std::sqrt(Kallen(s,m2,M2))/(2*root),k0=(s+m2)/(2*root),k=a/(2*root);
    double A=2+2*(m2-M2)/s+16*(M2+2*m2)*s/(a*a);
    double B=2-4*(M2+2*m2)/a-4*(4*m2*m2-M2*M2)/(a*a);
    return pi*alpha*alpha/(2*s)*(p*A/k+B*root/k*std::log((2*p0*k0+2*p*k-M2)/(2*p0*k0-2*p*k-M2)))*hbarc*hbarc;
}

struct Kinematics {Vector p,k,q,pprime;double outgoing;};
inline Kinematics CenterOfMass(double energy,double mass,double cosine){
    // Construct an independent on-shell event from a CM angle, then boost to lab.
    double s=electron*electron+2*electron*energy,root=std::sqrt(s);
    double momentum=std::sqrt(Kallen(s,electron*electron,mass*mass))/(2*root);
    double q0=(s+mass*mass-electron*electron)/(2*root),beta=energy/(energy+electron),boost=(energy+electron)/root;
    Kinematics v;v.p={electron,0,0,0};v.k={energy,0,0,energy};
    v.q={boost*(q0+beta*momentum*cosine),momentum*std::sqrt(1-cosine*cosine),0,boost*(momentum*cosine+beta*q0)};
    v.pprime=Sum(Sum(v.p,v.k),v.q,-1);v.outgoing=v.q[0];return v;
}
inline Matrix Vertex(const Kinematics&v,int nu,int mu){
    auto rs=Sum(v.p,v.k),ru=Sum(v.p,v.q,-1);
    return Complex(1/(Dot(rs,rs)-electron*electron))*(Gamma()[nu]*(Slash(rs)+Complex(electron)*Identity())*Gamma()[mu])+
           Complex(1/(Dot(ru,ru)-electron*electron))*(Gamma()[mu]*(Slash(ru)+Complex(electron)*Identity())*Gamma()[nu]);
}
inline double MatrixTrace(const Kinematics&v,bool physical=false){
    Matrix in=Slash(v.p)+Complex(electron)*Identity(),out=Slash(v.pprime)+Complex(electron)*Identity();
    Complex total=0;
    if(!physical){for(int nu=0;nu<4;++nu)for(int mu=0;mu<4;++mu){auto a=Vertex(v,nu,mu);total+=double((nu?1:-1)*(mu?1:-1))*Trace(out*a*in*Bar(a));}}
    else {
        // Two incident polarizations and three final massive-vector polarizations.
        double P=std::sqrt(v.q[0]*v.q[0]-Dot(v.q,v.q)),mass=std::sqrt(Dot(v.q,v.q));
        Vector ex={0,1,0,0},ey={0,0,1,0};
        Vector transverse={0,v.q[3]/P,0,-v.q[1]/P},longitudinal={P/mass,v.q[0]*v.q[1]/(mass*P),0,v.q[0]*v.q[3]/(mass*P)};
        for(auto incoming:{ex,ey})for(auto outgoing:{transverse,ey,longitudinal}){
            Matrix a;for(int nu=0;nu<4;++nu)for(int mu=0;mu<4;++mu)a=a+Complex((nu?-1:1)*(mu?-1:1)*outgoing[nu]*incoming[mu])*Vertex(v,nu,mu);
            total+=Trace(out*a*in*Bar(a));
        }
    }
    Require(std::abs(total.imag())<1e-8,"Nonreal matrix trace");return total.real()/4;
}
inline void CheckMatrices(){
    for(int mu=0;mu<4;++mu)for(int nu=0;nu<4;++nu){
        Matrix difference=Gamma()[mu]*Gamma()[nu]+Gamma()[nu]*Gamma()[mu]+Complex(mu==nu?-2*(mu?-1:1):0)*Identity();
        Require(Norm(difference)==0,"Clifford algebra failed");
    }
    for(double mass:{.1,.5,1.})for(double energy:{2.2,4.,8.})for(double z:{-.8,.1,.9}){
        auto v=CenterOfMass(energy,mass,z);
        Close(Dot(v.pprime,v.pprime),electron*electron,2e-10,"Electron shell failed");
        Close(MatrixTrace(v),TracePolynomial(energy,v.outgoing,mass),2e-10,"Matrix/polynomial mismatch");
        Close(MatrixTrace(v,true),MatrixTrace(v),2e-9,"Physical polarization mismatch");
        Matrix in=Slash(v.p)+Complex(electron)*Identity(),out=Slash(v.pprime)+Complex(electron)*Identity();
        for(int index=0;index<4;++index){Matrix photon,vector;
            for(int j=0;j<4;++j){photon=photon+Complex((j?-1:1)*v.k[j])*Vertex(v,index,j);vector=vector+Complex((j?-1:1)*v.q[j])*Vertex(v,j,index);}
            Require(Norm(out*photon*in)<1e-8&&Norm(out*vector*in)<1e-8,"Ward identity failed");
        }
    }
}

struct Source {double powerMW=1000,gammaMin=1,gammaMax=40,normalization=.58e18,scale=.91;int nodes=96;};
inline double PhotonSpectrum(double energy,const Source&s){return s.normalization*s.powerMW*std::exp(-energy/s.scale);}
inline double Spectrum(double outgoing,double mass,const Source&s,const Quadrature&q){
    // Exact lab-angle support; no fitted normalization or imported rate routine.
    if(outgoing<=mass)return 0;
    double momentum=std::sqrt(outgoing*outgoing-mass*mass),numerator=electron*outgoing-mass*mass/2;
    double positive=electron-outgoing+momentum,negative=electron-outgoing-momentum;
    if(positive<=0)return 0;
    double lo=std::max({s.gammaMin,Threshold(mass),numerator/positive});
    double hi=negative>0?std::min(s.gammaMax,numerator/negative):s.gammaMax;
    return q.Integrate([&](double energy){return PhotonSpectrum(energy,s)*Differential(energy,outgoing,mass)/KleinNishina(energy);},lo,hi);
}
inline double ForwardBin(double low,double high,double mass,const Source&s,const Quadrature&q){
    // Independent integration order: start from an incident photon, then clip
    // the CM-derived outgoing endpoints to the requested bin. Split at the
    // support transitions so quadrature does not straddle derivative kinks.
    std::vector<double> edges={s.gammaMin,s.gammaMax,Threshold(mass)};
    for(double x:{low,high})if(x>mass){double Q=std::sqrt(x*x-mass*mass),n=electron*x-mass*mass/2;
        for(double denominator:{electron-x+Q,electron-x-Q})if(denominator>0)edges.push_back(n/denominator);}
    std::sort(edges.begin(),edges.end());double integral=0;
    for(size_t i=1;i<edges.size();++i)integral+=q.Integrate([&](double E){
        auto [a,b]=Endpoints(E,mass);
        return PhotonSpectrum(E,s)/KleinNishina(E)*q.Integrate([&](double x){return Differential(E,x,mass);},std::max(a,low),std::min(b,high));
    },std::max(edges[i-1],s.gammaMin),std::min(edges[i],s.gammaMax));
    return integral;
}
}
