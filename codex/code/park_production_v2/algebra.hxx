// Exact dyadic polynomial Dirac traces, written for this task (no imported formula).
// Variables: a=s-m^2, b=u-m^2, c=m^2, d=M^2. All coefficients are dyadic
// rationals small enough to be represented exactly by binary double arithmetic.
#pragma once
#include <array>
#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace Algebra {
using Powers=std::array<int,4>;
struct Poly {
    std::map<Powers,double> terms;
    Poly(double x=0){if(x)terms[{0,0,0,0}]=x;}
    static Poly Var(int i){Poly p;Powers e{};e[i]=1;p.terms[e]=1;return p;}
    static Poly Inverse(int i){Poly p;Powers e{};e[i]=-1;p.terms[e]=1;return p;}
    Poly Derivative(int i)const{Poly p;for(auto&[e,k]:terms)if(e[i]){auto f=e;--f[i];p.terms[f]+=k*e[i];}p.Clean();return p;}
    void Clean(){for(auto i=terms.begin();i!=terms.end();)if(i->second==0)i=terms.erase(i);else ++i;}
    double Value(double a,double b,double c,double d)const{
        double sum=0;std::array<double,4>x={a,b,c,d};
        for(auto& [e,k]:terms){double v=k;for(int i=0;i<4;++i)v*=std::pow(x[i],e[i]);sum+=v;}return sum;
    }
    std::string Tex()const{
        std::ostringstream out;bool first=true;const char* names[]={"a","b","c","d"};
        for(auto& [e,k]:terms){if(!first&&k>0)out<<"+";out<<k;for(int i=0;i<4;++i)if(e[i]){out<<names[i];if(e[i]>1)out<<"^{"<<e[i]<<"}";}first=false;}
        return first?"0":out.str();
    }
};
inline Poly operator+(Poly a,const Poly& b){for(auto&[e,k]:b.terms)a.terms[e]+=k;a.Clean();return a;}
inline Poly operator*(const Poly&a,const Poly&b){Poly p;for(auto&[e,k]:a.terms)for(auto&[f,l]:b.terms){Powers g;for(int i=0;i<4;++i)g[i]=e[i]+f[i];p.terms[g]+=k*l;}p.Clean();return p;}
inline Poly operator-(Poly a,const Poly&b){return a+Poly(-1)*b;}

inline Poly Dot(int i,int j){
    // Vector ids 0=p, 1=p'=p+k-q, 2=r_s=p+k, 3=r_u=p-q.
    const int v[4][3]={{1,0,0},{1,1,-1},{1,1,0},{1,0,-1}};
    Poly a=Poly::Var(0),b=Poly::Var(1),c=Poly::Var(2),d=Poly::Var(3);
    Poly gram[3][3]={{c,Poly(.5)*a,Poly(.5)*(d-b)},
        {Poly(.5)*a,Poly(0),Poly(.5)*(a+b)},
        {Poly(.5)*(d-b),Poly(.5)*(a+b),d}};
    Poly result;
    for(int x=0;x<3;++x)for(int y=0;y<3;++y)result=result+Poly(v[i][x]*v[j][y])*gram[x][y];
    return result;
}

inline Poly PairingProduct(const std::vector<int>& labels,const std::vector<std::pair<int,int>>& pairs){
    // Contract repeated Lorentz dummy indices. Each closed metric loop is D=4;
    // each open chain terminates on two external vectors and gives their dot.
    std::vector<int> parent(labels.size());std::iota(parent.begin(),parent.end(),0);
    std::function<int(int)> find=[&](int i){return parent[i]==i?i:parent[i]=find(parent[i]);};
    auto join=[&](int i,int j){parent[find(i)]=find(j);};
    for(size_t i=0;i<labels.size();++i)for(size_t j=0;j<i;++j)
        if(labels[i]<0&&labels[i]==labels[j])join(i,j);
    for(auto [i,j]:pairs)join(i,j);
    std::map<int,std::vector<int>> groups;
    for(size_t i=0;i<labels.size();++i){auto& g=groups[find(i)];if(labels[i]>=0)g.push_back(labels[i]);}
    Poly value(1);
    for(auto&[key,ext]:groups){
        if(ext.empty())value=value*Poly(4);
        else if(ext.size()==2)value=value*Dot(ext[0],ext[1]);
        else throw std::runtime_error("Invalid Dirac contraction graph");
    }return value;
}

inline Poly Trace(const std::vector<int>& labels){
    // Tr(v1/...v2n/) = 4 sum_pairings (-1)^crossings product(dot products).
    if(labels.size()%2)return Poly(0);
    std::vector<int> remaining(labels.size());std::iota(remaining.begin(),remaining.end(),0);
    std::vector<std::pair<int,int>> pairs;Poly result;
    std::function<void(std::vector<int>,int)> visit=[&](std::vector<int> rest,int sign){
        if(rest.empty()){result=result+Poly(4*sign)*PairingProduct(labels,pairs);return;}
        for(size_t j=1;j<rest.size();++j){
            pairs.emplace_back(rest[0],rest[j]);std::vector<int> next;
            for(size_t k=1;k<rest.size();++k)if(k!=j)next.push_back(rest[k]);
            visit(next,sign*(j%2?1:-1));pairs.pop_back();
        }
    };visit(remaining,1);return result;
}

inline Poly DiagramTrace(const std::array<int,8>& pattern){
    // Expand the four (slash+mass) factors. Negative labels are gamma indices;
    // -1=mu and -2=nu. Odd numbers of scalar masses give odd traces and vanish.
    Poly result,c=Poly::Var(2);
    for(int mask=0;mask<16;++mask){
        int slot=0,masses=0;std::vector<int> labels;
        for(int id:pattern)if(id<0)labels.push_back(id);else{if(mask&(1<<slot))labels.push_back(id);else ++masses;++slot;}
        if(masses%2)continue;Poly factor(1);for(int i=0;i<masses/2;++i)factor=factor*c;
        result=result+factor*Trace(labels);
    }return result;
}

struct Result {Poly ss,uu,su,us,numerator;};
inline Result Derive(){
    Result r;
    r.ss=DiagramTrace({1,-2,2,-1,0,-1,2,-2});
    r.uu=DiagramTrace({1,-1,3,-2,0,-2,3,-1});
    r.su=DiagramTrace({1,-2,2,-1,0,-2,3,-1});
    r.us=DiagramTrace({1,-1,3,-2,0,-1,2,-2});
    Poly a=Poly::Var(0),b=Poly::Var(1),c=Poly::Var(2),d=Poly::Var(3);
    r.numerator=Poly(.25)*(r.ss*b*b+r.uu*a*a+(r.su+r.us)*a*b);
    Poly compact=Poly(-2)*(a*a*a*b+a*b*b*b)+Poly(4)*(Poly(2)*c+d)*((a+b)*a*b+c*(a+b)*(a+b)-d*a*b);
    if(!(compact-r.numerator).terms.empty())throw std::runtime_error("Compact formula failed exact trace identity");
    Poly ia=Poly::Inverse(0),ib=Poly::Inverse(1),H=Poly(2)*c+d;
    Poly primitive=Poly(-1)*b*b*ia+Poly(4)*H*(ia+c*ia*ia)*b-Poly(4)*H*c*ib;
    Poly logCoefficient=Poly(-2)*a+Poly(4)*H*(Poly(1)+(Poly(2)*c-d)*ia);
    Poly differentiated=(primitive.Derivative(1)+logCoefficient*ib)*a*a*b*b;
    if(!(differentiated-r.numerator).terms.empty())throw std::runtime_error("Analytic primitive derivative failed exact identity");
    if(!(r.su-r.us).terms.empty())throw std::runtime_error("Interference traces disagree");
    Poly zeroMass;
    for(auto&[e,k]:r.numerator.terms)if(e[3]==0)zeroMass.terms[e]=k;
    Poly kleinNishina=Poly(-2)*(a*a*a*b+a*b*b*b)+Poly(8)*c*(a+b)*a*b+Poly(8)*c*c*(a+b)*(a+b);
    if(!(zeroMass-kleinNishina).terms.empty())throw std::runtime_error("Symbolic photon limit failed");
    Poly exchanged;
    for(auto&[e,k]:r.numerator.terms){auto f=e;std::swap(f[0],f[1]);exchanged.terms[f]+=k;}
    if(!(exchanged-r.numerator).terms.empty())throw std::runtime_error("Symbolic crossing symmetry failed");
    return r;
}
}
