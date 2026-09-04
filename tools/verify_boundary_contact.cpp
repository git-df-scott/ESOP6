// Standalone finite part of BOUNDARY_CONTACT_6.md. C++17, no dependencies.
// Independently uses the U,V norm and explicit midpoint coefficients;
// tools/boundary_contact.py uses A,B powers and recursive fifth roots.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
using Poly = std::vector<long long>;
using Point = std::array<int,9>;
using Vec = std::array<int,24>;
using Matrix = std::array<std::array<int,9>,24>;
const Point base = {2,0,1,0,1,0,2,0,1};
long long rem(long long n, int m) { n %= m; return n < 0 ? n+m : n; }
int inverse(int a, int m) {
    for (int b=1;b<m;++b) if (a*b%m==1) return b;
    assert(false); return 0;
}
Poly multiply(const Poly& a,const Poly& b,int mod) {
    Poly c(a.size()+b.size()-1);
    for(size_t i=0;i<a.size();++i) for(size_t j=0;j<b.size();++j)
        c[i+j]=(c[i+j]+a[i]*b[j])%mod;
    return c;
}
Poly power(Poly a,int e,int mod) {
    Poly r={1};
    while(e) { if(e&1) r=multiply(r,a,mod); e>>=1; if(e) a=multiply(a,a,mod); }
    return r;
}
Vec residual(const Point& q,int mod) {
    long long u2=q[0],u3=q[1],u4=q[2],u5=q[3];
    long long v1=q[4],v2=q[5],v3=q[6],v4=q[7],v5=q[8];
    // All visited coefficient representatives lie in [0,80]. The displayed
    // integer formulas remain below 2^63 before modular reduction.
    for(int x:q) assert(0<=x && x<=80);
    long long v12=v1*v1,v13=v12*v1,v14=v12*v12,v16=v13*v13;
    Poly u={1,0,u2,u3,u4,u5},v={0,v1,v2,v3,v4,v5};
    Poly d=power(u,6,mod), vv=power(v,6,mod);
    Poly uv1=multiply(power(u,4,mod),power(v,2,mod),mod);
    Poly uv2=multiply(power(u,2,mod),power(v,4,mod),mod);
    for(size_t i=0;i<d.size();++i) d[i]=(d[i]+15*uv1[i]+15*uv2[i]+vv[i])%mod;
    Poly m(7); m[0]=1;
    m[2]=rem(3*(2*u2+5*v12),mod)*inverse(5,mod)%mod;
    m[3]=rem(6*(u3+5*v1*v2),mod)*inverse(5,mod)%mod;
    m[4]=rem(3*(u2*u2-20*u2*v12+10*u4-125*v14+50*v1*v3+25*v2*v2),mod)*inverse(25%mod,mod)%mod;
    m[5]=rem(6*(u2*u3-20*u2*v1*v2-10*u3*v12+5*u5-250*v13*v2+25*v1*v4+25*v2*v3),mod)*inverse(25%mod,mod)%mod;
    m[6]=rem(-(4*u2*u2*u2-270*u2*u2*v12-30*u2*u4-5250*u2*v14
        +600*u2*v1*v3+300*u2*v2*v2-15*u3*u3+600*u3*v1*v2+300*u4*v12
        -15775*v16+7500*v13*v3+11250*v12*v2*v2-750*v1*v5-750*v2*v4-375*v3*v3),mod)*inverse(125%mod,mod)%mod;
    Poly fifth=power(m,5,mod);
    for(int i=0;i<=6;++i) assert(rem(d[i]-fifth[i],mod)==0);
    Vec r{};
    for(int i=7;i<=30;++i) r[i-7]=static_cast<int>(rem(d[i]-fifth[i],mod));
    return r;
}
struct Space {
    Matrix L{};
    std::array<int,4> selected{};
    std::array<std::array<int,4>,4> inv{};
    std::array<Point,5> kernel{};
    Space() {
        Vec r=residual(base,27);
        for(int j=0;j<9;++j) {
            Point q=base; q[j]+=3; Vec d=residual(q,27);
            for(int i=0;i<24;++i) {
                int delta=static_cast<int>(rem(d[i]-r[i],27)); assert(delta%9==0);
                L[i][j]=delta/9;
            }
        }
        auto work=L; std::array<int,24> ids{};
        for(int i=0;i<24;++i) ids[i]=i;
        for(int col=0;col<4;++col) {
            int p=col; while(p<24 && work[p][col]==0) ++p; assert(p<24);
            std::swap(work[col],work[p]); std::swap(ids[col],ids[p]); selected[col]=ids[col];
            int v=inverse(work[col][col],3);
            for(int j=col;j<4;++j) work[col][j]=work[col][j]*v%3;
            for(int i=col+1;i<24;++i) {
                int c=work[i][col];
                for(int j=col;j<4;++j) work[i][j]=static_cast<int>(rem(work[i][j]-c*work[col][j],3));
            }
        }
        std::array<std::array<int,8>,4> square{};
        for(int i=0;i<4;++i) for(int j=0;j<4;++j) {
            square[i][j]=L[selected[i]][j]; square[i][j+4]=(i==j);
        }
        for(int col=0;col<4;++col) {
            int p=col; while(p<4 && !square[p][col]) ++p; assert(p<4);
            std::swap(square[p],square[col]); int v=inverse(square[col][col],3);
            for(int j=0;j<8;++j) square[col][j]=square[col][j]*v%3;
            for(int i=0;i<4;++i) if(i!=col) {
                int c=square[i][col];
                for(int j=0;j<8;++j) square[i][j]=static_cast<int>(rem(square[i][j]-c*square[col][j],3));
            }
        }
        for(int i=0;i<4;++i) for(int j=0;j<4;++j) inv[i][j]=square[i][j+4];
        for(int j=4;j<9;++j) {
            Vec b{}; for(int i=0;i<24;++i) b[i]=static_cast<int>(rem(-L[i][j],3));
            Point p{}; assert(solve(b,p)); p[j]=1; kernel[j-4]=p;
        }
    }
    bool solve(const Vec& b,Point& p) const {
        p.fill(0);
        for(int i=0;i<4;++i) for(int j=0;j<4;++j) p[i]=(p[i]+inv[i][j]*b[selected[j]])%3;
        for(int i=0;i<24;++i) {
            int value=0; for(int j=0;j<9;++j) value=(value+L[i][j]*p[j])%3;
            if(value!=b[i]) return false;
        }
        return true;
    }
};
std::string digest(const std::vector<Point>& points) {
    uint64_t h=UINT64_C(14695981039346656037);
    for(const Point& q:points) {
        std::string line;
        for(int j=0;j<9;++j) { if(j) line+=','; line+=std::to_string(q[j]); } line+='\n';
        for(unsigned char c:line) h=(h^c)*UINT64_C(1099511628211);
    }
    std::ostringstream s; s<<std::hex<<std::setw(16)<<std::setfill('0')<<h; return s.str();
}
int main() {
    Space space;
    std::vector<Point> prefixes={base};
    const int expected[4][3]={{1,1,243},{243,81,19683},{19683,162,39366},{39366,0,0}};
    std::cout<<"{\n  \"implementation\": \"C++17 explicit midpoint coefficients\",\n  \"levels\": [\n";
    int index=0;
    for(int place:{3,9,27,81}) {
        int modulus=9*place, liftable=0; std::sort(prefixes.begin(),prefixes.end());
        std::vector<Point> children;
        for(const Point& q:prefixes) {
            Vec r=residual(q,modulus), b{}; Point p{};
            for(int i=0;i<24;++i) { assert(r[i]%(3*place)==0); b[i]=static_cast<int>(rem(-r[i]/(3*place),3)); }
            if(!space.solve(b,p)) continue;
            ++liftable;
            for(int code=0;code<243;++code) {
                int k=code; Point h=p;
                for(int j=0;j<5;++j) { int c=k%3; k/=3;
                    for(int z=0;z<9;++z) h[z]=(h[z]+c*space.kernel[j][z])%3; }
                Point next{}; for(int j=0;j<9;++j) next[j]=q[j]+place*h[j];
                children.push_back(next);
            }
        }
        assert(static_cast<int>(prefixes.size())==expected[index][0]);
        assert(liftable==expected[index][1]);
        assert(static_cast<int>(children.size())==expected[index][2]);
        if(index) std::cout<<",\n";
        std::cout<<"    {\"modulus\": "<<modulus<<", \"prefixes\": "<<prefixes.size()
                 <<", \"liftable\": "<<liftable<<", \"children\": "<<children.size()
                 <<", \"prefix_fnv1a64\": \""<<digest(prefixes)<<"\"}";
        prefixes=std::move(children); ++index;
    }
    assert(prefixes.empty());
    std::cout<<"\n  ],\n  \"result\": \"NO NORMALIZED LIFT MODULO 729\"\n}\n";
}
