// Complete bounded repeated-summand search; no floating point and no Bloom filter.
// See MULTIPLICITY_SQUARE_CONICS_2026_09_05.md for completeness and overflow proofs.
#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>
using U = unsigned __int128;
using I = uint64_t;
U sixth(I x) { U y=U(x)*x; return y*y*y; }
std::string dec(U x) { if(!x)return "0"; std::string s; while(x){s+=char('0'+x%10);x/=10;} std::reverse(s.begin(),s.end());return s; }
I root(U x,I hi) { I lo=0; while(lo<hi){I m=lo+(hi-lo+1)/2;if(sixth(m)<=x)lo=m;else hi=m-1;}return lo; }
struct Entry { U value; uint32_t b,c; bool operator<(const Entry& r)const{return value<r.value;} };
struct Sieve {
 std::vector<int> mods{64,729,13,19,31,37,43};
 std::vector<std::vector<bool>> masks;
 Sieve(int terms,int scale) {
  for(int m:mods){std::vector<bool> s(m),cur(m);cur[0]=true;
   for(int x=0;x<m;++x)s[int(sixth(I(scale)*x)%m)]=true;
   for(int j=0;j<terms;++j){std::vector<bool> n(m);for(int a=0;a<m;++a)if(cur[a])for(int b=0;b<m;++b)if(s[b])n[(a+b)%m]=true;cur=n;}
   masks.push_back(cur);
  }
 }
 bool pass(U x)const {for(size_t i=0;i<mods.size();++i)if(!masks[i][size_t(x%mods[i])])return false;return true;}
};
bool quotient(I f,I t,I denom,U& out) {
 std::array<I,4> fac{f-t,f+t,f*f+f*t+t*t,f*f-f*t+t*t};
 for(I& x:fac){I g=std::gcd(x,denom);x/=g;denom/=g;}
 if(denom!=1)return false;
 out=1;for(I x:fac)out*=x;return true;
}
std::vector<I> roots7(){std::vector<I> r;for(I x=1;x<117649;++x)if(sixth(x)%117649==1)r.push_back(x);assert(r.size()==6);return r;}
std::vector<I> roots42(){std::vector<I> r;for(I z:roots7())for(I x=z;x<5489031744ULL;x+=117649)if(x%2&&x%3&&sixth(x%64)%64==1&&sixth(x%729)%729==1)r.push_back(x);assert(r.size()==144);return r;}
std::vector<std::array<I,2>> pair_solve(U target,I bound) {
 I a=1,b=root(target,bound);std::vector<std::array<I,2>> out;
 while(a<=b){U sum=sixth(a)+sixth(b);if(sum<target)++a;else if(sum>target)--b;else{out.push_back({a,b});++a;--b;}}
 return out;
}
void selftest(){
 // Positive controls include equal pairs and arbitrary targets: negative ESOP6
 // runs alone would never validate the pair oracle.
 I controls=0;
 for(I a=1;a<=75;++a)for(I b=a;b<=75;++b){auto found=pair_solve(sixth(a)+sixth(b),75);assert(std::find(found.begin(),found.end(),std::array<I,2>{a,b})!=found.end());++controls;}
 for(I t=0;t<150000;++t){std::vector<std::array<I,2>> want;for(I a=1;a<=10;++a)for(I b=a;b<=10;++b)if(sixth(a)+sixth(b)==t)want.push_back({a,b});assert(pair_solve(t,10)==want);}
 for(int terms:{1,2}){Sieve s(terms,7);for(I a=1;a<=100;++a)for(I b=1;b<=100;++b)assert(s.pass(sixth(7*a)+(terms==2?sixth(7*b):0)));}
 auto r7=roots7(),r42=roots42();
 for(I mod:{117649ULL,5489031744ULL}){const auto&r=mod==117649?r7:r42;for(I t=1;t<=300;++t)if(std::gcd(t,mod)==1)for(I f=1;f<=300;++f)if(std::gcd(f,mod)==1){bool actual=(sixth(t)%mod==sixth(f)%mod),got=false;for(I a:r)got|=(a*t)%mod==f;assert(actual==got);}}
 std::cout<<"{\"status\":\"PASS\",\"positive_pair_controls\":"<<controls<<",\"exhaustive_target_controls\":150000,\"sieve_controls\":20000,\"root_join_controls\":\"all unit f,t <= 300 at 7^6 and 42^6\"}\n";
}
int main(int argc,char**argv){
 if(argc==2&&std::string(argv[1])=="--roots") { for(I x:roots42())std::cout<<x<<"\n";return 0; }
 if(argc==2&&std::string(argv[1])=="--quotients") { I f,t,d;while(std::cin>>f>>t>>d){if(f<=t||f>100000000||d<5489031744ULL)throw std::runtime_error("unsafe diagnostic quotient");U q;if(quotient(f,t,d,q))std::cout<<dec(q)<<"\n";else std::cout<<"NONE\n";}return 0; }
 if(argc==2&&std::string(argv[1])=="--pairs") { std::string s;I bound;while(std::cin>>s>>bound){if(s.size()>39||bound>2000000)throw std::runtime_error("unsafe diagnostic pair");U q=0;for(char c:s){if(c<'0'||c>'9'||q>(~U(0)-I(c-'0'))/10)throw std::runtime_error("invalid 128-bit integer");q=q*10+(c-'0');}auto ps=pair_solve(q,bound);std::cout<<"[";for(size_t i=0;i<ps.size();++i){if(i)std::cout<<",";std::cout<<"["<<ps[i][0]<<","<<ps[i][1]<<"]";}std::cout<<"]\n";}return 0; }
 if(argc==2&&std::string(argv[1])=="--r-seeds"){
  I count=0;for(I b=1;b<=1000000;++b){if(std::gcd(b,5489031744ULL)!=1)continue;static auto rays=roots42();for(I ray:rays){I a=ray*b%5489031744ULL;if(a>b&&a<2*b&&std::gcd(a,b)==1){std::cout<<a<<","<<b<<"\n";if(++count==24)return 0;}}}return 0;
 }
 if(argc==2&&std::string(argv[1])=="--self-test"){selftest();return 0;}
 if(argc!=3)throw std::runtime_error("usage: repeated_campaign 221|311|2111 bound; or --self-test");
 std::string family=argv[1];I H=std::stoull(argv[2]);
 if(H<1||H>100000000)throw std::runtime_error("bound must be 1..100000000 (factored-quotient safety limit)");
 if(family!="221"&&family!="311"&&family!="2111")throw std::runtime_error("unknown family");
 if(family!="221" && H>1000000)throw std::runtime_error("311/2111 bound limited to 1000000 for 128-bit sum arithmetic");
 auto start=std::chrono::steady_clock::now();I pairs=0,queries=0,sieve_pass=0,rawhits=0,primitive=0;
 std::vector<std::array<I,6>> hits;std::array<std::vector<Entry>,4> tables;
 auto record=[&](std::array<I,6> v){if(family=="221"){U q;assert(quotient(v[5],v[4],5489031744ULL,q));assert(q==2*sixth(v[0]/42)+2*sixth(v[2]/42));}else assert(sixth(v[0])+sixth(v[1])+sixth(v[2])+sixth(v[3])+sixth(v[4])==sixth(v[5]));++rawhits;I g=0;for(I x:v)g=std::gcd(g,x);if(g==1){++primitive;std::sort(v.begin(),v.begin()+5);hits.push_back(v);}};
 if(family=="2111"){
  for(I b=7;b<H;b+=7)for(I c=b;c<H;c+=7){int n2=(b%2!=0)+(c%2!=0),n3=(b%3!=0)+(c%3!=0);if(n2>1||n3>1)continue;
   U v=sixth(b)+sixth(c);if(v+2*sixth(42)+1>sixth(H))break;
   tables[n2+2*n3].push_back({v,uint32_t(b),uint32_t(c)});
  }
  for(auto&tab:tables)std::sort(tab.begin(),tab.end());
  size_t n=0;for(auto&tab:tables)n+=tab.size();std::cerr<<"pair table entries "<<n<<"\n";
 }
 I mod=family=="221"?5489031744ULL:117649ULL;auto rays=family=="221"?roots42():roots7();
 Sieve sieve(family=="311"?1:2,family=="221"?1:7);
 for(I t=1;t<H;++t){if(std::gcd(t,mod)!=1)continue;
  for(I ray:rays){I f0=(ray*t)%mod;for(I f=f0;f<=H;f+=mod){if(f<=t||f%2==0||f%3==0)continue;++pairs;
   if(family=="221"){
    U target;if(!quotient(f,t,2*mod,target))continue;++queries;
    if(target<2||!sieve.pass(target))continue;
    ++sieve_pass;
    for(auto ab:pair_solve(target,H/42))record({42*ab[0],42*ab[0],42*ab[1],42*ab[1],t,f});
   }else{
    U delta=sixth(f)-sixth(t);
    I mult=family=="311"?3:2;int need2=1-(t%2!=0),need3=1-(t%3!=0);int key=need2+2*need3;
    U minrest=family=="311"?sixth(7):2*sixth(7);
    for(I a=42;a<H;a+=42){U used=mult*sixth(a);if(used+minrest>delta)break;U target=delta-used;++queries;
     if(!sieve.pass(target))continue;
     ++sieve_pass;
     if(family=="311"){I b=root(target,H);if(b&&b%7==0&&int(b%2!=0)==need2&&int(b%3!=0)==need3&&sixth(b)==target)record({a,a,a,b,t,f});}
     else{const auto&tab=tables[key];auto it=std::lower_bound(tab.begin(),tab.end(),Entry{target,0,0});for(;it!=tab.end()&&it->value==target;++it)record({a,a,it->b,it->c,t,f});}
    }
   }
  }}
 }
 std::sort(hits.begin(),hits.end());hits.erase(std::unique(hits.begin(),hits.end()),hits.end());
 double sec=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
 std::cout<<"{\"family\":\""<<family<<"\",\"max_f\":"<<H<<",\"complete\":true,\"rhs_singleton_pairs\":"<<pairs<<",\"queries\":"<<queries<<",\"local_pass\":"<<sieve_pass<<",\"raw_hits\":"<<rawhits<<",\"primitive_occurrences\":"<<primitive<<",\"unique_primitive_hits\":"<<hits.size()<<",\"pair_table_sizes\":[";
 for(int i=0;i<4;++i){if(i)std::cout<<',';std::cout<<tables[i].size();}std::cout<<"],\"seconds\":"<<sec<<",\"solutions\":[";
 for(size_t j=0;j<hits.size();++j){if(j)std::cout<<',';std::cout<<'[';for(int i=0;i<6;++i){if(i)std::cout<<',';std::cout<<hits[j][i];}std::cout<<']';}
 std::cout<<"]}\n";
}
