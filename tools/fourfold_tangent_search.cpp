#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <vector>
using I=__int128_t;
using V=std::array<I,6>;
I ab(I x){return x<0?-x:x;}
I gcd(I a,I b){a=ab(a);b=ab(b);while(b){I r=a%b;a=b;b=r;}return a;}
std::string str(I n){if(!n)return "0";bool neg=n<0;n=ab(n);std::string s;while(n){s.push_back('0'+n%10);n/=10;}if(neg)s.push_back('-');std::reverse(s.begin(),s.end());return s;}
I cube(I x){return x*x*x;}
I form(const V&x){I s=0;for(int i=0;i<6;i++)s+=(i==5?-1:1)*cube(x[i]);return s;}
I root(I n){I r=(I)sqrtl((long double)n);while(r*r>n)--r;while((r+1)*(r+1)<=n)++r;return r;}
bool square(I n){if(n<0)return false;I r=root(n);return r*r==n;}
int main(int argc,char**argv){
 int H=argc>1?atoi(argv[1]):24, maxseed=argc>2?atoi(argv[2]):4;
 // These bounds imply |C_i| <= 24*M*(4*M*M*H)^3 < 3.48e11,
 // so 6*max|C_i|^3 < 2.6e35 < 2^127. All cubic checks are exact.
 if(H<1||H>24||maxseed<1||maxseed>4){std::cerr<<"Supported bounds: 1<=H<=24, 1<=maxseed<=4\n";return 2;}
 const int primes[]={16,63,65,11,13,17,19,23,31};
 std::vector<std::vector<bool>> masks;
 for(int p:primes){std::vector<bool> m(p);for(int j=0;j<p;j++)m[(j*j)%p]=true;masks.push_back(m);}
 uint64_t directions=0,positive=0,checked=0,modpass=0,hits=0;I largest=0;std::array<int,3> bestseed{};V best{};int bestscore=-1;
 for(int a=1;a<=maxseed;a++)for(int b=a;b<=maxseed;b++)for(int c=1;c<=maxseed;c++){
  if(std::gcd(std::gcd(a,b),c)>1||(a==c&&b==c))continue;
  V B={(I)a,-a,(I)b,-b,(I)c,(I)c};uint64_t dp=0,pp=0;
  for(int x1=-H;x1<=H;x1++)for(int x2=-H;x2<=H;x2++)for(int x3=-H;x3<=H;x3++)for(int x4=-H;x4<=H;x4++){
   // D5=0 fixes addition of B; quotient out sign and integral scaling.
   int first=x1?x1:x2?x2:x3?x3:x4;if(first<=0)continue;
   int g=std::gcd(std::gcd(abs(x1),abs(x2)),std::gcd(abs(x3),abs(x4)));if(g!=1)continue;
   I L=(I)a*a*(x1+x2)+(I)b*b*(x3+x4),cc=c*c;
   V D={cc*x1,cc*x2,cc*x3,cc*x4,0,L};I gd=0;for(I v:D)gd=gcd(gd,v);for(I&v:D)v/=gd;
   directions++;dp++;
   I f=form(D),k=0;for(int i=0;i<6;i++)k+=(i==5?-1:1)*B[i]*D[i]*D[i];
   if(!f||!k)continue;
   V C;for(int i=0;i<6;i++)C[i]=f*B[i]-3*k*D[i];
   if(C[5]<0)for(I&v:C)v=-v;
   bool pos=true;for(I v:C)if(v<=0)pos=false;if(!pos)continue;
   positive++;pp++;I gc=0;for(I v:C)gc=gcd(gc,v);for(I&v:C)v/=gc;
   largest=std::max(largest,C[5]);
   // Check cubic identity before trusting any square sieve.
   if(form(C)!=0){std::cerr<<"CUBIC IDENTITY FAILURE\n";return 2;}checked++;
   int score=0;for(I v:C)if(square(v))score++;
   if(score>bestscore){bestscore=score;best=C;bestseed={a,b,c};}
   bool pass=true;for(I v:C){for(size_t j=0;j<masks.size();j++)if(!masks[j][(int)(v%primes[j])]){pass=false;break;}if(!pass)break;}
   if(!pass)continue;
   modpass++;
   if(score==6){hits++;std::cout<<"HIT [";for(int i=0;i<6;i++)std::cout<<(i?",":"")<<str(root(C[i]));std::cout<<"]\n";return 0;}
  }
  std::cout<<"seed "<<a<<","<<b<<","<<c<<" directions="<<dp<<" positive="<<pp<<std::endl;
 }
 std::cout<<"{\"bound\":"<<H<<",\"max_seed\":"<<maxseed<<",\"directions\":"<<directions<<",\"positive\":"<<positive<<",\"exact_cubic_checks\":"<<checked<<",\"square_sieve_pass\":"<<modpass<<",\"hits\":"<<hits<<",\"largest_primitive_cube_base\":\""<<str(largest)<<"\",\"best_square_count\":"<<bestscore<<",\"best_seed\":["<<bestseed[0]<<","<<bestseed[1]<<","<<bestseed[2]<<"],\"best_primitive_cube_bases\":[";
 for(int i=0;i<6;i++)std::cout<<(i?",":"")<<"\""<<str(best[i])<<"\"";
 std::cout<<"]}\n";
}
