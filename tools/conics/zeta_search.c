/* Exhaustive search for rational conics in the design stratum of the Fermat sextic fourfold:
   Gaussian integers zeta_1..zeta_m (m = 5 for (6,1,5), m = 6 for (6,1,6)) with
       S3 = sum zeta_i^6 = 0,   S2 = sum N_i zeta_i^4 = 0,   S1 = sum N_i^2 zeta_i^2 = 0,   N_i = |zeta_i|^2.
   Enumerate zeta_1..zeta_{m-1} up to sign (b>0 or b==0,a>0), with |zeta|^2 <= H, sorted;
   the last one is forced: N_m^3 = |S1'| where S1' = -sum_{i<m} N_i^2 zeta_i^2, zeta_m^2 = S1'/N_m^2.
   Also reports near-misses: 4-tuples where S1' has sixth-power norm (rare filter).
   usage: zeta_search m H  */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
typedef long long ll; typedef __int128 i128;
static ll isqrt_ll(i128 n){ if(n<0) return -1; ll r=(ll)sqrtl((long double)n); while((i128)r*r>n) r--; while((i128)(r+1)*(r+1)<=n) r++; return r; }
static ll icbrt_ll(ll n){ if(n<0) return -1; ll r=(ll)cbrtl((long double)n); while(r*r*r>n) r--; while((r+1)*(r+1)*(r+1)<=n) r++; return r; }
typedef struct { ll a,b; } G;
static G gmul(G x,G y){ G r={x.a*y.a-x.b*y.b, x.a*y.b+x.b*y.a}; return r; }
static ll gnorm(G x){ return x.a*x.a+x.b*x.b; }
int main(int argc,char**argv){
    int m = atoi(argv[1]); ll H = atoll(argv[2]);
    int lim = (int)sqrt((double)H)+1;
    G *Z = malloc(sizeof(G)*4*lim*lim); int n=0;
    for(ll a=-lim;a<=lim;a++) for(ll b=0;b<=lim;b++){ if(a*a+b*b==0||a*a+b*b>H) continue; if(b==0&&a<0) continue; Z[n].a=a; Z[n].b=b; n++; }
    fprintf(stderr,"m=%d H=%lld gaussian integers up to sign: %d\n",m,H,n);
    /* precompute powers */
    ll *N=malloc(sizeof(ll)*n); G *z2=malloc(sizeof(G)*n),*z4=malloc(sizeof(G)*n),*z6=malloc(sizeof(G)*n);
    for(int i=0;i<n;i++){ N[i]=gnorm(Z[i]); z2[i]=gmul(Z[i],Z[i]); z4[i]=gmul(z2[i],z2[i]); z6[i]=gmul(z4[i],z2[i]); }
    ll found=0, nearmiss=0; ll count=0;
    int k = m-1; int idx[8];
    /* iterate sorted k-subsets (with repetition allowed) of indices */
    for(int i=0;i<k;i++) idx[i]=0;
    while(1){
        /* compute S1' */
        i128 s1a=0,s1b=0,s2a=0,s2b=0,s3a=0,s3b=0;
        for(int i=0;i<k;i++){ int t=idx[i]; i128 N2=(i128)N[t]*N[t];
            s1a+=N2*z2[t].a; s1b+=N2*z2[t].b; s2a+=(i128)N[t]*z4[t].a; s2b+=(i128)N[t]*z4[t].b; s3a+=z6[t].a; s3b+=z6[t].b; }
        count++;
        /* need zeta_m^2 * N_m^2 = -(s1) ; N_m^3 = |s1| => |s1|^2 = N_m^6 */
        i128 ns1 = s1a*s1a+s1b*s1b;
        if(ns1!=0){
            ll r = isqrt_ll(ns1);
            if((i128)r*r==ns1){
                ll c = icbrt_ll(r);
                if(c*c*c==r){
                    nearmiss++;
                    /* zeta_m^2 = -s1 / c^2 */
                    i128 c2=(i128)c*c;
                    if(s1a % c2==0 && s1b % c2==0){
                        ll wa=(ll)(-s1a/c2), wb=(ll)(-s1b/c2);
                        /* square root in Z[i]: solve (x+iy)^2 = wa + i wb */
                        ll nw = isqrt_ll((i128)wa*wa+(i128)wb*wb);
                        if((i128)nw*nw==(i128)wa*wa+(i128)wb*wb){
                            ll x2 = (nw+wa); /* 2x^2 = nw + wa */
                            if(x2%2==0){ ll x=isqrt_ll(x2/2); if(x*x==x2/2){
                                ll y2=(nw-wa); if(y2%2==0){ ll y=isqrt_ll(y2/2); if(y*y==y2/2){
                                    for(int sy=-1;sy<=1;sy+=2){ if(2*x*(sy*y)!=wb) continue;
                                        G zm={x,sy*y}; G q2=gmul(zm,zm),q4=gmul(q2,q2),q6=gmul(q4,q2); ll Nm=gnorm(zm);
                                        i128 t2a=s2a+(i128)Nm*q4.a, t2b=s2b+(i128)Nm*q4.b, t3a=s3a+q6.a, t3b=s3b+q6.b;
                                        if(t2a==0&&t2b==0&&t3a==0&&t3b==0){
                                            found++; printf("SOLUTION m=%d:",m);
                                            for(int i=0;i<k;i++) printf(" (%lld%+lldi)",Z[idx[i]].a,Z[idx[i]].b);
                                            printf(" (%lld%+lldi)\n",zm.a,zm.b); fflush(stdout);
                                        }
                                        break; }
                                }}}}
                        }
                    }
                }
            }
        }
        /* next multiset */
        int p=k-1; while(p>=0 && idx[p]==n-1) p--; if(p<0) break; idx[p]++; for(int i=p+1;i<k;i++) idx[i]=idx[p];
    }
    printf("done m=%d H=%lld: tuples=%lld sixth-power-norm passes=%lld SOLUTIONS=%lld\n",m,H,count,nearmiss,found);
    return 0;
}
