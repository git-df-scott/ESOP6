/* Cheap independent replay of every concentrated-case (f,t) candidate set. */
#include <stdint.h>
#include <stdio.h>

typedef unsigned __int128 u128;
typedef uint64_t u64;
static const u64 M=5489031744ULL;
static u64 roots[160]; static int nroots;

static u64 pow6m(u64 x,u64 m){u64 r=x%m;u128 s=(u128)r*r%m;s=s*s%m;return (u64)((u128)s*((u128)r*r%m)%m);}
static u64 mix(u64 h){h^=h>>33;h*=0xff51afd7ed558ccdULL;h^=h>>33;h*=0xc4ceb9fe1a85ec53ULL;h^=h>>33;return h;}
static void build_roots(void){
    u64 r2[8],r3[8],r7[8]; int n2=0,n3=0,n7=0;
    for(u64 x=1;x<64;x+=2)if(pow6m(x,64)==1)r2[n2++]=x;
    for(u64 x=1;x<729;x++)if(x%3&&pow6m(x,729)==1)r3[n3++]=x;
    for(u64 x=1;x<117649;x++)if(x%7&&pow6m(x,117649)==1)r7[n7++]=x;
    for(int i=0;i<n2;i++)for(int j=0;j<n3;j++)for(int k=0;k<n7;k++){
        u64 x12=0,x=0;
        for(u64 t=r2[i];;t+=64)if(t%729==r3[j]){x12=t;break;}
        for(u64 t=x12;;t+=64*729ULL)if(t%117649==r7[k]){x=t;break;}
        roots[nroots++]=x;
    }
}
static u64 audit(u64 lo,u64 hi,u64 *digest_xor,u64 *digest_sum){
    u64 count=0,dx=0,ds=0;
    #pragma omp parallel for reduction(+:count,ds) reduction(^:dx) schedule(static)
    for(u64 f=lo;f<=hi;f++){
        if(f%2==0||f%3==0||f%7==0)continue;
        for(int i=0;i<nroots;i++){
            u64 t=(u64)((u128)roots[i]*f%M);
            if(t==0||t>=f)continue;
            u64 h=mix(f^(0x9e3779b97f4a7c15ULL*t));
            count++; dx^=h; ds+=h;
        }
    }
    *digest_xor=dx;*digest_sum=ds;return count;
}
int main(void){
    static const u64 bands[][3]={
        {700000,730000,124},{730000,1000000,1314},{1000000,1500000,3523},
        {1500000,2000000,5129},{2000000,2500000,7222},{2500000,3200000,12443},
        {3200000,4000000,18003},{4000000,4300000,8050}
    };
    build_roots(); if(nroots!=144){fprintf(stderr,"root count %d\n",nroots);return 1;}
    u64 production_total=0;
    for(size_t i=0;i<sizeof(bands)/sizeof(bands[0]);i++){
        u64 dx,ds,n=audit(bands[i][0],bands[i][1],&dx,&ds);
        printf("[%llu,%llu] candidates=%llu expected=%llu xor=%016llx sum=%016llx %s\n",
            (unsigned long long)bands[i][0],(unsigned long long)bands[i][1],
            (unsigned long long)n,(unsigned long long)bands[i][2],
            (unsigned long long)dx,(unsigned long long)ds,n==bands[i][2]?"PASS":"FAIL");
        if(i)production_total+=n;
        if(n!=bands[i][2])return 1;
    }
    printf("production_total=%llu expected=55684 %s\n",(unsigned long long)production_total,
        production_total==55684?"PASS":"FAIL");
    return production_total==55684?0:1;
}
