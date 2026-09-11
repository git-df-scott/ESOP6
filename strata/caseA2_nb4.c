/* caseA2_nb4.c: strata/caseA2_timed.c with the nb=3,2,1 passes removed, so that its
 * j2_nodes counter is the nb=4-only leaf count that strata/v3engine.c reports.
 * Nothing else is changed. */
/*
 * caseA2_timed.c — INSTRUMENTED COPY of src/caseA2.c (strata/, TASK 2).
 * Semantics are IDENTICAL to caseA2.c: no filter, bound, stride, mask or
 * traversal order is changed. The only additions are (a) clock_gettime
 * wall/core timers around the Bloom build, the candidate enumeration and
 * the DFS traversal, and (b) per-thread event counters (j=2 leaves, Bloom
 * queries/positives, pair_verify calls). Counters are thread-local and
 * merged once at the end, so no atomic is added to the hot path.
 *
 * caseA2.c — deep search of the "concentrated exemption" case of
 * a^6+b^6+c^6+d^6+e^6 = f^6, Bloom-filter edition.
 *
 * Background (see caseA.c/search.c): in a primitive solution exactly one
 * term is odd, one coprime to 3, one coprime to 7. When all three land on
 * one term t, the remaining four terms are divisible by 42 and
 *      f^6 ≡ t^6 (mod 42^6),  f,t coprime to 42,
 * so t/f is one of the 144 sixth roots of unity mod M = 42^6. Surviving
 * pairs (f,t) need m = (f^6-t^6)/42^6 to be a sum of four sixth powers
 * with bases <= (f-1)/42.
 *
 * This version answers "is R a sum of two sixth powers?" with a Bloom
 * filter over ALL pair sums x^6+y^6 (x>=y>=1, x<=Bmax), built once and
 * shared: the 4-part test enumerates the two largest parts (with exact
 * class budgets from m mod 8/9/7 forcing strides) and probes the filter.
 * Positives are verified exactly. Residue masks mod 64,27,49,13,43 prune.
 *
 * Build:  gcc -O3 -march=native -fopenmp -o caseA2 caseA2.c -lm
 * Run:    ./caseA2 <fmin> <fmax>   (fmax <= 1e8 hard; RAM limits ~4e6)
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <omp.h>

/* ---- instrumentation (additive only; no semantic effect) ---- */
#define MAXTH 256
typedef struct __attribute__((aligned(128))) { double t_dfs; long long j2_nodes, bloom_q, bloom_pos, pv_calls; char pad[128]; } tstat_t;
static tstat_t TS[MAXTH];
static int th_of(void){ int i=omp_get_thread_num(); return (i<MAXTH)?i:MAXTH-1; }
static double wnow(void){ struct timespec ts; clock_gettime(CLOCK_MONOTONIC,&ts); return ts.tv_sec+1e-9*ts.tv_nsec; }
static double cnow(void){ struct timespec ts; clock_gettime(CLOCK_THREAD_CPUTIME_ID,&ts); return ts.tv_sec+1e-9*ts.tv_nsec; }

typedef unsigned __int128 u128;
typedef uint64_t u64;
typedef uint32_t u32;

static const u64 M = 5489031744ULL;   /* 42^6 */
static u64 roots[160]; static int nroots=0;

static u64 Bmax;                       /* max base in the pair table */
static u128 *P6;                       /* P6[x]=x^6, x<=Bmax */

static u128 ipow6(u64 x){ u128 s=(u128)x*x; return s*s*s; }
static u64 pow6m(u64 x,u64 m){ u64 r=x%m; u128 s=(u128)r*r%m; s=s*s%m; return (u64)((u128)s*((u128)r*r%m)%m); }

static void build_roots(void){
    u64 r2[8]; int n2=0; for(u64 x=1;x<64;x+=2) if(pow6m(x,64)==1) r2[n2++]=x;
    u64 r3[8]; int n3=0; for(u64 x=1;x<729;x++) if(x%3&&pow6m(x,729)==1) r3[n3++]=x;
    u64 r7[8]; int n7=0; for(u64 x=1;x<117649;x++) if(x%7&&pow6m(x,117649)==1) r7[n7++]=x;
    for(int i=0;i<n2;i++)for(int j=0;j<n3;j++)for(int k=0;k<n7;k++){
        u64 m12=64*729ULL,x12=0,x=0;
        for(u64 t=r2[i];;t+=64) if(t%729==r3[j]){x12=t;break;}
        for(u64 t=x12;;t+=m12) if(t%117649==r7[k]){x=t;break;}
        roots[nroots++]=x;
    }
}

/* residue masks */
#define NM 4
static const u32 MODS[NM]={27,49,13,43};
static u64 sumres[NM][5], sumres64[5];
static u32 powmod6[NM][64];
static void init_sieves(void){
    for(int m=0;m<NM;m++){
        u32 md=MODS[m]; u64 single=0;
        for(u32 x=0;x<md;x++){u64 r=pow6m(x,md);powmod6[m][x]=(u32)r;single|=1ULL<<r;}
        sumres[m][0]=1;
        for(int j=1;j<5;j++){
            u64 prev=sumres[m][j-1],cur=0;
            for(u32 r=0;r<md;r++) if(prev>>r&1)
                for(u32 s=0;s<md;s++) if(single>>s&1)
                    cur|=1ULL<<((r+s)%md);
            sumres[m][j]=cur;
        }
    }
    u64 single=0;
    for(u32 x=0;x<64;x++) single|=1ULL<<(u32)pow6m(x,64);
    sumres64[0]=1;
    for(int j=1;j<5;j++){
        u64 prev=sumres64[j-1],cur=0;
        for(u32 r=0;r<64;r++) if(prev>>r&1)
            for(u32 s=0;s<64;s++) if(single>>s&1)
                cur|=1ULL<<((r+s)&63);
        sumres64[j]=cur;
    }
}

static u64 iroot6(u128 R){
    if(!R) return 0;
    u64 r=(u64)powl((long double)R,1.0L/6);
    while(ipow6(r+1)<=R) r++;
    while(r&&ipow6(r)>R) r--;
    return r;
}

/* ---- blocked Bloom filter over pair sums: 8 bits in one 64-byte line ---- */
static u64 *bloom; static u64 nlines;      /* 8 u64 per cache line */
static inline u64 mix(u64 h){ h^=h>>33; h*=0xff51afd7ed558ccdULL; h^=h>>33; h*=0xc4ceb9fe1a85ec53ULL; h^=h>>33; return h; }
#define BLOOMK 8
static inline void bloom_hashes(u128 s,u64 *h1,u64 *h2){
    u64 lo=(u64)s, hi=(u64)(s>>64);
    *h1=mix(lo^(0x9e3779b97f4a7c15ULL*hi));
    *h2=mix(hi^(0xbf58476d1ce4e5b9ULL+lo))|1ULL;
}
/* Lemire range reduction: map h1 uniformly into [0,nlines) without requiring
 * nlines to be a power of two. Rounding the line count up to a power of two
 * (the previous scheme) silently allocated up to 2x the requested
 * bits-per-pair -- e.g. 17.2 GB for a 12 bpp request at fmax=5e6. */
static inline u64 line_of(u64 h1){ return (u64)(((u128)h1*(u128)nlines)>>64); }
static inline void bloom_add(u128 s){
    u64 h1,h2; bloom_hashes(s,&h1,&h2);
    u64 *line=bloom+(line_of(h1)<<3);
    for(int i=0;i<BLOOMK;i++){ u32 b=(u32)((h2>>(8*i))&511);
        __atomic_or_fetch(&line[b>>6],1ULL<<(b&63),__ATOMIC_RELAXED); }
}
static inline int bloom_query(u128 s){
    TS[th_of()].bloom_q++;
    u64 h1,h2; bloom_hashes(s,&h1,&h2);
    const u64 *line=bloom+(line_of(h1)<<3);
    for(int i=0;i<BLOOMK;i++){ u32 b=(u32)((h2>>(8*i))&511);
        if(!(line[b>>6]>>(b&63)&1)) return 0; }
    TS[th_of()].bloom_pos++;
    return 1;
}

/* exact: R = x^6+y^6, x>=y>=1, x<=maxv, class budgets (o2,o3,o7 in 0..2) */
static int pair_verify(u128 R,u64 maxv,int o2,int o3,int o7,u64 *ox,u64 *oy){
    u64 x=iroot6(R); if(x>maxv) x=maxv;
    u64 xmin=iroot6((R+1)/2); if(ipow6(xmin)*2<R) xmin++;
    for(;x>=xmin&&x>0;x--){
        int x2=(int)(x&1),x3=(x%3!=0),x7=(x%7!=0);
        if(x2>o2||x3>o3||x7>o7) continue;
        u128 rest=R-P6[x];
        u64 y=iroot6(rest);
        if(y==0||y>x||ipow6(y)!=rest) continue;
        int y2=(int)(y&1),y3=(y%3!=0),y7=(y%7!=0);
        if(x2+y2!=o2||x3+y3!=o3||x7+y7!=o7) continue;
        *ox=x;*oy=y; return 1;
    }
    return 0;
}

static long long nfound=0;
static int use_valuation_prune=0;
static long long valuation_rejects=0;

/* If S is a sum of at most four sixth powers, repeated use of
 * x^6 in {0,1} modulo 8, 9, and 7 forces
 *   v2(S) mod 6 in {0,1,2}, v3(S) mod 6 in {0,1}, v7(S) mod 6 = 0.
 * This is optional so the historical node-count control remains unchanged. */
static int valuation_possible(u128 s){
    int v2=0,v3=0,v7=0;
    while(s%2==0){ s/=2; v2++; }
    while(s%3==0){ s/=3; v3++; }
    while(s%7==0){ s/=7; v7++; }
    return v2%6<=2 && v3%6<=1 && v7%6==0;
}

static void report(u64 f,u64 t,u64 *b,int nb){
    const char *kind=(nb==4)?"SOLUTION":"DEGENERATE";
    printf("%s f=%llu t=%llu parts42:",kind,(unsigned long long)f,(unsigned long long)t);
    for(int i=0;i<nb;i++) printf(" %llu",(unsigned long long)(b[i]*42));
    printf("\n"); fflush(stdout);
    if(nb==4){
        #pragma omp atomic
        nfound++;
    }
}

/* R as sum of j sixth powers (j=2,3,4), bases<=maxv, exact budgets. */
static int dfs(u128 R,int j,u64 maxv,int o2,int o3,int o7,
               u32 r27,u32 r49,u32 r13,u32 r43,u64 *out){
    if(o2>j||o3>j||o7>j) return 0;
    if(!(sumres64[j]>>((u32)R&63)&1)) return 0;
    if(!(sumres[0][j]>>r27&1)) return 0;
    if(!(sumres[1][j]>>r49&1)) return 0;
    if(!(sumres[2][j]>>r13&1)) return 0;
    if(!(sumres[3][j]>>r43&1)) return 0;
    if(j==1){
        u64 e=iroot6(R);
        if(e==0||e>maxv||ipow6(e)!=R) return 0;
        if((int)(e&1)!=o2||(int)(e%3!=0)!=o3||(int)(e%7!=0)!=o7) return 0;
        out[0]=e; return 1;
    }
    if(j==2){
        TS[th_of()].j2_nodes++;
        if(R>0 && !bloom_query(R)) return 0;      /* covers all pairs with x<=Bmax>=maxv */
        TS[th_of()].pv_calls++;
        return pair_verify(R,maxv,o2,o3,o7,out,out+1);
    }
    u64 hi=iroot6(R); if(hi>maxv) hi=maxv;
    u64 lo=iroot6((R+j-1)/j);
    if(ipow6(lo)*(u128)j<R) lo++;
    int d2=(o2==0)?2:1,d3=(o3==0)?3:1,d7=(o7==0)?7:1;
    int stride=d2*d3*d7;
    u64 x=hi; if(stride>1) x-=x%stride;
    for(;x>=lo&&x>0;x-=(u64)stride){
        int x2=(int)(x&1),x3=(x%3!=0),x7=(x%7!=0);
        if(o2==j&&!x2) continue;
        if(o2==0&&x2) continue;
        if(o3==j&&!x3) continue;
        if(o3==0&&x3) continue;
        if(o7==j&&!x7) continue;
        if(o7==0&&x7) continue;
        out[0]=x;
        u32 p27=powmod6[0][x%27],p49=powmod6[1][x%49],p13=powmod6[2][x%13],p43=powmod6[3][x%43];
        if(dfs(R-P6[x],j-1,x,o2-x2,o3-x3,o7-x7,
               r27+((r27<p27)?27:0)-p27,r49+((r49<p49)?49:0)-p49,
               r13+((r13<p13)?13:0)-p13,r43+((r43<p43)?43:0)-p43,out+1)) return 1;
        if(x<(u64)stride) break;
    }
    return 0;
}

static int try_decompose(u128 m,int nb,u64 maxb,u64 *out){
    if(use_valuation_prune && !valuation_possible(m)){
        #pragma omp atomic
        valuation_rejects++;
        return 0;
    }
    int k2=(int)(m%8),k3=(int)(m%9),k7=(int)(m%7);
    if(k2>nb||k3>nb||k7>nb) return 0;
    return dfs(m,nb,maxb,k2,k3,k7,(u32)(m%27),(u32)(m%49),(u32)(m%13),(u32)(m%43),out);
}

static long long ncand=0,nsurv=0;

int main(int argc,char **argv){
    if(argc<3){ fprintf(stderr,"usage: %s <fmin> <fmax> [bits_per_pair] [--dump-candidates] [--valuation-prune]\n",argv[0]); return 2; }
    u64 fmin=strtoull(argv[1],0,10),fmax=strtoull(argv[2],0,10);
    u64 bpp=16; int dump_candidates=0;
    for(int i=3;i<argc;i++){
        if(!strcmp(argv[i],"--dump-candidates")) dump_candidates=1;
        else if(!strcmp(argv[i],"--valuation-prune")) use_valuation_prune=1;
        else if(argv[i][0]!='-') bpp=strtoull(argv[i],0,10);
    }
    if(fmax>100000000ULL){ fprintf(stderr,"fmax capped at 1e8\n"); return 2; }
    Bmax=(fmax-1)/42+1;
    build_roots(); init_sieves();
    if(nroots!=144){ fprintf(stderr,"root count %d != 144\n",nroots); return 1; }
    for(int i=0;i<nroots;i++) if(pow6m(roots[i],M)!=1){ fprintf(stderr,"bad root\n"); return 1; }

    double p6_w0=wnow();
    P6=malloc(sizeof(u128)*(Bmax+1));
    for(u64 x=0;x<=Bmax;x++) P6[x]=ipow6(x);
    double p6_w=wnow()-p6_w0;

    /* blocked Bloom: ~16 bits per pair (argv[3] overrides), 512-bit lines.
     * Sized exactly -- no power-of-two rounding. */
    u64 npairs=Bmax*(Bmax+1)/2;
    nlines = (npairs*bpp + 511)/512; if(!nlines) nlines=1;
    bloom=calloc(nlines*64,1);
    if(!bloom){ fprintf(stderr,"bloom alloc failed (%llu bytes)\n",(unsigned long long)(nlines*64)); return 1; }
    fprintf(stderr,"building bloom: Bmax=%llu pairs=%llu lines=%llu (%.1f GB)\n",
        (unsigned long long)Bmax,(unsigned long long)npairs,
        (unsigned long long)nlines,nlines*64.0/1e9);
    double build_w0=wnow(); double build_core=0;
    #pragma omp parallel reduction(+:build_core)
    {   double c0=cnow();
        #pragma omp for schedule(dynamic,64)
        for(u64 x=1;x<=Bmax;x++)
            for(u64 y=1;y<=x;y++)
                bloom_add(P6[x]+P6[y]);
        build_core=cnow()-c0;
    }
    double build_w=wnow()-build_w0;
    fprintf(stderr,"bloom built\n");

    {   /* self-tests; dedicated powers keep tiny fmax runs in bounds */
        u128 *bigP6=P6;
        const u64 TB=200;
        u128 *testP6=malloc(sizeof(u128)*(TB+1));
        if(!testP6){ fprintf(stderr,"self-test powers alloc failed\n"); return 1; }
        for(u64 x=0;x<=TB;x++) testP6[x]=ipow6(x);
        P6=testP6;
        u64 out[4];
        if(!bloom_query(P6[5]+P6[9])){ fprintf(stderr,"BLOOM SELF-TEST FAILED\n"); return 1; }
        if(!try_decompose(ipow6(5)+ipow6(9)+ipow6(11)+ipow6(14),4,100,out))
            { fprintf(stderr,"SELF-TEST 1 FAILED\n"); return 1; }
        if(!try_decompose(ipow6(42)+ipow6(84)+ipow6(126)+ipow6(168),4,200,out))
            { fprintf(stderr,"SELF-TEST 2 FAILED\n"); return 1; }
        if(try_decompose((u128)12345677,4,100,out))
            { fprintf(stderr,"SELF-TEST 3 FAILED\n"); return 1; }
        free(testP6);
        P6=bigP6;
    }
    memset(TS,0,sizeof(TS));   /* exclude self-test events */

    double loop_w0=wnow(); double loop_core=0;
    #pragma omp parallel reduction(+:ncand,nsurv,loop_core)
    { double lc0=cnow();
    #pragma omp for schedule(dynamic,4096)
    for(u64 f=fmin;f<=fmax;f++){
        if(f%2==0||f%3==0||f%7==0) continue;
        for(int i=0;i<nroots;i++){
            u64 t=(u64)((u128)roots[i]*f%M);
            if(t==0||t>=f) continue;
            ncand++;
            if(dump_candidates){
                #pragma omp critical(candidate_dump)
                printf("CANDIDATE %llu %llu\n",(unsigned long long)f,(unsigned long long)t);
            }
            u64 g[2]={f-t,f+t};
            u128 qq[2]={(u128)f*f-(u128)f*t+(u128)t*t,(u128)f*f+(u128)f*t+(u128)t*t};
            u64 rem=M; u128 m=1;
            for(int s=0;s<2;s++){
                u64 a=g[s],b=rem; while(b){u64 w=a%b;a=b;b=w;}
                rem/=a; m*=g[s]/a;
            }
            for(int s=0;s<2;s++){
                u64 a=(u64)(qq[s]%rem);
                if(a){u64 x=a,y=rem;while(y){u64 w=x%y;x=y;y=w;}a=x;} else a=rem;
                rem/=a; m*=qq[s]/a;
            }
            if(rem!=1){ fprintf(stderr,"WARN rem!=1 f=%llu\n",(unsigned long long)f); continue; }
            nsurv++;
            u64 out[4],maxb=(f-1)/42;
            double d0=cnow();
            if(try_decompose(m,4,maxb,out)) report(f,t,out,4);
            TS[th_of()].t_dfs += cnow()-d0;
        }
    }
    loop_core=cnow()-lc0;
    }
    double loop_w=wnow()-loop_w0;
    fprintf(stderr,"done: f in [%llu,%llu] candidates=%lld processed=%lld found=%lld\n",
        (unsigned long long)fmin,(unsigned long long)fmax,ncand,nsurv,nfound);
    if(use_valuation_prune) fprintf(stderr,"valuation_rejects=%lld\n",valuation_rejects);
    {   double dfs_core=0; long long j2=0,bq=0,bp=0,pv=0;
        for(int i=0;i<MAXTH;i++){ dfs_core+=TS[i].t_dfs; j2+=TS[i].j2_nodes;
            bq+=TS[i].bloom_q; bp+=TS[i].bloom_pos; pv+=TS[i].pv_calls; }
        double enum_core = loop_core - dfs_core; if(enum_core<0) enum_core=0;
        int nth = omp_get_max_threads();
        fprintf(stderr,
          "TIMING threads=%d\n"
          "TIMING p6_table_wall=%.3f\n"
          "TIMING bloom_build_wall=%.3f bloom_build_core=%.3f\n"
          "TIMING search_loop_wall=%.3f search_loop_core=%.3f\n"
          "TIMING   enumeration_core=%.3f\n"
          "TIMING   dfs_traversal_core=%.3f\n"
          "COUNT j2_nodes=%lld bloom_queries=%lld bloom_positives=%lld pair_verify_calls=%lld\n"
          "RATE j2_leaves_per_core_sec=%.4g bloom_queries_per_core_sec=%.4g\n"
          "RATE bloom_inserts=%.6g bloom_inserts_per_core_sec=%.4g\n"
          "SPLIT build_wall=%.3f traversal_wall_est=%.3f build_frac=%.4f\n",
          nth, p6_w, build_w, build_core, loop_w, loop_core,
          enum_core, dfs_core, j2, bq, bp, pv,
          dfs_core>0? j2/dfs_core : 0.0, dfs_core>0? bq/dfs_core : 0.0,
          (double)((double)Bmax*(Bmax+1)/2),
          build_core>0? ((double)Bmax*(Bmax+1)/2)/build_core : 0.0,
          build_w+p6_w, loop_w,
          (build_w+p6_w)/((build_w+p6_w+loop_w)>0?(build_w+p6_w+loop_w):1));
    }
    return 0;
}
