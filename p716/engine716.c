/* engine716.c -- search for (7,1,6): f^7 = a^7+b^7+c^7+d^7+e^7+g^7
 *
 * Algorithm (SPEC_716.md): 4+3 meet in the middle with nested "largest part"
 * windows.  The four enumerated parts are f (outer), a, b, c; the three table
 * parts are d >= e >= g.  A blocked Bloom filter holds every 3-sum
 * d^7+e^7+g^7 with TB >= d >= e >= g >= 1; every Bloom positive is settled by
 * an exact verifier in unsigned __int128.
 *
 * Ordering / completeness: a >= b >= c >= d >= e >= g, so every unordered
 * solution is enumerated exactly once.  Equal bases are allowed.  The
 * verifier caps d at c, which is the only place where the boundary c >= d is
 * enforced exactly (the window only guarantees R3 <= 3 c^7).
 *
 * Bucket passes (--nb NB): bucket(v) is a deterministic function of the full
 * 128-bit value v.  Pass k inserts only 3-sums with bucket k and evaluates
 * only leaves whose R3 has bucket k.  Each leaf is therefore evaluated in
 * exactly one pass, and in that pass the generating triple (if any) is in the
 * table, because the triple's sum is exactly R3.  Table memory / NB.
 *
 * Arithmetic: everything is unsigned __int128.  f^7 < 2^128 for f < 318000.
 *
 * Bloom / hashing / iroot machinery follows strata/k7engine.c.
 */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <sys/mman.h>
#ifdef _OPENMP
#include <omp.h>
#endif

typedef unsigned __int128 u128;
typedef uint64_t u64;
typedef uint32_t u32;
typedef long long i64;

#define DIE(...) do{ fprintf(stderr,"FATAL " __VA_ARGS__); fputc('\n',stderr); exit(1);}while(0)
#define ASSERT(c,...) do{ if(!(c)){ fprintf(stderr,"ASSERT FAILED %s:%d: %s : ",__FILE__,__LINE__,#c); fprintf(stderr,__VA_ARGS__); fputc('\n',stderr); exit(1);} }while(0)

/* ------------------------------------------------------------- basics */
static u64 IR7MAX;                 /* largest r with r^7 < 2^128 */
static inline u128 ipow7(u64 x){ u128 a=x, s=a*a; s=s*s; return s*(a*a)*a; } /* x^4*x^2*x */

static u64 iroot7(u128 R){
    if(R==0) return 0;
    long double d=(long double)R;
    long double e=powl(d,1.0L/7.0L);
    if(!(e>=0)) e=0;
    if(e>(long double)IR7MAX) e=(long double)IR7MAX;
    u64 r=(u64)e;
    if(r>IR7MAX) r=IR7MAX;
    while(r>0 && ipow7(r)>R) r--;
    while(r<IR7MAX && ipow7(r+1)<=R) r++;
    return r;
}
static u128 KLIM[8];               /* KLIM[k] = floor((2^128-1)/k) */
static inline int kge(u128 p,u64 k,u128 X){ if(p>KLIM[k]) return 1; return p*(u128)k>=X; }
/* smallest n >= 1 with k*n^7 >= X */
static u64 ceil_root7_div(u128 X,u64 k){
    if(X<=(u128)k) return 1;
    long double e=powl((long double)X/(long double)k,1.0L/7.0L);
    if(!(e>=0)) e=0;
    if(e>(long double)IR7MAX) e=(long double)IR7MAX;
    u64 n=(u64)e; if(n<1) n=1; if(n>IR7MAX) n=IR7MAX;
    while(n>1 && kge(ipow7(n-1),k,X)) n--;
    while(n<IR7MAX && !kge(ipow7(n),k,X)) n++;
    return n;
}
static u64 pow7m(u64 x,u64 m){ u64 r=x%m; u64 s=(u64)((u128)r*r%m); s=(u64)((u128)s*s%m); return (u64)((u128)s*((u128)r*r%m)%m*r%m); }

static void u128_to_str(u128 v,char *buf){ char t[48]; int n=0;
    if(!v){ buf[0]='0'; buf[1]=0; return; }
    while(v){ t[n++]=(char)('0'+(int)(v%10)); v/=10; }
    for(int i=0;i<n;i++) buf[i]=t[n-1-i]; buf[n]=0; }
static u128 str_to_u128(const char *s){ u128 v=0; for(;*s;s++){ if(*s<'0'||*s>'9') DIE("bad u128 literal '%s'",s); v=v*10+(u64)(*s-'0'); } return v; }

/* --------------------------------------------------- residue masks 49,29,43 */
/* sumres[i][j] = bitset of residues attainable by a sum of j seventh powers */
#define NMOD 3
static const u32 MODS[NMOD]={49,29,43};
static u64 R64MOD[NMOD];
static u64 sumres[NMOD][4];
static double passrate[NMOD][4];
static void init_masks(void){
    for(int i=0;i<NMOD;i++){
        u32 md=MODS[i];
        R64MOD[i]=(u64)(((u128)1<<64)%md);
        u64 single=0;
        for(u32 x=0;x<md;x++) single|=1ULL<<pow7m(x,md);
        sumres[i][0]=1;
        for(int j=1;j<4;j++){
            u64 prev=sumres[i][j-1],cur=0;
            for(u32 r=0;r<md;r++) if(prev>>r&1)
                for(u32 s=0;s<md;s++) if(single>>s&1) cur|=1ULL<<((r+s)%md);
            sumres[i][j]=cur;
        }
        for(int j=0;j<4;j++) passrate[i][j]=(double)__builtin_popcountll(sumres[i][j])/(double)md;
    }
}
#define MODFN(NAME,MD,IDX) static inline u32 NAME(u64 hi,u64 lo){ \
    return (u32)((((hi%(MD))*R64MOD[IDX]) + (lo%(MD))) % (MD)); }
MODFN(md49_,49,0) MODFN(md29_,29,1) MODFN(md43_,43,2)
static inline int mask3_ok(u128 v){
    u64 hi=(u64)(v>>64), lo=(u64)v;
    if(!((sumres[0][3]>>md49_(hi,lo))&1)) return 0;
    if(!((sumres[1][3]>>md29_(hi,lo))&1)) return 0;
    if(!((sumres[2][3]>>md43_(hi,lo))&1)) return 0;
    return 1;
}
static inline int mask2_ok(u128 v){
    u64 hi=(u64)(v>>64), lo=(u64)v;
    if(!((sumres[0][2]>>md49_(hi,lo))&1)) return 0;
    if(!((sumres[1][2]>>md29_(hi,lo))&1)) return 0;
    if(!((sumres[2][2]>>md43_(hi,lo))&1)) return 0;
    return 1;
}

/* -------------------------------------------------- blocked Bloom filters */
static inline u64 mixA(u64 h){ h^=h>>30; h*=0xbf58476d1ce4e5b9ULL; h^=h>>27; h*=0x94d049bb133111ebULL; h^=h>>31; return h; }
static inline u64 mixB(u64 h){ h^=h>>33; h*=0xff51afd7ed558ccdULL; h^=h>>33; h*=0xc4ceb9fe1a85ec53ULL; h^=h>>33; return h; }
static inline u64 mixC(u64 h){ h^=h>>29; h*=0xd6e8feb86659fd93ULL; h^=h>>32; h*=0xa5cb3b1f2f5a1d7bULL; h^=h>>29; return h; }
static inline u64 mixD(u64 h){ h^=h>>31; h*=0x9e3779b97f4a7c15ULL; h^=h>>30; h*=0xc2b2ae3d27d4eb4fULL; h^=h>>32; return h; }

typedef struct { u64 *w; u64 nlines; int k; u64 seed; } Filt;

/* a,b,d give 4+4+4 = 12 independent 9-bit probe slices; c selects the line
 * (from its HIGH bits, which are not used for probes). */
static inline void filt_h(const Filt *F,u128 s,u64 *a,u64 *b,u64 *c,u64 *d){
    u64 lo=(u64)s ^ F->seed, hi=(u64)(s>>64);
    *a=mixA(lo ^ (0x9e3779b97f4a7c15ULL*(hi+1)));
    *b=mixB(hi ^ (0xbf58476d1ce4e5b9ULL*(lo|1)));
    *c=mixC((lo+0x165667b19e3779f9ULL) ^ __builtin_bswap64(hi));
    *d=mixD((hi+0x27d4eb2f165667c5ULL) ^ (lo*0xff51afd7ed558ccdULL));
}
static inline void filt_probes(u64 a,u64 b,u64 d,int k,u32 *p){
    u64 src[3]={a,b,d};
    for(int i=0;i<k;i++) p[i]=(u32)((src[i>>2]>>(9*(i&3)))&511);
}
static inline void filt_add(const Filt *F,u128 s){
    u64 a,b,c,d; filt_h(F,s,&a,&b,&c,&d);
    u64 *line=F->w+((u64)(((u128)c*(u128)F->nlines)>>64)<<3);
    u32 p[12]; filt_probes(a,b,d,F->k,p);
    for(int i=0;i<F->k;i++) __atomic_or_fetch(&line[p[i]>>6],1ULL<<(p[i]&63),__ATOMIC_RELAXED);
}
static inline int filt_query(const Filt *F,u128 s){
    u64 a,b,c,d; filt_h(F,s,&a,&b,&c,&d);
    const u64 *line=F->w+((u64)(((u128)c*(u128)F->nlines)>>64)<<3);
    u32 p[12]; filt_probes(a,b,d,F->k,p);
    for(int i=0;i<F->k;i++) if(!((line[p[i]>>6]>>(p[i]&63))&1)) return 0;
    return 1;
}
static Filt F3, F2;

/* big zero-filled allocation, backed by transparent huge pages when possible:
 * the 3-sum filter is many GB and every query is a random access, so 4 KB
 * pages cost a TLB miss on top of the cache miss. */
static void *alloc_big(size_t bytes){
    void *p=mmap(NULL,bytes,PROT_READ|PROT_WRITE,MAP_PRIVATE|MAP_ANONYMOUS|MAP_NORESERVE,-1,0);
    if(p==MAP_FAILED) return NULL;
#ifdef MADV_HUGEPAGE
    madvise(p,bytes,MADV_HUGEPAGE);
#endif
    return p;
}                 /* 3-sum table, 2-sum helper table */

/* ------------------------------------------------------------- buckets */
static inline u64 bucket_of(u128 s,u64 nb){
    if(nb<=1) return 0;
    u64 lo=(u64)s, hi=(u64)(s>>64);
    u64 h=mixB(lo*0x9e3779b97f4a7c15ULL ^ mixA(hi^0x51ed270bULL));
    return (u64)(((u128)h*(u128)nb)>>64);
}

/* ------------------------------------------------------------ counters */
typedef struct { i64 leaves,masked,bqueries,positives,pairq,verified,sols; char pad[128-7*8]; } Ctr;
static Ctr *ctrs=NULL; static int nthreads=1;
static inline int tid(void){
#ifdef _OPENMP
    return omp_get_thread_num();
#else
    return 0;
#endif
}

/* --------------------------------------------------------- power table */
static u128 *P7=NULL; static u64 TB=0;   /* P7[x]=x^7, x <= TB */

/* ------------------------------------------------------ exact verifier */
/* Is R = d^7+e^7+g^7 with cmax >= d >= e >= g >= 1 ?  Exact, u128 only. */
static int verify3(u128 R,u64 cmax,Ctr *C,u64 *od,u64 *oe,u64 *og){
    if(R<3) return 0;
    u64 dhi=iroot7(R); if(dhi>cmax) dhi=cmax; if(dhi>TB) dhi=TB;
    u64 dlo=ceil_root7_div(R,3);
    for(u64 d=dhi;d>=dlo && d>=1;d--){
        u128 Rp=R-P7[d];
        if(Rp<2) continue;
        if(!mask2_ok(Rp)) continue;
        C->pairq++;
        if(!filt_query(&F2,Rp)) continue;
        C->verified++;
        u64 ehi=iroot7(Rp); if(ehi>d) ehi=d;
        u64 elo=ceil_root7_div(Rp,2);
        for(u64 e=ehi;e>=elo && e>=1;e--){
            u128 rest=Rp-P7[e];
            if(rest<1) continue;
            u64 g=iroot7(rest);
            if(g>=1 && g<=e && ipow7(g)==rest){ *od=d;*oe=e;*og=g; return 1; }
        }
    }
    return 0;
}

/* -------------------------------------------------------- solution slot */
static volatile int sol_found=0;
static u64 sol_f,sol_a,sol_b,sol_c,sol_d,sol_e,sol_g;

/* ------------------------------------------------- the 4+3 enumeration */
/* F7 = f^7 (or the planted S); amax = f-1 (or floor(S^(1/7))). */
static void enumerate(u128 F7,u64 fval,u64 amax,u64 nb,u64 pass,Ctr *C){
    if(F7<6) return;
    u64 amin=ceil_root7_div(F7,6);
    {   u64 r=iroot7(F7); if(amax>r) amax=r; }
    if(amax>TB) amax=TB;
    for(u64 a=amax;a>=amin && a>=1;a--){
        u128 R1=F7-P7[a];
        if(R1<5) continue;
        u64 bhi=iroot7(R1); if(bhi>a) bhi=a;
        u64 blo=ceil_root7_div(R1,5);
        for(u64 b=bhi;b>=blo && b>=1;b--){
            u128 R2=R1-P7[b];
            if(R2<4) continue;
            u64 chi=iroot7(R2); if(chi>b) chi=b;
            u64 clo=ceil_root7_div(R2,4);
            for(u64 c=chi;c>=clo && c>=1;c--){
                u128 R3=R2-P7[c];
                if(R3<3) continue;
                if(R3>3*P7[c]) continue;          /* d <= c must be possible */
                if(bucket_of(R3,nb)!=pass) continue;
                C->leaves++;
                if(!mask3_ok(R3)) continue;
                C->masked++;
                C->bqueries++;
                if(!filt_query(&F3,R3)) continue;
                C->positives++;
                u64 d,e,g;
                if(verify3(R3,c,C,&d,&e,&g)){
                    C->sols++;
                    if(__sync_bool_compare_and_swap(&sol_found,0,1)){
                        sol_f=fval; sol_a=a; sol_b=b; sol_c=c; sol_d=d; sol_e=e; sol_g=g;
                    }
                    return;
                }
            }
        }
    }
}

/* --------------------------------------------------------------- table */
static void build_table(u64 nb,u64 pass,u64 *n_inserted){
    i64 cnt=0;
#pragma omp parallel for schedule(dynamic,8) reduction(+:cnt)
    for(u64 d=1;d<=TB;d++){
        u128 pd=P7[d];
        for(u64 e=1;e<=d;e++){
            u128 pde=pd+P7[e];
            for(u64 g=1;g<=e;g++){
                u128 s=pde+P7[g];
                if(bucket_of(s,nb)!=pass) continue;
                filt_add(&F3,s); cnt++;
            }
        }
    }
    *n_inserted=(u64)cnt;
}
static void build_pairs(void){
#pragma omp parallel for schedule(dynamic,8)
    for(u64 e=1;e<=TB;e++){
        u128 pe=P7[e];
        for(u64 g=1;g<=e;g++) filt_add(&F2,pe+P7[g]);
    }
}

static double now(void){ struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t); return t.tv_sec+1e-9*t.tv_nsec; }
static void ctr_total(Ctr *T){ memset(T,0,sizeof(*T));
    for(int i=0;i<nthreads;i++){ T->leaves+=ctrs[i].leaves; T->masked+=ctrs[i].masked;
        T->bqueries+=ctrs[i].bqueries; T->positives+=ctrs[i].positives; T->pairq+=ctrs[i].pairq;
        T->verified+=ctrs[i].verified; T->sols+=ctrs[i].sols; } }

static u64 FMIN,FMAX;
static double T0;
static void band_line(u64 lo,u64 hi){
    Ctr T; ctr_total(&T);
    printf("BAND %llu %llu %lld %lld %lld %lld %lld %lld %.1f\n",
        (unsigned long long)lo,(unsigned long long)hi,(long long)T.leaves,(long long)T.masked,
        (long long)T.bqueries,(long long)T.positives,(long long)T.verified,(long long)T.sols,now()-T0);
    fflush(stdout);
}
static void emit_solution(void){
    printf("SOLUTION %llu %llu %llu %llu %llu %llu %llu\n",
        (unsigned long long)sol_f,(unsigned long long)sol_a,(unsigned long long)sol_b,
        (unsigned long long)sol_c,(unsigned long long)sol_d,(unsigned long long)sol_e,
        (unsigned long long)sol_g);
    fflush(stdout);
    exit(3);
}

static void selftest(void){
    for(u64 x=1;x<=20000;x+=7){
        ASSERT(iroot7(ipow7(x))==x,"iroot7 exact %llu",(unsigned long long)x);
        ASSERT(iroot7(ipow7(x)-1)==x-1,"iroot7 -1 %llu",(unsigned long long)x);
        ASSERT(iroot7(ipow7(x)+1)==x,"iroot7 +1 %llu",(unsigned long long)x);
    }
    for(u64 k=2;k<=6;k++) for(u64 x=1;x<=3000;x+=13){
        u128 X=ipow7(x)*(u128)k;
        ASSERT(ceil_root7_div(X,k)==x,"ceil_root7_div exact k=%llu x=%llu",(unsigned long long)k,(unsigned long long)x);
        ASSERT(ceil_root7_div(X-1,k)==x,"ceil_root7_div -1 k=%llu x=%llu",(unsigned long long)k,(unsigned long long)x);
        ASSERT(ceil_root7_div(X+1,k)==x+1,"ceil_root7_div +1 k=%llu x=%llu",(unsigned long long)k,(unsigned long long)x);
    }
    { u128 v=ipow7(123456); ASSERT(iroot7(v)==123456,"big iroot7"); }
}

int main(int argc,char **argv){
    if(argc<3){ fprintf(stderr,"usage: %s FMIN FMAX [--nb NB] [--plant S] [--bpp N] [--query3 V] [--bloomtest N] [--bloomstat N] [--threads T]\n",argv[0]); return 1; }
    FMIN=strtoull(argv[1],0,10); FMAX=strtoull(argv[2],0,10);
    u64 nb=1,bpp=16,bloomtest=0,bloomstat=0; int plant=0,do_query=0; u128 plantS=0,queryV=0; int thr=0;
    for(int i=3;i<argc;i++){
        if(!strcmp(argv[i],"--nb")&&i+1<argc) nb=strtoull(argv[++i],0,10);
        else if(!strcmp(argv[i],"--bpp")&&i+1<argc) bpp=strtoull(argv[++i],0,10);
        else if(!strcmp(argv[i],"--plant")&&i+1<argc){ plant=1; plantS=str_to_u128(argv[++i]); }
        else if(!strcmp(argv[i],"--query3")&&i+1<argc){ do_query=1; queryV=str_to_u128(argv[++i]); }
        else if(!strcmp(argv[i],"--bloomtest")&&i+1<argc) bloomtest=strtoull(argv[++i],0,10);
        else if(!strcmp(argv[i],"--bloomstat")&&i+1<argc) bloomstat=strtoull(argv[++i],0,10);
        else if(!strcmp(argv[i],"--threads")&&i+1<argc) thr=atoi(argv[++i]);
        else DIE("unknown argument %s",argv[i]);
    }
    if(FMAX<=FMIN&&!plant&&!do_query&&!bloomtest&&!bloomstat) DIE("need FMAX > FMIN");
    if(FMAX>200000) DIE("FMAX too large for 128-bit headroom");
    if(nb<1) nb=1;
#ifdef _OPENMP
    if(thr>0) omp_set_num_threads(thr);
#pragma omp parallel
    { if(omp_get_thread_num()==0) nthreads=omp_get_num_threads(); }
#endif
    T0=now();
    { u64 lo=1,hi=1ULL<<20; while(lo<hi){ u64 m=(lo+hi+1)/2; u128 a=m,s=a*a; s=s*s;
        long double est=powl((long double)m,7.0L);
        if(est<3.4e38L) lo=m; else hi=m-1; (void)s; } IR7MAX=lo; }
    while(IR7MAX>1){ long double e=powl((long double)IR7MAX,7.0L); if(e<3.402e38L) break; IR7MAX--; }
    for(u64 k=1;k<8;k++) KLIM[k]=(~(u128)0)/k;
    init_masks();
    selftest();
    fprintf(stderr,"IR7MAX=%llu threads=%d\n",(unsigned long long)IR7MAX,nthreads);
    for(int i=0;i<NMOD;i++)
        fprintf(stderr,"mask mod %u: 1-sum %.4f 2-sum %.4f 3-sum %.4f\n",MODS[i],passrate[i][1],passrate[i][2],passrate[i][3]);
    fprintf(stderr,"combined 3-sum mask pass rate %.4f, 2-sum %.4f\n",
        passrate[0][3]*passrate[1][3]*passrate[2][3], passrate[0][2]*passrate[1][2]*passrate[2][2]);

    TB=FMAX;
    if(plant){ u64 r=iroot7(plantS); if(r>TB) TB=r; }
    if(do_query){ /* table bound stays FMAX */ }
    if(TB<1) TB=1;
    P7=malloc(sizeof(u128)*(TB+2)); if(!P7) DIE("alloc P7");
    for(u64 x=0;x<=TB+1;x++) P7[x]=ipow7(x);

    /* filters */
    u128 ntri=(u128)TB*(TB+1)*(TB+2)/6;
    u64 per=(u64)(ntri/nb)+1;
    F3.k=(int)(bpp*0.693+0.5); if(F3.k<4) F3.k=4; if(F3.k>12) F3.k=12; F3.seed=0x716a5ULL;
    F3.nlines=(u64)(((u128)per*bpp*103/100+511)/512); if(F3.nlines<1) F3.nlines=1;
    F2.k=8; F2.seed=0xb2b2b2ULL;
    { u128 npair=(u128)TB*(TB+1)/2; F2.nlines=(u64)(((u128)npair*32+511)/512); if(F2.nlines<1) F2.nlines=1; }
    fprintf(stderr,"table bound TB=%llu triples=%.4g nb=%llu per-pass=%.4g bpp=%llu probes=%d bloom=%.3f GB; pair bloom %.3f GB\n",
        (unsigned long long)TB,(double)ntri,(unsigned long long)nb,(double)per,(unsigned long long)bpp,F3.k,
        F3.nlines*64.0/1e9,F2.nlines*64.0/1e9);
    F3.w=alloc_big((size_t)F3.nlines*64); if(!F3.w) DIE("bloom alloc failed (%.2f GB)",F3.nlines*64.0/1e9);
    F2.w=alloc_big((size_t)F2.nlines*64); if(!F2.w) DIE("pair bloom alloc failed");
    build_pairs();
    { u64 x=TB>266?266:TB, y=TB>258?258:1; ASSERT(filt_query(&F2,P7[x]+P7[y]),"pair bloom self-test"); }
    ASSERT(filt_query(&F2,P7[1]+P7[1]),"pair bloom self-test 1+1");
    fprintf(stderr,"pair bloom built (%.1f s)\n",now()-T0);

    ctrs=calloc(nthreads?nthreads:1,sizeof(Ctr)); if(!ctrs) DIE("alloc ctrs");

    for(u64 pass=0;pass<nb;pass++){
        if(pass>0) memset(F3.w,0,(size_t)F3.nlines*64);
        u64 nins=0; double tb0=now();
        build_table(nb,pass,&nins);
        fprintf(stderr,"pass %llu/%llu: inserted %llu 3-sums in %.1f s (load %.3f bits/entry)\n",
            (unsigned long long)pass+1,(unsigned long long)nb,(unsigned long long)nins,now()-tb0,
            nins?F3.nlines*512.0/nins:0.0);
        /* self-test: a few inserted values must query positive */
        for(u64 d=1;d<=TB && d<=97;d+=17) for(u64 e=1;e<=d;e+=13) for(u64 g=1;g<=e;g+=11){
            u128 s=P7[d]+P7[e]+P7[g];
            if(bucket_of(s,nb)==pass) ASSERT(filt_query(&F3,s),"BLOOM FALSE NEGATIVE %llu %llu %llu",
                (unsigned long long)d,(unsigned long long)e,(unsigned long long)g);
        }

        if(bloomtest){
            u64 fn=0; unsigned seed=12345;
            for(u64 i=0;i<bloomtest;i++){
                u64 d=1+(u64)(rand_r(&seed)%(unsigned)TB);
                u64 e=1+(u64)(rand_r(&seed)%(unsigned)d);
                u64 g=1+(u64)(rand_r(&seed)%(unsigned)e);
                u128 s=P7[d]+P7[e]+P7[g];
                if(bucket_of(s,nb)!=pass) continue;
                if(!filt_query(&F3,s)) fn++;
            }
            printf("BLOOMTEST n=%llu false_negatives=%llu\n",(unsigned long long)bloomtest,(unsigned long long)fn);
            fflush(stdout);
        }
        if(bloomstat){
            u64 fp=0,q=0; unsigned seed=999; u128 base=P7[TB]*3;
            for(u64 i=0;i<bloomstat;i++){
                u128 v=base+(u128)rand_r(&seed)*1000003u+i*7919u;  /* far above any table value */
                if(bucket_of(v,nb)!=pass) continue;
                q++; if(filt_query(&F3,v)) fp++;
            }
            printf("BLOOMSTAT queries=%llu positives=%llu fp_rate=%.3g\n",
                (unsigned long long)q,(unsigned long long)fp,q?(double)fp/q:0.0);
            fflush(stdout);
        }
        if(do_query){
            u64 d=0,e=0,g=0; Ctr C; memset(&C,0,sizeof C);
            int bq=filt_query(&F3,queryV);
            int ex=verify3(queryV,TB,&C,&d,&e,&g);
            char buf[48]; u128_to_str(queryV,buf);
            printf("QUERY3 %s bucket=%llu pass=%llu bloom=%d exact=%d %llu %llu %llu\n",buf,
                (unsigned long long)bucket_of(queryV,nb),(unsigned long long)pass,bq,ex,
                (unsigned long long)d,(unsigned long long)e,(unsigned long long)g);
            fflush(stdout);
        }
        if(plant){
            Ctr *C=&ctrs[0];
            enumerate(plantS,0,TB,nb,pass,C);
            if(sol_found) emit_solution();
            continue;
        }
        if(bloomtest||bloomstat||do_query) continue;

        for(u64 lo=FMIN+1;lo<=FMAX;lo+=100){
            u64 hi=lo+99; if(hi>FMAX) hi=FMAX;
#pragma omp parallel for schedule(dynamic,1)
            for(u64 f=lo;f<=hi;f++){
                if(sol_found) continue;
                enumerate(P7[f],f,f-1,nb,pass,&ctrs[tid()]);
            }
            band_line(FMIN,hi);
            if(sol_found) emit_solution();
        }
    }
    if(!plant&&!bloomtest&&!bloomstat&&!do_query) band_line(FMIN,FMAX);
    else { band_line(FMIN,FMAX); }
    if(sol_found) emit_solution();
    return 0;
}
