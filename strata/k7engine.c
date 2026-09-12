/*
 * k7engine.c -- stratified class-1 engine for the k7 in {0,1,2} stratum of
 *               a^6 + b^6 + c^6 + d^6 + e^6 = f^6.
 *
 * Implements strata/SPEC_K7.md exactly.  See strata/K7_LOG.md for runs.
 *
 * SETTING.  Class-1 (Meyrignac) primitive solution: one term t carries all three
 * exemptions (odd, coprime to 3, coprime to 7); the other four are 42*b_i.  Then
 *      f^6 - t^6 = 42^6 * m,   m = b1^6+b2^6+b3^6+b4^6,   1 <= b_i <= B=(f-1)/42,
 * and t = zeta*f mod 42^6 for one of the 144 sixth roots of unity zeta mod 42^6,
 * with 0 < t < f, gcd(f,42)=1.  Exact counting budgets:
 *      k2 = m mod 8 = #{i: b_i odd},  k3 = m mod 9 = #{i: 3 !| b_i},
 *      k7 = m mod 7 = #{i: 7 !| b_i}.        Each must be <= 4.
 *
 * BAND CONVENTION.  This engine processes every f with FMIN < f <= FMAX, i.e. the
 * half-open interval (FMIN, FMAX].  src/caseA2.c uses the CLOSED interval
 * [fmin,fmax].  The two agree whenever gcd(FMIN,42) > 1 (then f=FMIN contributes
 * no candidate anyway), which holds for every band endpoint used in this project
 * (2, 700000, 730000, 1000000, 4300000, 4400000, ...).  The half-open convention
 * is the correct one for chaining bands without double-counting the joint.
 *
 * STRATA.  Let R7 = 7^6 = 117649, B' = B/7.
 *   k7=0 : all four b_i divisible by 7.  Requires R7 | m (valuation lemma);
 *          m' = m/R7 is solved as a generic 4-sum with bases <= B'.
 *   k7=1 : exactly one base b4 coprime to 7.  m mod R7 must be a unit sixth-power
 *          residue; b4 runs over r + j*R7 for r in ROOTS[m mod R7].  The residue
 *          R' = (m-b4^6)/R7 is solved as a 3-sum with bases <= B'.
 *   k7=2 : exactly two bases b3 <= b4 coprime to 7.  b3 is enumerated, b4 comes
 *          from ROOTS[(m-b3^6) mod R7]; R' = (m-b3^6-b4^6)/R7 is a 2-sum.
 *   k7 in {3,4} : NOT COVERED, counted and reported.
 *   k7 in {5,6} : impossible (budget), eliminated.
 *
 * ENGINEERING (lessons from strata/BASELINE_LOG.md):
 *   * The 2-sum residue masks (mod 64,27,49,13,19,31,37,43) are applied BEFORE
 *     any Bloom query, so the filter is only touched by ~2% of leaves.
 *   * Blocked Bloom, 512-bit lines, 16 bits per pair, 8 bit positions taken as
 *     8 INDEPENDENT 9-bit slices of two independent 64-bit mixes (the old
 *     family used 8-bit strides of a 9-bit mask -> overlapping, correlated
 *     slices, measured 1.67x the ideal false-positive rate).  A third mix picks
 *     the line, so line choice and bit choice are independent.
 *   * OpenMP schedule(dynamic,1) over f (schedule(dynamic,4096) serialised any
 *     band narrower than 4096).
 *   * Exact division by R7 in the hot path is done by multiplying with the
 *     2-adic inverse of R7 mod 2^128 (R7 is odd), never by __udivti3.
 *
 * ARITHMETIC.  Everything is unsigned __int128.  f^6 does NOT fit in 128 bits for
 * f > 2.64e6, so m = (f^6-t^6)/42^6 is computed from the cyclotomic factorisation
 * (f-t)(f+t)(f^2-ft+t^2)(f^2+ft+t^2) with 42^6 stripped by gcd, exactly as
 * src/caseA2.c:282-294 does.  Every candidate is additionally checked with
 * m*42^6 + t^6 == f^6 modulo the prime 2^61-1.
 *
 * Build: gcc -O3 -march=native -fopenmp -std=c11 -Wall -o k7engine k7engine.c -lm
 * Run:   k7engine FMIN FMAX [--plant M B] [--candidates FILE] [--bpp N]
 */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif

typedef unsigned __int128 u128;
typedef uint64_t u64;
typedef uint32_t u32;
typedef long long i64;

#define DIE(...) do{ fprintf(stderr,"FATAL " __VA_ARGS__); fputc('\n',stderr); exit(1);}while(0)
#define ASSERT(c,...) do{ if(!(c)){ fprintf(stderr,"ASSERT FAILED %s:%d: %s : ",__FILE__,__LINE__,#c); fprintf(stderr,__VA_ARGS__); fputc('\n',stderr); exit(1);} }while(0)

static const u64 M42 = 5489031744ULL;   /* 42^6 */
static const u64 R7   = 117649ULL;      /* 7^6  */

/* ------------------------------------------------------------------ basics */
static u128 IR7INV;          /* inverse of R7 mod 2^128 (R7 odd) -> exact division */
static u64  IR6MAX = 2642245;/* largest r with r^6 < 2^128 */

static inline u128 ipow6(u64 x){ u128 s=(u128)x*x; return s*s*s; }

static int pow6_ovf(u64 x,u128 *out){
    u128 a=x,s,t;
    if(__builtin_mul_overflow(a,a,&s)) return 1;
    if(__builtin_mul_overflow(s,s,&t)) return 1;
    if(__builtin_mul_overflow(t,s,&t)) return 1;
    *out=t; return 0;
}

static u64 iroot6(u128 R){
    if(R==0) return 0;
    long double d=(long double)R;
    long double e=powl(d,1.0L/6.0L);
    if(!(e>=0)) e=0;
    if(e> (long double)IR6MAX) e=(long double)IR6MAX;
    u64 r=(u64)e;
    if(r>IR6MAX) r=IR6MAX;
    while(r>0 && ipow6(r)>R) r--;
    while(r<IR6MAX && ipow6(r+1)<=R) r++;
    return r;
}

static u64 pow6m(u64 x,u64 m){ u64 r=x%m; u128 s=(u128)r*r%m; s=s*s%m; return (u64)((u128)s*((u128)r*r%m)%m); }

/* 2-adic inverse of an odd 128-bit value, by Newton iteration */
static u128 inv_mod_2_128(u128 a){
    u128 x=a;                 /* a*a == 1 (mod 8) for odd a */
    for(int i=0;i<6;i++) x = x*(2-a*x);
    return x;
}

/* ------------------------------------------------------- u128 printing/parsing */
static void u128_to_str(u128 v,char *buf){   /* buf >= 41 bytes */
    char tmp[48]; int n=0;
    if(v==0){ buf[0]='0'; buf[1]=0; return; }
    while(v){ tmp[n++]=(char)('0'+(int)(v%10)); v/=10; }
    for(int i=0;i<n;i++) buf[i]=tmp[n-1-i];
    buf[n]=0;
}
static u128 str_to_u128(const char *s){
    u128 v=0;
    for(;*s;s++){ if(*s<'0'||*s>'9') DIE("bad u128 literal"); v=v*10+(u64)(*s-'0'); }
    return v;
}

/* -------------------------------------------------- 144 sixth roots mod 42^6 */
static u64 roots144[160]; static int nroots=0;
static void build_roots(void){
    u64 r2[8]; int n2=0; for(u64 x=1;x<64;x+=2) if(pow6m(x,64)==1) r2[n2++]=x;
    u64 r3[8]; int n3=0; for(u64 x=1;x<729;x++) if(x%3&&pow6m(x,729)==1) r3[n3++]=x;
    u64 r7[8]; int n7=0; for(u64 x=1;x<117649;x++) if(x%7&&pow6m(x,117649)==1) r7[n7++]=x;
    ASSERT(n2==4&&n3==6&&n7==6,"root factor counts %d %d %d",n2,n3,n7);
    for(int i=0;i<n2;i++)for(int j=0;j<n3;j++)for(int k=0;k<n7;k++){
        u64 m12=64*729ULL,x12=0,x=0;
        for(u64 t=r2[i];;t+=64) if(t%729==r3[j]){x12=t;break;}
        for(u64 t=x12;;t+=m12) if(t%117649==r7[k]){x=t;break;}
        roots144[nroots++]=x;
    }
    ASSERT(nroots==144,"nroots=%d",nroots);
    for(int i=0;i<nroots;i++) ASSERT(pow6m(roots144[i],M42)==1,"bad root %llu",(unsigned long long)roots144[i]);
}

/* -------------------------------------------- sixth roots mod 7^6 (ROOTS[u]) */
static u32 *P6R7;        /* x^6 mod R7 for x in [0,R7) */
static uint8_t *ROOTCNT; /* ROOTCNT[u] in {0,6} */
static u32 *ROOTLIST;    /* 6 slots per u, ascending */
static void build_r7_roots(void){
    P6R7   = malloc(sizeof(u32)*R7);
    ROOTCNT= calloc(R7,1);
    ROOTLIST=malloc(sizeof(u32)*R7*6);
    if(!P6R7||!ROOTCNT||!ROOTLIST) DIE("alloc r7 roots");
    for(u64 x=0;x<R7;x++) P6R7[x]=(u32)pow6m(x,R7);
    for(u64 x=1;x<R7;x++){
        if(x%7==0) continue;
        u32 u=P6R7[x];
        ASSERT(ROOTCNT[u]<6,"more than 6 sixth roots mod 7^6 for u=%u",u);
        ROOTLIST[(size_t)u*6+ROOTCNT[u]++]=(u32)x;
    }
    long nonempty=0; for(u64 u=0;u<R7;u++){ if(ROOTCNT[u]){ ASSERT(ROOTCNT[u]==6,"root count %d",ROOTCNT[u]); nonempty++; } }
    ASSERT(nonempty==(long)(R7/7),"unit sixth-power residues mod 7^6 = %ld (want %ld)",nonempty,(long)(R7/7));
}

/* --------------------------------------------------------- residue masks */
/* sumres[i][j] = bitset of achievable residues of a sum of j sixth powers mod MODS[i] */
#define NMOD 8
static const u32 MODS[NMOD]={64,27,49,13,19,31,37,43};
static u64 R64MOD[NMOD];
static u64 sumres[NMOD][5];
static u32 pw6tab[NMOD][64];

static void init_masks(void){
    for(int i=0;i<NMOD;i++){
        u32 md=MODS[i];
        R64MOD[i]=(u64)(((u128)1<<64)%md);
        u64 single=0;
        for(u32 x=0;x<md;x++){ u32 r=(u32)pow6m(x,md); pw6tab[i][x]=r; single|=1ULL<<r; }
        sumres[i][0]=1;
        for(int j=1;j<5;j++){
            u64 prev=sumres[i][j-1],cur=0;
            for(u32 r=0;r<md;r++) if(prev>>r&1)
                for(u32 s=0;s<md;s++) if(single>>s&1)
                    cur|=1ULL<<((r+s)%md);
            sumres[i][j]=cur;
        }
    }
}

/* literal-modulus reductions of a 128-bit value: hi*(2^64 mod md) + lo, mod md */
#define MODFN(NAME,MD,IDX) static inline u32 NAME(u64 hi,u64 lo){ \
    return (u32)((((hi%(MD))*R64MOD[IDX]) + (lo%(MD))) % (MD)); }
MODFN(md27_,27,1) MODFN(md49_,49,2) MODFN(md13_,13,3)
MODFN(md19_,19,4) MODFN(md31_,31,5) MODFN(md37_,37,6) MODFN(md43_,43,7)
static inline u32 md9_(u64 hi,u64 lo){ return (u32)((((hi%9)*7u)+(lo%9))%9); }   /* 2^64 = 7 mod 9 */
static inline u32 md7_(u64 hi,u64 lo){ return (u32)((((hi%7)*2u)+(lo%7))%7); }   /* 2^64 = 2 mod 7 */

/* masks ordered most-selective-first for j=2 */
static inline int mask_ok(u128 v,int j){
    u64 hi=(u64)(v>>64), lo=(u64)v;
    if(!((sumres[1][j]>>md27_(hi,lo))&1)) return 0;
    if(!((sumres[3][j]>>md13_(hi,lo))&1)) return 0;
    if(!((sumres[4][j]>>md19_(hi,lo))&1)) return 0;
    if(!((sumres[0][j]>>(lo&63))&1))      return 0;
    if(!((sumres[2][j]>>md49_(hi,lo))&1)) return 0;
    if(!((sumres[5][j]>>md31_(hi,lo))&1)) return 0;
    if(!((sumres[6][j]>>md37_(hi,lo))&1)) return 0;
    if(!((sumres[7][j]>>md43_(hi,lo))&1)) return 0;
    return 1;
}

/* ----------------------------------------------------------- power table */
static u64  Bmax_global;   /* largest original base = (FMAX-1)/42 */
static u64  Bp_global;     /* largest primed base = Bmax_global/7 */
static u128 *P6;           /* P6[x]=x^6 for x <= Bmax_global */

/* ------------------------------------------------- blocked Bloom filter */
static u64 *bloom=NULL; static u64 nlines=0;
static inline u64 mixA(u64 h){ h^=h>>30; h*=0xbf58476d1ce4e5b9ULL; h^=h>>27; h*=0x94d049bb133111ebULL; h^=h>>31; return h; }
static inline u64 mixB(u64 h){ h^=h>>33; h*=0xff51afd7ed558ccdULL; h^=h>>33; h*=0xc4ceb9fe1a85ec53ULL; h^=h>>33; return h; }
static inline u64 mixC(u64 h){ h^=h>>29; h*=0xd6e8feb86659fd93ULL; h^=h>>32; h*=0xa5cb3b1f2f5a1d7bULL; h^=h>>29; return h; }

/* Two independent 64-bit mixes give 8 independent 9-bit slices (bits 0-8, 9-17,
 * 18-26, 27-35 of each).  A third, independent mix selects the 512-bit line. */
static inline void bloom_h(u128 s,u64 *a,u64 *b,u64 *c){
    u64 lo=(u64)s, hi=(u64)(s>>64);
    *a=mixA(lo ^ (0x9e3779b97f4a7c15ULL*(hi+1)));
    *b=mixB(hi ^ (0xbf58476d1ce4e5b9ULL*(lo|1)));
    *c=mixC((lo+0x165667b19e3779f9ULL) ^ __builtin_bswap64(hi));
}
static inline u64 line_of(u64 c){ return (u64)(((u128)c*(u128)nlines)>>64); }
static inline void bloom_add(u128 s){
    u64 a,b,c; bloom_h(s,&a,&b,&c);
    u64 *line=bloom+(line_of(c)<<3);
    u32 p0=(u32)(a&511),p1=(u32)((a>>9)&511),p2=(u32)((a>>18)&511),p3=(u32)((a>>27)&511);
    u32 p4=(u32)(b&511),p5=(u32)((b>>9)&511),p6=(u32)((b>>18)&511),p7=(u32)((b>>27)&511);
    u32 ps[8]={p0,p1,p2,p3,p4,p5,p6,p7};
    for(int i=0;i<8;i++) __atomic_or_fetch(&line[ps[i]>>6],1ULL<<(ps[i]&63),__ATOMIC_RELAXED);
}
static inline int bloom_query(u128 s){
    u64 a,b,c; bloom_h(s,&a,&b,&c);
    const u64 *line=bloom+(line_of(c)<<3);
    u32 ps[8]={(u32)(a&511),(u32)((a>>9)&511),(u32)((a>>18)&511),(u32)((a>>27)&511),
               (u32)(b&511),(u32)((b>>9)&511),(u32)((b>>18)&511),(u32)((b>>27)&511)};
    for(int i=0;i<8;i++) if(!((line[ps[i]>>6]>>(ps[i]&63))&1)) return 0;
    return 1;
}

/* ------------------------------------------------------------- counters */
typedef struct {
    i64 leaves, bqueries, bpositives, exactver, dfsnodes;
    char pad[128-5*8];
} Ctr;
static Ctr *ctrs=NULL; static int nthreads=1;

/* ------------------------------------------------------ exact 2-sum search */
/* R = x^6 + y^6 with maxv >= x >= y >= 1 and exact budgets (o2,o3,o7 in 0..2) */
static int pair_verify(u128 R,u64 maxv,int o2,int o3,int o7,u64 *ox,u64 *oy){
    if(R==0) return 0;
    u64 x=iroot6(R); if(x>maxv) x=maxv;
    u64 xmin=iroot6((R+1)/2); if(ipow6(xmin)*2<R) xmin++;
    if(xmin<1) xmin=1;
    for(;x>=xmin&&x>0;x--){
        int x2=(int)(x&1),x3=(x%3!=0),x7=(x%7!=0);
        if(x2>o2||x3>o3||x7>o7) continue;
        u128 rest=R-P6[x];
        if(rest==0) continue;
        u64 y=iroot6(rest);
        if(y==0||y>x||ipow6(y)!=rest) continue;
        int y2=(int)(y&1),y3=(y%3!=0),y7=(y%7!=0);
        if(x2+y2!=o2||x3+y3!=o3||x7+y7!=o7) continue;
        *ox=x;*oy=y; return 1;
    }
    return 0;
}

/* R = sum of exactly j positive sixth powers, bases <= maxv, exact budgets */
static int dfs(u128 R,int j,u64 maxv,int o2,int o3,int o7,u64 *out,Ctr *C){
    if(o2<0||o3<0||o7<0) return 0;
    if(o2>j||o3>j||o7>j) return 0;
    if(R==0) return 0;
    if(j==2){ C->leaves++; }
    if(!mask_ok(R,j)) return 0;
    if(j==1){
        u64 e=iroot6(R);
        if(e==0||e>maxv||ipow6(e)!=R) return 0;
        if((int)(e&1)!=o2||(int)(e%3!=0)!=o3||(int)(e%7!=0)!=o7) return 0;
        out[0]=e; return 1;
    }
    if(j==2){
        C->bqueries++;
        if(!bloom_query(R)) return 0;
        C->bpositives++;
        C->exactver++;
        return pair_verify(R,maxv,o2,o3,o7,out,out+1);
    }
    u64 hi=iroot6(R); if(hi>maxv) hi=maxv;
    u64 lo=iroot6((R+j-1)/j);
    if(ipow6(lo)*(u128)j<R) lo++;
    if(lo<1) lo=1;
    int d2=(o2==0)?2:1,d3=(o3==0)?3:1,d7=(o7==0)?7:1;
    u64 stride=(u64)(d2*d3*d7);
    u64 x=hi; if(stride>1) x-=x%stride;
    for(;x>=lo&&x>0;){
        C->dfsnodes++;
        int x2=(int)(x&1),x3=(x%3!=0),x7=(x%7!=0);
        int skip=0;
        if(o2==j&&!x2) skip=1;
        if(o2==0&&x2) skip=1;
        if(o3==j&&!x3) skip=1;
        if(o3==0&&x3) skip=1;
        if(o7==j&&!x7) skip=1;
        if(o7==0&&x7) skip=1;
        if(!skip){
            out[0]=x;
            if(dfs(R-P6[x],j-1,x,o2-x2,o3-x3,o7-x7,out+1,C)) return 1;
        }
        if(x<stride) break;
        x-=stride;
    }
    return 0;
}

/* exact division by R7 with a divisibility assertion */
/* R*IR7INV mod 2^128 equals R/R7 exactly when R7 | R, and lands above
 * (2^128-1)/R7 otherwise (x -> x*R7 mod 2^128 is a bijection and does not
 * overflow precisely for x <= (2^128-1)/R7).  So the magnitude test IS the
 * divisibility test -- no 128-bit division anywhere in the hot path. */
static u128 R7_QLIMIT;
static inline int div_r7(u128 R,u128 *q){ u128 t=R*IR7INV; if(t>R7_QLIMIT) return 0; *q=t; return 1; }

/* --------------------------------------------------------------- strata */
/* Each returns 1 and writes the four ORIGINAL bases into b[0..3]. */

static int solve_k70(u128 m,u64 B,u64 *b,Ctr *C,i64 *val_elim){
    u128 mp;
    if(!div_r7(m,&mp)){ (*val_elim)++; return 0; }  /* valuation lemma: 7^6 | m required */
    u64 Bp=B/7;
    if(Bp==0) return 0;
    int k2=(int)((u64)mp&7), k3=(int)md9_((u64)(mp>>64),(u64)mp), k7=(int)md7_((u64)(mp>>64),(u64)mp);
    if(k2>4||k3>4||k7>4) return 0;
    u64 out[4];
    if(dfs(mp,4,Bp,k2,k3,k7,out,C)){
        for(int i=0;i<4;i++) b[i]=7*out[i];
        return 1;
    }
    return 0;
}

static int solve_k71(u128 m,u64 B,u64 *b,Ctr *C){
    u64 Bp=B/7; if(Bp==0) return 0;
    u32 u=(u32)(m%R7);
    if(u==0) return 0;                       /* would force 7 | b4 */
    if(!ROOTCNT[u]) return 0;                /* m mod 7^6 not a sixth power */
    u64 b4hi=iroot6(m); if(b4hi>B) b4hi=B;
    u64 out[3];
    for(int i=0;i<6;i++){
        u64 r=ROOTLIST[(size_t)u*6+i];
        for(u64 b4=r; b4<=b4hi; b4+=R7){
            if(b4==0) continue;
            u128 p=P6[b4];
            if(p>=m) break;
            u128 R=m-p,Rp;
            ASSERT(div_r7(R,&Rp),"k7=1: R7 does not divide m-b4^6 (b4=%llu)",(unsigned long long)b4);
            int k2=(int)((u64)Rp&7),k3=(int)md9_((u64)(Rp>>64),(u64)Rp),k7=(int)md7_((u64)(Rp>>64),(u64)Rp);
            if(k2>3||k3>3||k7>3) continue;
            if(dfs(Rp,3,Bp,k2,k3,k7,out,C)){
                b[0]=7*out[0]; b[1]=7*out[1]; b[2]=7*out[2]; b[3]=b4;
                return 1;
            }
        }
    }
    return 0;
}

static int solve_k72(u128 m,u64 B,int k2,int k3,u64 *b,Ctr *C){
    u64 Bp=B/7; if(Bp==0) return 0;
    u32 mr=(u32)(m%R7);
    /* residues mod 42 that b3 and b4 may take */
    uint8_t allow[42];
    for(int r=0;r<42;r++){
        int ok = (r%7!=0);
        if(k2==0 && (r&1)) ok=0;
        if(k2==4 && !(r&1)) ok=0;
        if(k3==0 && (r%3)) ok=0;
        if(k3==4 && (r%3==0)) ok=0;
        allow[r]=(uint8_t)ok;
    }
    u64 b3hi=iroot6(m/2); if(b3hi>B) b3hi=B;
    u64 out[2];
    for(u64 b3=1;b3<=b3hi;b3++){
        if(!allow[b3%42]) continue;
        u32 p=P6R7[b3%R7];
        u32 u = (mr>=p)? (mr-p) : (mr+(u32)R7-p);
        if(u==0||!ROOTCNT[u]) continue;
        u128 b36=P6[b3];
        for(int i=0;i<6;i++){
            u64 r=ROOTLIST[(size_t)u*6+i];
            u64 b4=r;
            if(b4<b3){ u64 k=(b3-b4+R7-1)/R7; b4+=k*R7; }
            for(;b4<=B;b4+=R7){
                if(!allow[b4%42]) continue;
                u128 s=b36+P6[b4];
                if(s>=m) break;
                u128 R=m-s,Rp;
                ASSERT(div_r7(R,&Rp),"k7=2: R7 does not divide m-b3^6-b4^6 (b3=%llu b4=%llu)",
                       (unsigned long long)b3,(unsigned long long)b4);
                int o2=(int)((u64)Rp&7),o3=(int)md9_((u64)(Rp>>64),(u64)Rp),o7=(int)md7_((u64)(Rp>>64),(u64)Rp);
                if(o2>2||o3>2||o7>2) continue;
                if(dfs(Rp,2,Bp,o2,o3,o7,out,C)){
                    b[0]=7*out[0]; b[1]=7*out[1]; b[2]=b3; b[3]=b4;
                    return 1;
                }
            }
        }
    }
    return 0;
}

/* ------------------------------------------------------- table construction */
static void build_tables(u64 Bmax,u64 bpp){
    Bmax_global=Bmax; Bp_global=Bmax/7;
    P6=malloc(sizeof(u128)*(Bmax+2));
    if(!P6) DIE("alloc P6 (%llu entries)",(unsigned long long)Bmax+2);
    for(u64 x=0;x<=Bmax+1;x++) P6[x]=ipow6(x);

    u64 Bp=Bp_global;
    u64 npairs = Bp? Bp*(Bp+1)/2 : 1;
    nlines=(npairs*bpp+511)/512; if(!nlines) nlines=1;
    bloom=calloc(nlines*64,1);
    if(!bloom) DIE("bloom alloc failed (%.2f GB)",nlines*64.0/1e9);
    fprintf(stderr,"bloom: Bp=%llu pairs=%llu lines=%llu (%.3f GB, %llu bits/pair)\n",
            (unsigned long long)Bp,(unsigned long long)npairs,(unsigned long long)nlines,
            nlines*64.0/1e9,(unsigned long long)bpp);
    #pragma omp parallel for schedule(dynamic,64)
    for(u64 x=1;x<=Bp;x++)
        for(u64 y=1;y<=x;y++)
            bloom_add(P6[x]+P6[y]);
    fprintf(stderr,"bloom built\n");
}

static void self_test(void){
    /* Bloom must never produce a false negative on an inserted pair */
    if(Bp_global>=9){
        if(!bloom_query(P6[5]+P6[9])) DIE("BLOOM SELF-TEST (false negative)");
        if(!bloom_query(P6[1]+P6[1])) DIE("BLOOM SELF-TEST (false negative 1+1)");
    }
    /* exact division by 7^6 */
    u128 q; ASSERT(div_r7((u128)R7*123456789u,&q) && q==123456789u,"div_r7 broken");
    for(u64 z=1;z<R7;z+=13) ASSERT(!div_r7((u128)R7*1000003u+z,&q),"div_r7 accepted a non-multiple (+%llu)",(unsigned long long)z);
    ASSERT(div_r7(R7_QLIMIT*(u128)R7,&q)&&q==R7_QLIMIT,"div_r7 at limit");
    /* iroot6 */
    for(u64 x=1;x<100000;x+=997){
        ASSERT(iroot6(ipow6(x))==x,"iroot6 exact %llu",(unsigned long long)x);
        ASSERT(iroot6(ipow6(x)-1)==x-1,"iroot6 minus one %llu",(unsigned long long)x);
    }
}

/* ------------------------------------------------------------ candidates */
typedef struct { u64 f,t; u32 k2,k3,k7; u128 m; } Cand;

static int cand_cmp(const void *A,const void *B){
    const Cand *a=A,*b=B;
    if(a->f!=b->f) return a->f<b->f?-1:1;
    if(a->t!=b->t) return a->t<b->t?-1:1;
    return 0;
}

/* m*42^6 + t^6 == f^6 checked mod p = 2^61-1 */
#define MP61 2305843009213693951ULL
static inline u64 mulmod61(u64 a,u64 b){ return (u64)((u128)a*b % MP61); }
static inline u64 pow6mod61(u64 x){ u64 r=x%MP61,s=mulmod61(r,r); return mulmod61(mulmod61(s,s),s); }
static int check_m_mod_p(u64 f,u64 t,u128 m){
    u64 mm=(u64)(m%MP61);
    u64 lhs=(mulmod61(mm,M42%MP61) + pow6mod61(t))%MP61;
    return lhs==pow6mod61(f);
}

/* ------------------------------------------------------------------- main */
static void report_solution(u64 f,u64 t,const u64 *b,int k7){
    char buf[64];
    printf("SOLUTION f=%llu t=%llu k7=%d b1=%llu b2=%llu b3=%llu b4=%llu\n",
           (unsigned long long)f,(unsigned long long)t,k7,
           (unsigned long long)b[0],(unsigned long long)b[1],
           (unsigned long long)b[2],(unsigned long long)b[3]);
    printf("SOLUTION-SIX %llu %llu %llu %llu %llu %llu\n",
           (unsigned long long)(42*b[0]),(unsigned long long)(42*b[1]),
           (unsigned long long)(42*b[2]),(unsigned long long)(42*b[3]),
           (unsigned long long)t,(unsigned long long)f);
    (void)buf;
    fflush(stdout);
}

int main(int argc,char **argv){
    if(argc<3){
        fprintf(stderr,"usage: %s FMIN FMAX [--plant M B] [--candidates FILE] [--bpp N]\n",argv[0]);
        fprintf(stderr,"  processes every f with FMIN < f <= FMAX  (half-open, (FMIN,FMAX])\n");
        return 2;
    }
    u64 fmin=strtoull(argv[1],0,10), fmax=strtoull(argv[2],0,10);
    const char *candfile=NULL;
    int plant=0; u128 plant_m=0; u64 plant_B=0; u64 bpp=16; const char *isqrt6file=NULL, *plantfile=NULL; u64 bs_B=0,bs_Q=0; int bloomstat=0;
    for(int i=3;i<argc;i++){
        if(!strcmp(argv[i],"--plant")){ if(i+2>=argc) DIE("--plant needs M and B"); plant=1; plant_m=str_to_u128(argv[i+1]); plant_B=strtoull(argv[i+2],0,10); i+=2; }
        else if(!strcmp(argv[i],"--candidates")){ if(i+1>=argc) DIE("--candidates needs FILE"); candfile=argv[i+1]; i++; }
        else if(!strcmp(argv[i],"--bloomstat")){ if(i+2>=argc) DIE("--bloomstat needs B NQ"); bloomstat=1; bs_B=strtoull(argv[i+1],0,10); bs_Q=strtoull(argv[i+2],0,10); i+=2; }
        else if(!strcmp(argv[i],"--plantfile")){ if(i+1>=argc) DIE("--plantfile needs FILE"); plantfile=argv[i+1]; i++; }
        else if(!strcmp(argv[i],"--isqrt6")){ if(i+1>=argc) DIE("--isqrt6 needs FILE"); isqrt6file=argv[i+1]; i++; }
        else if(!strcmp(argv[i],"--bpp")){ if(i+1>=argc) DIE("--bpp needs N"); bpp=strtoull(argv[i+1],0,10); i++; }
        else DIE("unknown argument %s",argv[i]);
    }

    /* largest r with r^6 < 2^128 */
    { u64 lo=1,hi=1ULL<<22,v; u128 tmp;
      while(lo<hi){ u64 mid=lo+(hi-lo+1)/2; if(!pow6_ovf(mid,&tmp)) lo=mid; else hi=mid-1; }
      IR6MAX=lo; v=lo; ASSERT(!pow6_ovf(v,&tmp),"IR6MAX"); ASSERT(pow6_ovf(v+1,&tmp),"IR6MAX+1"); }
    IR7INV=inv_mod_2_128((u128)R7);
    ASSERT(IR7INV*(u128)R7==1,"2-adic inverse of 7^6");
    R7_QLIMIT=(~(u128)0)/R7;

    if(isqrt6file){
        FILE *fp=fopen(isqrt6file,"r"); if(!fp) DIE("cannot open %s",isqrt6file);
        char line[64];
        while(fgets(line,sizeof line,fp)){
            size_t L=strlen(line); while(L&&(line[L-1]=='\n'||line[L-1]=='\r')) line[--L]=0;
            if(!L) continue;
            u128 v=str_to_u128(line);
            printf("%llu\n",(unsigned long long)iroot6(v));
        }
        fclose(fp); return 0;
    }

    build_roots(); init_masks(); build_r7_roots();

#ifdef _OPENMP
    nthreads=omp_get_max_threads();
#endif
    ctrs=calloc(nthreads,sizeof(Ctr));
    if(!ctrs) DIE("alloc counters");

    double t0=0;
#ifdef _OPENMP
    t0=omp_get_wtime();
#else
    { struct timespec ts; clock_gettime(CLOCK_MONOTONIC,&ts); t0=ts.tv_sec+1e-9*ts.tv_nsec; }
#endif
    if(bloomstat){
        build_tables(bs_B*7,bpp);
        u64 Bp=Bp_global;
        /* (a) no false negative on any inserted pair */
        i64 fn=0;
        #pragma omp parallel for schedule(dynamic,64) reduction(+:fn)
        for(u64 x=1;x<=Bp;x++) for(u64 y=1;y<=x;y++) if(!bloom_query(P6[x]+P6[y])) fn++;
        /* (b) false positives on random 128-bit non-members.  Every pair sum is
         * <= 2*Bp^6; drawn values are forced above that bound, so no draw can be
         * a member and every positive is a genuine false positive. */
        u128 bound=2*P6[Bp];
        i64 fp=0;
        #pragma omp parallel for schedule(static) reduction(+:fp)
        for(u64 q=0;q<bs_Q;q++){
            u64 a=mixA(q*0x9e3779b97f4a7c15ULL+12345), b=mixB(q*0xc2b2ae3d27d4eb4fULL+6789);
            u128 v=((u128)a<<64)|b;
            if(v<=bound) v+=bound+1;
            if(bloom_query(v)) fp++;
        }
        double ideal=pow(1.0-exp(-8.0/(double)bpp),8.0);
        printf("BLOOMSTAT Bp=%llu pairs=%llu bpp=%llu false_negatives=%lld "
               "queries=%llu false_positives=%lld fp_rate=%.6g ideal=%.6g ratio=%.4f\n",
               (unsigned long long)Bp,(unsigned long long)(Bp*(Bp+1)/2),(unsigned long long)bpp,
               fn,(unsigned long long)bs_Q,fp,(double)fp/(double)bs_Q,ideal,
               ((double)fp/(double)bs_Q)/ideal);
        return fn?1:0;
    }

    if(plantfile){
        /* batch plant mode: each line "M B"; tables sized to the largest B */
        FILE *fp=fopen(plantfile,"r"); if(!fp) DIE("cannot open %s",plantfile);
        size_t n=0,cap=1024; u128 *ms=malloc(cap*sizeof(u128)); u64 *bs=malloc(cap*sizeof(u64));
        char line[128]; u64 Bmx=1;
        while(fgets(line,sizeof line,fp)){
            char *sp=strchr(line,' '); if(!sp) continue; *sp=0;
            if(n==cap){ cap*=2; ms=realloc(ms,cap*sizeof(u128)); bs=realloc(bs,cap*sizeof(u64)); }
            ms[n]=str_to_u128(line); bs[n]=strtoull(sp+1,0,10);
            if(bs[n]>Bmx) Bmx=bs[n];
            n++;
        }
        fclose(fp);
        build_tables(Bmx,bpp); self_test();
        Ctr *C=&ctrs[0];
        for(size_t i=0;i<n;i++){
            u128 m=ms[i]; u64 B=bs[i];
            int k2=(int)((u64)m&7),k3=(int)md9_((u64)(m>>64),(u64)m),k7=(int)md7_((u64)(m>>64),(u64)m);
            u64 b[4]; int found=0; i64 ve=0;
            if(k2<=4&&k3<=4&&k7<=2){
                if(k7==0) found=solve_k70(m,B,b,C,&ve);
                else if(k7==1) found=solve_k71(m,B,b,C);
                else found=solve_k72(m,B,k2,k3,b,C);
            }
            if(found){
                u128 ss=0; for(int q=0;q<4;q++) ss+=ipow6(b[q]);
                ASSERT(ss==m,"batch plant: decomposition does not sum to m (line %zu)",i);
                u64 sb[4]={b[0],b[1],b[2],b[3]};
                for(int a2=0;a2<4;a2++) for(int b2=a2+1;b2<4;b2++) if(sb[b2]<sb[a2]){u64 w=sb[a2];sb[a2]=sb[b2];sb[b2]=w;}
                printf("%zu FOUND %llu %llu %llu %llu\n",i,(unsigned long long)sb[0],(unsigned long long)sb[1],
                       (unsigned long long)sb[2],(unsigned long long)sb[3]);
            } else printf("%zu NONE\n",i);
        }
        return 0;
    }

    if(plant){
        build_tables(plant_B,bpp);
        self_test();
        Ctr *C=&ctrs[0];
        u128 m=plant_m;
        int k2=(int)((u64)m&7),k3=(int)md9_((u64)(m>>64),(u64)m),k7=(int)md7_((u64)(m>>64),(u64)m);
        char mb[48]; u128_to_str(m,mb);
        printf("PLANT m=%s B=%llu k2=%d k3=%d k7=%d\n",mb,(unsigned long long)plant_B,k2,k3,k7);
        if(k2>4||k3>4){ printf("PLANTNONE reason=budget\n"); return 0; }
        if(k7>=3){ printf("PLANTNONE reason=k7-not-in-stratum\n"); return 0; }
        u64 b[4]; int found=0; i64 ve=0;
        if(k7==0) found=solve_k70(m,plant_B,b,C,&ve);
        else if(k7==1) found=solve_k71(m,plant_B,b,C);
        else found=solve_k72(m,plant_B,k2,k3,b,C);
        if(found){
            /* exact re-verification in 128-bit */
            u128 s=0; for(int i=0;i<4;i++) s+=ipow6(b[i]);
            ASSERT(s==m,"planted decomposition does not sum to m");
            u64 sb[4]={b[0],b[1],b[2],b[3]};
            for(int i=0;i<4;i++) for(int j=i+1;j<4;j++) if(sb[j]<sb[i]){u64 w=sb[i];sb[i]=sb[j];sb[j]=w;}
            printf("PLANTFOUND %llu %llu %llu %llu\n",(unsigned long long)sb[0],(unsigned long long)sb[1],
                   (unsigned long long)sb[2],(unsigned long long)sb[3]);
        } else printf("PLANTNONE reason=no-decomposition\n");
        return 0;
    }

    if(fmax<=fmin) DIE("need FMAX > FMIN");
    if(fmax>100000000ULL) DIE("FMAX capped at 1e8");
    u64 Bmax=(fmax-1)/42;
    build_tables(Bmax,bpp);
    double tbuild;
#ifdef _OPENMP
    tbuild=omp_get_wtime()-t0;
#else
    { struct timespec ts; clock_gettime(CLOCK_MONOTONIC,&ts); tbuild=ts.tv_sec+1e-9*ts.tv_nsec-t0; }
#endif
    self_test();

    /* per-thread candidate buffers */
    Cand **cbuf=calloc(nthreads,sizeof(Cand*));
    size_t *cn=calloc(nthreads,sizeof(size_t)),*ccap=calloc(nthreads,sizeof(size_t));
    if(!cbuf||!cn||!ccap) DIE("alloc candidate buffers");

    i64 cand_total=0, cnt_k7[7]={0,0,0,0,0,0,0};
    i64 elim_budget=0, elim_k7=0, notcovered=0, proc0=0,proc1=0,proc2=0, val_elim=0;
    i64 nsol=0;

    double s0;
#ifdef _OPENMP
    s0=omp_get_wtime();
#else
    { struct timespec ts; clock_gettime(CLOCK_MONOTONIC,&ts); s0=ts.tv_sec+1e-9*ts.tv_nsec; }
#endif

    #pragma omp parallel for schedule(dynamic,1) \
        reduction(+:cand_total,elim_budget,elim_k7,notcovered,proc0,proc1,proc2,val_elim,nsol) \
        reduction(+:cnt_k7[:7])
    for(u64 f=fmin+1; f<=fmax; f++){
        if(f%2==0||f%3==0||f%7==0) continue;
        int tid=0;
#ifdef _OPENMP
        tid=omp_get_thread_num();
#endif
        Ctr *C=&ctrs[tid];
        u64 B=(f-1)/42;
        for(int i=0;i<nroots;i++){
            u64 t=(u64)((u128)roots144[i]*f%M42);
            if(t==0||t>=f) continue;
            cand_total++;

            /* m = (f^6-t^6)/42^6 via the cyclotomic factorisation, gcd-peeling 42^6
             * exactly as src/caseA2.c:282-294 (f^6 itself does not fit in 128 bits). */
            u64 g[2]={f-t,f+t};
            u128 qq[2]={(u128)f*f-(u128)f*t+(u128)t*t,(u128)f*f+(u128)f*t+(u128)t*t};
            u64 rem=M42; u128 m=1;
            for(int s=0;s<2;s++){
                u64 a=g[s],b2=rem; while(b2){u64 w=a%b2;a=b2;b2=w;}
                rem/=a; m*=g[s]/a;
            }
            for(int s=0;s<2;s++){
                u64 a=(u64)(qq[s]%rem);
                if(a){u64 x=a,y=rem;while(y){u64 w=x%y;x=y;y=w;}a=x;} else a=rem;
                rem/=a; m*=qq[s]/a;
            }
            ASSERT(rem==1,"42^6 did not divide f^6-t^6 (f=%llu t=%llu rem=%llu)",
                   (unsigned long long)f,(unsigned long long)t,(unsigned long long)rem);
            ASSERT(check_m_mod_p(f,t,m),"m*42^6+t^6 != f^6 mod 2^61-1 (f=%llu t=%llu)",
                   (unsigned long long)f,(unsigned long long)t);

            int k2=(int)((u64)m&7), k3=(int)md9_((u64)(m>>64),(u64)m), k7=(int)md7_((u64)(m>>64),(u64)m);
            cnt_k7[k7]++;

            if(candfile){
                if(cn[tid]==ccap[tid]){ ccap[tid]=ccap[tid]?ccap[tid]*2:1024; cbuf[tid]=realloc(cbuf[tid],ccap[tid]*sizeof(Cand)); if(!cbuf[tid]) DIE("realloc candidates"); }
                Cand *cc=&cbuf[tid][cn[tid]++];
                cc->f=f; cc->t=t; cc->k2=(u32)k2; cc->k3=(u32)k3; cc->k7=(u32)k7; cc->m=m;
            }

            if(k2>4||k3>4){ elim_budget++; continue; }   /* proven impossible */
            if(k7>4){ elim_k7++; continue; }             /* proven impossible */
            if(k7>=3){ notcovered++; continue; }         /* NOT COVERED by this engine */

            u64 b[4]; int found=0;
            if(k7==0){ proc0++; found=solve_k70(m,B,b,C,&val_elim); }
            else if(k7==1){ proc1++; found=solve_k71(m,B,b,C); }
            else { proc2++; found=solve_k72(m,B,k2,k3,b,C); }

            if(found){
                u128 s=0; for(int q=0;q<4;q++) s+=ipow6(b[q]);
                ASSERT(s==m,"reported decomposition does not sum to m (f=%llu)",(unsigned long long)f);
                nsol++;
                #pragma omp critical(sol)
                report_solution(f,t,b,k7);
            }
        }
    }

    double s1;
#ifdef _OPENMP
    s1=omp_get_wtime();
#else
    { struct timespec ts; clock_gettime(CLOCK_MONOTONIC,&ts); s1=ts.tv_sec+1e-9*ts.tv_nsec; }
#endif

    i64 leaves=0,bq=0,bp=0,ev=0,dn=0;
    for(int i=0;i<nthreads;i++){ leaves+=ctrs[i].leaves; bq+=ctrs[i].bqueries; bp+=ctrs[i].bpositives; ev+=ctrs[i].exactver; dn+=ctrs[i].dfsnodes; }

    if(candfile){
        size_t tot=0; for(int i=0;i<nthreads;i++) tot+=cn[i];
        Cand *all=malloc(tot*sizeof(Cand)); if(!all&&tot) DIE("alloc merged candidates");
        size_t k=0; for(int i=0;i<nthreads;i++){ memcpy(all+k,cbuf[i],cn[i]*sizeof(Cand)); k+=cn[i]; }
        qsort(all,tot,sizeof(Cand),cand_cmp);
        FILE *fp=fopen(candfile,"w"); if(!fp) DIE("cannot open %s",candfile);
        fprintf(fp,"# CANDIDATE f t k2 k3 k7 m   (band (%llu,%llu], class 1)\n",
                (unsigned long long)fmin,(unsigned long long)fmax);
        char mb[48];
        for(size_t i=0;i<tot;i++){
            u128_to_str(all[i].m,mb);
            fprintf(fp,"CANDIDATE %llu %llu %u %u %u %s\n",(unsigned long long)all[i].f,
                    (unsigned long long)all[i].t,all[i].k2,all[i].k3,all[i].k7,mb);
        }
        fclose(fp); free(all);
    }

    printf("BAND fmin=%llu fmax=%llu convention=(fmin,fmax] candidates=%lld "
           "k7_0=%lld k7_1=%lld k7_2=%lld k7_3=%lld k7_4=%lld k7_5=%lld k7_6=%lld "
           "elim_k2k3=%lld elim_k7=%lld notcovered_k7_3_4=%lld "
           "processed_k7_0=%lld processed_k7_1=%lld processed_k7_2=%lld val_elim_k7_0=%lld "
           "leaves=%lld bloom_queries=%lld bloom_positives=%lld exact_verifications=%lld dfs_nodes=%lld "
           "solutions=%lld bloom_build_s=%.3f search_s=%.3f elapsed_s=%.3f threads=%d\n",
           (unsigned long long)fmin,(unsigned long long)fmax,cand_total,
           cnt_k7[0],cnt_k7[1],cnt_k7[2],cnt_k7[3],cnt_k7[4],cnt_k7[5],cnt_k7[6],
           elim_budget,elim_k7,notcovered,proc0,proc1,proc2,val_elim,
           leaves,bq,bp,ev,dn,nsol,tbuild,s1-s0,s1-t0,nthreads);
    fflush(stdout);
    return nsol?3:0;
}
