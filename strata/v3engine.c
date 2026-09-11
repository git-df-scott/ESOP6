/*
 * v3engine.c -- table-free-leaf class-1 engine for a^6+b^6+c^6+d^6+e^6 = f^6.
 *
 * Implements strata/SPEC_V3.md.  See strata/V3_LOG.md for runs, and for the
 * ONE place where this code deliberately departs from the spec (section
 * "COMPLETENESS" below; the spec's step-1 completeness claim is false for the
 * primes p with p^(6*v_p(b)) <= Bmax, and the ROUGH table has to be enlarged
 * accordingly).
 *
 * SETTING (identical to k7engine.c / src/caseA2.c).  Class-1 primitive
 * solution: one term t carries all three exemptions; the other four are 42*b_i,
 *      f^6 - t^6 = 42^6 * m,   m = b1^6+b2^6+b3^6+b4^6,  1 <= b_i <= B=(f-1)/42,
 * t = zeta*f mod 42^6 for one of the 144 sixth roots of unity zeta, 0<t<f,
 * gcd(f,42)=1.  Exact budgets k2 = m mod 8, k3 = m mod 9, k7 = m mod 7.
 *
 * BAND CONVENTION: (FMIN, FMAX], half-open, as in k7engine.c.
 *
 * OUTER STRUCTURE: src/caseA2.c:dfs() replicated *exactly* (same nested
 * windows, same strides, same skip rules, same five residue masks 64/27/49/13/43
 * with the same incremental residues, same traversal order), with nb = 4 only.
 * The j == 2 leaf is the only thing that changes: instead of a Bloom filter over
 * all pair sums it calls is_two_sum() below.  Leaves are counted at exactly the
 * point caseA2_timed.c counts j2_nodes (after the five masks), so the counts are
 * directly comparable.
 *
 * LEAF: is_two_sum(R, bound) -- is R = b1^6 + b2^6 with 1 <= b1,b2 <= bound?
 *   Step 0  two composite 2-sum residue masks, Q1 = 64*27*49 = 84672 and
 *           Q2 = 13*19*31*37*43 = 12182287, applied incrementally in the DFS
 *           (tables PQ1[b], PQ2[b]); two bitmap lookups per leaf, no division.
 *   Step 1  for every prime p <= 700, with k_p minimal such that p^k_p > Bmax:
 *           (i)  DETERMINATION.  r = R mod p^k_p (incremental from T_p[b]).  If
 *                r is a unit sixth-power residue, each sixth root rho of r mod
 *                p^k_p is a candidate for the base NOT divisible by p (unique
 *                representative below the bound because p^k_p > Bmax >= bound).
 *                Test R - rho^6 for being a sixth power.
 *                Roots: mod p by table, then Hensel-lifted to p^k_p (p >= 5);
 *                p = 2 and p = 3 divide 6 so Hensel does not apply and one root
 *                per residue is tabulated directly mod 2^k_2 / 3^k_3, the other
 *                roots coming from the sixth roots of unity.
 *           (ii) RECURSION.  If p^6 | R (exact test by 2-adic inverse, no
 *                division), recurse on (R/p^6, bound/p): that is the case p |
 *                both bases.
 *           The prime loop is NEVER broken out of on a failed hypothesis.
 *   Step 2  RESIDUAL TABLE over V x V (see COMPLETENESS), blocked Bloom filter,
 *           exact verification restricted to V x V on a hit.
 *
 * COMPLETENESS (and the departure from SPEC_V3.md).
 *   Determination at p is valid only when p^k_p divides b^6 for the base b that
 *   p divides, i.e. when 6*v_p(b) >= k_p, i.e. when p^(6*v_p(b)) > Bmax.  For
 *   p >= 11 (with Bmax < 11^6) k_p <= 6 and v_p >= 1 suffices, so every prime
 *   factor in [11,700] determines.  For p in {2,3,5,7} it does NOT: e.g. b1 =
 *   2*q, b2 = q' leaves R mod 2^k_2 unrelated to b2^6.  SPEC_V3.md's proof
 *   ("if p divides exactly one of them, step 1(i) at p finds it") is false in
 *   exactly that case, and the residual class is not {1} u {primes > 700}.
 *   Correct statement.  Let lim = floor(Bmax^(1/6)) and
 *       E   = { e : every prime power p^v || e has p^v <= lim }  (a set of
 *             7-smooth numbers, |E| = 12..32 for Bmax = 1e5..5e5),
 *       V   = { e*q <= Bmax : e in E, q = 1 or q prime > 700 }.
 *   Induction on R.  If some prime p <= 700 divides both bases, p^6 | R and
 *   step 1(ii) recurses on a strictly smaller instance.  Otherwise the bases are
 *   coprime except possibly for a common prime q > 700; in that case both bases
 *   are q*x, q*y with x,y <= Bmax/700 < 701 and gcd(x,y) = 1, so x,y have no
 *   prime factor > 700 either.  If either base has a prime power p^v with
 *   p <= 700 and p^(6v) > Bmax then p divides that base only, and step 1(i) at p
 *   determines the other base exactly (it is < p^k_p, so it is the unique
 *   representative of its class).  Otherwise every prime power p^v || b with
 *   p <= 700 satisfies p^v <= lim, so the <=700-part of each base lies in E, and
 *   since 701^2 > Bmax each base has at most one prime factor > 700: both bases
 *   lie in V, and step 2 finds the pair.  QED.
 *   Cost of the fix: |V| ~ 2.8 * pi(Bmax) instead of pi(Bmax), i.e. ~7.8x the
 *   pairs of the spec's ROUGH table.  RAM therefore forces a blocked Bloom
 *   filter (8 bits/pair by default) rather than the spec's sorted u64
 *   fingerprint array; every hit is verified exactly against V x V.
 *
 * Build: gcc -O3 -march=native -fopenmp -std=gnu11 -Wall -Wextra -o v3engine v3engine.c -lm
 * Run:   v3engine FMIN FMAX [--k7 LIST] [--candidates FILE] [--plant M B]
 *                 [--twosum R BOUND] [--bmax N] [--bpp N] [--nores]
 *                 [--nostep1] [--notable] [--selftest]
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

/* ------------------------------------------------------------------ basics */
static u64 IR6MAX = 2642245;            /* largest r with r^6 < 2^128 */
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
    double e=pow((double)R,1.0/6.0);
    if(!(e>=0)) e=0;
    if(e>(double)IR6MAX) e=(double)IR6MAX;
    u64 r=(u64)e; if(r>IR6MAX) r=IR6MAX;
    while(r>0 && ipow6(r)>R) r--;
    while(r<IR6MAX && ipow6(r+1)<=R) r++;
    return r;
}
static u64 pow6m(u64 x,u64 m){ u64 r=x%m; u128 s=(u128)r*r%m; s=s*s%m; return (u64)((u128)s*((u128)r*r%m)%m); }
static u128 inv_mod_2_128(u128 a){ u128 x=a; for(int i=0;i<7;i++) x=x*(2-a*x); return x; }

static void u128_to_str(u128 v,char *buf){
    char tmp[48]; int n=0;
    if(v==0){ buf[0]='0'; buf[1]=0; return; }
    while(v){ tmp[n++]=(char)('0'+(int)(v%10)); v/=10; }
    for(int i=0;i<n;i++) buf[i]=tmp[n-1-i];
    buf[n]=0;
}
static u128 str_to_u128(const char *s){
    u128 v=0; for(;*s;s++){ if(*s<'0'||*s>'9') DIE("bad u128 literal '%s'",s); v=v*10+(u64)(*s-'0'); } return v;
}

/* fast reduction of a u64 by a modulus md < 2^32 using magic = floor(2^64/md) */
static inline u32 redu(u64 t,u32 md,u64 magic){
    u64 q=(u64)(((u128)t*magic)>>64);
    u64 r=t-q*(u64)md;
    while(r>=md) r-=md;
    return (u32)r;
}
static inline u32 mulmod32(u32 a,u32 b,u32 md,u64 magic){ return redu((u64)a*b,md,magic); }
static inline u32 submod32(u32 a,u32 b,u32 md){ return a>=b? a-b : a+md-b; }

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
    for(int i=0;i<nroots;i++) ASSERT(pow6m(roots144[i],M42)==1,"bad root");
}

/* ----------------------- caseA2's residue masks (identical set and order) --- */
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
                for(u32 s=0;s<md;s++) if(single>>s&1) cur|=1ULL<<((r+s)%md);
            sumres[m][j]=cur;
        }
    }
    u64 single=0;
    for(u32 x=0;x<64;x++) single|=1ULL<<(u32)pow6m(x,64);
    sumres64[0]=1;
    for(int j=1;j<5;j++){
        u64 prev=sumres64[j-1],cur=0;
        for(u32 r=0;r<64;r++) if(prev>>r&1)
            for(u32 s=0;s<64;s++) if(single>>s&1) cur|=1ULL<<((r+s)&63);
        sumres64[j]=cur;
    }
}

/* --------------------------- composite 2-sum masks Q1, Q2 (spec step 0) ----- */
#define Q1 84672u          /* 64*27*49   */
#define Q2 12182287u       /* 13*19*31*37*43  (SPEC_V3.md prints 12190001, which is not that product) */
static u64 *TWO1=NULL,*TWO2=NULL;      /* bitmaps of {x^6+y^6 mod Q} */
static u32 *PQ1=NULL,*PQ2=NULL;        /* b^6 mod Q1 / Q2 for b <= Bmax */
static u64 MAGQ1,MAGQ2;
static u64 *SIXQ1=NULL;    /* bitmap of sixth-power residues mod Q1 (288 of 84672) */
static double mask1_frac=0,mask2_frac=0;

static void build_two_sum_mask(u32 Q,u64 **bm,double *frac,u64 **sixthout){
    u64 nw=(Q+63)/64;
    u64 *sixth=calloc(nw,8); u64 *bits=calloc(nw,8);
    if(!sixth||!bits) DIE("alloc two-sum mask");
    u32 *vals=malloc(sizeof(u32)*Q); u32 nv=0;
    for(u32 x=0;x<Q;x++){
        u32 r=(u32)pow6m(x,Q);
        if(!(sixth[r>>6]>>(r&63)&1)){ sixth[r>>6]|=1ULL<<(r&63); vals[nv++]=r; }
    }
    for(u32 i=0;i<nv;i++) for(u32 j=0;j<nv;j++){
        u64 s=(u64)vals[i]+vals[j]; if(s>=Q) s-=Q;
        bits[s>>6]|=1ULL<<(s&63);
    }
    u64 cnt=0; for(u32 r=0;r<Q;r++) if(bits[r>>6]>>(r&63)&1) cnt++;
    *frac=(double)cnt/(double)Q;
    free(vals);
    if(sixthout) *sixthout=sixth; else free(sixth);
    *bm=bits;
}

/* ----------------------------- primes <= 700 and their determination data --- */
#define NPMAX 130
static u32 PCUT=700;   /* least prime with PCUT^2 > Bmax: guarantees every base <= Bmax
                          has at most one prime factor > PCUT (see COMPLETENESS) */
static int NP=0;
static u32 PR[NPMAX];        /* p */
static u32 PKS[NPMAX];       /* p^k_p > Bmax */
static int PKE[NPMAX];       /* k_p */
static u64 PKMAG[NPMAX];     /* floor(2^64/PKS) */
static u64 P6V[NPMAX];       /* p^6 */
static u128 P6INV[NPMAX];    /* (p^6)^{-1} mod 2^128 */
static u128 P6LIM[NPMAX];    /* (2^128-1)/p^6 */
/* roots of x^6 = r mod p, for p >= 5 (flat) */
static u32 *RP_off=NULL;     /* offset into RP_list per prime */
static uint8_t *RP_cnt=NULL; /* per (prime,residue) count */
static u32 *RP_list=NULL;    /* 6 slots per (prime,residue) */
static u32 *RP_w=NULL;       /* (6*x0^5)^{-1} mod p per (prime, x0) */
/* direct root tables for p = 2 and p = 3 (Hensel does not apply) */
static u32 *ROOT2=NULL,*ROOT3=NULL;    /* one root mod p^k, 0xffffffff if none */
static u32 U2[8]; static int NU2=0;    /* sixth roots of unity mod 2^k_2 */
static u32 U3[8]; static int NU3=0;
static int IDX2=-1,IDX3=-1;

static u64 Bmax_global=0;
static u128 *P6=NULL;                  /* x^6 for x <= Bmax */
static u32 *TP=NULL;                   /* TP[b*NP+i] = b^6 mod PKS[i] */
static int use_tp=1;


static u64 modinv_small(u64 a,u64 m){
    long long t=0,nt=1; long long r=(long long)m,nr=(long long)(a%m);
    while(nr){ long long q=r/nr; long long tmp=t-q*nt; t=nt; nt=tmp; tmp=r-q*nr; r=nr; nr=tmp; }
    if(r>1) return 0;
    if(t<0) t+=(long long)m;
    return (u64)t;
}

static void build_primes(u64 Bmax){
    /* cutoff: least prime whose square exceeds Bmax (never above 700 for FMAX<=2e7) */
    {   u32 c=2;
        for(;;){ int isp=1;
                 for(u32 d=2;d*d<=c;d++) if(c%d==0){isp=0;break;}
                 if(isp&&(u64)c*c>Bmax) break;
                 c++; }
        PCUT=c;
        ASSERT((u64)PCUT*PCUT>Bmax,"PCUT^2 <= Bmax");
        ASSERT(PCUT<=700,"PCUT=%u > 700 (FMAX too large)",PCUT); }
    char *sv=calloc(PCUT+1,1);
    for(u32 i=2;i<=PCUT;i++) sv[i]=1;
    for(u32 i=2;(u64)i*i<=PCUT;i++) if(sv[i]) for(u32 j=i*i;j<=PCUT;j+=i) sv[j]=0;
    NP=0;
    for(u32 p=2;p<=PCUT;p++) if(sv[p]){
        u64 pk=p; int k=1;
        while(pk<=Bmax){ pk*=p; k++; }
        ASSERT(pk<(1ULL<<32),"p^k_p too large: p=%d k=%d",p,k);
        PR[NP]=(u32)p; PKS[NP]=(u32)pk; PKE[NP]=k;
        PKMAG[NP]=(u64)(((u128)1<<64)/pk);       /* floor(2^64/p^k_p) */
        u64 p6=1; for(int i=0;i<6;i++) p6*=(u64)p;
        P6V[NP]=p6;
        if(p==2){ IDX2=NP; } else if(p==3){ IDX3=NP; }
        else { P6INV[NP]=inv_mod_2_128((u128)p6); P6LIM[NP]=(~(u128)0)/p6; }
        NP++;
    }
    free(sv);
    ASSERT(NP<=NPMAX,"too many primes");
    /* p=3 exact division also by inverse (3^6 odd) */
    if(IDX3>=0){ P6INV[IDX3]=inv_mod_2_128((u128)P6V[IDX3]); P6LIM[IDX3]=(~(u128)0)/P6V[IDX3]; }

    /* root tables mod p (p >= 5) */
    u64 tot=0; RP_off=malloc(sizeof(u32)*NP);
    for(int i=0;i<NP;i++){ RP_off[i]=(u32)tot; tot+=PR[i]; }
    RP_cnt=calloc(tot,1); RP_list=calloc(tot*6,sizeof(u32)); RP_w=calloc(tot,sizeof(u32));
    if(!RP_cnt||!RP_list||!RP_w) DIE("alloc root tables");
    for(int i=0;i<NP;i++){
        u32 p=PR[i]; if(p<5) continue;
        u32 off=RP_off[i];
        for(u32 x=1;x<p;x++){
            u32 r=(u32)pow6m(x,p);
            ASSERT(RP_cnt[off+r]<6,"more than 6 sixth roots mod p=%u",p);
            RP_list[(size_t)(off+r)*6+RP_cnt[off+r]++]=x;
            /* w = (6 x^5)^{-1} mod p */
            u64 x5=1; for(int e=0;e<5;e++) x5=x5*x%p;
            u64 d=(6*x5)%p;
            u64 w=modinv_small(d,p);
            ASSERT(w!=0||d==0,"no inverse of 6x^5 mod %u",p);
            RP_w[off+x]=(u32)w;
        }
    }
    /* p = 2 */
    if(IDX2>=0){
        u32 pk=PKS[IDX2];
        ROOT2=malloc(sizeof(u32)*pk);
        for(u32 r=0;r<pk;r++) ROOT2[r]=0xffffffffu;
        for(u32 x=1;x<pk;x+=2){ u32 r=(u32)pow6m(x,pk); if(ROOT2[r]==0xffffffffu) ROOT2[r]=x; }
        NU2=0; for(u32 x=1;x<pk;x+=2) if(pow6m(x,pk)==1){ ASSERT(NU2<8,"too many sixth roots of unity mod 2^k"); U2[NU2++]=x; }
    }
    /* p = 3 */
    if(IDX3>=0){
        u32 pk=PKS[IDX3];
        ROOT3=malloc(sizeof(u32)*pk);
        for(u32 r=0;r<pk;r++) ROOT3[r]=0xffffffffu;
        for(u32 x=1;x<pk;x++){ if(x%3==0) continue; u32 r=(u32)pow6m(x,pk); if(ROOT3[r]==0xffffffffu) ROOT3[r]=x; }
        NU3=0; for(u32 x=1;x<pk;x++){ if(x%3==0) continue; if(pow6m(x,pk)==1){ ASSERT(NU3<8,"too many sixth roots of unity mod 3^k"); U3[NU3++]=x; } }
    }
}

/* Hensel lift of a root mod p to a root mod p^k (p >= 5, p does not divide 6) */
static inline u32 hensel_lift(int i,u32 x0,u32 r){
    u32 pk=PKS[i]; u64 mg=PKMAG[i];
    u32 w=RP_w[RP_off[i]+x0];
    u32 x=x0;
    for(int j=1;j<PKE[i];j++){
        u32 x2=mulmod32(x,x,pk,mg);
        u32 x3=mulmod32(x2,x,pk,mg);
        u32 x6=mulmod32(x3,x3,pk,mg);
        u32 f=submod32(x6,r,pk);
        if(!f) break;
        x=submod32(x,mulmod32(f,w,pk,mg),pk);
    }
    return x;
}

/* ------------------------------------- residual set V and its pair filter --- */
static u64 *Vlist=NULL; static u64 nV=0;
static u64 *Vbit=NULL;             /* membership bitmap over [0,Bmax] */
static u64 *bloom=NULL; static u64 nlines=0; static u64 bpp_global=8;

static inline u64 mixA(u64 h){ h^=h>>30; h*=0xbf58476d1ce4e5b9ULL; h^=h>>27; h*=0x94d049bb133111ebULL; h^=h>>31; return h; }
static inline u64 mixB(u64 h){ h^=h>>33; h*=0xff51afd7ed558ccdULL; h^=h>>33; h*=0xc4ceb9fe1a85ec53ULL; h^=h>>33; return h; }
static inline u64 mixC(u64 h){ h^=h>>29; h*=0xd6e8feb86659fd93ULL; h^=h>>32; h*=0xa5cb3b1f2f5a1d7bULL; h^=h>>29; return h; }
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
    u32 ps[8]={(u32)(a&511),(u32)((a>>9)&511),(u32)((a>>18)&511),(u32)((a>>27)&511),
               (u32)(b&511),(u32)((b>>9)&511),(u32)((b>>18)&511),(u32)((b>>27)&511)};
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

static int u64cmp(const void *a,const void *b){ u64 x=*(const u64*)a,y=*(const u64*)b; return x<y?-1:(x>y?1:0); }

static double vbuild_s=0, tbuild_s=0;
static u64 nEset=0; static u64 Eset[64];

static void build_V(u64 Bmax){
    /* lim = floor(Bmax^(1/6)) */
    u64 lim=1; while(ipow6(lim+1)<=(u128)Bmax) lim++;
    /* allowed prime powers p^v <= lim, one per prime */
    u64 tmp[64]; u64 nt=1; tmp[0]=1;
    for(u64 p=2;p<=lim;p++){
        int isp=1; for(u64 d=2;d*d<=p;d++) if(p%d==0){isp=0;break;}
        if(!isp) continue;
        u64 cur=nt;
        for(u64 v=p;v<=lim;v*=p)
            for(u64 i=0;i<cur;i++){ u64 e=tmp[i]*v; if(e<=Bmax){ ASSERT(nt<64,"E too large"); tmp[nt++]=e; } }
    }
    qsort(tmp,nt,8,u64cmp);
    nEset=nt; for(u64 i=0;i<nt;i++) Eset[i]=tmp[i];

    /* primes in (700, Bmax] */
    char *sv=calloc(Bmax+1,1);
    for(u64 i=2;i<=Bmax;i++) sv[i]=1;
    for(u64 i=2;i*i<=Bmax;i++) if(sv[i]) for(u64 j=i*i;j<=Bmax;j+=i) sv[j]=0;
    /* V */
    u64 cap=1024; Vlist=malloc(cap*8); nV=0;
    for(u64 i=0;i<nt;i++){
        u64 e=tmp[i];
        if(nV==cap){cap*=2;Vlist=realloc(Vlist,cap*8);}
        Vlist[nV++]=e;                                 /* q = 1 */
        for(u64 q=PCUT+1;q<=Bmax/e;q++) if(sv[q]){
            if(nV==cap){cap*=2;Vlist=realloc(Vlist,cap*8);}
            Vlist[nV++]=e*q;
        }
    }
    qsort(Vlist,nV,8,u64cmp);
    u64 k=0; for(u64 i=0;i<nV;i++){ if(i==0||Vlist[i]!=Vlist[i-1]) Vlist[k++]=Vlist[i]; } nV=k;
    Vbit=calloc((Bmax+64)/64,8);
    for(u64 i=0;i<nV;i++) Vbit[Vlist[i]>>6]|=1ULL<<(Vlist[i]&63);
    free(sv);
}
static inline int inV(u64 y){ return (Vbit[y>>6]>>(y&63))&1; }

static void build_pairfilter(void){
    u64 npairs=nV*(nV+1)/2; if(!npairs) npairs=1;
    nlines=(npairs*bpp_global+511)/512; if(!nlines) nlines=1;
    bloom=calloc(nlines*64,1);
    if(!bloom) DIE("pair filter alloc failed (%.2f GB)",nlines*64.0/1e9);
    fprintf(stderr,"vtable: |V|=%llu pairs=%llu bits/pair=%llu size=%.3f GB\n",
        (unsigned long long)nV,(unsigned long long)npairs,(unsigned long long)bpp_global,nlines*64.0/1e9);
    #pragma omp parallel for schedule(dynamic,64)
    for(u64 i=0;i<nV;i++){
        u128 a=P6[Vlist[i]];
        for(u64 j=0;j<=i;j++) bloom_add(a+P6[Vlist[j]]);
    }
}

/* ------------------------------------------------------------- counters ---- */
typedef struct {
    i64 leaves, step1, primecand, iroot_calls, recursions, tblq, tblpos, tblscan, dfsnodes, sols;
    char pad[128];
} Ctr;
static Ctr *ctrs=NULL; static int nthreads=1;
static int opt_nostep1=0, opt_notable=0;

static inline u32 modu128(u128 R,u32 md,u64 mag);
/* ------------------------------------------------- exact V x V verification */
static int verify_V(u128 R,u64 bound,u64 *ox,u64 *oy,Ctr *C){
    u64 xhi=iroot6(R); if(xhi>bound) xhi=bound;
    u64 xlo=iroot6((R+1)/2); if(ipow6(xlo)*2<R) xlo++;
    if(xlo<1) xlo=1;
    if(xhi<xlo) return 0;
    /* lower_bound in V */
    u64 lo=0,hi=nV;
    while(lo<hi){ u64 mid=(lo+hi)/2; if(Vlist[mid]<xlo) lo=mid+1; else hi=mid; }
    u32 rq1=modu128(R,Q1,MAGQ1);
    for(u64 i=lo;i<nV&&Vlist[i]<=xhi;i++){
        C->tblscan++;
        u64 x=Vlist[i];
        u32 s=submod32(rq1,PQ1[x],Q1);
        if(!((SIXQ1[s>>6]>>(s&63))&1)) continue;
        u128 S=R-P6[x];
        if(S==0) continue;
        u64 y=iroot6(S);
        if(y>=1&&y<=x&&y<=bound&&ipow6(y)==S&&inV(y)){ *ox=x; *oy=y; return 1; }
    }
    return 0;
}

/* ---------------------------------------------------------- the leaf test -- */
static inline u32 modu128(u128 R,u32 md,u64 mag){
    u64 hi=(u64)(R>>64), lo=(u64)R;
    if(!hi) return redu(lo,md,mag);
    u64 h=redu(hi,md,mag);
    u64 c=redu(~(u64)0,md,mag)+1; if(c>=md) c-=md;      /* 2^64 mod md */
    u64 t=(u64)h*c; u32 a=redu(t,md,mag);
    u32 b=redu(lo,md,mag);
    u32 s=a+b; if(s>=md) s-=md;
    return s;
}

static int is_two_sum(u128 R,u64 bound,const u32 *res,u64 *ox,u64 *oy,Ctr *C,int depth)
{
    if(bound==0||R<2) return 0;
    u64 lo=(u64)R;
    if((lo&7)>2) return 0;
    { u64 hi=(u64)(R>>64);
      u32 m9=(u32)((((hi%9)*7u)+(lo%9))%9); if(m9>2) return 0;
      u32 m7=(u32)((((hi%7)*2u)+(lo%7))%7); if(m7>2) return 0; }
    if(res==NULL){                         /* recursion level: apply masks here */
        u32 r1=modu128(R,Q1,MAGQ1), r2=modu128(R,Q2,MAGQ2);
        if(!((TWO1[r1>>6]>>(r1&63))&1)) return 0;
        if(!((TWO2[r2>>6]>>(r2&63))&1)) return 0;
    }
    C->step1++;
    u32 rq1=modu128(R,Q1,MAGQ1);

    u32 resbuf[NPMAX];
    if(res==NULL){ for(int i=0;i<NP;i++) resbuf[i]=modu128(R,PKS[i],PKMAG[i]); res=resbuf; }

    if(!opt_nostep1)
    for(int i=0;i<NP;i++){
        u32 p=PR[i], pk=PKS[i], r=res[i];
        if(r){
            u32 cand[8]; int nc=0;
            if(i==IDX2){
                if(r&1){ u32 x0=ROOT2[r];
                    if(x0!=0xffffffffu) for(int u=0;u<NU2;u++) cand[nc++]=mulmod32(x0,U2[u],pk,PKMAG[i]); }
            } else if(i==IDX3){
                if(r%3){ u32 x0=ROOT3[r];
                    if(x0!=0xffffffffu) for(int u=0;u<NU3;u++) cand[nc++]=mulmod32(x0,U3[u],pk,PKMAG[i]); }
            } else {
                u32 r0=r%p;
                if(r0){
                    u32 off=RP_off[i]; int cnt=RP_cnt[off+r0];
                    for(int t=0;t<cnt;t++) cand[nc++]=hensel_lift(i,RP_list[(size_t)(off+r0)*6+t],r);
                }
            }
            for(int t=0;t<nc;t++){
                u64 b2=cand[t];
                if(b2<1||b2>bound) continue;
                C->primecand++;
                u32 s=submod32(rq1,PQ1[b2],Q1);
                if(!((SIXQ1[s>>6]>>(s&63))&1)) continue;
                u128 P=P6[b2];
                if(P>=R) continue;
                u128 S=R-P;
                C->iroot_calls++;
                u64 y=iroot6(S);
                if(y>=1&&y<=bound&&ipow6(y)==S){ *ox=b2; *oy=y; return 1; }
            }
        }
        /* (ii) recursion: p^6 | R  ->  both bases divisible by p */
        {
            u128 q;
            int div=0;
            if(i==IDX2){ if((lo&63)==0){ q=R>>6; div=1; } }
            else { u128 tq=R*P6INV[i]; if(tq<=P6LIM[i]){ q=tq; div=1; } }
            if(div&&q>=2){
                C->recursions++;
                u64 x,y;
                if(is_two_sum(q,bound/p,NULL,&x,&y,C,depth+1)){ *ox=p*x; *oy=p*y; return 1; }
            }
        }
    }

    if(opt_notable) return 0;
    C->tblq++;
    if(!bloom_query(R)) return 0;
    C->tblpos++;
    return verify_V(R,bound,ox,oy,C);
}

/* ------------------------------ caseA2 dfs replica, nb = 4, v3 leaf --------- */
static const u32 *MRES_tls=NULL;   /* set per candidate (thread-private) */
#pragma omp threadprivate(MRES_tls)

static int dfs(u128 R,int j,u64 maxv,int o2,int o3,int o7,
               u32 r27,u32 r49,u32 r13,u32 r43,u32 rq1,u32 rq2,
               u64 *out,u64 *base,Ctr *C){
    if(o2>j||o3>j||o7>j) return 0;
    if(!(sumres64[j]>>((u32)R&63)&1)) return 0;
    if(!(sumres[0][j]>>r27&1)) return 0;
    if(!(sumres[1][j]>>r49&1)) return 0;
    if(!(sumres[2][j]>>r13&1)) return 0;
    if(!(sumres[3][j]>>r43&1)) return 0;
    if(j==2){
        C->leaves++;                       /* exactly caseA2_timed.c's j2_nodes */
        if(!((TWO1[rq1>>6]>>(rq1&63))&1)) return 0;
        if(!((TWO2[rq2>>6]>>(rq2&63))&1)) return 0;
        u32 resbuf[NPMAX]; const u32 *res=NULL;
        if(use_tp){
            u64 b4=base[0],b3=base[1];
            const u32 *t4=TP+(size_t)b4*NP,*t3=TP+(size_t)b3*NP;
            for(int i=0;i<NP;i++){
                u32 v=MRES_tls[i], pk=PKS[i];
                v=submod32(v,t4[i],pk); v=submod32(v,t3[i],pk);
                resbuf[i]=v;
            }
            res=resbuf;
        }
        return is_two_sum(R,maxv,res,out,out+1,C,0);
    }
    u64 hi=iroot6(R); if(hi>maxv) hi=maxv;
    u64 lo=iroot6((R+j-1)/j);
    if(ipow6(lo)*(u128)j<R) lo++;
    int d2=(o2==0)?2:1,d3=(o3==0)?3:1,d7=(o7==0)?7:1;
    int stride=d2*d3*d7;
    u64 x=hi; if(stride>1) x-=x%stride;
    for(;x>=lo&&x>0;x-=(u64)stride){
        C->dfsnodes++;
        int x2=(int)(x&1),x3=(x%3!=0),x7=(x%7!=0);
        if(o2==j&&!x2) continue;
        if(o2==0&&x2) continue;
        if(o3==j&&!x3) continue;
        if(o3==0&&x3) continue;
        if(o7==j&&!x7) continue;
        if(o7==0&&x7) continue;
        out[0]=x;
        u32 p27=powmod6[0][x%27],p49=powmod6[1][x%49],p13=powmod6[2][x%13],p43=powmod6[3][x%43];
        u32 q1=PQ1[x],q2=PQ2[x];
        if(dfs(R-P6[x],j-1,x,o2-x2,o3-x3,o7-x7,
               r27+((r27<p27)?27:0)-p27,r49+((r49<p49)?49:0)-p49,
               r13+((r13<p13)?13:0)-p13,r43+((r43<p43)?43:0)-p43,
               rq1+((rq1<q1)?Q1:0)-q1, rq2+((rq2<q2)?Q2:0)-q2,
               out+1,base,C)) return 1;
        if(x<(u64)stride) break;
    }
    return 0;
}

/* decompose m into four sixth powers with bases <= maxb */
static int decompose4(u128 m,u64 maxb,const u32 *mres,u64 *out,Ctr *C){
    int k2=(int)(m%8),k3=(int)(m%9),k7=(int)(m%7);
    if(k2>4||k3>4||k7>4) return 0;
    MRES_tls=mres;
    return dfs(m,4,maxb,k2,k3,k7,(u32)(m%27),(u32)(m%49),(u32)(m%13),(u32)(m%43),
               (u32)(m%Q1),(u32)(m%Q2),out,out,C);
}

/* ------------------------------------------------------- table construction */
static void build_tables(u64 Bmax){
    double t0,t1;
#ifdef _OPENMP
    t0=omp_get_wtime();
#else
    { struct timespec ts; clock_gettime(CLOCK_MONOTONIC,&ts); t0=ts.tv_sec+1e-9*ts.tv_nsec; }
#endif
    Bmax_global=Bmax;
    ASSERT((u128)Bmax < (u128)701*701,"Bmax=%llu >= 701^2: a base could have two prime factors > 700",
           (unsigned long long)Bmax);
    P6=malloc(sizeof(u128)*(Bmax+2)); if(!P6) DIE("alloc P6");
    for(u64 x=0;x<=Bmax+1;x++) P6[x]=ipow6(x);
    MAGQ1=(u64)(((u128)1<<64)/Q1); MAGQ2=(u64)(((u128)1<<64)/Q2);
    build_two_sum_mask(Q1,&TWO1,&mask1_frac,&SIXQ1);
    build_two_sum_mask(Q2,&TWO2,&mask2_frac,NULL);
    PQ1=malloc(sizeof(u32)*(Bmax+2)); PQ2=malloc(sizeof(u32)*(Bmax+2));
    if(!PQ1||!PQ2) DIE("alloc PQ");
    for(u64 x=0;x<=Bmax+1;x++){ PQ1[x]=(u32)pow6m(x,Q1); PQ2[x]=(u32)pow6m(x,Q2); }
    build_primes(Bmax);
    if(use_tp){
        size_t nbytes=(size_t)(Bmax+2)*NP*4;
        TP=malloc(nbytes); if(!TP) DIE("alloc TP (%.2f GB)",nbytes/1e9);
        #pragma omp parallel for schedule(static)
        for(u64 x=0;x<=Bmax+1;x++)
            for(int i=0;i<NP;i++) TP[(size_t)x*NP+i]=(u32)pow6m(x,PKS[i]);
    }
#ifdef _OPENMP
    t1=omp_get_wtime();
#else
    { struct timespec ts; clock_gettime(CLOCK_MONOTONIC,&ts); t1=ts.tv_sec+1e-9*ts.tv_nsec; }
#endif
    tbuild_s=t1-t0;
    build_V(Bmax);
    build_pairfilter();
#ifdef _OPENMP
    vbuild_s=omp_get_wtime()-t1;
#else
    { struct timespec ts; clock_gettime(CLOCK_MONOTONIC,&ts); vbuild_s=ts.tv_sec+1e-9*ts.tv_nsec-t1; }
#endif
    /* memory report */
    double mb=1e6;
    double m_p6=(double)(Bmax+2)*16, m_tp=use_tp?(double)(Bmax+2)*NP*4:0,
           m_pq=(double)(Bmax+2)*8, m_masks=(double)(Q1+Q2)/8,
           m_roots=(double)((IDX2>=0?PKS[IDX2]:0)+(IDX3>=0?PKS[IDX3]:0))*4,
           m_v=(double)nV*8+(double)Bmax/8, m_bl=(double)nlines*64;
    fprintf(stderr,"memory: P6=%.1fMB T_p=%.1fMB PQ=%.1fMB masks=%.1fMB root2,3=%.1fMB V=%.1fMB pairfilter=%.1fMB TOTAL=%.2fGB\n",
        m_p6/mb,m_tp/mb,m_pq/mb,m_masks/mb,m_roots/mb,m_v/mb,m_bl/mb,
        (m_p6+m_tp+m_pq+m_masks+m_roots+m_v+m_bl)/1e9);
    fprintf(stderr,"masks: Q1 pass=%.4f Q2 pass=%.4f combined=%.5f ; pcut=%u primes=%d |E|=%llu maxE=%llu\n",
        mask1_frac,mask2_frac,mask1_frac*mask2_frac,PCUT,NP,(unsigned long long)nEset,(unsigned long long)(nEset?Eset[nEset-1]:0));
}

/* --------------------------------------------------------------- self test */
static void self_test(void){
    Ctr C; memset(&C,0,sizeof C);
    /* iroot6 */
    for(u64 x=1;x<100000;x+=997){
        ASSERT(iroot6(ipow6(x))==x,"iroot6 exact %llu",(unsigned long long)x);
        ASSERT(iroot6(ipow6(x)-1)==x-1,"iroot6-1 %llu",(unsigned long long)x);
    }
    /* Hensel lift correctness for every prime >= 5 on a few residues */
    for(int i=0;i<NP;i++){
        u32 p=PR[i],pk=PKS[i]; if(p<5) continue;
        for(u32 x=1;x<pk&&x<50;x++){
            if(x%p==0) continue;
            u32 r=(u32)pow6m(x,pk);
            u32 r0=r%p; u32 off=RP_off[i]; int cnt=RP_cnt[off+r0];
            int seen=0;
            for(int t=0;t<cnt;t++){
                u32 z=hensel_lift(i,RP_list[(size_t)(off+r0)*6+t],r);
                ASSERT(pow6m(z,pk)==r,"hensel: p=%u k=%d z=%u r=%u",p,PKE[i],z,r);
                if(z==x) seen=1;
            }
            ASSERT(seen,"hensel missed the root x=%u mod p^k (p=%u)",x,p);
        }
    }
    /* p=2 and p=3 root tables */
    if(IDX2>=0){ u32 pk=PKS[IDX2];
        for(u32 x=1;x<pk&&x<200;x+=2){ u32 r=(u32)pow6m(x,pk); u32 x0=ROOT2[r];
            ASSERT(x0!=0xffffffffu,"ROOT2 missing"); int seen=0;
            for(int u=0;u<NU2;u++){ u32 z=mulmod32(x0,U2[u],pk,PKMAG[IDX2]); ASSERT(pow6m(z,pk)==r,"ROOT2 bad"); if(z==x) seen=1; }
            ASSERT(seen,"ROOT2 missed x=%u",x); } }
    if(IDX3>=0){ u32 pk=PKS[IDX3];
        for(u32 x=1;x<pk&&x<200;x++){ if(x%3==0) continue; u32 r=(u32)pow6m(x,pk); u32 x0=ROOT3[r];
            ASSERT(x0!=0xffffffffu,"ROOT3 missing"); int seen=0;
            for(int u=0;u<NU3;u++){ u32 z=mulmod32(x0,U3[u],pk,PKMAG[IDX3]); ASSERT(pow6m(z,pk)==r,"ROOT3 bad"); if(z==x) seen=1; }
            ASSERT(seen,"ROOT3 missed x=%u",x); } }
    /* exact division by p^6 */
    for(int i=0;i<NP;i++){
        if(i==IDX2) continue;
        u128 v=(u128)P6V[i]*123457; u128 t=v*P6INV[i];
        ASSERT(t<=P6LIM[i]&&t==123457,"p^6 exact division broken for p=%u",PR[i]);
        for(u64 z=1;z<P6V[i]&&z<50;z++){ u128 t2=((u128)P6V[i]*1000003+z)*P6INV[i]; ASSERT(t2>P6LIM[i],"p^6 divisibility false positive p=%u",PR[i]); }
    }
    /* pair filter: no false negatives on inserted pairs */
    for(u64 i=0;i<nV&&i<50;i++) for(u64 j=0;j<=i;j++)
        ASSERT(bloom_query(P6[Vlist[i]]+P6[Vlist[j]]),"pair filter false negative");
    /* modu128 */
    for(int i=0;i<NP;i++){
        u128 v=((u128)0x123456789abcdefULL<<64)|0xfedcba9876543210ULL;
        ASSERT(modu128(v,PKS[i],PKMAG[i])==(u32)(v%PKS[i]),"modu128 wrong for p=%u",PR[i]);
    }
    ASSERT(modu128((u128)12345,Q1,MAGQ1)==12345%Q1,"modu128 Q1");
    /* is_two_sum on a handful of hand-made cases */
    if(Bmax_global>=200){
        u64 x,y;
        ASSERT(is_two_sum(ipow6(3)+ipow6(5),200,NULL,&x,&y,&C,0),"is_two_sum 3,5");
        ASSERT(is_two_sum(ipow6(64)+ipow6(81),200,NULL,&x,&y,&C,0),"is_two_sum 64,81");
        ASSERT(!is_two_sum(ipow6(3)+ipow6(5)+1,200,NULL,&x,&y,&C,0),"is_two_sum false positive");
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
#define MP61 2305843009213693951ULL
static inline u64 mulmod61(u64 a,u64 b){ return (u64)((u128)a*b % MP61); }
static inline u64 pow6mod61(u64 x){ u64 r=x%MP61,s=mulmod61(r,r); return mulmod61(mulmod61(s,s),s); }
static int check_m_mod_p(u64 f,u64 t,u128 m){
    u64 mm=(u64)(m%MP61);
    u64 lhs=(mulmod61(mm,M42%MP61) + pow6mod61(t))%MP61;
    return lhs==pow6mod61(f);
}

static void report_solution(u64 f,u64 t,const u64 *b,int k7){
    printf("SOLUTION f=%llu t=%llu k7=%d b1=%llu b2=%llu b3=%llu b4=%llu\n",
           (unsigned long long)f,(unsigned long long)t,k7,
           (unsigned long long)b[0],(unsigned long long)b[1],
           (unsigned long long)b[2],(unsigned long long)b[3]);
    printf("SOLUTION-SIX %llu %llu %llu %llu %llu %llu\n",
           (unsigned long long)(42*b[0]),(unsigned long long)(42*b[1]),
           (unsigned long long)(42*b[2]),(unsigned long long)(42*b[3]),
           (unsigned long long)t,(unsigned long long)f);
    fflush(stdout);
}

int main(int argc,char **argv){
    if(argc<3){
        fprintf(stderr,"usage: %s FMIN FMAX [--k7 LIST] [--candidates FILE] [--plant M B] [--twosum R BOUND]\n",argv[0]);
        fprintf(stderr,"       [--bmax N] [--bpp N] [--nores] [--nostep1] [--notable] [--selftest]\n");
        fprintf(stderr,"  band convention: every f with FMIN < f <= FMAX\n");
        return 2;
    }
    u64 fmin=0,fmax=0;
    int have_band=0;
    if(argv[1][0]!='-'){ fmin=strtoull(argv[1],0,10); fmax=strtoull(argv[2],0,10); have_band=1; }
    const char *candfile=NULL,*twosumfile=NULL,*plantfile=NULL;
    int plant=0,twosum=0,selftest=0;
    u128 plant_m=0,ts_R=0; u64 plant_B=0,ts_bound=0,bmax_override=0;
    int k7sel[7]={0,0,0,1,1,0,0};
    for(int i=have_band?3:1;i<argc;i++){
        if(!strcmp(argv[i],"--k7")){ if(i+1>=argc) DIE("--k7 needs LIST");
            for(int q=0;q<7;q++) k7sel[q]=0;
            for(const char *s=argv[i+1];*s;s++) if(*s>='0'&&*s<='6') k7sel[*s-'0']=1;
            i++; }
        else if(!strcmp(argv[i],"--candidates")){ if(i+1>=argc) DIE("--candidates needs FILE"); candfile=argv[i+1]; i++; }
        else if(!strcmp(argv[i],"--plant")){ if(i+2>=argc) DIE("--plant needs M B"); plant=1; plant_m=str_to_u128(argv[i+1]); plant_B=strtoull(argv[i+2],0,10); i+=2; }
        else if(!strcmp(argv[i],"--plantfile")){ if(i+1>=argc) DIE("--plantfile needs FILE"); plantfile=argv[i+1]; i++; }
        else if(!strcmp(argv[i],"--twosum")){ if(i+2>=argc) DIE("--twosum needs R BOUND"); twosum=1; ts_R=str_to_u128(argv[i+1]); ts_bound=strtoull(argv[i+2],0,10); i+=2; }
        else if(!strcmp(argv[i],"--twosumfile")){ if(i+1>=argc) DIE("--twosumfile needs FILE"); twosumfile=argv[i+1]; i++; }
        else if(!strcmp(argv[i],"--bmax")){ if(i+1>=argc) DIE("--bmax needs N"); bmax_override=strtoull(argv[i+1],0,10); i++; }
        else if(!strcmp(argv[i],"--bpp")){ if(i+1>=argc) DIE("--bpp needs N"); bpp_global=strtoull(argv[i+1],0,10); i++; }
        else if(!strcmp(argv[i],"--nores")) use_tp=0;
        else if(!strcmp(argv[i],"--nostep1")) opt_nostep1=1;
        else if(!strcmp(argv[i],"--notable")) opt_notable=1;
        else if(!strcmp(argv[i],"--selftest")) selftest=1;
        else DIE("unknown argument %s",argv[i]);
    }
    { u64 lo=1,hi=1ULL<<22; u128 tmp;
      while(lo<hi){ u64 mid=lo+(hi-lo+1)/2; if(!pow6_ovf(mid,&tmp)) lo=mid; else hi=mid-1; }
      IR6MAX=lo; }
    build_roots(); init_sieves();
#ifdef _OPENMP
    nthreads=omp_get_max_threads();
#endif
    ctrs=calloc(nthreads,sizeof(Ctr)); if(!ctrs) DIE("alloc counters");

    double t0;
#ifdef _OPENMP
    t0=omp_get_wtime();
#else
    { struct timespec ts; clock_gettime(CLOCK_MONOTONIC,&ts); t0=ts.tv_sec+1e-9*ts.tv_nsec; }
#endif

    /* ---- table-only modes -------------------------------------------------- */
    if(twosumfile){
        /* each line: "R BOUND"; tables sized by --bmax (required) */
        if(!bmax_override) DIE("--twosumfile needs --bmax");
        build_tables(bmax_override); self_test();
        FILE *fp=fopen(twosumfile,"r"); if(!fp) DIE("cannot open %s",twosumfile);
        char line[128]; Ctr *C=&ctrs[0]; u64 n=0;
        while(fgets(line,sizeof line,fp)){
            char *sp=strchr(line,' '); if(!sp) continue; *sp=0;
            u128 R=str_to_u128(line); u64 bd=strtoull(sp+1,0,10);
            u64 x=0,y=0;
            int f=is_two_sum(R,bd,NULL,&x,&y,C,0);
            if(f){ u128 s=ipow6(x)+ipow6(y);
                   ASSERT(s==R,"twosum: reported pair does not sum to R (line %llu)",(unsigned long long)n);
                   ASSERT(x<=bd&&y<=bd&&x>=1&&y>=1,"twosum: pair out of range"); }
            printf("%llu %s %llu %llu\n",(unsigned long long)n,f?"YES":"NO",
                   (unsigned long long)x,(unsigned long long)y);
            n++;
        }
        fclose(fp); return 0;
    }
    if(twosum){
        build_tables(bmax_override?bmax_override:ts_bound); self_test();
        u64 x=0,y=0; Ctr *C=&ctrs[0];
        int f=is_two_sum(ts_R,ts_bound,NULL,&x,&y,C,0);
        if(f){ u128 s=ipow6(x)+ipow6(y); ASSERT(s==ts_R,"twosum pair does not sum to R"); }
        printf("TWOSUM %s %llu %llu\n",f?"YES":"NO",(unsigned long long)x,(unsigned long long)y);
        return 0;
    }
    if(plantfile){
        if(!bmax_override) DIE("--plantfile needs --bmax");
        build_tables(bmax_override); self_test();
        FILE *fp=fopen(plantfile,"r"); if(!fp) DIE("cannot open %s",plantfile);
        char line[128]; Ctr *C=&ctrs[0]; u64 n=0;
        u32 mres[NPMAX];
        while(fgets(line,sizeof line,fp)){
            char *sp=strchr(line,' '); if(!sp) continue; *sp=0;
            u128 m=str_to_u128(line); u64 B=strtoull(sp+1,0,10);
            for(int i=0;i<NP;i++) mres[i]=modu128(m,PKS[i],PKMAG[i]);
            u64 b[4];
            int f=decompose4(m,B,mres,b,C);
            if(f){ u128 s=0; for(int q=0;q<4;q++) s+=ipow6(b[q]);
                   ASSERT(s==m,"plant: decomposition does not sum to m (line %llu)",(unsigned long long)n);
                   u64 sb[4]={b[0],b[1],b[2],b[3]};
                   for(int a=0;a<4;a++) for(int c=a+1;c<4;c++) if(sb[c]<sb[a]){u64 w=sb[a];sb[a]=sb[c];sb[c]=w;}
                   for(int a=0;a<4;a++) ASSERT(sb[a]>=1&&sb[a]<=B,"plant: base out of range");
                   printf("%llu FOUND %llu %llu %llu %llu\n",(unsigned long long)n,
                       (unsigned long long)sb[0],(unsigned long long)sb[1],
                       (unsigned long long)sb[2],(unsigned long long)sb[3]);
            } else printf("%llu NONE\n",(unsigned long long)n);
            n++;
        }
        fclose(fp); return 0;
    }
    if(plant){
        build_tables(bmax_override?bmax_override:plant_B); self_test();
        Ctr *C=&ctrs[0]; u128 m=plant_m;
        u32 mres[NPMAX]; for(int i=0;i<NP;i++) mres[i]=modu128(m,PKS[i],PKMAG[i]);
        char mb[48]; u128_to_str(m,mb);
        int k2=(int)(m%8),k3=(int)(m%9),k7=(int)(m%7);
        printf("PLANT m=%s B=%llu k2=%d k3=%d k7=%d\n",mb,(unsigned long long)plant_B,k2,k3,k7);
        u64 b[4];
        if(decompose4(m,plant_B,mres,b,C)){
            u128 s=0; for(int i=0;i<4;i++) s+=ipow6(b[i]);
            ASSERT(s==m,"planted decomposition does not sum to m");
            u64 sb[4]={b[0],b[1],b[2],b[3]};
            for(int i=0;i<4;i++) for(int j=i+1;j<4;j++) if(sb[j]<sb[i]){u64 w=sb[i];sb[i]=sb[j];sb[j]=w;}
            printf("PLANTFOUND %llu %llu %llu %llu\n",(unsigned long long)sb[0],(unsigned long long)sb[1],
                   (unsigned long long)sb[2],(unsigned long long)sb[3]);
        } else printf("PLANTNONE\n");
        return 0;
    }
    if(!have_band) DIE("need FMIN FMAX");
    if(fmax<=fmin) DIE("need FMAX > FMIN");
    if(fmax>20000000ULL) DIE("FMAX capped at 2e7 (701^2 assertion on Bmax)");

    u64 Bmax=bmax_override?bmax_override:(fmax-1)/42;
    build_tables(Bmax);
    if(selftest) self_test();
    double tbuild;
#ifdef _OPENMP
    tbuild=omp_get_wtime()-t0;
#else
    { struct timespec ts; clock_gettime(CLOCK_MONOTONIC,&ts); tbuild=ts.tv_sec+1e-9*ts.tv_nsec-t0; }
#endif

    Cand **cbuf=calloc(nthreads,sizeof(Cand*));
    size_t *cn=calloc(nthreads,sizeof(size_t)),*ccap=calloc(nthreads,sizeof(size_t));
    if(!cbuf||!cn||!ccap) DIE("alloc candidate buffers");

    i64 cand_total=0,cnt_k7[7]={0,0,0,0,0,0,0};
    i64 elim_budget=0,elim_k7=0,skipped_k7=0,processed=0,nsol=0;

    double s0;
#ifdef _OPENMP
    s0=omp_get_wtime();
#else
    { struct timespec ts; clock_gettime(CLOCK_MONOTONIC,&ts); s0=ts.tv_sec+1e-9*ts.tv_nsec; }
#endif

    #pragma omp parallel for schedule(dynamic,1) \
        reduction(+:cand_total,elim_budget,elim_k7,skipped_k7,processed,nsol) \
        reduction(+:cnt_k7[:7])
    for(u64 f=fmin+1;f<=fmax;f++){
        if(f%2==0||f%3==0||f%7==0) continue;
        int tid=0;
#ifdef _OPENMP
        tid=omp_get_thread_num();
#endif
        Ctr *C=&ctrs[tid];
        u64 B=(f-1)/42;
        u32 mres[NPMAX];
        for(int i=0;i<nroots;i++){
            u64 t=(u64)((u128)roots144[i]*f%M42);
            if(t==0||t>=f) continue;
            cand_total++;
            u64 g[2]={f-t,f+t};
            u128 qq[2]={(u128)f*f-(u128)f*t+(u128)t*t,(u128)f*f+(u128)f*t+(u128)t*t};
            u64 rem=M42; u128 m=1;
            for(int s=0;s<2;s++){ u64 a=g[s],b2=rem; while(b2){u64 w=a%b2;a=b2;b2=w;} rem/=a; m*=g[s]/a; }
            for(int s=0;s<2;s++){
                u64 a=(u64)(qq[s]%rem);
                if(a){u64 x=a,y=rem;while(y){u64 w=x%y;x=y;y=w;}a=x;} else a=rem;
                rem/=a; m*=qq[s]/a;
            }
            ASSERT(rem==1,"42^6 did not divide f^6-t^6 (f=%llu t=%llu)",(unsigned long long)f,(unsigned long long)t);
            ASSERT(check_m_mod_p(f,t,m),"m*42^6+t^6 != f^6 mod 2^61-1 (f=%llu)",(unsigned long long)f);
            int k2=(int)(m%8),k3=(int)(m%9),k7=(int)(m%7);
            cnt_k7[k7]++;
            if(candfile){
                if(cn[tid]==ccap[tid]){ ccap[tid]=ccap[tid]?ccap[tid]*2:1024; cbuf[tid]=realloc(cbuf[tid],ccap[tid]*sizeof(Cand)); if(!cbuf[tid]) DIE("realloc candidates"); }
                Cand *cc=&cbuf[tid][cn[tid]++];
                cc->f=f; cc->t=t; cc->k2=(u32)k2; cc->k3=(u32)k3; cc->k7=(u32)k7; cc->m=m;
            }
            if(k2>4||k3>4){ elim_budget++; continue; }
            if(k7>4){ elim_k7++; continue; }
            if(!k7sel[k7]){ skipped_k7++; continue; }
            processed++;
            for(int q=0;q<NP;q++) mres[q]=modu128(m,PKS[q],PKMAG[q]);
            u64 b[4];
            if(decompose4(m,B,mres,b,C)){
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

    i64 leaves=0,step1=0,pc=0,ir=0,rec=0,tq=0,tp_=0,tscan=0,dn=0;
    for(int i=0;i<nthreads;i++){
        leaves+=ctrs[i].leaves; step1+=ctrs[i].step1; pc+=ctrs[i].primecand;
        ir+=ctrs[i].iroot_calls; rec+=ctrs[i].recursions; tq+=ctrs[i].tblq;
        tp_+=ctrs[i].tblpos; tscan+=ctrs[i].tblscan; dn+=ctrs[i].dfsnodes;
    }
    if(candfile){
        size_t tot=0; for(int i=0;i<nthreads;i++) tot+=cn[i];
        Cand *all=malloc((tot?tot:1)*sizeof(Cand)); if(!all) DIE("alloc merged candidates");
        size_t k=0; for(int i=0;i<nthreads;i++){ if(cn[i]) memcpy(all+k,cbuf[i],cn[i]*sizeof(Cand)); k+=cn[i]; }
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
           "elim_k2k3=%lld elim_k7=%lld skipped_k7_not_selected=%lld processed=%lld "
           "leaves=%lld step1=%lld prime_cands=%lld iroot_calls=%lld recursions=%lld "
           "table_queries=%lld table_positives=%lld table_scanned=%lld dfs_nodes=%lld "
           "solutions=%lld build_s=%.3f search_s=%.3f elapsed_s=%.3f threads=%d\n",
           (unsigned long long)fmin,(unsigned long long)fmax,cand_total,
           cnt_k7[0],cnt_k7[1],cnt_k7[2],cnt_k7[3],cnt_k7[4],cnt_k7[5],cnt_k7[6],
           elim_budget,elim_k7,skipped_k7,processed,
           leaves,step1,pc,ir,rec,tq,tp_,tscan,dn,nsol,tbuild,s1-s0,s1-t0,nthreads);
    fflush(stdout);
    return nsol?3:0;
}
