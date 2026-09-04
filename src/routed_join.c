/*
 * routed_join.c -- exact external bucket-routing oracle prototype.
 *
 * Pair sums and queries are each enumerated once and routed to the same hash
 * bucket.  One bucket is loaded, sorted, and joined at a time.  This is an
 * exact proof-of-concept for removing caseA3's repeated enumeration, with
 * disk O(N^2+Q) and RAM O(N^2/NB+Q/NB).
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

typedef unsigned __int128 u128;
typedef uint64_t u64;

typedef struct { u128 value; unsigned char expected; } Query;

static u128 *p6;
static u128 ipow6(u64 x){ u128 s=(u128)x*x; return s*s*s; }
static u64 mix(u64 h){ h^=h>>33; h*=0xff51afd7ed558ccdULL; h^=h>>33; h*=0xc4ceb9fe1a85ec53ULL; h^=h>>33; return h; }
static u64 bucket_of(u128 s,u64 nb){
    u64 h=mix((u64)s^(0xd6e8feb86659fd93ULL*((u64)(s>>64)+1)));
    return (u64)(((u128)h*nb)>>64);
}
static int cmp128(const void *aa,const void *bb){
    u128 a=*(const u128*)aa,b=*(const u128*)bb;
    return (a>b)-(a<b);
}
static int contains(const u128 *a,size_t n,u128 key){
    size_t lo=0,hi=n;
    while(lo<hi){ size_t m=lo+(hi-lo)/2; if(a[m]<key) lo=m+1; else hi=m; }
    return lo<n&&a[lo]==key;
}
static int exact_pair(u128 target,u64 n){
    for(u64 x=1;x<=n;x++){
        if(p6[x]>=target) continue;
        u128 rest=target-p6[x];
        size_t lo=1,hi=n+1;
        while(lo<hi){ size_t m=lo+(hi-lo)/2; if(p6[m]<rest) lo=m+1; else hi=m; }
        if(lo<=n&&p6[lo]==rest) return 1;
    }
    return 0;
}
static u64 rnd(u64 *state){
    *state^=*state>>12; *state^=*state<<25; *state^=*state>>27;
    return *state*2685821657736338717ULL;
}

int main(int argc,char **argv){
    u64 n=argc>1?strtoull(argv[1],0,10):2000;
    u64 nb=argc>2?strtoull(argv[2],0,10):8;
    u64 nq=argc>3?strtoull(argv[3],0,10):20000;
    if(!n||!nb||nb>128){ fprintf(stderr,"usage: %s [N>0] [1<=NB<=128] [Q]\n",argv[0]); return 2; }
    FILE **pf=calloc(nb,sizeof(*pf)),**qf=calloc(nb,sizeof(*qf));
    u64 *pc=calloc(nb,sizeof(*pc)),*qc=calloc(nb,sizeof(*qc));
    p6=malloc(sizeof(*p6)*(n+1));
    if(!pf||!qf||!pc||!qc||!p6){ fprintf(stderr,"allocation failed\n"); return 1; }
    for(u64 b=0;b<nb;b++) if(!(pf[b]=tmpfile())||!(qf[b]=tmpfile())){ fprintf(stderr,"tmpfile failed\n"); return 1; }
    for(u64 x=0;x<=n;x++) p6[x]=ipow6(x);

    u64 pairs=0;
    for(u64 x=1;x<=n;x++) for(u64 y=1;y<=x;y++){
        u128 sum=p6[x]+p6[y]; u64 b=bucket_of(sum,nb);
        if(fwrite(&sum,sizeof(sum),1,pf[b])!=1){ fprintf(stderr,"pair write failed\n"); return 1; }
        pc[b]++; pairs++;
    }

    u64 state=0x615e50f6d00dULL;
    for(u64 i=0;i<nq;i++){
        Query q;
        if((i&1)==0){
            u64 x=1+rnd(&state)%n,y=1+rnd(&state)%n;
            q.value=p6[x]+p6[y]; q.expected=1;
        }else{
            q.value=(u128)(1+rnd(&state)%nq)+p6[1+rnd(&state)%n];
            q.expected=(unsigned char)exact_pair(q.value,n);
        }
        u64 b=bucket_of(q.value,nb);
        if(fwrite(&q,sizeof(q),1,qf[b])!=1){ fprintf(stderr,"query write failed\n"); return 1; }
        qc[b]++;
    }

    u64 peak=0,mismatches=0,hits=0;
    for(u64 b=0;b<nb;b++){
        if(pc[b]>peak) peak=pc[b];
        u128 *table=malloc(sizeof(*table)*pc[b]);
        Query *queries=malloc(sizeof(*queries)*qc[b]);
        if((pc[b]&&!table)||(qc[b]&&!queries)){ fprintf(stderr,"bucket allocation failed\n"); return 1; }
        rewind(pf[b]); rewind(qf[b]);
        if(fread(table,sizeof(*table),pc[b],pf[b])!=pc[b]){ fprintf(stderr,"pair read failed\n"); return 1; }
        if(fread(queries,sizeof(*queries),qc[b],qf[b])!=qc[b]){ fprintf(stderr,"query read failed\n"); return 1; }
        qsort(table,pc[b],sizeof(*table),cmp128);
        for(u64 i=0;i<qc[b];i++){
            int got=contains(table,pc[b],queries[i].value); hits+=got;
            if(got!=queries[i].expected) mismatches++;
        }
        free(table); free(queries); fclose(pf[b]); fclose(qf[b]);
    }
    printf("routed_join N=%llu buckets=%llu pairs=%llu pair_passes=1 queries=%llu query_passes=1 peak_bucket_entries=%llu disk_pair_bytes=%llu hits=%llu mismatches=%llu\n",
        (unsigned long long)n,(unsigned long long)nb,(unsigned long long)pairs,
        (unsigned long long)nq,(unsigned long long)peak,(unsigned long long)(pairs*sizeof(u128)),
        (unsigned long long)hits,(unsigned long long)mismatches);
    return mismatches?1:0;
}
