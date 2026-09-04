/*
 * mmap_bloom.c -- disk-backed monolithic Bloom oracle prototype.
 *
 * It enumerates pairs once and answers queries once.  mmap lets the kernel
 * page the filter, but does not guarantee a hard RSS cap; random page faults
 * are the expected failure mode when the file is much larger than RAM.
 */
#define _POSIX_C_SOURCE 200809L
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <unistd.h>

typedef unsigned __int128 u128;
typedef uint64_t u64;
typedef uint32_t u32;

static u64 *bits,nlines;
static u128 *p6;
static u128 ipow6(u64 x){u128 s=(u128)x*x;return s*s*s;}
static u64 mix(u64 h){h^=h>>33;h*=0xff51afd7ed558ccdULL;h^=h>>33;h*=0xc4ceb9fe1a85ec53ULL;h^=h>>33;return h;}
static void hashes(u128 s,u64 *h1,u64 *h2){u64 lo=(u64)s,hi=(u64)(s>>64);*h1=mix(lo^(0x9e3779b97f4a7c15ULL*hi));*h2=mix(hi^(0xbf58476d1ce4e5b9ULL+lo))|1ULL;}
static u64 line_of(u64 h){return (u64)(((u128)h*nlines)>>64);}
static void add(u128 s){u64 h1,h2;hashes(s,&h1,&h2);u64 *line=bits+(line_of(h1)<<3);for(int i=0;i<8;i++){u32 b=(u32)((h2>>(8*i))&511);line[b>>6]|=1ULL<<(b&63);}}
static int maybe(u128 s){u64 h1,h2;hashes(s,&h1,&h2);u64 *line=bits+(line_of(h1)<<3);for(int i=0;i<8;i++){u32 b=(u32)((h2>>(8*i))&511);if(!(line[b>>6]>>(b&63)&1))return 0;}return 1;}
static int exact_pair(u128 target,u64 n){for(u64 x=1;x<=n;x++){if(p6[x]>=target)continue;u128 r=target-p6[x];size_t lo=1,hi=n+1;while(lo<hi){size_t m=lo+(hi-lo)/2;if(p6[m]<r)lo=m+1;else hi=m;}if(lo<=n&&p6[lo]==r)return 1;}return 0;}
static u64 rnd(u64 *s){*s^=*s>>12;*s^=*s<<25;*s^=*s>>27;return *s*2685821657736338717ULL;}

int main(int argc,char **argv){
    u64 n=argc>1?strtoull(argv[1],0,10):5000,bpp=argc>2?strtoull(argv[2],0,10):12,nq=argc>3?strtoull(argv[3],0,10):20000;
    u64 pairs=n*(n+1)/2; nlines=(pairs*bpp+511)/512; size_t bytes=(size_t)nlines*64;
    FILE *file=tmpfile(); if(!file||ftruncate(fileno(file),(off_t)bytes)){perror("backing file");return 1;}
    bits=mmap(0,bytes,PROT_READ|PROT_WRITE,MAP_SHARED,fileno(file),0); if(bits==MAP_FAILED){perror("mmap");return 1;}
    p6=malloc(sizeof(*p6)*(n+1)); if(!p6){fprintf(stderr,"powers allocation failed\n");return 1;}
    for(u64 x=0;x<=n;x++)p6[x]=ipow6(x);
    for(u64 x=1;x<=n;x++)for(u64 y=1;y<=x;y++)add(p6[x]+p6[y]);

    u64 state=0x615e50f6d00dULL,positives=0,false_positives=0,false_negatives=0;
    for(u64 i=0;i<nq;i++){
        u128 q;
        if((i&1)==0){u64 x=1+rnd(&state)%n,y=1+rnd(&state)%n;q=p6[x]+p6[y];}
        else q=(u128)(1+rnd(&state)%nq)+p6[1+rnd(&state)%n];
        int truth=exact_pair(q,n),filter=maybe(q); positives+=filter;
        if(filter&&!truth) false_positives++;
        if(!filter&&truth) false_negatives++;
    }
    printf("mmap_bloom N=%llu pairs=%llu pair_passes=1 queries=%llu query_passes=1 file_bytes=%zu filter_positives=%llu false_positives=%llu false_negatives=%llu\n",
        (unsigned long long)n,(unsigned long long)pairs,(unsigned long long)nq,bytes,
        (unsigned long long)positives,(unsigned long long)false_positives,(unsigned long long)false_negatives);
    munmap(bits,bytes); fclose(file); free(p6);
    return false_negatives?1:0;
}
