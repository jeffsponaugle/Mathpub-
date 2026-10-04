// Enumerate (n-3)-prefixes of nondecreasing unit-fraction representations of 1 with n terms,
// and count "heads" = (n-2)-prefixes, bucketed by the size of the x_{n-2} range.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64; typedef __int128 i128; typedef unsigned __int128 u128;
static u64 gcd64(u64 a,u64 b){while(b){u64 t=a%b;a=b;b=t;}return a;}
static int N;
static u64 nheads_bin[40], npref_bin[40]; static u64 total_heads=0,total_pref=0;
static u64 maxq=0;
// remainder p/q (reduced), k terms placed, last denominator last
static void rec(int k,u64 p,u64 q,u64 last){
  // need: terms k+1..N each >= last; remainder p/q > 0
  int rem=N-k; // terms remaining
  if(rem==3){ // this is an (N-3)-prefix; x_{N-2} ranges [lo,hi]
    u64 lo = q/p+1; if(lo<last) lo=last;
    u64 hi = (u64)(((u128)3*q)/p);
    total_pref++;
    if(hi>=lo){ u64 cnt=hi-lo+1; total_heads+=cnt; int b=0; u64 c=cnt; while(c>=10){c/=10;b++;} nheads_bin[b]+=cnt; npref_bin[b]++; }
    if(q>maxq) maxq=q;
    return;
  }
  u64 lo = q/p+1; if(lo<last) lo=last;
  u64 hi = (u64)(((u128)rem*q)/p);
  for(u64 x=lo;x<=hi;x++){
    // new remainder p/q - 1/x = (p*x - q)/(q*x)
    u128 np=(u128)p*x-q; u128 nq=(u128)q*x;
    u64 g=gcd64((u64)np,(u64)(nq%np)); // gcd(np,nq)=gcd(np, nq mod np)
    if(nq>>64){ // reduce carefully
      // compute gcd via mod
      u128 a=np,b=nq; while(b){u128 t=a%b;a=b;b=t;} g=(u64)a;
    }
    u128 rp=np/g, rq=nq/g;
    if(rq>>64){ fprintf(stderr,"overflow q\n"); exit(1);}
    rec(k+1,(u64)rp,(u64)rq,x);
  }
}
int main(int argc,char**argv){
  N=atoi(argv[1]);
  rec(0,1,1,2);
  printf("N=%d: (N-3)-prefixes=%llu heads=%llu maxq=%llu\n",N,total_pref,total_heads,maxq);
  for(int b=0;b<40;b++) if(npref_bin[b]) printf("range 1e%d..1e%d: prefixes %llu heads %llu\n",b,b+1,npref_bin[b],nheads_bin[b]);
  return 0;
}
