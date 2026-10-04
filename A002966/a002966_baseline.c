// Baseline: count nondecreasing n-term unit fraction representations of 1.
// Brute-force first n-2 terms; count last two by divisor enumeration of q^2 in a residue class.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned long long u64; typedef unsigned __int128 u128;
static int N; static u64 LIM; static uint32_t *spf;  // smallest prime factor sieve
static u64 gcd64(u64 a,u64 b){while(b){u64 t=a%b;a=b;b=t;}return a;}
static u128 gcd128(u128 a,u128 b){while(b){u128 t=a%b;a=b;b=t;}return a;}
// prime factorization of current q maintained as (prime, exponent) list; q = prod
#define MAXP 24
static u64 pr[MAXP]; static int ex[MAXP]; static int np=0;
static void addfac(u64 x,int sign){ // multiply/divide tracked factorization by x (x <= LIM)
  while(x>1){ u64 p=spf[x]; int e=0; while(x%p==0){x/=p;e++;}
    int i; for(i=0;i<np;i++) if(pr[i]==p) break;
    if(i==np){ pr[np]=p; ex[np]=0; np++; }
    ex[i]+=sign*e;
    if(ex[i]==0){ pr[i]=pr[np-1]; ex[i]=ex[np-1]; np--; }
  }
}
// count divisors d of q^2 with d <= q, d ≡ -q (mod p), d >= lo  (lo = (p*s - q) possibly negative -> 0)
// q = prod pr[i]^qe[i] where qe is exponent of p_i in q (reduced q may have smaller exponents than tracked lcm!)
static u64 Q,P,LO; static u64 ans;
static u64 qe_pr[MAXP]; static int qe_ex[MAXP]; static int qnp;
static u64 negq_modp;
static void divrec(int i,u128 d){
  if(i==qnp){ if(d<=Q && d>=LO){ if((u64)(d%P)==negq_modp) ans++; } return; }
  u128 dd=d; for(int e=0;e<=2*qe_ex[i];e++){ if(dd>Q) break; divrec(i+1,dd); dd*=qe_pr[i]; }
}
static u64 count2(u64 p,u64 q,u64 s){
  // factor q using tracked factorization (q | lcm tracked)
  qnp=0; u64 r=q; for(int i=0;i<np;i++){ int e=0; while(r%pr[i]==0){r/=pr[i];e++;} if(e){qe_pr[qnp]=pr[i];qe_ex[qnp]=e;qnp++;} }
  if(r!=1){ fprintf(stderr,"factor fail q=%llu r=%llu\n",q,r); exit(1);}
  Q=q;P=p; u128 lo=(u128)p*s; LO = lo>q ? (u64)(lo-q) : 0; negq_modp=(p-q%p)%p; ans=0; divrec(0,1); return ans;
}
static u64 total=0;
static void rec(int k,u64 p,u64 q,u64 last){
  int rem=N-k;
  if(rem==2){ total+=count2(p,q,last); return; }
  u64 lo=q/p+1; if(lo<last) lo=last; u64 hi=(u64)(((u128)rem*q)/p);
  for(u64 x=lo;x<=hi;x++){
    u128 np_=(u128)p*x-q, nq=(u128)q*x; u128 g=gcd128(np_,nq); np_/=g; nq/=g;
    if(nq>>64){fprintf(stderr,"q overflow\n");exit(1);}
    addfac(x,+1); rec(k+1,(u64)np_,(u64)nq,x); addfac(x,-1);
  }
}
int main(int argc,char**argv){
  N=atoi(argv[1]); LIM=(argc>2)?strtoull(argv[2],0,10):20000000ULL;
  spf=calloc(LIM+1,4); for(u64 i=2;i<=LIM;i++) if(!spf[i]) for(u64 j=i;j<=LIM;j+=i) if(!spf[j]) spf[j]=i;
  rec(0,1,1,2);
  printf("a(%d) = %llu\n",N,total);
  return 0;
}
