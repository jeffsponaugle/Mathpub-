// Instrument 3-term tails: count s<=x<=y<=z with 1/x+1/y+1/z = p/q, by iterating x and
// enumerating divisors d of (q x)^2 with d ≡ -q x (mod a), a = p x - q, (a-q)x <= d <= q x.
// Report histogram of solutions by log2(a), and by log2 of the "x-part" w of d relative to a/p.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef unsigned long long u64; typedef unsigned __int128 u128;
static u64 P,Q,S; static uint32_t *spf; static u64 LIM;
static u64 fp[40]; static int fe[40]; static int nf;
static void addfac(u64 x){ while(x>1){u64 p=(x<=LIM)?spf[x]:x; if(x>LIM){ // trial divide
      p=0; for(u64 t=2;t*t<=x;t++) if(x%t==0){p=t;break;} if(!p)p=x; }
    int e=0; while(x%p==0){x/=p;e++;} int i; for(i=0;i<nf;i++) if(fp[i]==p)break; if(i==nf){fp[nf]=p;fe[nf]=0;nf++;} fe[i]+=e; } }
static u64 hist_a[70], hist_small[70], hist_big[70]; static u64 nsol=0, nx_with=0;
static u128 QX, LO; static u64 A; static u64 X; static u64 sol_x;
static void divrec(int i,u128 d){
  if(i==nf){ if(d<=QX && d>=LO && (u64)(d%A)==(u64)((A-(QX%A))%A)){ nsol++; sol_x++; int b=0; u64 t=A; while(t>1){t>>=1;b++;} hist_a[b]++;
      // x-part w of d: gcd(d, x^2)-ish: compute w = product over primes of x
      u128 w=1; u64 xx=X; while(xx>1){u64 p=(xx<=LIM)?spf[xx]:xx; if(xx>LIM){p=0;for(u64 t=2;t*t<=xx;t++)if(xx%t==0){p=t;break;} if(!p)p=xx;} while(xx%p==0)xx/=p; u128 dd=d; while(dd%p==0){dd/=p;w*=p;} }
      // small if p*w < a
      if((u128)P*w < A) hist_small[b]++; else hist_big[b]++; }
    return; }
  u128 dd=d; for(int e=0;e<=2*fe[i];e++){ if(dd>QX) break; divrec(i+1,dd); dd*=fp[i]; }
}
int main(int argc,char**argv){
  P=strtoull(argv[1],0,10); Q=strtoull(argv[2],0,10); S=strtoull(argv[3],0,10);
  u64 lo=Q/P+1; if(lo<S) lo=S; u64 hi=(u64)(((u128)3*Q)/P); LIM=hi+1; 
  spf=calloc(LIM+1,4); for(u64 i=2;i<=LIM;i++) if(!spf[i]) for(u64 j=i;j<=LIM;j+=i) if(!spf[j]) spf[j]=i;
  for(u64 x=lo;x<=hi;x++){
    nf=0; addfac(Q); addfac(x); X=x; A=P*x-Q; QX=(u128)Q*x; u128 l=(u128)A*x; LO = l>QX ? l-QX : 0; // (a-q)x = a x - q x
    sol_x=0; divrec(0,1); if(sol_x) nx_with++;
  }
  printf("p=%llu q=%llu s=%llu range=%llu solutions=%llu x_with_solution=%llu\n",P,Q,S,hi-lo+1,nsol,nx_with);
  for(int b=0;b<70;b++) if(hist_a[b]) printf("a in [2^%d,2^%d): %llu  (small-w %llu, big-w %llu)\n",b,b+1,hist_a[b],hist_small[b],hist_big[b]);
  return 0;
}
