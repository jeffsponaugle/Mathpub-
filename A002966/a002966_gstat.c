// For all solutions of 1/x+1/y+1/z=p/q (s<=x<=y<=z): histogram of x1 = x/gcd(x,y) and of D/a where D=ab-q^2.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef unsigned long long u64; typedef unsigned __int128 u128;
static u64 P,Q,S; static uint32_t *spf; static u64 LIM;
static u64 fp[40]; static int fe[40]; static int nf;
static void addfac(u64 x){ while(x>1){u64 p=spf[x]; int e=0; while(x%p==0){x/=p;e++;} int i; for(i=0;i<nf;i++) if(fp[i]==p)break; if(i==nf){fp[nf]=p;fe[nf]=0;nf++;} fe[i]+=e; } }
static u128 gcd128(u128 a,u128 b){while(b){u128 t=a%b;a=b;b=t;}return a;}
static u128 QX, LO; static u64 A, X; static u64 nsol; static u64 hx1[80], hDa[80], hx1_biga[80]; static u64 cnt_big_a;
static void divrec(int i,u128 d){
  if(i==nf){ if(d<=QX && d>=LO && (u64)(d%A)==(u64)((A-(QX%A))%A)){ nsol++;
      u128 y=(d+QX)/A; u128 g=gcd128(X,y); u128 x1=X/g; int b=0; while(x1>1){x1>>=1;b++;} hx1[b]++; if(A>Q){hx1_biga[b]++;cnt_big_a++;}
      u128 D=(u128)P*d; u128 r=D/A; int c=0; while(r){r>>=1;c++;} hDa[c]++; }
    return; }
  u128 dd=d; for(int e=0;e<=2*fe[i];e++){ if(dd>QX) break; divrec(i+1,dd); dd*=fp[i]; }
}
int main(int argc,char**argv){
  P=strtoull(argv[1],0,10); Q=strtoull(argv[2],0,10); S=strtoull(argv[3],0,10);
  u64 lo=Q/P+1; if(lo<S) lo=S; u64 hi=(u64)(((u128)3*Q)/P); LIM=hi+1;
  spf=calloc(LIM+1,4); for(u64 i=2;i<=LIM;i++) if(!spf[i]) for(u64 j=i;j<=LIM;j+=i) if(!spf[j]) spf[j]=i;
  for(u64 x=lo;x<=hi;x++){ A=P*x-Q; nf=0; addfac(Q); addfac(x); X=x; QX=(u128)Q*x; u128 l=(u128)A*x; LO = l>QX ? l-QX : 0; divrec(0,1); }
  printf("solutions=%llu (with a>q: %llu)\n",nsol,cnt_big_a);
  printf("x1=x/gcd(x,y) histogram (log2 bins): all | a>q\n"); for(int b=0;b<80;b++) if(hx1[b]) printf("  x1 in [2^%d,2^%d): %llu | %llu\n",b,b+1,hx1[b],hx1_biga[b]);
  printf("D/a histogram (bits):\n"); for(int c=0;c<80;c++) if(hDa[c]) printf("  D/a < 2^%d: %llu\n",c,hDa[c]);
  return 0; }
