// egy2: count nondecreasing N-term unit-fraction representations of 1.
// Enumerate (N-3)-prefixes with remainder p/q; scan x=x_{N-2} over its range with a segmented
// sieve for factorizations; for each x count pairs (y,z), x<=y<=z, 1/y+1/z = (px-q)/(qx) via
// meet-in-the-middle over divisors of (qx)^2 in the residue class  d ≡ -qx (mod a), a=px-q,
// with (a-q)x <= d <= qx.  Also supports timing a single head over a sub-range.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
typedef unsigned long long u64; typedef unsigned __int128 u128; typedef uint32_t u32;
static u64 gcd64(u64 a,u64 b){while(b){u64 t=a%b;a=b;b=t;}return a;}
static u128 gcd128(u128 a,u128 b){while(b){u128 t=a%b;a=b;b=t;}return a;}
static inline u64 mulmod(u64 a,u64 b,u64 m){return (u64)(((u128)a*b)%m);}
static u64 invmod(u64 a,u64 m){ // a,m coprime; extended euclid
  __int128 t=0,nt=1; u64 r=m,nr=a%m; while(nr){u64 qq=r/nr; __int128 tmp=t-(__int128)qq*nt; t=nt; nt=tmp; u64 tr=r-qq*nr; r=nr; nr=tr;} if(t<0)t+=m; return (u64)t; }
// ---------- primes for sieving ----------
static u32 *primes; static int nprimes;
static void init_primes(u64 lim){ char*s=calloc(lim+1,1); primes=malloc(sizeof(u32)*(lim/2+10)); nprimes=0;
  for(u64 i=2;i<=lim;i++){ if(!s[i]){primes[nprimes++]=i; for(u64 j=i*i;j<=lim;j+=i)s[j]=1;} } free(s); }
// ---------- per-head state ----------
#define MAXQP 16
static u64 qpr[MAXQP]; static int qex[MAXQP]; static int nqp; // factorization of q
static u64 P,Q,S; static u64 total_count; static u64 nx_scanned;
// segmented sieve storage
#define SEG 262144
static u64 segbase; static u64 rem_[SEG]; // remaining cofactor after removing small primes
static u32 segfac[SEG][16]; static unsigned char segnf[SEG]; // small prime factors (distinct) per x, exponents recomputed
static void sieve_segment(u64 lo,u64 len){ // factor lo..lo+len-1 by primes up to sqrt
  segbase=lo; for(u64 i=0;i<len;i++){rem_[i]=lo+i; segnf[i]=0;}
  for(int k=0;k<nprimes;k++){ u64 pp=primes[k]; if(pp*pp>lo+len) break; u64 st=(lo+pp-1)/pp*pp; for(u64 j=st;j<lo+len;j+=pp){ u64 i=j-lo; if(segnf[i]<16) segfac[i][segnf[i]++]=pp; u64 r=rem_[i]; while(r%pp==0) r/=pp; rem_[i]=r; } }
}
// MITM count for one x
static u64 dbg_direct=0;
static u64 count_x(u64 x, u32*fac, int nf, u64 cof){
  // split x's primes into q-primes (exponents merge) and fresh primes
  u64 fresh_p[24]; int fresh_e[24]; int nfr=0; int xv[MAXQP]; for(int i=0;i<nqp;i++) xv[i]=0;
  u64 xx=x;
  for(int k=0;k<nf;k++){ u64 pp=fac[k]; int e=0; while(xx%pp==0){xx/=pp;e++;} int j; for(j=0;j<nqp;j++) if(qpr[j]==pp) break; if(j<nqp) xv[j]=e; else {fresh_p[nfr]=pp; fresh_e[nfr]=e; nfr++;} }
  if(cof>1){ // cof is prime (or 1) -- remaining cofactor after sieving by primes up to sqrt
    int j; for(j=0;j<nqp;j++) if(qpr[j]==cof) break; if(j<nqp) xv[j]+=1; else {fresh_p[nfr]=cof; fresh_e[nfr]=1; nfr++;} }
  // h = gcd(x,q) = prod qpr^min(xv,qex)
  u64 h=1; int hv[MAXQP]; for(int j=0;j<nqp;j++){ int m=xv[j]<qex[j]?xv[j]:qex[j]; hv[j]=m; for(int e=0;e<m;e++) h*=qpr[j]; }
  u64 q1=Q/h, x1=x/h; u64 a=P*x-Q; u64 a1=a/h; // a divisible by h
  // reduce a1/(h q1 x1): g2 = gcd(a1,h)
  int g2v[MAXQP]; u64 g2=1; for(int j=0;j<nqp;j++){ g2v[j]=0; if(hv[j]){ u64 t=a1; int e=0; while(e<hv[j] && t%qpr[j]==0){t/=qpr[j];e++;} g2v[j]=e; for(int k=0;k<e;k++) g2*=qpr[j]; } }
  u64 ap=a1/g2; u128 M=((u128)h*q1*x1)/g2;
  // interval for d: [max(1, ap*x - M), M]
  u128 hi=M; u128 lo=1; { u128 t=(u128)ap*x; if(t>M){ lo=t-M; if(lo<1) lo=1; } }
  if(lo>hi) return 0;
  u64 rho = (u64)(M % ap); rho = (ap-rho)%ap;
  u64 a1_unused=a1; (void)a1_unused;
  // exponents of M^2 for q-primes: 2*(max(qex,xv) - v(g2))
  int E[MAXQP]; for(int j=0;j<nqp;j++){ int Mx=qex[j]>xv[j]?qex[j]:xv[j]; E[j]=2*(Mx-g2v[j]); }
#define a1 ap
  // fresh divisors list (values and residues)
  static u128 vval[4096]; static u64 vres[4096]; int nv=1; vval[0]=1; vres[0]=1%a1;
  for(int k=0;k<nfr;k++){ int cur=nv; u128 pw=1; u64 pr=fresh_p[k]%a1, prw=1%a1; for(int e=1;e<=2*fresh_e[k];e++){ pw*=fresh_p[k]; prw=mulmod(prw,pr,a1); for(int i=0;i<cur;i++){ u128 nvv=vval[i]*pw; if(nvv>hi) continue; vval[nv]=nvv; vres[nv]=mulmod(vres[i],prw,a1); nv++; } } }
  // divisors of Qp: split primes into A (first part) and B
  // compute total tau
  u64 tauQ=1; for(int j=0;j<nqp;j++) tauQ*=(E[j]+1);
  // choose split: accumulate primes into A until |A| >= sqrt(tauQ*nv)
  double target = __builtin_sqrt((double)tauQ*(double)nv);
  static u128 aval[1<<16]; static u64 ares[1<<16]; int na=1; aval[0]=1; ares[0]=1%a1;
  static u128 bval[1<<12]; static u64 bres[1<<12]; int nb=1; bval[0]=1; bres[0]=1%a1;
  int j=0;
  for(;j<nqp;j++){ if(na>=target) break; int cur=na; u128 pw=1; u64 pr=qpr[j]%a1, prw=1%a1; for(int e=1;e<=E[j];e++){ pw*=qpr[j]; prw=mulmod(prw,pr,a1); for(int i=0;i<cur;i++){ u128 nvv=aval[i]*pw; if(nvv>hi) continue; aval[na]=nvv; ares[na]=mulmod(ares[i],prw,a1); na++; } } }
  for(;j<nqp;j++){ int cur=nb; u128 pw=1; u64 pr=qpr[j]%a1, prw=1%a1; for(int e=1;e<=E[j];e++){ pw*=qpr[j]; prw=mulmod(prw,pr,a1); for(int i=0;i<cur;i++){ u128 nvv=bval[i]*pw; if(nvv>hi) continue; bval[nb]=nvv; bres[nb]=mulmod(bres[i],prw,a1); nb++; } } }
  u64 cnt=0;
  if(a1 < 4096 || (u64)na*nb*nv < 64){ // direct enumeration
    dbg_direct++;
    for(int ia=0;ia<na;ia++) for(int ib=0;ib<nb;ib++){ u128 ab=aval[ia]*bval[ib]; if(ab>hi) continue; u64 rab=mulmod(ares[ia],bres[ib],a1); for(int iv=0;iv<nv;iv++){ u128 d=ab*vval[iv]; if(d>hi||d<lo) continue; if(mulmod(rab,vres[iv],a1)==rho) cnt++; } }
    return cnt;
  }
  // hash of A residues -> indices (open addressing)
  static int htab[1<<18]; static int hnext[1<<16]; int hbits=1; while((1<<hbits) < 2*na) hbits++; int hmask=(1<<hbits)-1;
  for(int i=0;i<=hmask;i++) htab[i]=-1;
  for(int ia=0;ia<na;ia++){ u64 key=ares[ia]; int slot=(int)((key*0x9E3779B97F4A7C15ULL)>>(64-hbits)); hnext[ia]=htab[slot]; htab[slot]=ia; }
  // for each (b,v): target = rho * inv(bres*vres); batch inversion
  int npairs=nb*nv; static u64 prod[1<<14], pref[1<<14];
  int idx=0; u64 acc=1;
  for(int ib=0;ib<nb;ib++) for(int iv=0;iv<nv;iv++){ u64 r=mulmod(bres[ib],vres[iv],a1); if(r==0) r=1; /* shouldn't happen (coprime) */ prod[idx]=r; acc=mulmod(acc,r,a1); pref[idx]=acc; idx++; }
  u64 inv=invmod(acc,a1);
  for(int k=npairs-1;k>=0;k--){ u64 invk = (k>0)? mulmod(inv,pref[k-1],a1) : inv; inv=mulmod(inv,prod[k],a1);
    int ib=k/nv, iv=k%nv; u64 tgt=mulmod(rho,invk,a1);
    int slot=(int)((tgt*0x9E3779B97F4A7C15ULL)>>(64-hbits));
    for(int ia=htab[slot]; ia!=-1; ia=hnext[ia]){ if(ares[ia]!=tgt) continue; u128 d=aval[ia]*bval[ib]; if(d>hi) continue; d*=vval[iv]; if(d>=lo && d<=hi) cnt++; }
  }
  return cnt;
}

#undef a1
static u64 brute_x(u64 x){ u128 qx=(u128)Q*x; u64 a=P*x-Q; u128 g=gcd128(a,qx); u64 ap=(u64)(a/g); u128 M=qx/g; u128 lo=(u128)ap*x; lo = lo>M? lo-M:1; u64 c=0; if(M>20000000) return (u64)-1;
  for(u128 d=1; d<=M; d++){ if(d<lo) continue; if((M*M)%d) continue; if((d+M)%ap==0) c++; } return c; }
static int DBG=0;
// scan a head over x in [xlo,xhi]
static u64 scan_head(u64 p,u64 q,u64 xlo,u64 xhi){
  P=p;Q=q; // factor q
  nqp=0; u64 r=q; for(int k=0;k<nprimes && (u64)primes[k]*primes[k]<=r;k++){ u64 pp=primes[k]; if(r%pp==0){ int e=0; while(r%pp==0){r/=pp;e++;} qpr[nqp]=pp;qex[nqp]=e;nqp++; } } if(r>1){qpr[nqp]=r;qex[nqp]=1;nqp++;}
  u64 cnt=0;
  for(u64 lo=xlo; lo<=xhi; lo+=SEG){ u64 len=xhi-lo+1; if(len>SEG) len=SEG; sieve_segment(lo,len);
    for(u64 i=0;i<len;i++){ u64 x=lo+i; u64 c=count_x(x,segfac[i],segnf[i],rem_[i]); if(DBG){u64 b=brute_x(x); if(b!=(u64)-1 && b!=c) printf("MISMATCH p=%llu q=%llu x=%llu mitm=%llu brute=%llu\n",P,Q,x,c,b);} cnt+=c; nx_scanned++; } }
  return cnt;
}
// prefix enumeration
static int N;
static void rec(int k,u64 p,u64 q,u64 last){
  int rem=N-k;
  if(rem==3){ u64 lo=q/p+1; if(lo<last) lo=last; u64 hi=(u64)(((u128)3*q)/p); if(hi<lo) return; total_count+=scan_head(p,q,lo,hi); return; }
  u64 lo=q/p+1; if(lo<last) lo=last; u64 hi=(u64)(((u128)rem*q)/p);
  for(u64 x=lo;x<=hi;x++){ u128 np_=(u128)p*x-q, nq=(u128)q*x; u128 g=gcd128(np_,nq); np_/=g; nq/=g; if(nq>>64){fprintf(stderr,"q overflow\n");exit(1);} rec(k+1,(u64)np_,(u64)nq,x); }
}
int main(int argc,char**argv){
  if(argc>=2 && strcmp(argv[1],"head")==0){ // head p q xlo xhi
    u64 p=strtoull(argv[2],0,10), q=strtoull(argv[3],0,10), xlo=strtoull(argv[4],0,10), xhi=strtoull(argv[5],0,10);
    init_primes(6000000);
    struct timespec t0,t1; clock_gettime(CLOCK_MONOTONIC,&t0);
    u64 c=scan_head(p,q,xlo,xhi);
    clock_gettime(CLOCK_MONOTONIC,&t1); double dt=(t1.tv_sec-t0.tv_sec)+1e-9*(t1.tv_nsec-t0.tv_nsec);
    printf("p=%llu q=%llu x in [%llu,%llu]: count=%llu  xs=%llu  time=%.3fs  %.3f us/x  (direct=%llu)\n",p,q,xlo,xhi,c,nx_scanned,dt,1e6*dt/nx_scanned,dbg_direct);
    return 0; }
  N=atoi(argv[1]); if(argc>2) DBG=1; init_primes(6000000);
  struct timespec t0,t1; clock_gettime(CLOCK_MONOTONIC,&t0);
  rec(0,1,1,2);
  clock_gettime(CLOCK_MONOTONIC,&t1); double dt=(t1.tv_sec-t0.tv_sec)+1e-9*(t1.tv_nsec-t0.tv_nsec);
  printf("a(%d) = %llu   heads=%llu  time=%.2fs  %.3f us/head\n",N,total_count,nx_scanned,dt,1e6*dt/nx_scanned);
  return 0;
}
