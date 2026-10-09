/*
 * a171740c.cu — CUDA production worker for the A171740 MITM search.
 * Direct port of a171740g (Metal): same CLI, stage semantics, banners,
 * certificate lines, exact -H resume (slab counters are contiguous).
 *
 * build: nvcc -O3 -arch=native -o a171740c a171740c.cu -lpthread
 * usage: ./a171740c -n N [-B b2start] [-E b2end] [-L limit]
 *                   [-H done [-P b1]] [-T tablecap] [-S slab] [-q]
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <ctime>
#include <unistd.h>
#include <pthread.h>
#include <atomic>
#include <thread>
#include <vector>
#include <cuda_runtime.h>
#include <algorithm>

typedef unsigned __int128 u128;
typedef uint64_t u64;
typedef uint32_t u32;

/* ---------------- host w192 ---------------- */
struct w192 { u64 d[3]; };
static w192 w_from_u128(u128 x){ w192 r; r.d[0]=(u64)x; r.d[1]=(u64)(x>>64); r.d[2]=0; return r; }
static w192 w_max_(){ w192 r; r.d[0]=r.d[1]=r.d[2]=~0ULL; return r; }
static int w_is_max(w192 a){ return a.d[0]==~0ULL&&a.d[1]==~0ULL&&a.d[2]==~0ULL; }
static int w_cmp(w192 a, w192 b){
    for (int i=2;i>=0;i--){ if(a.d[i]<b.d[i])return -1; if(a.d[i]>b.d[i])return 1; }
    return 0;
}
static int w_is_zero(w192 a){ return !(a.d[0]|a.d[1]|a.d[2]); }
static w192 w_add(w192 a, w192 b){
    w192 r; unsigned c=0;
    for (int i=0;i<3;i++){ u64 s=a.d[i]+b.d[i]; u64 s2=s+c; c=(s<a.d[i])||(s2<s); r.d[i]=s2; }
    return r;
}
static w192 w_addu(w192 a, u64 b){ w192 r=a; r.d[0]+=b; if(r.d[0]<b){ if(++r.d[1]==0) ++r.d[2]; } return r; }
static w192 w_mulu_ov(w192 a, u64 m, int *ov){
    u128 t0=(u128)a.d[0]*m, t1=(u128)a.d[1]*m+(u64)(t0>>64), t2=(u128)a.d[2]*m+(u64)(t1>>64);
    if (t2>>64) *ov=1;
    w192 r; r.d[0]=(u64)t0; r.d[1]=(u64)t1; r.d[2]=(u64)t2; return r;
}
static w192 w_mulu(w192 a, u64 m){ int ov=0; return w_mulu_ov(a,m,&ov); }
static w192 w_pow_sat(u64 b, int e){
    w192 r=w_from_u128(1);
    while (e-->0){ int ov=0; r=w_mulu_ov(r,b,&ov); if(ov) return w_max_(); }
    return r;
}
static unsigned w_divmod_small(w192 *v, unsigned d){
    u64 r=0;
    for (int i=2;i>=0;i--){
        u64 limb=v->d[i];
        u64 a=(r<<32)|(limb>>32); u64 q1=a/d; r=a%d;
        u64 b=(r<<32)|(limb&0xffffffffu); u64 q0=b/d; r=b%d;
        v->d[i]=(q1<<32)|q0;
    }
    return (unsigned)r;
}
static u64 w_div_big(w192 num, w192 den, w192 *rem){
    w192 r; r.d[0]=r.d[1]=r.d[2]=0; u64 q=0;
    for (int i=191;i>=0;i--){
        r.d[2]=(r.d[2]<<1)|(r.d[1]>>63); r.d[1]=(r.d[1]<<1)|(r.d[0]>>63);
        r.d[0]=(r.d[0]<<1)|((num.d[i>>6]>>(i&63))&1);
        if (w_cmp(r,den)>=0){
            unsigned br=0; w192 nr;
            for(int k=0;k<3;k++){ u64 s=r.d[k]-den.d[k]; u64 s2=s-br;
                br=(r.d[k]<den.d[k])||(s<br); nr.d[k]=s2; }
            r=nr; if(i<64) q|=1ULL<<i;
        }
    }
    if (rem) *rem=r;
    return q;
}
static u64 w_mod_u64(w192 v, u64 d){
    u64 r=0;
    for (int i=191;i>=0;i--){
        int hib=(int)(r>>63);
        r=(r<<1)|((v.d[i>>6]>>(i&63))&1);
        if (hib||r>=d) r-=d;
    }
    return r;
}
static char *w_str(w192 v, char *buf){
    char *p=buf+63; *--p=0; if(w_is_zero(v)) *--p='0';
    while(!w_is_zero(v)){ unsigned r=w_divmod_small(&v,1000000000u);
        for(int i=0;i<9&&!(w_is_zero(v)&&r==0);i++){ *--p=(char)('0'+r%10); r/=10; } }
    return p;
}
static w192 w_parse(const char *s){
    w192 v; v.d[0]=v.d[1]=v.d[2]=0;
    if(!*s){ fprintf(stderr,"empty number\n"); exit(1); }
    for(;*s;s++){ if(*s<'0'||*s>'9'){ fprintf(stderr,"bad number '%s'\n",s); exit(1); }
        v=w_addu(w_mulu(v,10),(u64)(*s-'0')); }
    return v;
}
static u64 pw64(u64 b,int e){ u64 r=1; while(e-->0) r*=b; return r; }
static u64 revdig(u64 x, unsigned b, int d){ u64 r=0; for(int i=0;i<d;i++){ r=r*b+x%b; x/=b; } return r; }

/* ---------------- globals / options ---------------- */
static int g_n, g_kdig, g_odd, g_m;
static u64 g_tablecap = 300000000ULL;
static int g_quiet = 0;
static w192 g_limitW; static int g_limit_set = 0;
static u64 g_skip = 0; static unsigned g_skip_b1 = 0;
static int g_nbuild = 0;
static u64 g_slab = 268435456ULL;
static u64 g_qcap = 0;            /* -J: cap on Q (bucket count); 0 = table-size rule */
static w192 g_best; static int g_have = 0; static unsigned g_bb1, g_bb2;

struct MCfg {
    unsigned b1, b2;
    int s, t, j, je, K;
    u64 Q, QF, c1f, c2f, c3f;
    w192 PM, PS, PMT, PB1NJE, hi_p;
    u64 yh_start, yh_end;
};

static int mcfg_init(MCfg *M, unsigned b2, unsigned b1, w192 p_hiW, w192 stage_hi){
    memset(M,0,sizeof *M);
    M->b1=b1; M->b2=b2;
    int m=g_m;
    M->K = g_odd ? g_n/2 : m;
    int t=0; u64 e=1;
    while (t+1<=m/2 && e*b2<=g_tablecap){ e*=b2; t++; }
    M->t=t; M->s=m-t;
    int j=0; u64 q=1;
    u64 qlim = (g_qcap && g_qcap < e) ? g_qcap : e;
    while (q*b1<=qlim){ q*=b1; j++; }
    if (j<1){ j=1; q=b1; }
    M->j=j; M->Q=q;
    M->PMT=w_pow_sat(b2, M->K+t);
    if (w_is_max(M->PMT)) return -1;
    int je=j; u64 qf=q;
    while (je+1<=g_kdig && je+1<=g_n-1 &&
           (u128)qf*b1<=((u64)1<<40) && (qf/q)*(u64)b1<=((u64)1<<30)){
        w192 pnje=w_pow_sat(b1, g_n-(je+1));
        if (w_is_max(pnje)) break;
        if (w_div_big(M->PMT, pnje, NULL) > 8) break;
        je++; qf*=b1;
    }
    M->je=je; M->QF=qf;
    M->PM=w_pow_sat(b2, M->K);
    M->PS=w_pow_sat(b2, M->s);
    M->PB1NJE=w_pow_sat(b1, g_n-M->je);
    if (w_is_max(M->PB1NJE)) return -1;
    M->hi_p=p_hiW; if (w_cmp(stage_hi,M->hi_p)<0) M->hi_p=stage_hi;
    M->c1f=w_mod_u64(M->PM, M->QF);
    M->c2f=w_mod_u64(M->PS, M->QF);
    M->c3f=w_mod_u64(M->PMT, M->QF);
    M->yh_start=pw64(b2, M->s-1);
    M->yh_end=w_div_big(M->hi_p, M->PMT, NULL)+1;
    u64 cap=pw64(b2, M->s);
    if (M->yh_end>cap) M->yh_end=cap;
    return 0;
}

/* ---------------- parallel table build (into managed memory) ---------------- */
static u64 *g_off; static u64 *g_ent; static u64 g_Nyl;
/* entry packing: yl in high 34 bits, second-level residue in low 30
 * (je selection caps bhi < 2^30); off is u64 (counts exceed u32 at t=7) */
/* memory-lean build: counts and fill cursors live in the (managed) off
 * array itself; after the fill pass off[k] holds each bucket's END, so
 * the kernel reads bucket k as [k ? off[k-1] : 0, off[k]). */
struct BArg { const MCfg *M; u64 lo, hi; std::atomic<u64> *offA; int pass; };
static void *build_worker(void *ap){
    BArg *A=(BArg*)ap; const MCfg *M=A->M;
    for (u64 yl=A->lo; yl<A->hi; yl++){
        u64 r = g_odd ? revdig(yl/M->b2, M->b2, M->t-1) : revdig(yl, M->b2, M->t);
        u64 Bf = (u64)(((u128)yl*M->c1f + (u128)r*M->c2f) % M->QF);
        u64 key = Bf % M->Q;
        if (A->pass==1){
            A->offA[key].fetch_add(1,std::memory_order_relaxed);
        } else {
            u64 slot=A->offA[key].fetch_add(1,std::memory_order_relaxed);
            g_ent[slot]=(yl<<30)|(Bf/M->Q);
        }
    }
    return NULL;
}
static int build_table(const MCfg *M){
    g_Nyl=pw64(M->b2, M->t);
    if (cudaMallocManaged(&g_off,(M->Q+1)*8)!=cudaSuccess) return -1;
    memset(g_off,0,(M->Q+1)*8);
    if (cudaMallocManaged(&g_ent,g_Nyl*8)!=cudaSuccess) return -1;
    std::atomic<u64> *offA=(std::atomic<u64>*)g_off;
    int nth=g_nbuild;
    pthread_t tids[64]; BArg args[64];
    for (int pass=1; pass<=2; pass++){
        if (pass==2){
            /* turn counts into exclusive-prefix cursors in place */
            u64 run=0;
            for (u64 q=0;q<M->Q;q++){
                u64 c=g_off[q];
                g_off[q]=run;
                run+=c;
            }
            g_off[M->Q]=run;
        }
        u64 chunk=(g_Nyl+nth-1)/nth;
        for (int i=0;i<nth;i++){
            args[i].M=M; args[i].lo=i*chunk;
            args[i].hi=(u64)(i+1)*chunk>g_Nyl?g_Nyl:(u64)(i+1)*chunk;
            args[i].offA=offA; args[i].pass=pass;
            pthread_create(&tids[i],NULL,build_worker,&args[i]);
        }
        for (int i=0;i<nth;i++) pthread_join(tids[i],NULL);
    }
    /* off now holds bucket ENDS (cursor final positions) */
    /* pass 3: sort each bucket by its low-30-bit residue so the sweep
     * can binary-search needhi instead of scanning — this defuses the
     * degenerate mega-buckets of shared-factor base pairs */
    {
        u64 chunk=(M->Q+nth-1)/nth;
        std::atomic<int> tix(0);
        auto sorter=[&](){
            int i=tix.fetch_add(1);
            u64 qlo=(u64)i*chunk, qhi=(u64)(i+1)*chunk; if(qhi>M->Q) qhi=M->Q;
            for (u64 q=qlo;q<qhi;q++){
                u64 s=q?g_off[q-1]:0, e=g_off[q];
                if (e>s+1) std::sort(g_ent+s, g_ent+e,
                    [](u64 a,u64 b){ return (a&0x3fffffffu)<(b&0x3fffffffu); });
            }
        };
        std::vector<std::thread> th;
        for (int i=0;i<nth;i++) th.emplace_back(sorter);
        for (auto &t:th) t.join();
    }
    return 0;
}

/* ---------------- exact host verification ---------------- */
static w192 vof(const MCfg *M, u64 yh, u64 yl){
    u64 rs=revdig(yh,M->b2,M->s);
    u64 rt=g_odd?revdig(yl/M->b2,M->b2,M->t-1):revdig(yl,M->b2,M->t);
    w192 v=w_mulu(M->PMT, yh);
    v=w_add(v, w_mulu(M->PM, yl));
    v=w_add(v, w_mulu(M->PS, rt));
    return w_addu(v, rs);
}
static int fullpal_host(unsigned b1, w192 v){
    w192 w=v; u128 rr=0;
    for (int i=0;i<g_kdig;i++) rr=rr*b1 + w_divmod_small(&w, b1);
    if (g_odd) w_divmod_small(&w, b1);
    return w.d[2]==0 && ((((u128)w.d[1]<<64)|w.d[0])==rr);
}

/* ---------------- device ---------------- */
struct GpuParams {
    u32 b1,b2,s,t,j,je,kdig,odd;
    u64 Q,QF,c3f;
    u64 Mb1,Mb2,MQ,MQF;
    u64 PM0,PM1,PM2, PS0,PS1,PS2, PT0,PT1,PT2, PB0,PB1x,PB2, ST0,ST1,ST2;
    u64 yh0,yhEnd,G; u32 outCap; u32 pad;
};

struct W { u64 d0,d1,d2; };
__device__ __forceinline__ u64 divqd(u64 x, u64 M, u64 d, u64 &r){
    u64 q=__umul64hi(x,M); r=x-q*d; if(r>=d){r-=d;q++;} return q;
}
__device__ __forceinline__ W mkw(u64 a,u64 b,u64 c){ W r; r.d0=a; r.d1=b; r.d2=c; return r; }
__device__ __forceinline__ int wcmpd(W a, W b){
    if (a.d2!=b.d2) return a.d2<b.d2?-1:1;
    if (a.d1!=b.d1) return a.d1<b.d1?-1:1;
    if (a.d0!=b.d0) return a.d0<b.d0?-1:1; return 0;
}
__device__ __forceinline__ W waddd(W a, W b){
    W r; r.d0=a.d0+b.d0; u64 c0=r.d0<a.d0?1:0;
    u64 s1=a.d1+b.d1; u64 c1=s1<a.d1?1:0; r.d1=s1+c0; if(r.d1<s1) c1++;
    r.d2=a.d2+b.d2+c1; return r;
}
__device__ __forceinline__ W wsubd(W a, W b){
    W r; r.d0=a.d0-b.d0; u64 br0=a.d0<b.d0?1:0;
    u64 s1=a.d1-b.d1; u64 br1=a.d1<b.d1?1:0; r.d1=s1-br0; if(s1<br0) br1++;
    r.d2=a.d2-b.d2-br1; return r;
}
__device__ __forceinline__ W wsub1d(W a){ W r=a; if(r.d0--==0){ if(r.d1--==0) r.d2--; } return r; }
__device__ __forceinline__ W waddud(W a, u64 b){ W r=a; r.d0+=b; if(r.d0<b){ if(++r.d1==0) ++r.d2; } return r; }
__device__ __forceinline__ W wmulud(W a, u64 m){
    u64 l0=a.d0*m, h0=__umul64hi(a.d0,m);
    u64 l1=a.d1*m, h1=__umul64hi(a.d1,m);
    u64 l2=a.d2*m;
    W r; r.d0=l0; r.d1=l1+h0; u64 c=(r.d1<l1)?1:0; r.d2=l2+h1+c; return r;
}
__device__ u64 wdivbigd(W num, W den, W &rem){
    W r=mkw(0,0,0); u64 q=0;
    for (int i=191;i>=0;i--){
        r.d2=(r.d2<<1)|(r.d1>>63); r.d1=(r.d1<<1)|(r.d0>>63);
        u64 bit=(i>=128)?((num.d2>>(i-128))&1):(i>=64)?((num.d1>>(i-64))&1):((num.d0>>i)&1);
        r.d0=(r.d0<<1)|bit;
        if (wcmpd(r,den)>=0){ r=wsubd(r,den); if(i<64) q|=(u64)1<<i; }
    }
    rem=r; return q;
}
__device__ u64 mod128d(u64 hi, u64 lo, u64 m){
    u64 r=0;
    for (int i=127;i>=0;i--){
        u64 top=r>>63; r=(r<<1)|((i>=64)?((hi>>(i-64))&1):((lo>>i)&1));
        if (top||r>=m) r-=m;
    }
    return r;
}
__device__ __forceinline__ u64 revmd(u64 x, u32 b, u32 d, u64 Mb){
    u64 r=0;
    for(u32 i=0;i<d;i++){ u64 rem; x=divqd(x,Mb,b,rem); r=r*b+rem; }
    return r;
}
__device__ __forceinline__ void wdm3(W &w, u32 b1, u64 Mb1, u64 &dg){
    u64 r64=0, limb;
    limb=w.d2; { u64 a=(r64<<32)|(limb>>32); u64 q1=divqd(a,Mb1,b1,r64); u64 bl=(r64<<32)|(limb&0xffffffffu); u64 q0=divqd(bl,Mb1,b1,r64); w.d2=(q1<<32)|q0; }
    limb=w.d1; { u64 a=(r64<<32)|(limb>>32); u64 q1=divqd(a,Mb1,b1,r64); u64 bl=(r64<<32)|(limb&0xffffffffu); u64 q0=divqd(bl,Mb1,b1,r64); w.d1=(q1<<32)|q0; }
    limb=w.d0; { u64 a=(r64<<32)|(limb>>32); u64 q1=divqd(a,Mb1,b1,r64); u64 bl=(r64<<32)|(limb&0xffffffffu); u64 q0=divqd(bl,Mb1,b1,r64); w.d0=(q1<<32)|q0; }
    dg=r64;
}

__global__ void sweep(const u64 *__restrict__ off, const u64 *__restrict__ ent,
                      GpuParams C, u32 *outN, u64 *outP, u32 *stats){
    u64 tid = (u64)blockIdx.x*blockDim.x + threadIdx.x;
    u64 yhb = C.yh0 + tid*C.G;
    if (yhb >= C.yhEnd) return;
    u64 yhe = yhb+C.G; if (yhe>C.yhEnd) yhe=C.yhEnd;
    W PMT=mkw(C.PT0,C.PT1,C.PT2), DEN=mkw(C.PB0,C.PB1x,C.PB2);
    W PMw=mkw(C.PM0,C.PM1,C.PM2), PSw=mkw(C.PS0,C.PS1,C.PS2), STOP=mkw(C.ST0,C.ST1,C.ST2);
    W vbase=wmulud(PMT,yhb);
    if (wcmpd(vbase,STOP)>=0) return;
    W prem; u64 p=wdivbigd(vbase,DEN,prem);
    u64 a0; divqd(yhb,C.MQF,C.QF,a0);
    u64 term3=mod128d(__umul64hi(a0,C.c3f), a0*C.c3f, C.QF);
    u32 recon=0, fulls=0;
    for (u64 yh=yhb; yh<yhe; yh++){
        u64 rs=revmd(yh,C.b2,C.s,C.Mb2);
        u64 rsm; divqd(rs,C.MQF,C.QF,rsm);
        u64 Af=term3+rsm; if(Af>=C.QF) Af-=C.QF;
        W vtop=wsub1d(waddd(vbase,PMT));
        if (wcmpd(vtop,STOP)>=0) vtop=wsub1d(STOP);
        for (u64 pp=p;;pp++){
            W plo=wmulud(DEN,pp);
            if (wcmpd(plo,vtop)>0) break;
            W phi=waddd(plo,DEN);
            u64 rjf=revmd(pp,C.b1,C.je,C.Mb1);
            u64 R=rjf>=Af?rjf-Af:rjf+C.QF-Af;
            u64 keyr; u64 hiq=divqd(R,C.MQ,C.Q,keyr);
            u64 needhi=hiq;
            u64 i0=keyr?off[keyr-1]:0, i1=off[keyr];
            /* buckets are sorted by low-30 residue: binary search needhi */
            { u64 blo=i0, bhi2=i1;
              while (blo<bhi2){ u64 mid=(blo+bhi2)>>1;
                  if ((ent[mid]&0x3fffffffu) < needhi) blo=mid+1; else bhi2=mid; }
              i0=blo; }
            for (u64 idx=i0; idx<i1; idx++){
                u64 en=ent[idx];
                if ((en&0x3fffffffu)!=needhi) break;
                u64 yl=en>>30;
                recon++;
                u64 xr=yl; if(C.odd){ u64 rr0; xr=divqd(yl,C.Mb2,C.b2,rr0); }
                u32 td=C.odd? C.t-1 : C.t;
                u64 rt=revmd(xr,C.b2,td,C.Mb2);
                W v=waddd(waddd(vbase,wmulud(PMw,yl)),wmulud(PSw,rt)); v=waddud(v,rs);
                if (wcmpd(v,plo)<0 || wcmpd(v,phi)>=0) continue;
                if (wcmpd(v,STOP)>=0) continue;
                fulls++;
                W w=v; u64 rhi=0, rlo=0;
                for (u32 i=0;i<C.kdig;i++){
                    u64 dg; wdm3(w,C.b1,C.Mb1,dg);
                    u64 nl=rlo*C.b1+dg; u64 nh=rhi*C.b1+__umul64hi(rlo,(u64)C.b1);
                    if (nl<dg) nh++;
                    rlo=nl; rhi=nh;
                }
                if (C.odd){ u64 dg; wdm3(w,C.b1,C.Mb1,dg); }
                if ((w.d2==0)&&(w.d1==rhi)&&(w.d0==rlo)){
                    u32 slot=atomicAdd(outN,1u);
                    if (slot<C.outCap){ outP[2*slot]=yh; outP[2*slot+1]=yl; }
                }
            }
        }
        vbase=waddd(vbase,PMT);
        prem=waddd(prem,PMT);
        while (wcmpd(prem,DEN)>=0){ prem=wsubd(prem,DEN); p++; }
        term3+=C.c3f; if(term3>=C.QF) term3-=C.QF;
        if (wcmpd(vbase,STOP)>=0) break;
    }
    if (recon) atomicAdd(&stats[0],recon);
    if (fulls) atomicAdd(&stats[1],fulls);
}

static double nowsec(void){ struct timespec ts; clock_gettime(CLOCK_MONOTONIC,&ts);
    return (double)ts.tv_sec+1e-9*ts.tv_nsec; }

static void report_match(w192 v, unsigned b1, unsigned b2){
    if (!g_have || w_cmp(v,g_best)<0){
        g_best=v; g_have=1; g_bb1=b1; g_bb2=b2;
        char buf[64];
        fprintf(stderr,"\n*** MATCH: %s is an %d-digit palindrome in bases %u and %u\n",
                w_str(v,buf), g_n, b1, b2);
    }
}

static const char *DIGCH="0123456789abcdefghijklmnopqrstuvwxyz";
static void print_all_bases(w192 v){
    unsigned d[200];
    for (unsigned b=2;b<=128;b++){
        w192 plo=w_pow_sat(b,g_n-1), phi=w_pow_sat(b,g_n);
        if (w_cmp(v,plo)<0||w_cmp(v,phi)>=0) continue;
        int nd=0; w192 x=v;
        while(!w_is_zero(x)&&nd<200) d[nd++]=w_divmod_small(&x,b);
        if (nd!=g_n) continue;
        int pal=1;
        for(int i=0,j2=nd-1;i<j2;i++,j2--) if(d[i]!=d[j2]){pal=0;break;}
        if(!pal) continue;
        printf("  base %-3u: ",b);
        for(int i=nd-1;i>=0;i--){ if(b<=36) putchar(DIGCH[d[i]]); else printf("%u%s",d[i],i?".":""); }
        printf("\n");
    }
}

#define CUCHECK(x) do{ cudaError_t e=(x); if(e!=cudaSuccess){ \
    fprintf(stderr,"CUDA error %s at line %d\n",cudaGetErrorString(e),__LINE__); exit(1);} }while(0)

int main(int argc, char **argv){
    int opt; unsigned b2start=3, b2end=0;
    while ((opt=getopt(argc,argv,"n:B:E:L:H:P:T:S:J:qh"))!=-1){
        switch(opt){
        case 'n': g_n=atoi(optarg); break;
        case 'B': b2start=(unsigned)atoi(optarg); break;
        case 'E': b2end=(unsigned)atoi(optarg); break;
        case 'L': g_limitW=w_parse(optarg); g_limit_set=1; break;
        case 'H': g_skip=strtoull(optarg,0,10); break;
        case 'P': g_skip_b1=(unsigned)atoi(optarg); break;
        case 'T': g_tablecap=strtoull(optarg,0,10); break;
        case 'S': g_slab=strtoull(optarg,0,10); break;
        case 'J': g_qcap=strtoull(optarg,0,10); break;
        case 'q': g_quiet=1; break;
        default:
            fprintf(stderr,"usage: %s -n digits [-B b2start] [-E b2end] [-L limit] "
                    "[-H done [-P b1]] [-T tablecap] [-S slab] [-q]\n", argv[0]);
            return 2;
        }
    }
    if (g_n<8||g_n>40){ fprintf(stderr,"need -n in 8..40\n"); return 2; }
    if (g_skip && b2start==3){ fprintf(stderr,"-H requires -B\n"); return 2; }
    g_kdig=g_n/2; g_odd=g_n&1; g_m=(g_n+1)/2;
    g_nbuild=(int)sysconf(_SC_NPROCESSORS_ONLN); if(g_nbuild<1)g_nbuild=1; if(g_nbuild>64)g_nbuild=64;

    cudaDeviceProp prop; CUCHECK(cudaGetDeviceProperties(&prop,0));
    char lb[64];
    fprintf(stderr,"A171740 CUDA search: n=%d, device=%s (%d SMs), build threads=%d, limit=%s\n",
            g_n, prop.name, prop.multiProcessorCount, g_nbuild,
            g_limit_set?w_str(g_limitW,lb):"none");

    u32 *dN; u64 *dP; u32 *dS;
    u32 outCap=1<<16;
    CUCHECK(cudaMallocManaged(&dN,4));
    CUCHECK(cudaMallocManaged(&dP,(size_t)outCap*16));
    CUCHECK(cudaMallocManaged(&dS,8));
    dS[0]=dS[1]=0;

    double t_start=nowsec();
    u64 total_done=0;

    for (unsigned b2=b2start;; b2++){
        if (b2end && b2>=b2end){
            fprintf(stderr,"reached -E %u: stopping before stage b2=%u\n",b2end,b2);
            if (!g_limit_set){ g_limitW=w_pow_sat(b2,g_n-1); g_limit_set=1; }
            break;
        }
        w192 bound = g_have ? g_best : w_max_();
        if (g_limit_set && w_cmp(g_limitW,bound)<0) bound=g_limitW;
        w192 loW=w_pow_sat(b2,g_n-1);
        if (w_is_max(loW)||w_cmp(loW,bound)>=0) break;
        unsigned b1s[24]; w192 phis[24]; int nb1=0;
        w192 himax; himax.d[0]=himax.d[1]=himax.d[2]=0;
        for (unsigned b1=2;b1<b2;b1++){
            if (!(g_n&1) && b1==b2-1) continue;
            w192 h=w_pow_sat(b1,g_n);
            if (w_cmp(h,loW)>0 && nb1<24){
                b1s[nb1]=b1; phis[nb1]=h; nb1++;
                if (w_cmp(h,himax)>0) himax=h;
            }
        }
        if (!nb1) continue;
        w192 hiW=himax;
        if (g_limit_set && w_cmp(g_limitW,hiW)<0) hiW=g_limitW;
        {
            char a[64],b[64]; char pl[128]; size_t o=0; pl[0]=0;
            for (int i=0;i<nb1;i++) o+=(size_t)snprintf(pl+o,sizeof pl-o,"%s%u",i?",":"",b1s[i]);
            fprintf(stderr,"stage b2=%u vs b1={%s}: values [%s, %s) (CUDA MITM)\n",
                    b2, pl, w_str(loW,a), w_str(hiW,b));
        }
        if (g_skip && g_skip_b1){
            int found=0;
            for (int i=0;i<nb1;i++) if (b1s[i]==g_skip_b1) found=1;
            if (!found){ fprintf(stderr,"warning: -P %u not a partner of b2=%u; ignoring -H\n",g_skip_b1,b2); g_skip=0; }
        }
        for (int pi=0; pi<nb1; pi++){
            MCfg M;
            if (g_skip){
                unsigned target=g_skip_b1?g_skip_b1:b1s[0];
                if (b1s[pi]<target){
                    fprintf(stderr,"  MITM b1=%u: skipped (resume asserts completed)\n",b1s[pi]);
                    continue;
                }
            }
            if (mcfg_init(&M,b2,b1s[pi],phis[pi],hiW)){
                fprintf(stderr,"  MITM b1=%u: exceeds w192 range, skipping\n",b1s[pi]);
                continue;
            }
            if (M.yh_end<=M.yh_start) continue;
            w192 stop0 = g_have ? g_best : w_max_();
            if (w_cmp(M.hi_p,stop0)<0) stop0=M.hi_p;
            if (w_cmp(w_mulu(M.PMT,M.yh_start),stop0)>=0) continue;

            u64 skip0=0;
            if (g_skip){
                skip0=g_skip; g_skip=0;
                u64 span=M.yh_end-M.yh_start;
                if (skip0>span) skip0=span;
                fprintf(stderr,"  resuming partner b1=%u at +%llu of %llu yh (exact)\n",
                        M.b1,(unsigned long long)skip0,(unsigned long long)span);
            }
            fprintf(stderr,"  MITM b1=%u: s=%d t=%d j=%d je=%d Q=%llu table=%llu yh=%llu\n",
                    M.b1,M.s,M.t,M.j,M.je,(unsigned long long)M.Q,
                    (unsigned long long)pw64(M.b2,M.t),
                    (unsigned long long)(M.yh_end-M.yh_start));
            double tb=nowsec();
            if (build_table(&M)) return 1;
            fprintf(stderr,"  table built in %.1fs (%d threads)\n",nowsec()-tb,g_nbuild);

            GpuParams P; memset(&P,0,sizeof P);
            P.b1=M.b1; P.b2=M.b2; P.s=(u32)M.s; P.t=(u32)M.t; P.j=(u32)M.j; P.je=(u32)M.je;
            P.kdig=(u32)g_kdig; P.odd=(u32)g_odd;
            P.Q=M.Q; P.QF=M.QF; P.c3f=M.c3f;
            P.Mb1=(u64)(((u128)1<<64)/M.b1); P.Mb2=(u64)(((u128)1<<64)/M.b2);
            P.MQ=(u64)(((u128)1<<64)/M.Q);  P.MQF=(u64)(((u128)1<<64)/M.QF);
            P.PM0=M.PM.d[0];P.PM1=M.PM.d[1];P.PM2=M.PM.d[2];
            P.PS0=M.PS.d[0];P.PS1=M.PS.d[1];P.PS2=M.PS.d[2];
            P.PT0=M.PMT.d[0];P.PT1=M.PMT.d[1];P.PT2=M.PMT.d[2];
            P.PB0=M.PB1NJE.d[0];P.PB1x=M.PB1NJE.d[1];P.PB2=M.PB1NJE.d[2];
            P.G=256; P.outCap=outCap;

            u64 span=M.yh_end-M.yh_start;
            u64 done=skip0;
            double lastprint=0, rate=0;
            int partner_done=0;
            for (u64 base=M.yh_start+skip0; base<M.yh_end && !partner_done; base+=g_slab){
                w192 stopv = g_have ? g_best : w_max_();
                if (w_cmp(M.hi_p,stopv)<0) stopv=M.hi_p;
                if (w_cmp(w_mulu(M.PMT,base),stopv)>=0) break;
                u64 end=base+g_slab; if (end>M.yh_end) end=M.yh_end;
                P.yh0=base; P.yhEnd=end;
                P.ST0=stopv.d[0]; P.ST1=stopv.d[1]; P.ST2=stopv.d[2];
                *dN=0;
                u64 nthreads=(end-base+P.G-1)/P.G;
                u32 tpb=256;
                u64 blocks=(nthreads+tpb-1)/tpb;
                double s0=nowsec();
                sweep<<<(unsigned)blocks,tpb>>>(g_off,g_ent,P,dN,dP,dS);
                CUCHECK(cudaGetLastError());
                CUCHECK(cudaDeviceSynchronize());
                double dt=nowsec()-s0;
                done+=end-base; total_done+=end-base;
                double inst=(double)(end-base)/(dt>0?dt:1e-9);
                rate = rate>0 ? 0.7*rate+0.3*inst : inst;
                u32 nm=*dN; if(nm>outCap)nm=outCap;
                if (nm){
                    for (u32 i2=0;i2<nm;i2++){
                        u64 yh=dP[2*i2], yl=dP[2*i2+1];
                        w192 v=vof(&M,yh,yl);
                        if (w_cmp(v,M.hi_p)<0 && fullpal_host(M.b1,v))
                            report_match(v,M.b1,M.b2);
                        else
                            fprintf(stderr,"  (GPU candidate failed host verify — ignored)\n");
                    }
                    partner_done=1;
                }
                double tn=nowsec();
                if (!g_quiet && tn-lastprint>2.0){
                    lastprint=tn;
                    double eta=rate>0?(double)(span-done)/rate:0;
                    char bb[64];
                    fprintf(stderr,"[b2=%u b1=%u CUDA] stage %5.1f%% (%llu/%llu) | %.1fM/s | ETA %ld:%02ld:%02ld | best %s\n",
                            b2,M.b1,100.0*(double)done/(double)span,
                            (unsigned long long)done,(unsigned long long)span,
                            rate/1e6,(long)(eta/3600),((long)eta/60)%60,(long)eta%60,
                            g_have?w_str(g_best,bb):"-");
                }
            }
            cudaFree(g_off); cudaFree(g_ent); g_off=NULL; g_ent=NULL;
        }
    }

    fprintf(stderr,"done: %.2fs, %llu yh swept\n",nowsec()-t_start,(unsigned long long)total_done);
    if (g_have){
        char buf[64];
        printf("a(%d) = %s   (palindromic with %d digits in bases %u and %u)\n",
               g_n,w_str(g_best,buf),g_n,g_bb1,g_bb2);
        print_all_bases(g_best);
        return 0;
    }
    char b[64];
    printf("no doubly %d-digit palindrome found below %s\n",
           g_n, g_limit_set?w_str(g_limitW,b):"(unbounded?!)");
    return 1;
}
