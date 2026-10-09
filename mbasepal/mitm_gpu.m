/*
 * mitm_gpu.m — Metal GPU prototype of the A171740 MITM sweep.
 *
 * Benchmarks the yh-sweep of one (n, b2, b1) stage window on CPU
 * (pthreads, same algorithm as a171740.c) and on the GPU, cross-checking
 * filter statistics and matches. Table build stays on the CPU; the GPU
 * reads it via unified memory. Survivors passing every in-kernel filter
 * are re-verified exactly on the host.
 *
 * build: clang -O3 -fobjc-arc -framework Metal -framework Foundation \
 *          -o mitm_gpu mitm_gpu.m
 * usage: ./mitm_gpu <n> <b2> <b1> <yh_offset> <yh_count> [tablecap]
 *   yh_offset is relative to the stage's yh_start.
 */

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#include <pthread.h>
#include <stdatomic.h>

typedef unsigned __int128 u128;
typedef uint64_t u64;
typedef uint32_t u32;

/* ---------------- host w192 (subset of a171740.c) ---------------- */
typedef struct { u64 d[3]; } w192;
static w192 w_from_u128(u128 x){ w192 r={{(u64)x,(u64)(x>>64),0}}; return r; }
static int w_cmp(w192 a, w192 b){
    for (int i=2;i>=0;i--){ if(a.d[i]<b.d[i])return -1; if(a.d[i]>b.d[i])return 1; }
    return 0;
}
static w192 w_add(w192 a, w192 b){
    w192 r; unsigned c=0;
    for (int i=0;i<3;i++){ u64 s=a.d[i]+b.d[i]; u64 s2=s+c; c=(s<a.d[i])||(s2<s); r.d[i]=s2; }
    return r;
}
static w192 w_addu(w192 a, u64 b){ w192 r=a; r.d[0]+=b; if(r.d[0]<b){ if(++r.d[1]==0) ++r.d[2]; } return r; }
static w192 w_mulu(w192 a, u64 m){
    u128 t0=(u128)a.d[0]*m, t1=(u128)a.d[1]*m+(u64)(t0>>64), t2=(u128)a.d[2]*m+(u64)(t1>>64);
    w192 r={{(u64)t0,(u64)t1,(u64)t2}}; return r;
}
static w192 w_pow(u64 b, int e){ w192 r=w_from_u128(1); while(e-->0) r=w_mulu(r,b); return r; }
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
static int w_is_zero(w192 a){ return !(a.d[0]|a.d[1]|a.d[2]); }
static char *w_str(w192 v, char *buf){
    char *p=buf+63; *--p=0; if(w_is_zero(v)) *--p='0';
    while(!w_is_zero(v)){ unsigned r=w_divmod_small(&v,1000000000u);
        for(int i=0;i<9&&!(w_is_zero(v)&&r==0);i++){ *--p=(char)('0'+r%10); r/=10; } }
    return p;
}
static u64 pw64(u64 b,int e){ u64 r=1; while(e-->0) r*=b; return r; }
static u64 revdig(u64 x, unsigned b, int d){ u64 r=0; for(int i=0;i<d;i++){ r=r*b+x%b; x/=b; } return r; }

/* ---------------- stage parameters (mirrors run_mitm_stage) ---------------- */
typedef struct {
    int n, K, kdig, odd, s, t, j, je;
    unsigned b1, b2;
    u64 Q, QF, c1f, c2f, c3f;
    u128 PM, PS, PMT, PB1NJE;   /* all fit u128 for n <= ~40 */
    w192 STOP;                  /* b1^n */
    u64 yh_start, yh_end;
} Cfg;

static void cfg_init(Cfg *C, int n, unsigned b2, unsigned b1, u64 tablecap){
    C->n=n; C->b1=b1; C->b2=b2;
    C->kdig = n/2; C->odd = n&1;
    int m = (n+1)/2;
    C->K = C->odd ? n/2 : m;
    int t=0; u64 e=1;
    while (t+1 <= m/2 && e*b2 <= tablecap){ e*=b2; t++; }
    C->t=t; C->s=m-t;
    int j=0; u64 q=1;
    while (q*b1 <= e){ q*=b1; j++; }
    if (j<1){ j=1; q=b1; }
    C->j=j; C->Q=q;
    int je=j; u64 qf=q;
    u128 PMT = 1; for (int i=0;i<C->K+t;i++) PMT*=b2;
    while (je+1 <= C->kdig && je+1 <= n-1 &&
           (u128)qf*b1 <= ((u64)1<<40) && (qf/q)*(u64)b1 <= ((u64)1<<30)){
        u128 pnje=1; int ok=1;
        for (int i=0;i<n-(je+1);i++){ pnje*=b1; if(!pnje){ok=0;break;} }
        if(!ok) break;
        if (PMT/pnje > 8) break;
        je++; qf*=b1;
    }
    C->je=je; C->QF=qf;
    C->PM=1; for(int i=0;i<C->K;i++) C->PM*=b2;
    C->PS=1; for(int i=0;i<C->s;i++) C->PS*=b2;
    C->PMT=PMT;
    C->PB1NJE=1; for(int i=0;i<n-je;i++) C->PB1NJE*=b1;
    C->c1f=(u64)(C->PM%qf); C->c2f=(u64)(C->PS%qf); C->c3f=(u64)(C->PMT%qf);
    C->STOP = w_pow(b1, n);
    C->yh_start = pw64(b2, C->s-1);
    {
        /* yh_end = floor(b1^n / PMT) + 1: via u128 when b1^n fits,
         * else binary search on w_mulu (exact) */
        w192 hi=C->STOP;
        if (hi.d[2]==0){ u128 h=((u128)hi.d[1]<<64)|hi.d[0]; C->yh_end=(u64)(h/C->PMT)+1; }
        else {
            u64 lo=1, hi2=(u64)1<<62;
            while (lo<hi2){ u64 mid=lo+(hi2-lo)/2;
                w192 prod=w_mulu(w_from_u128(C->PMT), mid);
                if (w_cmp(prod, C->STOP)<=0) lo=mid+1; else hi2=mid; }
            C->yh_end=lo; /* first yh with yh*PMT > STOP  ==  floor+1 */
        }
        u64 cap=pw64(b2, C->s);
        if (C->yh_end>cap) C->yh_end=cap;
    }
}

/* ---------------- table build (CPU) ---------------- */
static u32 *g_off; static u64 *g_ent; static u64 g_Nyl;
static void build_table(const Cfg *C){
    g_Nyl = pw64(C->b2, C->t);
    u32 *tmpB = malloc(g_Nyl*4);
    g_off = calloc(C->Q+1, 4);
    g_ent = malloc(g_Nyl*8);
    if(!tmpB||!g_off||!g_ent){ fprintf(stderr,"alloc fail\n"); exit(1); }
    for (u64 yl=0; yl<g_Nyl; yl++){
        u64 r = C->odd ? revdig(yl/C->b2, C->b2, C->t-1) : revdig(yl, C->b2, C->t);
        u64 Bf = (u64)(((u128)yl*C->c1f + (u128)r*C->c2f) % C->QF);
        u32 key=(u32)(Bf%C->Q); tmpB[yl]=key; g_off[key+1]++;
    }
    for (u64 q=0;q<C->Q;q++) g_off[q+1]+=g_off[q];
    u32 *fill=malloc(C->Q*4); memcpy(fill,g_off,C->Q*4);
    for (u64 yl=0; yl<g_Nyl; yl++){
        u64 r = C->odd ? revdig(yl/C->b2, C->b2, C->t-1) : revdig(yl, C->b2, C->t);
        u64 Bf = (u64)(((u128)yl*C->c1f + (u128)r*C->c2f) % C->QF);
        g_ent[fill[tmpB[yl]]++]=((u64)yl<<32)|(u32)(Bf/C->Q);
    }
    free(fill); free(tmpB);
}

/* ---------------- exact host check ---------------- */
static int fullpal_host(const Cfg *C, w192 v){
    w192 w=v; u128 rr=0;
    for (int i=0;i<C->kdig;i++) rr=rr*C->b1 + w_divmod_small(&w, C->b1);
    if (C->odd) w_divmod_small(&w, C->b1);
    return w.d[2]==0 && ((((u128)w.d[1]<<64)|w.d[0])==rr);
}
static w192 vof(const Cfg *C, u64 yh, u64 yl){
    u64 rs=revdig(yh,C->b2,C->s);
    u64 rt=C->odd?revdig(yl/C->b2,C->b2,C->t-1):revdig(yl,C->b2,C->t);
    w192 v=w_mulu(w_from_u128(C->PMT), yh);
    v=w_add(v, w_mulu(w_from_u128(C->PM), yl));
    v=w_add(v, w_mulu(w_from_u128(C->PS), rt));
    return w_addu(v, rs);
}

/* ---------------- CPU reference sweep ---------------- */
typedef struct { const Cfg *C; u64 y0, y1; _Atomic(u64) *cursor;
                 u64 recon, fulls; u64 matches[16]; int nmatch; } CpuArg;
#define CPBLOCK 65536ULL
static void *cpu_worker(void *ap){
    CpuArg *A=ap; const Cfg *C=A->C;
    for(;;){
        u64 y0=atomic_fetch_add(A->cursor, CPBLOCK);
        if (y0>=A->y1) break;
        u64 y1=y0+CPBLOCK; if(y1>A->y1) y1=A->y1;
        w192 vbase=w_mulu(w_from_u128(C->PMT), y0);
        w192 den=w_from_u128(C->PB1NJE);
        /* seed p, prem */
        w192 r={{0,0,0}}; u64 p=0;
        for (int i=191;i>=0;i--){
            r.d[2]=(r.d[2]<<1)|(r.d[1]>>63); r.d[1]=(r.d[1]<<1)|(r.d[0]>>63);
            r.d[0]=(r.d[0]<<1)|((vbase.d[i>>6]>>(i&63))&1);
            if (w_cmp(r,den)>=0){ unsigned br=0; w192 nr;
                for(int k2=0;k2<3;k2++){ u64 s=r.d[k2]-den.d[k2]; u64 s2=s-br;
                    br=(r.d[k2]<den.d[k2])||(s<br); nr.d[k2]=s2; }
                r=nr; if(i<64) p|=1ULL<<i; }
        }
        w192 prem=r;
        u64 term3=(u64)(((u128)(y0%C->QF)*C->c3f)%C->QF);
        for (u64 yh=y0; yh<y1; yh++){
            u64 rs=revdig(yh,C->b2,C->s);
            u64 Af=term3+rs%C->QF; if(Af>=C->QF) Af-=C->QF;
            w192 vtop=w_add(vbase, w_from_u128(C->PMT)); /* -1 */
            if (vtop.d[0]--==0){ if(vtop.d[1]--==0) vtop.d[2]--; }
            if (w_cmp(vtop, C->STOP)>=0){ vtop=C->STOP;
                if (vtop.d[0]--==0){ if(vtop.d[1]--==0) vtop.d[2]--; } }
            for (u64 pp=p;;pp++){
                w192 plo=w_mulu(den, pp);
                if (w_cmp(plo,vtop)>0) break;
                w192 phi=w_add(plo,den);
                u64 rjf=revdig(pp,C->b1,C->je);
                u64 R=rjf>=Af?rjf-Af:rjf+C->QF-Af;
                u64 key=R%C->Q; u32 needhi=(u32)(R/C->Q);
                for (u32 idx=g_off[key]; idx<g_off[key+1]; idx++){
                    u64 en=g_ent[idx];
                    if ((u32)en!=needhi) continue;
                    u64 yl=en>>32;
                    A->recon++;
                    w192 v=vof(C,yh,yl);
                    if (w_cmp(v,plo)<0||w_cmp(v,phi)>=0) continue;
                    if (w_cmp(v,C->STOP)>=0) continue;
                    A->fulls++;
                    if (fullpal_host(C,v) && A->nmatch<16){
                        A->matches[A->nmatch*2]=yh; A->matches[A->nmatch*2+1]=yl; A->nmatch++;
                    }
                }
            }
            vbase=w_add(vbase, w_from_u128(C->PMT));
            prem=w_add(prem, w_from_u128(C->PMT));
            while (w_cmp(prem,den)>=0){ unsigned br=0; w192 nr;
                for(int k2=0;k2<3;k2++){ u64 s=prem.d[k2]-den.d[k2]; u64 s2=s-br;
                    br=(prem.d[k2]<den.d[k2])||(s<br); nr.d[k2]=s2; }
                prem=nr; p++; }
            term3+=C->c3f; if(term3>=C->QF) term3-=C->QF;
        }
    }
    return NULL;
}

/* ---------------- Metal kernel source ---------------- */
static NSString *kSrc = @""
"#include <metal_stdlib>\n"
"using namespace metal;\n"
"struct W { ulong d0,d1,d2; };\n"
"struct P {\n"
"  uint b1,b2,s,t,j,je,kdig,odd;\n"
"  ulong Q,QF,c3f;\n"
"  ulong Mb1,Mb2,MQ,MQF;\n"
"  ulong PM0,PM1,PM2, PS0,PS1,PS2, PT0,PT1,PT2, PB0,PB1x,PB2, ST0,ST1,ST2;\n"
"  ulong yh0,yhEnd,G; uint outCap; uint pad;\n"
"};\n"
"static inline ulong divq(ulong x, ulong M, ulong d, thread ulong &r){\n"
"  ulong q=mulhi(x,M); r=x-q*d; if(r>=d){r-=d;q++;} return q; }\n"
"static inline W mkw(ulong a, ulong b, ulong c){ W r; r.d0=a; r.d1=b; r.d2=c; return r; }\n"
"static inline int wcmp(W a, W b){\n"
"  if (a.d2!=b.d2) return a.d2<b.d2?-1:1;\n"
"  if (a.d1!=b.d1) return a.d1<b.d1?-1:1;\n"
"  if (a.d0!=b.d0) return a.d0<b.d0?-1:1; return 0; }\n"
"static inline W wadd(W a, W b){ W r; r.d0=a.d0+b.d0; ulong c0=r.d0<a.d0?1:0;\n"
"  ulong s1=a.d1+b.d1; ulong c1=s1<a.d1?1:0; r.d1=s1+c0; if(r.d1<s1) c1++;\n"
"  r.d2=a.d2+b.d2+c1; return r; }\n"
"static inline W wsub(W a, W b){ W r; r.d0=a.d0-b.d0; ulong br0=a.d0<b.d0?1:0;\n"
"  ulong s1=a.d1-b.d1; ulong br1=a.d1<b.d1?1:0; r.d1=s1-br0; if(s1<br0) br1++;\n"
"  r.d2=a.d2-b.d2-br1; return r; }\n"
"static inline W wsub1(W a){ W r=a; if(r.d0--==0){ if(r.d1--==0) r.d2--; } return r; }\n"
"static inline W waddu(W a, ulong b){ W r=a; r.d0+=b; if(r.d0<b){ if(++r.d1==0) ++r.d2; } return r; }\n"
"static inline W wmulu(W a, ulong m){\n"
"  ulong l0=a.d0*m, h0=mulhi(a.d0,m);\n"
"  ulong l1=a.d1*m, h1=mulhi(a.d1,m);\n"
"  ulong l2=a.d2*m;\n"
"  W r; r.d0=l0; r.d1=l1+h0; ulong c=(r.d1<l1)?1:0; r.d2=l2+h1+c; return r; }\n"
"static inline ulong wdivbig(W num, W den, thread W &rem){\n"
"  W r=mkw(0,0,0); ulong q=0;\n"
"  for (int i=191;i>=0;i--){\n"
"    r.d2=(r.d2<<1)|(r.d1>>63); r.d1=(r.d1<<1)|(r.d0>>63);\n"
"    ulong bit = (i>=128)?((num.d2>>(i-128))&1):(i>=64)?((num.d1>>(i-64))&1):((num.d0>>i)&1);\n"
"    r.d0=(r.d0<<1)|bit;\n"
"    if (wcmp(r,den)>=0){ r=wsub(r,den); if(i<64) q|=(ulong)1<<i; }\n"
"  }\n"
"  rem=r; return q; }\n"
"static inline ulong mod128(ulong hi, ulong lo, ulong m){\n"
"  ulong r=0;\n"
"  for (int i=127;i>=0;i--){\n"
"    ulong top = r>>63; r=(r<<1)|((i>=64)?((hi>>(i-64))&1):((lo>>i)&1));\n"
"    if (top || r>=m) r-=m; }\n"
"  return r; }\n"
"static inline ulong revm(ulong x, uint b, uint d, ulong Mb){\n"
"  ulong r=0;\n"
"  for(uint i=0;i<d;i++){ ulong rem; x=divq(x,Mb,b,rem); r=r*b+rem; }\n"
"  return r; }\n"
"kernel void sweep(device const uint *off [[buffer(0)]],\n"
"                  device const ulong *ent [[buffer(1)]],\n"
"                  constant P &C [[buffer(2)]],\n"
"                  device atomic_uint *outN [[buffer(3)]],\n"
"                  device ulong *outP [[buffer(4)]],\n"
"                  device atomic_uint *stats [[buffer(5)]],\n"
"                  uint tid [[thread_position_in_grid]]){\n"
"  ulong yhb = C.yh0 + (ulong)tid*C.G;\n"
"  if (yhb >= C.yhEnd) return;\n"
"  ulong yhe = yhb+C.G; if (yhe>C.yhEnd) yhe=C.yhEnd;\n"
"  W PMT=mkw(C.PT0,C.PT1,C.PT2), DEN=mkw(C.PB0,C.PB1x,C.PB2);\n"
"  W PMw=mkw(C.PM0,C.PM1,C.PM2), PSw=mkw(C.PS0,C.PS1,C.PS2), STOP=mkw(C.ST0,C.ST1,C.ST2);\n"
"  W vbase=wmulu(PMT,yhb);\n"
"  if (wcmp(vbase,STOP)>=0) return;\n"
"  W prem; ulong p=wdivbig(vbase,DEN,prem);\n"
"  ulong a0; divq(yhb,C.MQF,C.QF,a0);\n"
"  ulong term3=mod128(mulhi(a0,C.c3f), a0*C.c3f, C.QF);\n"
"  uint recon=0, fulls=0;\n"
"  for (ulong yh=yhb; yh<yhe; yh++){\n"
"    ulong rs=revm(yh,C.b2,C.s,C.Mb2);\n"
"    ulong rsm; divq(rs,C.MQF,C.QF,rsm);\n"
"    ulong Af=term3+rsm; if(Af>=C.QF) Af-=C.QF;\n"
"    W vtop=wsub1(wadd(vbase,PMT));\n"
"    if (wcmp(vtop,STOP)>=0) vtop=wsub1(STOP);\n"
"    for (ulong pp=p;;pp++){\n"
"      W plo=wmulu(DEN,pp);\n"
"      if (wcmp(plo,vtop)>0) break;\n"
"      W phi=wadd(plo,DEN);\n"
"      ulong rjf=revm(pp,C.b1,C.je,C.Mb1);\n"
"      ulong R=rjf>=Af?rjf-Af:rjf+C.QF-Af;\n"
"      ulong keyr; ulong hiq=divq(R,C.MQ,C.Q,keyr);\n"
"      ulong key=keyr; uint needhi=(uint)hiq;\n"
"      uint i0=off[key], i1=off[key+1];\n"
"      for (uint idx=i0; idx<i1; idx++){\n"
"        ulong en=ent[idx];\n"
"        if ((uint)(en&0xffffffff)!=needhi) continue;\n"
"        ulong yl=en>>32;\n"
"        recon++;\n"
"        ulong xr=yl; if(C.odd){ ulong rr0; xr=divq(yl,C.Mb2,C.b2,rr0); }\n"
"        uint td=C.odd? C.t-1 : C.t;\n"
"        ulong rt=revm(xr,C.b2,td,C.Mb2);\n"
"        W v=wadd(wadd(vbase,wmulu(PMw,yl)),wmulu(PSw,rt)); v=waddu(v,rs);\n"
"        if (wcmp(v,plo)<0 || wcmp(v,phi)>=0) continue;\n"
"        if (wcmp(v,STOP)>=0) continue;\n"
"        fulls++;\n"
"        W w=v; ulong rhi=0, rlo=0;\n"
"        for (uint i=0;i<C.kdig;i++){\n"
"          ulong r64=0; ulong limb;\n"
"          limb=w.d2; { ulong a=(r64<<32)|(limb>>32); ulong q1=divq(a,C.Mb1,C.b1,r64); ulong bl=(r64<<32)|(limb&0xffffffff); ulong q0=divq(bl,C.Mb1,C.b1,r64); w.d2=(q1<<32)|q0; }\n"
"          limb=w.d1; { ulong a=(r64<<32)|(limb>>32); ulong q1=divq(a,C.Mb1,C.b1,r64); ulong bl=(r64<<32)|(limb&0xffffffff); ulong q0=divq(bl,C.Mb1,C.b1,r64); w.d1=(q1<<32)|q0; }\n"
"          limb=w.d0; { ulong a=(r64<<32)|(limb>>32); ulong q1=divq(a,C.Mb1,C.b1,r64); ulong bl=(r64<<32)|(limb&0xffffffff); ulong q0=divq(bl,C.Mb1,C.b1,r64); w.d0=(q1<<32)|q0; }\n"
"          ulong dg=r64;\n"
"          ulong nl=rlo*C.b1+dg; ulong nh=rhi*C.b1+mulhi(rlo,(ulong)C.b1);\n"
"          if (nl<dg) nh++;\n"
"          rlo=nl; rhi=nh; }\n"
"        if (C.odd){ ulong r64=0; ulong limb;\n"
"          limb=w.d2; { ulong a=(r64<<32)|(limb>>32); ulong q1=divq(a,C.Mb1,C.b1,r64); ulong bl=(r64<<32)|(limb&0xffffffff); ulong q0=divq(bl,C.Mb1,C.b1,r64); w.d2=(q1<<32)|q0; }\n"
"          limb=w.d1; { ulong a=(r64<<32)|(limb>>32); ulong q1=divq(a,C.Mb1,C.b1,r64); ulong bl=(r64<<32)|(limb&0xffffffff); ulong q0=divq(bl,C.Mb1,C.b1,r64); w.d1=(q1<<32)|q0; }\n"
"          limb=w.d0; { ulong a=(r64<<32)|(limb>>32); ulong q1=divq(a,C.Mb1,C.b1,r64); ulong bl=(r64<<32)|(limb&0xffffffff); ulong q0=divq(bl,C.Mb1,C.b1,r64); w.d0=(q1<<32)|q0; }\n"
"        }\n"
"        if ((w.d2==0) && (w.d1==rhi) && (w.d0==rlo)){\n"
"          uint slot=atomic_fetch_add_explicit(outN,1u,memory_order_relaxed);\n"
"          if (slot<C.outCap){ outP[2*slot]=yh; outP[2*slot+1]=yl; } }\n"
"      }\n"
"    }\n"
"    vbase=wadd(vbase,PMT);\n"
"    prem=wadd(prem,PMT);\n"
"    while (wcmp(prem,DEN)>=0){ prem=wsub(prem,DEN); p++; }\n"
"    term3+=C.c3f; if(term3>=C.QF) term3-=C.QF;\n"
"    if (wcmp(vbase,STOP)>=0) break;\n"
"  }\n"
"  if (recon) atomic_fetch_add_explicit(&stats[0],recon,memory_order_relaxed);\n"
"  if (fulls) atomic_fetch_add_explicit(&stats[1],fulls,memory_order_relaxed);\n"
"}\n";

typedef struct {
    u32 b1,b2,s,t,j,je,kdig,odd;
    u64 Q,QF,c3f;
    u64 Mb1,Mb2,MQ,MQF;
    u64 PM0,PM1,PM2, PS0,PS1,PS2, PT0,PT1,PT2, PB0,PB1x,PB2, ST0,ST1,ST2;
    u64 yh0,yhEnd,G; u32 outCap; u32 pad;
} GpuParams;

static double nowsec(void){ return (double)clock_gettime_nsec_np(CLOCK_MONOTONIC)/1e9; }

int main(int argc, char **argv){ @autoreleasepool {
    if (argc < 6){ fprintf(stderr,"usage: %s n b2 b1 yh_off yh_cnt [tablecap]\n",argv[0]); return 2; }
    int n = atoi(argv[1]);
    unsigned b2=(unsigned)atoi(argv[2]), b1=(unsigned)atoi(argv[3]);
    u64 off0=strtoull(argv[4],0,10), cnt=strtoull(argv[5],0,10);
    u64 cap = argc>6? strtoull(argv[6],0,10) : 300000000ULL;

    Cfg C; cfg_init(&C,n,b2,b1,cap);
    u64 y0=C.yh_start+off0, y1=y0+cnt;
    if (y1>C.yh_end) y1=C.yh_end;
    if (y0>=y1){ fprintf(stderr,"window empty/out of range (yh_end=%llu)\n",(unsigned long long)C.yh_end); return 2; }
    fprintf(stderr,"cfg: n=%d b2=%u b1=%u s=%d t=%d j=%d je=%d Q=%llu QF=%llu table=%llu\n",
            n,b2,b1,C.s,C.t,C.j,C.je,(unsigned long long)C.Q,(unsigned long long)C.QF,
            (unsigned long long)pw64(b2,C.t));
    fprintf(stderr,"window: yh [%llu, %llu) = %llu steps\n",
            (unsigned long long)y0,(unsigned long long)y1,(unsigned long long)(y1-y0));

    double tb=nowsec(); build_table(&C);
    fprintf(stderr,"table built in %.1fs (%llu entries)\n",nowsec()-tb,(unsigned long long)g_Nyl);

    /* ---- CPU reference ---- */
    int nth=10;
    CpuArg args[64]; _Atomic(u64) cursor=y0;
    pthread_t tids[64];
    double t0=nowsec();
    for (int i=0;i<nth;i++){ memset(&args[i],0,sizeof args[i]);
        args[i].C=&C; args[i].y0=y0; args[i].y1=y1; args[i].cursor=&cursor;
        pthread_create(&tids[i],NULL,cpu_worker,&args[i]); }
    u64 crecon=0,cfulls=0; int cn=0; u64 cm[32];
    for (int i=0;i<nth;i++){ pthread_join(tids[i],NULL);
        crecon+=args[i].recon; cfulls+=args[i].fulls;
        for(int k=0;k<args[i].nmatch&&cn<16;k++){ cm[cn*2]=args[i].matches[k*2]; cm[cn*2+1]=args[i].matches[k*2+1]; cn++; } }
    double cpuT=nowsec()-t0;
    fprintf(stderr,"CPU (%d thr): %.2fs = %.1f M/s | recon=%llu fulls=%llu matches=%d\n",
            nth,cpuT,(double)(y1-y0)/cpuT/1e6,(unsigned long long)crecon,(unsigned long long)cfulls,cn);

    /* ---- GPU ---- */
    id<MTLDevice> dev = MTLCreateSystemDefaultDevice();
    if (!dev){ fprintf(stderr,"no Metal device\n"); return 1; }
    fprintf(stderr,"GPU: %s\n",dev.name.UTF8String);
    NSError *err=nil;
    id<MTLLibrary> lib = [dev newLibraryWithSource:kSrc options:nil error:&err];
    if (!lib){ fprintf(stderr,"kernel compile failed: %s\n",err.localizedDescription.UTF8String); return 1; }
    id<MTLFunction> fn = [lib newFunctionWithName:@"sweep"];
    id<MTLComputePipelineState> pso = [dev newComputePipelineStateWithFunction:fn error:&err];
    if (!pso){ fprintf(stderr,"pso failed: %s\n",err.localizedDescription.UTF8String); return 1; }
    id<MTLCommandQueue> q = [dev newCommandQueue];

    id<MTLBuffer> bOff = [dev newBufferWithBytesNoCopy:g_off length:((C.Q+1)*4+16383)&~16383ULL
                          options:MTLResourceStorageModeShared deallocator:nil];
    id<MTLBuffer> bEnt = [dev newBufferWithBytesNoCopy:g_ent length:(g_Nyl*8+16383)&~16383ULL
                          options:MTLResourceStorageModeShared deallocator:nil];
    if (!bOff || !bEnt){
        /* newBufferWithBytesNoCopy needs page alignment from malloc — fall back to copy */
        bOff = [dev newBufferWithBytes:g_off length:(C.Q+1)*4 options:MTLResourceStorageModeShared];
        bEnt = [dev newBufferWithBytes:g_ent length:g_Nyl*8 options:MTLResourceStorageModeShared];
    }
    u32 outCap=1<<16;
    id<MTLBuffer> bN = [dev newBufferWithLength:4 options:MTLResourceStorageModeShared];
    id<MTLBuffer> bP = [dev newBufferWithLength:outCap*16 options:MTLResourceStorageModeShared];
    id<MTLBuffer> bS = [dev newBufferWithLength:8 options:MTLResourceStorageModeShared];
    memset(bN.contents,0,4); memset(bS.contents,0,8);

    GpuParams P; memset(&P,0,sizeof P);
    P.b1=b1; P.b2=b2; P.s=(u32)C.s; P.t=(u32)C.t; P.j=(u32)C.j; P.je=(u32)C.je;
    P.kdig=(u32)C.kdig; P.odd=(u32)C.odd;
    P.Q=C.Q; P.QF=C.QF; P.c3f=C.c3f;
    P.Mb1=(u64)(((u128)1<<64)/C.b1); P.Mb2=(u64)(((u128)1<<64)/C.b2);
    P.MQ=(u64)(((u128)1<<64)/C.Q);  P.MQF=(u64)(((u128)1<<64)/C.QF);
    w192 wpm=w_from_u128(C.PM), wps=w_from_u128(C.PS), wpt=w_from_u128(C.PMT), wpb=w_from_u128(C.PB1NJE);
    P.PM0=wpm.d[0];P.PM1=wpm.d[1];P.PM2=wpm.d[2];
    P.PS0=wps.d[0];P.PS1=wps.d[1];P.PS2=wps.d[2];
    P.PT0=wpt.d[0];P.PT1=wpt.d[1];P.PT2=wpt.d[2];
    P.PB0=wpb.d[0];P.PB1x=wpb.d[1];P.PB2=wpb.d[2];
    P.ST0=C.STOP.d[0];P.ST1=C.STOP.d[1];P.ST2=C.STOP.d[2];
    P.G=256; P.outCap=outCap;

    /* warmup + timed run over the same window, dispatched in slabs */
    u64 slab = 128ULL*1024*1024;  /* yh per command buffer */
    double gT=0;
    for (u64 base=y0; base<y1; base+=slab){
        u64 end = base+slab<y1? base+slab : y1;
        P.yh0=base; P.yhEnd=end;
        u64 nthreads=(end-base+P.G-1)/P.G;
        id<MTLCommandBuffer> cb=[q commandBuffer];
        id<MTLComputeCommandEncoder> enc=[cb computeCommandEncoder];
        [enc setComputePipelineState:pso];
        [enc setBuffer:bOff offset:0 atIndex:0];
        [enc setBuffer:bEnt offset:0 atIndex:1];
        [enc setBytes:&P length:sizeof P atIndex:2];
        [enc setBuffer:bN offset:0 atIndex:3];
        [enc setBuffer:bP offset:0 atIndex:4];
        [enc setBuffer:bS offset:0 atIndex:5];
        NSUInteger tg = pso.maxTotalThreadsPerThreadgroup; if (tg>256) tg=256;
        [enc dispatchThreads:MTLSizeMake((NSUInteger)nthreads,1,1)
              threadsPerThreadgroup:MTLSizeMake(tg,1,1)];
        [enc endEncoding];
        [cb commit];
        [cb waitUntilCompleted];
        if (cb.error){ fprintf(stderr,"GPU error: %s\n",cb.error.localizedDescription.UTF8String); return 1; }
        gT += cb.GPUEndTime - cb.GPUStartTime;
    }
    u32 gN=*(u32*)bN.contents; if(gN>outCap) gN=outCap;
    u32 *st=(u32*)bS.contents;
    fprintf(stderr,"GPU: %.2fs = %.1f M/s | recon=%u fulls=%u matches=%u\n",
            gT,(double)(y1-y0)/gT/1e6,st[0],st[1],gN);
    fprintf(stderr,"speedup vs CPU(%d thr): %.2fx\n",nth,cpuT/gT);

    /* cross-check */
    int okstats = (st[0]==(u32)crecon)&&(st[1]==(u32)cfulls);
    fprintf(stderr,"stats match CPU: %s (cpu recon=%llu fulls=%llu)\n",
            okstats?"YES":"NO",(unsigned long long)crecon,(unsigned long long)cfulls);
    u64 *op=(u64*)bP.contents;
    for (u32 i=0;i<gN;i++){
        u64 yh=op[2*i], yl=op[2*i+1];
        w192 v=vof(&C,yh,yl);
        char buf[64];
        int exact = fullpal_host(&C,v);
        printf("GPU match: yh=%llu yl=%llu v=%s host-verify=%s\n",
               (unsigned long long)yh,(unsigned long long)yl,w_str(v,buf),exact?"PASS":"FAIL");
    }
    for (int i=0;i<cn;i++){
        w192 v=vof(&C,cm[i*2],cm[i*2+1]); char buf[64];
        printf("CPU match: yh=%llu yl=%llu v=%s\n",
               (unsigned long long)cm[i*2],(unsigned long long)cm[i*2+1],w_str(v,buf));
    }
    return 0;
}}
