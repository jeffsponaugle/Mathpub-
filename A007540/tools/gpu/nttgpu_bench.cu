// GB10 micro-benchmark for a GPU big-integer multiplier feasibility study:
//  (a) device memory bandwidth, (b) 64-bit Montgomery mulmod throughput,
//  (c) radix-2 Stockham NTT over a 62-bit prime, validated vs naive DFT (small N) and inverse round trip.
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cuda_runtime.h>
typedef uint64_t u64; typedef unsigned __int128 u128; typedef uint32_t u32;
#define CK(x) do{cudaError_t e=(x); if(e!=cudaSuccess){printf("CUDA error %s at %d\n",cudaGetErrorString(e),__LINE__);exit(1);} }while(0)
static const u64 P = 4611685941117976577ULL; // 0x3fffffee00000001, p-1 = 2^33 * 311 * 1726273, primitive root 3
static const u64 G = 3;
// host helpers
static u64 h_mulmod(u64 a,u64 b){ return (u64)(((u128)a*b)%P); }
static u64 h_pow(u64 a,u64 e){ u64 r=1; while(e){ if(e&1) r=h_mulmod(r,a); a=h_mulmod(a,a); e>>=1;} return r; }
static u64 h_inv(u64 a){ return h_pow(a,P-2); }
static u64 R2;      // R^2 mod p, R=2^64
static u64 PINV;    // -p^{-1} mod 2^64
static u64 to_mont(u64 a){ return h_mulmod(a, (u64)(((u128)1<<64)%P)); }
__constant__ u64 c_P, c_PINV;
__device__ __forceinline__ u64 mont_mul(u64 a,u64 b){
  u64 lo=a*b, hi=__umul64hi(a,b);
  u64 m=lo*c_PINV;
  u64 mp_hi=__umul64hi(m,c_P);
  u64 u=hi+mp_hi+(lo!=0);
  return u>=c_P?u-c_P:u;
}
__device__ __forceinline__ u64 add_mod(u64 a,u64 b){ u64 s=a+b; return s>=c_P?s-c_P:s; }
__device__ __forceinline__ u64 sub_mod(u64 a,u64 b){ return a>=b?a-b:a+c_P-b; }
__global__ void k_copy(const u64* __restrict__ a,u64* __restrict__ b,size_t n){ size_t i=blockIdx.x*(size_t)blockDim.x+threadIdx.x; size_t st=(size_t)gridDim.x*blockDim.x; for(;i<n;i+=st) b[i]=a[i]; }
__global__ void k_mulmod(u64* out,u64 seed,int iters){
  u64 x0=seed+threadIdx.x+blockIdx.x*977ULL, x1=x0+1, x2=x0+2, x3=x0+3; u64 y=seed|1;
  for(int i=0;i<iters;i++){ x0=mont_mul(x0,y); x1=mont_mul(x1,y); x2=mont_mul(x2,y); x3=mont_mul(x3,y); }
  out[blockIdx.x*blockDim.x+threadIdx.x]=x0^x1^x2^x3;
}
__global__ void k_pass2(const u64* __restrict__ in,u64* __restrict__ out,const u64* __restrict__ t1,const u64* __restrict__ t2,u32 half,u32 Ns,u32 tw_stride){
  u32 j=blockIdx.x*blockDim.x+threadIdx.x; if(j>=half) return;
  u32 k=j&(Ns-1);
  u64 a=in[j], b=in[j+half];
  u32 idx=k*tw_stride;
  u64 w=mont_mul(t1[idx>>13],t2[idx&8191]);
  u64 wb=mont_mul(w,b);
  u32 base=((j-k)<<1)+k;
  out[base]=add_mod(a,wb); out[base+Ns]=sub_mod(a,wb);
}
__global__ void k_pass(const u64* __restrict__ in,u64* __restrict__ out,const u64* __restrict__ tw,u32 half,u32 Ns,u32 tw_stride){
  u32 j=blockIdx.x*blockDim.x+threadIdx.x; if(j>=half) return;
  u32 k=j&(Ns-1);
  u64 a=in[j], b=in[j+half];
  u64 w=tw[(size_t)k*tw_stride];
  u64 wb=mont_mul(w,b);
  u32 base=((j-k)<<1)+k;
  out[base]=add_mod(a,wb); out[base+Ns]=sub_mod(a,wb);
}
__global__ void k_scale(u64* a,u64 s,u32 n){ u32 i=blockIdx.x*blockDim.x+threadIdx.x; if(i<n) a[i]=mont_mul(a[i],s); }
// forward NTT of length N (power of 2) with root 'root' (Montgomery-form twiddle table tw[i]=root^i, i<N/2). ping-pong buffers; result in *pa.
static float ntt(u64** pa,u64** pb,const u64* tw,u32 N){
  cudaEvent_t e0,e1; cudaEventCreate(&e0); cudaEventCreate(&e1); cudaEventRecord(e0);
  u32 half=N/2; u32 blocks=(half+255)/256;
  for(u32 Ns=1;Ns<N;Ns<<=1){ k_pass<<<blocks,256>>>(*pa,*pb,tw,half,Ns,N/(2*Ns)); u64* t=*pa; *pa=*pb; *pb=t; }
  cudaEventRecord(e1); cudaEventSynchronize(e1); float ms; cudaEventElapsedTime(&ms,e0,e1); return ms;
}
static float ntt2(u64** pa,u64** pb,const u64* t1,const u64* t2,u32 N){
  cudaEvent_t e0,e1; cudaEventCreate(&e0); cudaEventCreate(&e1); cudaEventRecord(e0);
  u32 half=N/2; u32 blocks=(half+255)/256;
  for(u32 Ns=1;Ns<N;Ns<<=1){ k_pass2<<<blocks,256>>>(*pa,*pb,t1,t2,half,Ns,N/(2*Ns)); u64* t=*pa; *pa=*pb; *pb=t; }
  cudaEventRecord(e1); cudaEventSynchronize(e1); float ms; cudaEventElapsedTime(&ms,e0,e1); return ms;
}
int main(int argc,char**argv){
  int logN = argc>1?atoi(argv[1]):26;
  // Montgomery constants
  { u64 inv=1; for(int i=0;i<6;i++) inv*=2-P*inv; PINV=(u64)(0-inv); }
  CK(cudaMemcpyToSymbol(c_P,&P,8)); CK(cudaMemcpyToSymbol(c_PINV,&PINV,8));
  cudaDeviceProp prop; cudaGetDeviceProperties(&prop,0); printf("GPU: %s, %d SMs, %.0f MB\n",prop.name,prop.multiProcessorCount,prop.totalGlobalMem/1048576.0);
  cudaEvent_t e0,e1; cudaEventCreate(&e0); cudaEventCreate(&e1); float ms;
  // (a) bandwidth: copy 2 GB
  { size_t n=(size_t)1<<28; u64 *a,*b; CK(cudaMalloc(&a,n*8)); CK(cudaMalloc(&b,n*8)); k_copy<<<4096,256>>>(a,b,n); CK(cudaDeviceSynchronize());
    cudaEventRecord(e0); for(int r=0;r<3;r++) k_copy<<<4096,256>>>(a,b,n); cudaEventRecord(e1); cudaEventSynchronize(e1); cudaEventElapsedTime(&ms,e0,e1);
    printf("(a) copy kernel bandwidth: %.1f GB/s (read+write)\n", 3*2*n*8/ (ms*1e-3)/1e9); cudaFree(a); cudaFree(b); }
  // (b) mulmod throughput
  { int blocks=prop.multiProcessorCount*32, threads=256, iters=2048; u64* out; CK(cudaMalloc(&out,(size_t)blocks*threads*8));
    k_mulmod<<<blocks,threads>>>(out,12345,16); CK(cudaDeviceSynchronize());
    cudaEventRecord(e0); k_mulmod<<<blocks,threads>>>(out,12345,iters); cudaEventRecord(e1); cudaEventSynchronize(e1); cudaEventElapsedTime(&ms,e0,e1);
    double ops=(double)blocks*threads*iters*4; printf("(b) 62-bit Montgomery mulmod throughput: %.1f Gop/s\n",ops/(ms*1e-3)/1e9); cudaFree(out); }
  // (c) NTT
  u32 N=1u<<logN; u64 root=h_pow(G,(P-1)/N), iroot=h_inv(root), Ninv=h_inv(N);
  u64 *tw,*itw,*a,*b; CK(cudaMallocManaged(&tw,(N/2)*8)); CK(cudaMallocManaged(&itw,(N/2)*8)); CK(cudaMalloc(&a,(size_t)N*8)); CK(cudaMalloc(&b,(size_t)N*8));
  { u64 w=1,iw=1; for(u32 i=0;i<N/2;i++){ tw[i]=to_mont(w); itw[i]=to_mont(iw); w=h_mulmod(w,root); iw=h_mulmod(iw,iroot);} }
  // validation 1: small N naive DFT
  { u32 n=1024; u64 r=h_pow(G,(P-1)/n); u64 *x=(u64*)malloc(n*8),*xm=(u64*)malloc(n*8),*y=(u64*)malloc(n*8);
    srand(7); for(u32 i=0;i<n;i++){ x[i]=(((u64)rand()<<40)^((u64)rand()<<20)^rand())%P; xm[i]=to_mont(x[i]); }
    u64 *stw; CK(cudaMallocManaged(&stw,(n/2)*8)); { u64 w=1; for(u32 i=0;i<n/2;i++){ stw[i]=to_mont(w); w=h_mulmod(w,r);} }
    u64 *da,*db; CK(cudaMalloc(&da,n*8)); CK(cudaMalloc(&db,n*8)); CK(cudaMemcpy(da,xm,n*8,cudaMemcpyHostToDevice));
    ntt(&da,&db,stw,n); CK(cudaMemcpy(y,da,n*8,cudaMemcpyDeviceToHost));
    int bad=0; for(u32 k=0;k<n;k++){ u64 s=0, wk=h_pow(r,k), wn=1; for(u32 i=0;i<n;i++){ s=(s+h_mulmod(x[i],wn))%P; wn=h_mulmod(wn,wk);} u64 got=h_mulmod(y[k],1); /* y in Montgomery form: y = Y*R */ got=h_mulmod(y[k],h_inv((u64)(((u128)1<<64)%P))); if(got!=s) bad++; }
    printf("(c) validation vs naive DFT, N=%u: %s\n",n,bad?"MISMATCH":"ok"); }
  // validation 2 + timing: forward, inverse, compare
  { u64 *x=(u64*)malloc((size_t)N*8); srand(11); for(u32 i=0;i<N;i++) x[i]=to_mont((((u64)rand()<<40)^((u64)rand()<<20)^rand())%P);
    CK(cudaMemcpy(a,x,(size_t)N*8,cudaMemcpyHostToDevice)); ntt(&a,&b,tw,N); CK(cudaDeviceSynchronize()); // warm-up
    CK(cudaMemcpy(a,x,(size_t)N*8,cudaMemcpyHostToDevice));
    float f=ntt(&a,&b,tw,N); float g=ntt(&a,&b,itw,N); k_scale<<<(N+255)/256,256>>>(a,to_mont(Ninv),N); CK(cudaDeviceSynchronize());
    u64 *y=(u64*)malloc((size_t)N*8); CK(cudaMemcpy(y,a,(size_t)N*8,cudaMemcpyDeviceToHost)); size_t bad=0; for(u32 i=0;i<N;i++) if(y[i]!=x[i]) bad++;
    printf("(c) N=2^%d radix-2 Stockham NTT: forward %.1f ms, inverse %.1f ms (%d passes each, %.1f ms/pass, %.0f GB/s effective); round trip %s\n",logN,f,g,logN,f/logN, (double)N*8*2*logN/(f*1e-3)/1e9, bad?"MISMATCH":"ok"); }
  if(logN>=14){ u64 *t1,*t2; CK(cudaMallocManaged(&t1,8192*8)); CK(cudaMallocManaged(&t2,8192*8));
    { u64 w=1,r13=h_pow(root,8192); for(u32 i=0;i<8192;i++){ t2[i]=to_mont(w); w=h_mulmod(w,root);} w=1; for(u32 i=0;i<8192;i++){ t1[i]=to_mont(w); w=h_mulmod(w,r13);} }
    u64 *x=(u64*)malloc((size_t)N*8); srand(11); for(u32 i=0;i<N;i++) x[i]=to_mont((((u64)rand()<<40)^((u64)rand()<<20)^rand())%P);
    CK(cudaMemcpy(a,x,(size_t)N*8,cudaMemcpyHostToDevice)); ntt2(&a,&b,t1,t2,N); CK(cudaDeviceSynchronize());
    CK(cudaMemcpy(a,x,(size_t)N*8,cudaMemcpyHostToDevice)); float f=ntt2(&a,&b,t1,t2,N);
    u64 *y=(u64*)malloc((size_t)N*8); CK(cudaMemcpy(y,a,(size_t)N*8,cudaMemcpyDeviceToHost));
    // compare with single-table forward result
    CK(cudaMemcpy(a,x,(size_t)N*8,cudaMemcpyHostToDevice)); ntt(&a,&b,tw,N); u64 *z=(u64*)malloc((size_t)N*8); CK(cudaMemcpy(z,a,(size_t)N*8,cudaMemcpyDeviceToHost));
    size_t bad=0; for(u32 i=0;i<N;i++) if(y[i]!=z[i]) bad++;
    printf("(c2) two-level twiddle table: forward %.1f ms (%.1f ms/pass, %.0f GB/s effective), matches single-table: %s\n",f,f/logN,(double)N*8*2*logN/(f*1e-3)/1e9,bad?"NO":"yes"); }
  return 0;
}
