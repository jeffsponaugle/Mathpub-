// kernels.cu - CUDA kernels, compiled at run time by NVRTC (gpu_cuda.cpp) with the SPEC_* #defines of
// the current search, exactly like kernels.metal.  The per-node logic is shared with the CPU (core.h).
// Grids are 1-D (2-D work is flattened: CUDA limits grid y to 65535), the deferred-candidate passes
// loop over the count left by the previous kernel of the same stream (no indirect dispatch needed).
#include "core.h"

__constant__ Params c_prm;  // set by the host before each run

// minimum resident blocks of 256 threads per SM (caps registers per thread: occupancy vs spills)
#ifndef NODES_MINB
#define NODES_MINB 1
#endif
#ifndef AGEN_MINB
#define AGEN_MINB 3  // measured on GB10: 3 (<= 85 registers) beats 2 by ~10%; 5 spills
#endif

struct KArgs {
  u32 nSuffix;
  u32 nPrefix;
  u32 hitCap;
  u32 candCap;
};

__device__ __forceinline__ u32 warp_sum(u32 v) { return __reduce_add_sync(0xffffffffu, v); }

// counters: [0] hits  [1] probes  [2] candidates  [3] overflow pieces  [4] deferred candidates
struct GpuSink {
  u32 probes;
  u32 cands;
  u32 x, y;
  u32* counters;
  u32* cbuf;
  u32 ccap;
  const u32* revtab;
  const NumQ* recHi;
  const NumQ* recLo;
  NumQ* hits;
  u32 hcap;
  __device__ void emit(NumQ n) {
    u32 i = atomicAdd(&counters[0], 1u);
    if (i < hcap) hits[i] = n;
  }
  // defer the full test to check_cands (no divergence here); test inline if the buffer is full
  __device__ void candidate(NumQ C, u32 idx) {
    u32 i = atomicAdd(&counters[4], 1u);
    if (i < ccap) {
      cbuf[3 * i] = x;
      cbuf[3 * i + 1] = y;
      cbuf[3 * i + 2] = idx;
    } else {
      check_candidate(C, idx, c_prm, revtab, recHi, recLo, *this);
    }
  }
};

// node_requests() hands every prepared piece to add(); on the GPU we probe immediately.
struct GpuOut {
  const u32* offs;
  const u64* ents;
  NumQ C;
  GpuSink* sink;
  u32 over;
  __device__ void add(ProbeReq rq) {
    if (rq.flag == 2) {
      over++;
      return;
    }
    run_probe(rq, C, c_prm, offs, ents, *sink);
  }
};

extern "C" __global__ void __launch_bounds__(256, NODES_MINB) search_nodes(const NumQ* prefixes, const NumQ* suffix, const u32* revtab, const u32* offs,
                                        const u64* ents, const NumQ* recHi, const NumQ* recLo, u32* counters,
                                        NumQ* hits, KArgs ka, u32* cbuf) {
  const Params& p = c_prm;
  u64 gid = (u64)blockIdx.x * blockDim.x + threadIdx.x;
  GpuSink sink;
  sink.probes = 0;
  sink.cands = 0;
  sink.x = 0;
  sink.y = 0;
  sink.counters = counters;
  sink.cbuf = cbuf;
  sink.ccap = ka.candCap;
  sink.revtab = revtab;
  sink.recHi = recHi;
  sink.recLo = recLo;
  sink.hits = hits;
  sink.hcap = ka.hitCap;
  u32 over = 0;
  if (gid < (u64)ka.nSuffix * ka.nPrefix) {
    u32 gy = (u32)(gid / ka.nSuffix);
    u32 gx = (u32)(gid - (u64)gy * ka.nSuffix);
    sink.x = gx;
    sink.y = gy;
    GpuOut out;
    out.offs = offs;
    out.ents = ents;
    out.C = nq_add(prefixes[gy], suffix[gx], PRM_RHO);
    out.sink = &sink;
    out.over = 0;
    node_requests(out.C, p, revtab, out);
    over = out.over;
  }
  u32 pr = warp_sum(sink.probes);
  u32 ca = warp_sum(sink.cands);
  u32 ov = warp_sum(over);
  if ((threadIdx.x & 31) == 0) {
    atomicAdd(&counters[1], pr);
    atomicAdd(&counters[2], ca);
    if (ov) atomicAdd(&counters[3], ov);
  }
}

struct CandSink {
  u32* counters;
  NumQ* hits;
  u32 hcap;
  __device__ void emit(NumQ n) {
    u32 i = atomicAdd(&counters[0], 1u);
    if (i < hcap) hits[i] = n;
  }
};

extern "C" __global__ void check_cands(const NumQ* prefixes, const NumQ* suffix, const u32* revtab, const NumQ* recHi,
                                       const NumQ* recLo, u32* counters, NumQ* hits, KArgs ka, const u32* cbuf) {
  const Params& p = c_prm;
  u32 n = UMIN(counters[4], ka.candCap);  // final: written by search_nodes earlier in this stream
  for (u32 i = blockIdx.x * blockDim.x + threadIdx.x; i < n; i += gridDim.x * blockDim.x) {
    u32 x = cbuf[3 * i], y = cbuf[3 * i + 1], idx = cbuf[3 * i + 2];
    NumQ C = nq_add(prefixes[y], suffix[x], PRM_RHO);
    CandSink cs;
    cs.counters = counters;
    cs.hits = hits;
    cs.hcap = ka.hitCap;
    check_candidate(C, idx, c_prm, revtab, recHi, recLo, cs);
  }
}

// Host/device struct layout check.
extern "C" __global__ void layout_check(u32* out) {
  out[0] = sizeof(Params);
  out[1] = sizeof(SegInfo);
  out[2] = sizeof(NumQ);
  out[3] = sizeof(FastDiv);
  out[4] = c_prm.Q;
  out[5] = c_prm.seg[1].LQ;
  out[6] = (u32)(c_prm.qpow[3] & 0xffffffffu);
  out[7] = c_prm.centerDig[1];
  out[8] = (u32)(c_prm.Wmax.l[1] & 0xffffffffu);
  out[9] = c_prm.recLoN;
}

// =============================================================================================
// v2 kernels (class-partitioned search).  Batch table layout: A[0..nBuckets+2), counts are added
// at A[key+2], an inclusive scan makes A[key+1] the start of bucket key, the scatter advances
// A[key+1] to the end of bucket key; afterwards bucket key = [A[key], A[key+1]).
// =============================================================================================
struct V2Args {
  u64 c0;          // first class of the batch
  u32 ww;          // classes in this batch
  u32 nBuckets;    // w * Q^tp
  u32 nB1, subB;   // B1 threads x sub-threads per b1
  u32 nPieces, subA;
  u32 hitCap, candCap;
  u32 nBlocks;     // scan blocks of 1024
  u32 pad;
};

extern "C" __global__ void v2_zero(u32* A, V2Args va, u32* counters) {
  u32 gid = blockIdx.x * blockDim.x + threadIdx.x;
  if (gid < va.nBuckets + 2) A[gid] = 0;
  if (gid == 0) counters[4] = 0;
}

// cyclic residue range [lo, lo+w) of a CSR keyed by x mod Qr
__device__ __forceinline__ void csr_range(const u32* offs, u64 Qr, u64 lo, u64 w, u32* r) {
  if (lo + w <= Qr) {
    r[0] = offs[lo];
    r[1] = offs[lo + w];
    r[2] = 0;
    r[3] = 0;
  } else {
    r[0] = offs[lo];
    r[1] = offs[Qr];
    r[2] = offs[0];
    r[3] = offs[lo + w - Qr];
  }
}

extern "C" __global__ void v2_bgen(const u64* b1V, const u64* b2V, const u32* b2Idx, const u32* b2Offs, u32* A,
                                   u64* ents, V2Args va, u32 scatter) {
  const Params& p = c_prm;
  u64 gid = (u64)blockIdx.x * blockDim.x + threadIdx.x;
  if (gid >= (u64)va.nB1 * va.subB) return;
  u32 gy = (u32)(gid / va.nB1);
  u32 gx = (u32)(gid - (u64)gy * va.nB1);
  u64 Qr = p.Qr;
  u64 V1 = b1V[gx];
  u64 lo = (va.c0 + Qr - QDIV(V1, PRM_R0) % Qr) % Qr;
  u32 rg[4];
  csr_range(b2Offs, Qr, lo, va.ww, rg);
  u64 idxBase = (u64)gx * p.nB2;
  for (int part = 0; part < 2; part++)
    for (u32 i = rg[2 * part] + gy; i < rg[2 * part + 1]; i += va.subB) {
      u64 V = V1 + b2V[i];
      if (V >= PRM_RHO) V -= PRM_RHO;
      u32 key = (u32)v2_key(V, va.c0, p);
      if (scatter) {
        u32 pos = atomicAdd(&A[key + 1], 1u);
        ents[pos] = v2_entry(V, idxBase + b2Idx[i], p);
      } else {
        atomicAdd(&A[key + 2], 1u);
      }
    }
}

// warp-inclusive scan
__device__ __forceinline__ u32 warp_incl(u32 v, u32 lane) {
  for (u32 o = 1; o < 32; o <<= 1) {
    u32 y = __shfl_up_sync(0xffffffffu, v, o);
    if (lane >= o) v += y;
  }
  return v;
}

// inclusive scan, 1024 elements per block of 256 threads
extern "C" __global__ void scan_local(u32* A, u32* blockSums, V2Args va) {
  __shared__ u32 warpTot[8];
  u32 n = va.nBuckets + 2;
  u32 tid = threadIdx.x, lane = tid & 31, wid = tid >> 5;
  u32 base = blockIdx.x * 1024 + tid * 4;
  u32 v[4];
  u32 s = 0;
  for (int j = 0; j < 4; j++) {
    u32 x = (base + j < n) ? A[base + j] : 0;
    s += x;
    v[j] = s;
  }
  u32 incl = warp_incl(s, lane);
  if (lane == 31) warpTot[wid] = incl;
  __syncthreads();
  u32 off = incl - s;
  for (u32 j = 0; j < wid; j++) off += warpTot[j];
  for (int j = 0; j < 4; j++)
    if (base + j < n) A[base + j] = v[j] + off;
  if (tid == 255) blockSums[blockIdx.x] = off + s;
}

// exclusive scan of the block sums (single block of 1024)
extern "C" __global__ void scan_blocks(u32* blockSums, V2Args va) {
  __shared__ u32 warpTot[32];
  u32 nb = va.nBlocks;
  u32 tid = threadIdx.x, lane = tid & 31, wid = tid >> 5;
  u32 per = (nb + 1023) / 1024;
  u32 start = tid * per;
  u32 s = 0;
  for (u32 j = 0; j < per; j++)
    if (start + j < nb) s += blockSums[start + j];
  u32 incl = warp_incl(s, lane);
  if (lane == 31) warpTot[wid] = incl;
  __syncthreads();
  u32 off = incl - s;
  for (u32 j = 0; j < wid; j++) off += warpTot[j];
  for (u32 j = 0; j < per; j++)
    if (start + j < nb) {
      u32 x = blockSums[start + j];
      blockSums[start + j] = off;
      off += x;
    }
}

extern "C" __global__ void scan_add(u32* A, const u32* blockSums, V2Args va) {
  u32 gid = blockIdx.x * blockDim.x + threadIdx.x;
  if (gid < va.nBuckets + 2 && gid >= 1024) A[gid] += blockSums[gid / 1024];
}

struct V2GpuSink {
  u32 probes, cands, pre, over;
  u32 pi, ai;
  u32* counters;
  uint4* cbuf;
  u32 ccap;
  // inline fallback when the candidate buffer is full
  const u32* revtab;
  const NumQ* recB1;
  const NumQ* midW;
  const NumQ* cenW;
  NumQ* hits;
  u32 hcap;
  __device__ void overflow() { over++; }
  __device__ void emit(NumQ n) {
    u32 i = atomicAdd(&counters[0], 1u);
    if (i < hcap) hits[i] = n;
  }
  __device__ void candidate(NumQ C, u64 idxB) {
    u32 i = atomicAdd(&counters[4], 1u);
    if (i < ccap) {
      cbuf[i] = make_uint4(pi, ai, (u32)(idxB & 0xffffffffu), (u32)(idxB >> 32));
    } else {
      NumQ N = v2_rebuild(C, idxB, c_prm, recB1, midW, cenW);
      if (nq_is_pal(N, c_prm, revtab)) emit(N);
    }
  }
};

extern "C" __global__ void __launch_bounds__(256, AGEN_MINB) v2_agen(const AhiPiece* ahi, const NumQ* aloC, const u32* aloOffs, const u32* A,
                                   const u64* ents, const u32* revtab, u32* counters, uint4* cbuf,
                                   const NumQ* recB1, const NumQ* midW, NumQ* hits, const NumQ* cenW,
                                   const u64* ahiRb, V2Args va) {
  const Params& p = c_prm;
  V2GpuSink sink;
  sink.probes = 0;
  sink.cands = 0;
  sink.pre = 0;
  sink.over = 0;
  sink.pi = 0;
  sink.ai = 0;
  sink.counters = counters;
  sink.cbuf = cbuf;
  sink.ccap = va.candCap;
  sink.revtab = revtab;
  sink.recB1 = recB1;
  sink.midW = midW;
  sink.cenW = cenW;
  sink.hits = hits;
  sink.hcap = va.hitCap;
  u32 nodes = 0;
  u64 gid = (u64)blockIdx.x * blockDim.x + threadIdx.x;
  if (gid < (u64)va.nPieces * va.subA) {
    u32 gy = (u32)(gid / va.nPieces);
    u32 gx = (u32)(gid - (u64)gy * va.nPieces);
    // the class base alone decides whether this piece has outer parts in the batch;
    // the full piece is loaded only when it does
    u64 Qr = p.Qr;
    u64 base = (ahiRb[gx] + Qr - va.c0) % Qr;
    u64 lo = (base + Qr - (va.ww - 1)) % Qr;
    u32 rg[4];
    csr_range(aloOffs, Qr, lo, va.ww, rg);
    if (rg[0] + gy < rg[1] || rg[2] + gy < rg[3]) {
      const AhiPiece& pc = ahi[gx];
      sink.pi = gx;
      NumQ pa = pc.a, pb = pc.b, pC = pc.C;
      u32 seg = pc.seg;
      for (int part = 0; part < 2; part++)
        for (u32 i = rg[2 * part] + gy; i < rg[2 * part + 1]; i += va.subA) {
          NumQ C = nq_add(pC, aloC[i], PRM_RHO);
          sink.ai = i;
          nodes++;
          v2_node(C, pa, pb, seg, va.c0, p, revtab, A, ents, sink);
        }
    }
  }
  u32 nd = warp_sum(nodes);
  u32 pr = warp_sum(sink.probes);
  u32 ca = warp_sum(sink.cands);
  u32 ov = warp_sum(sink.over);
  u32 pf = warp_sum(sink.pre);
  if ((threadIdx.x & 31) == 0) {
    atomicAdd(&counters[5], nd);
    atomicAdd(&counters[6], pf);
    atomicAdd(&counters[1], pr);
    atomicAdd(&counters[2], ca);
    if (ov) atomicAdd(&counters[3], ov);
  }
}

extern "C" __global__ void v2_check(const AhiPiece* ahi, const NumQ* aloC, const u32* revtab, u32* counters,
                                    const uint4* cbuf, const NumQ* recB1, const NumQ* midW, NumQ* hits,
                                    const NumQ* cenW, V2Args va) {
  const Params& p = c_prm;
  u32 n = UMIN(counters[4], va.candCap);
  for (u32 i = blockIdx.x * blockDim.x + threadIdx.x; i < n; i += gridDim.x * blockDim.x) {
    uint4 c = cbuf[i];
    u64 idxB = (u64)c.z | ((u64)c.w << 32);
    NumQ N = v2_rebuild(nq_add(ahi[c.x].C, aloC[c.y], PRM_RHO), idxB, c_prm, recB1, midW, cenW);
    if (nq_is_pal(N, c_prm, revtab)) {
      u32 j = atomicAdd(&counters[0], 1u);
      if (j < va.hitCap) hits[j] = N;
    }
  }
}
