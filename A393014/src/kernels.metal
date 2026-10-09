// kernels.metal - GPU node enumeration; the per-node logic is shared with the CPU (core.h).
// The host prepends SPEC_* #defines (constants of the current search) before compiling.
#include "core.h"

struct KArgs {
  uint nSuffix;
  uint nPrefix;
  uint hitCap;
  uint candCap;
};

// counters: [0] hits  [1] probes  [2] candidates  [3] overflow pieces  [4] deferred candidates
struct GpuSink {
  uint probes;
  uint cands;
  uint x, y;
  device atomic_uint* counters;
  device uint* cbuf;
  uint ccap;
  constant Params* p;
  device const uint* revtab;
  device const NumQ* recHi;
  device const NumQ* recLo;
  device NumQ* hits;
  uint hcap;
  void emit(NumQ n) {
    uint i = atomic_fetch_add_explicit(&counters[0], 1u, memory_order_relaxed);
    if (i < hcap) hits[i] = n;
  }
  // defer the full test to check_cands (no divergence here); test inline if the buffer is full
  void candidate(NumQ C, uint idx) {
    uint i = atomic_fetch_add_explicit(&counters[4], 1u, memory_order_relaxed);
    if (i < ccap) {
      cbuf[3 * i] = x;
      cbuf[3 * i + 1] = y;
      cbuf[3 * i + 2] = idx;
    } else {
      check_candidate(C, idx, *p, revtab, recHi, recLo, *this);
    }
  }
};

// node_requests() hands every prepared piece to add(); on the GPU we probe immediately.
struct GpuOut {
  constant Params* p;
  device const uint* offs;
  device const ulong* ents;
  NumQ C;
  thread GpuSink* sink;
  uint over;
  void add(ProbeReq rq) {
    if (rq.flag == 2) {
      over++;
      return;
    }
    run_probe(rq, C, *p, offs, ents, *sink);
  }
};

kernel void search_nodes(constant Params& p [[buffer(0)]],
                         device const NumQ* prefixes [[buffer(1)]],
                         device const NumQ* suffix [[buffer(2)]],
                         device const uint* revtab [[buffer(3)]],
                         device const uint* offs [[buffer(4)]],
                         device const ulong* ents [[buffer(5)]],
                         device const NumQ* recHi [[buffer(6)]],
                         device const NumQ* recLo [[buffer(7)]],
                         device atomic_uint* counters [[buffer(8)]],
                         device NumQ* hits [[buffer(9)]],
                         constant KArgs& ka [[buffer(10)]],
                         device uint* cbuf [[buffer(11)]],
                         uint2 gid [[thread_position_in_grid]],
                         uint lane [[thread_index_in_simdgroup]]) {
  GpuSink sink;
  sink.probes = 0;
  sink.cands = 0;
  sink.x = gid.x;
  sink.y = gid.y;
  sink.counters = counters;
  sink.cbuf = cbuf;
  sink.ccap = ka.candCap;
  sink.p = &p;
  sink.revtab = revtab;
  sink.recHi = recHi;
  sink.recLo = recLo;
  sink.hits = hits;
  sink.hcap = ka.hitCap;
  uint over = 0;
  if (gid.x < ka.nSuffix && gid.y < ka.nPrefix) {
    GpuOut out;
    out.p = &p;
    out.offs = offs;
    out.ents = ents;
    out.C = nq_add(prefixes[gid.y], suffix[gid.x], PRM_RHO);
    out.sink = &sink;
    out.over = 0;
    node_requests(out.C, p, revtab, out);
    over = out.over;
  }
  uint pr = simd_sum(sink.probes);
  uint ca = simd_sum(sink.cands);
  uint ov = simd_sum(over);
  if (lane == 0) {
    atomic_fetch_add_explicit(&counters[1], pr, memory_order_relaxed);
    atomic_fetch_add_explicit(&counters[2], ca, memory_order_relaxed);
    if (ov) atomic_fetch_add_explicit(&counters[3], ov, memory_order_relaxed);
  }
}

// indirect dispatch size for check_cands
kernel void prep_indirect(device const uint* counters [[buffer(0)]], device uint* ind [[buffer(1)]],
                          constant KArgs& ka [[buffer(2)]]) {
  uint n = min(counters[4], ka.candCap);
  ind[0] = (n + 255) / 256;
  ind[1] = 1;
  ind[2] = 1;
}

struct CandSink {
  device atomic_uint* counters;
  device NumQ* hits;
  uint hcap;
  void emit(NumQ n) {
    uint i = atomic_fetch_add_explicit(&counters[0], 1u, memory_order_relaxed);
    if (i < hcap) hits[i] = n;
  }
};

kernel void check_cands(constant Params& p [[buffer(0)]],
                        device const NumQ* prefixes [[buffer(1)]],
                        device const NumQ* suffix [[buffer(2)]],
                        device const uint* revtab [[buffer(3)]],
                        device const NumQ* recHi [[buffer(6)]],
                        device const NumQ* recLo [[buffer(7)]],
                        device atomic_uint* counters [[buffer(8)]],
                        device NumQ* hits [[buffer(9)]],
                        constant KArgs& ka [[buffer(10)]],
                        device const uint* cbuf [[buffer(11)]],
                        uint gid [[thread_position_in_grid]]) {
  uint n = min(atomic_load_explicit(&counters[4], memory_order_relaxed), ka.candCap);
  if (gid >= n) return;
  uint x = cbuf[3 * gid], y = cbuf[3 * gid + 1], idx = cbuf[3 * gid + 2];
  NumQ C = nq_add(prefixes[y], suffix[x], PRM_RHO);
  CandSink cs;
  cs.counters = counters;
  cs.hits = hits;
  cs.hcap = ka.hitCap;
  check_candidate(C, idx, p, revtab, recHi, recLo, cs);
}

// Host/device struct layout check.
kernel void layout_check(device uint* out [[buffer(0)]], constant Params& p [[buffer(1)]]) {
  out[0] = sizeof(Params);
  out[1] = sizeof(SegInfo);
  out[2] = sizeof(NumQ);
  out[3] = sizeof(FastDiv);
  out[4] = p.Q;
  out[5] = p.seg[1].LQ;
  out[6] = (uint)(p.qpow[3] & 0xffffffffu);
  out[7] = p.centerDig[1];
  out[8] = (uint)(p.Wmax.l[1] & 0xffffffffu);
  out[9] = p.recLoN;
}

// =============================================================================================
// v2 kernels (class-partitioned search).  Batch table layout: A[0..nBuckets+2), counts are added
// at A[key+2], an inclusive scan makes A[key+1] the start of bucket key, the scatter advances
// A[key+1] to the end of bucket key; afterwards bucket key = [A[key], A[key+1]).
// =============================================================================================
struct V2Args {
  ulong c0;        // first class of the batch
  uint ww;         // classes in this batch
  uint nBuckets;   // w * Q^tp
  uint nB1, subB;  // B1 threads x sub-threads per b1
  uint nPieces, subA;
  uint hitCap, candCap;
  uint nBlocks;    // scan blocks of 1024
  uint pad;
};

kernel void v2_zero(device uint* A [[buffer(0)]], constant V2Args& va [[buffer(1)]],
                    device atomic_uint* counters [[buffer(2)]], uint gid [[thread_position_in_grid]]) {
  if (gid < va.nBuckets + 2) A[gid] = 0;
  if (gid == 0) atomic_store_explicit(&counters[4], 0u, memory_order_relaxed);
}

// cyclic residue range [lo, lo+w) of a CSR keyed by x mod Qr
inline void csr_range(device const uint* offs, ulong Qr, ulong lo, ulong w, thread uint* r) {
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

kernel void v2_bgen(constant Params& p [[buffer(0)]], constant V2Args& va [[buffer(1)]],
                    device const ulong* b1V [[buffer(2)]], device const ulong* b2V [[buffer(3)]],
                    device const uint* b2Idx [[buffer(4)]], device const uint* b2Offs [[buffer(5)]],
                    device atomic_uint* A [[buffer(6)]], device ulong* ents [[buffer(7)]],
                    constant uint& scatter [[buffer(8)]], uint2 gid [[thread_position_in_grid]]) {
  if (gid.x >= va.nB1 || gid.y >= va.subB) return;
  ulong Qr = p.Qr;
  ulong V1 = b1V[gid.x];
  ulong lo = (va.c0 + Qr - QDIV(V1, PRM_R0) % Qr) % Qr;
  uint rg[4];
  csr_range(b2Offs, Qr, lo, va.ww, rg);
  ulong idxBase = (ulong)gid.x * p.nB2;
  for (int part = 0; part < 2; part++)
    for (uint i = rg[2 * part] + gid.y; i < rg[2 * part + 1]; i += va.subB) {
      ulong V = V1 + b2V[i];
      if (V >= PRM_RHO) V -= PRM_RHO;
      uint key = (uint)v2_key(V, va.c0, p);
      if (scatter) {
        uint pos = atomic_fetch_add_explicit(&A[key + 1], 1u, memory_order_relaxed);
        ents[pos] = v2_entry(V, idxBase + b2Idx[i], p);
      } else {
        atomic_fetch_add_explicit(&A[key + 2], 1u, memory_order_relaxed);
      }
    }
}

// inclusive scan, 1024 elements per threadgroup of 256 threads
kernel void scan_local(device uint* A [[buffer(0)]], device uint* blockSums [[buffer(1)]],
                       constant V2Args& va [[buffer(2)]], uint tid [[thread_index_in_threadgroup]],
                       uint tg [[threadgroup_position_in_grid]], uint lane [[thread_index_in_simdgroup]],
                       uint sg [[simdgroup_index_in_threadgroup]]) {
  threadgroup uint simdTot[8];
  uint n = va.nBuckets + 2;
  uint base = tg * 1024 + tid * 4;
  uint v[4];
  uint s = 0;
  for (int j = 0; j < 4; j++) {
    uint x = (base + j < n) ? A[base + j] : 0;
    s += x;
    v[j] = s;
  }
  uint ex = simd_prefix_exclusive_sum(s);
  if (lane == 31) simdTot[sg] = ex + s;
  threadgroup_barrier(mem_flags::mem_threadgroup);
  uint off = 0;
  for (uint j = 0; j < sg; j++) off += simdTot[j];
  off += ex;
  for (int j = 0; j < 4; j++)
    if (base + j < n) A[base + j] = v[j] + off;
  if (tid == 255) blockSums[tg] = off + s;
}

// exclusive scan of the block sums (single threadgroup of 1024)
kernel void scan_blocks(device uint* blockSums [[buffer(1)]], constant V2Args& va [[buffer(2)]],
                        uint tid [[thread_index_in_threadgroup]], uint lane [[thread_index_in_simdgroup]],
                        uint sg [[simdgroup_index_in_threadgroup]]) {
  threadgroup uint simdTot[32];
  uint nb = va.nBlocks;
  uint per = (nb + 1023) / 1024;
  uint start = tid * per;
  uint s = 0;
  for (uint j = 0; j < per; j++)
    if (start + j < nb) s += blockSums[start + j];
  uint ex = simd_prefix_exclusive_sum(s);
  if (lane == 31) simdTot[sg] = ex + s;
  threadgroup_barrier(mem_flags::mem_threadgroup);
  uint off = 0;
  for (uint j = 0; j < sg; j++) off += simdTot[j];
  off += ex;
  for (uint j = 0; j < per; j++)
    if (start + j < nb) {
      uint x = blockSums[start + j];
      blockSums[start + j] = off;
      off += x;
    }
}

kernel void scan_add(device uint* A [[buffer(0)]], device const uint* blockSums [[buffer(1)]],
                     constant V2Args& va [[buffer(2)]], uint gid [[thread_position_in_grid]]) {
  if (gid < va.nBuckets + 2 && gid >= 1024) A[gid] += blockSums[gid / 1024];
}

struct V2GpuSink {
  uint probes, cands, pre, over;
  uint pi, ai;
  device atomic_uint* counters;
  device uint4* cbuf;
  uint ccap;
  // inline fallback when the candidate buffer is full
  constant Params* p;
  device const uint* revtab;
  device const NumQ* recB1;
  device const NumQ* midW;
  device const NumQ* cenW;
  device NumQ* hits;
  uint hcap;
  void overflow() { over++; }
  void emit(NumQ n) {
    uint i = atomic_fetch_add_explicit(&counters[0], 1u, memory_order_relaxed);
    if (i < hcap) hits[i] = n;
  }
  void candidate(NumQ C, ulong idxB) {
    uint i = atomic_fetch_add_explicit(&counters[4], 1u, memory_order_relaxed);
    if (i < ccap) {
      cbuf[i] = uint4(pi, ai, (uint)(idxB & 0xffffffffu), (uint)(idxB >> 32));
    } else {
      NumQ N = v2_rebuild(C, idxB, *p, recB1, midW, cenW);
      if (nq_is_pal(N, *p, revtab)) emit(N);
    }
  }
};

kernel void v2_agen(constant Params& p [[buffer(0)]], constant V2Args& va [[buffer(1)]],
                    device const AhiPiece* ahi [[buffer(2)]], device const NumQ* aloC [[buffer(3)]],
                    device const uint* aloOffs [[buffer(4)]], device const uint* A [[buffer(5)]],
                    device const ulong* ents [[buffer(6)]], device const uint* revtab [[buffer(7)]],
                    device atomic_uint* counters [[buffer(8)]], device uint4* cbuf [[buffer(9)]],
                    device const NumQ* recB1 [[buffer(10)]], device const NumQ* midW [[buffer(11)]],
                    device NumQ* hits [[buffer(12)]], device const NumQ* cenW [[buffer(14)]],
                    device const ulong* ahiRb [[buffer(15)]],
                    uint2 gid [[thread_position_in_grid]], uint lane [[thread_index_in_simdgroup]]) {
  V2GpuSink sink;
  sink.probes = 0;
  sink.cands = 0;
  sink.pre = 0;
  sink.over = 0;
  sink.counters = counters;
  sink.cbuf = cbuf;
  sink.ccap = va.candCap;
  sink.p = &p;
  sink.revtab = revtab;
  sink.recB1 = recB1;
  sink.midW = midW;
  sink.cenW = cenW;
  sink.hits = hits;
  sink.hcap = va.hitCap;
  uint nodes = 0;
  if (gid.x < va.nPieces && gid.y < va.subA) {
    // the class base alone decides whether this piece has outer parts in the batch;
    // the full 112-byte piece is loaded only when it does
    ulong Qr = p.Qr;
    ulong base = (ahiRb[gid.x] + Qr - va.c0) % Qr;
    ulong lo = (base + Qr - (va.ww - 1)) % Qr;
    uint rg[4];
    csr_range(aloOffs, Qr, lo, va.ww, rg);
    if (rg[0] + gid.y < rg[1] || rg[2] + gid.y < rg[3]) {
      device const AhiPiece& pc = ahi[gid.x];
      sink.pi = gid.x;
      NumQ pa = pc.a, pb = pc.b, pC = pc.C;
      uint seg = pc.seg;
      for (int part = 0; part < 2; part++)
        for (uint i = rg[2 * part] + gid.y; i < rg[2 * part + 1]; i += va.subA) {
          NumQ C = nq_add(pC, aloC[i], PRM_RHO);
          sink.ai = i;
          nodes++;
          v2_node(C, pa, pb, seg, va.c0, p, revtab, A, ents, sink);
        }
    }
  }
  uint nd = simd_sum(nodes);
  uint pr = simd_sum(sink.probes);
  uint ca = simd_sum(sink.cands);
  uint ov = simd_sum(sink.over);
  uint pf = simd_sum(sink.pre);
  if (lane == 0) {
    atomic_fetch_add_explicit(&counters[5], nd, memory_order_relaxed);
    atomic_fetch_add_explicit(&counters[6], pf, memory_order_relaxed);
    atomic_fetch_add_explicit(&counters[1], pr, memory_order_relaxed);
    atomic_fetch_add_explicit(&counters[2], ca, memory_order_relaxed);
    if (ov) atomic_fetch_add_explicit(&counters[3], ov, memory_order_relaxed);
  }
}

kernel void v2_check(constant Params& p [[buffer(0)]], constant V2Args& va [[buffer(1)]],
                     device const AhiPiece* ahi [[buffer(2)]], device const NumQ* aloC [[buffer(3)]],
                     device const uint* revtab [[buffer(7)]], device atomic_uint* counters [[buffer(8)]],
                     device const uint4* cbuf [[buffer(9)]], device const NumQ* recB1 [[buffer(10)]],
                     device const NumQ* midW [[buffer(11)]], device NumQ* hits [[buffer(12)]],
                     device const NumQ* cenW [[buffer(14)]], uint gid [[thread_position_in_grid]]) {
  uint n = min(atomic_load_explicit(&counters[4], memory_order_relaxed), va.candCap);
  if (gid >= n) return;
  uint4 c = cbuf[gid];
  ulong idxB = (ulong)c.z | ((ulong)c.w << 32);
  NumQ N = v2_rebuild(nq_add(ahi[c.x].C, aloC[c.y], PRM_RHO), idxB, p, recB1, midW, cenW);
  if (nq_is_pal(N, p, revtab)) {
    uint i = atomic_fetch_add_explicit(&counters[0], 1u, memory_order_relaxed);
    if (i < va.hitCap) hits[i] = N;
  }
}

kernel void v2_indirect(device const uint* counters [[buffer(8)]], device uint* ind [[buffer(13)]],
                        constant V2Args& va [[buffer(1)]]) {
  uint n = min(counters[4], va.candCap);
  ind[0] = (n + 255) / 256;
  ind[1] = 1;
  ind[2] = 1;
}
