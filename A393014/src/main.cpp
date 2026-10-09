// main.cpp - dbpal: search for numbers that are palindromic in two bases (optionally prime).
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <map>
#include <set>
#include <string>
#include <thread>
#include <vector>

#include "engine.h"
#include "engine2.h"
#include "gpu.h"
#include "plan.h"
#include "table.h"

static void usage() {
  fprintf(stderr,
          "dbpal - numbers palindromic in two bases (OEIS A393014, A259385, A046472..A046484, ...)\n\n"
          "usage: dbpal <command> [options]\n"
          "commands:\n"
          "  search   exhaustive search, one base-P length at a time\n"
          "  plan     show the parameters chosen for each length\n"
          "  brute    exhaustive reference search (all base-P palindromes)\n"
          "  selftest compare the table method against brute force on small lengths\n"
          "  selftest2 compare the v2 (class-partitioned) method against brute force\n"
          "  check    test numbers given on the command line\n\n"
          "options:\n"
          "  -b, --bases B1,B2     the two bases (default 2,9)\n"
          "  --pal-base B          base whose palindromes are enumerated/tabulated (default: auto)\n"
          "  --Lmin N --Lmax N     range of lengths in the enumerated base\n"
          "  --below X             search all N < X  (X like 10^40, 2^130, 1e35)\n"
          "  --above X             start at the length containing X\n"
          "  --primes              only primes (skips even lengths; output is filtered)\n"
          "  --engine cpu|gpu      (default gpu if available)\n"
          "  --table-max N         max table entries (default: from RAM; 8 bytes each)\n"
          "  --mem GB              cap the planned memory (tables, batch tables), e.g. beside a GPU run\n"
          "  --k K --t T           force depth / key digits\n"
          "  --algo v1|v2|auto     v1: one big table; v2: class-partitioned (coprime bases); default auto\n"
          "  --khi N --h1 N --r N --tp N --batch N   v2 overrides\n"
          "  --part i/n            search only part i of n of every length (split work across machines;\n"
          "                        a length is complete once all n parts are logged)\n"
          "  --threads N           CPU threads (default all)\n"
          "  --out FILE            append verified results (TSV) to FILE\n"
          "  --log FILE            append one JSON line per finished length to FILE\n"
          "  -v                    verbose\n");
}

static mpz_class parse_num(const std::string& s) {
  // forms: 12345, a^b, aeb
  auto pos = s.find('^');
  if (pos != std::string::npos) {
    mpz_class b(s.substr(0, pos)), r;
    mpz_pow_ui(r.get_mpz_t(), b.get_mpz_t(), std::stoul(s.substr(pos + 1)));
    return r;
  }
  pos = s.find_first_of("eE");
  if (pos != std::string::npos) {
    mpz_class m(s.substr(0, pos)), r;
    mpz_ui_pow_ui(r.get_mpz_t(), 10, std::stoul(s.substr(pos + 1)));
    return m * r;
  }
  return mpz_class(s);
}

struct Hit {
  mpz_class n;
  bool prime;
};

static std::string now_iso() {
  time_t t = time(nullptr);
  char buf[64];
  strftime(buf, sizeof buf, "%Y-%m-%dT%H:%M:%S", localtime(&t));
  return buf;
}

static int real_main(int argc, char** argv);

int main(int argc, char** argv) {
  try {
    return real_main(argc, argv);
  } catch (std::exception& e) {
    fprintf(stderr, "dbpal: error: %s\n", e.what());
    return 3;
  }
}

static int real_main(int argc, char** argv) {
  if (argc < 2) {
    usage();
    return 1;
  }
  std::string cmd = argv[1];
  SearchOpts o;
  u32 B1 = 2, B2 = 9;
  int palBase = -1;
  int Lmin = -1, Lmax = -1;
  std::string below, above, outFile, logFile, engine = "auto", algo = "auto";
  double memCapGB = 0;  // --mem: cap on planned memory (GB), e.g. for a CPU campaign beside a GPU one
  V2Opts vo;
  std::vector<std::string> rest;
  for (int i = 2; i < argc; i++) {
    std::string a = argv[i];
    auto next = [&]() -> std::string {
      if (i + 1 >= argc) {
        fprintf(stderr, "missing value for %s\n", a.c_str());
        exit(1);
      }
      return argv[++i];
    };
    if (a == "-b" || a == "--bases") {
      std::string v = next();
      if (sscanf(v.c_str(), "%u,%u", &B1, &B2) != 2) {
        fprintf(stderr, "bad --bases\n");
        return 1;
      }
    } else if (a == "--pal-base")
      palBase = std::stoi(next());
    else if (a == "--Lmin")
      Lmin = std::stoi(next());
    else if (a == "--Lmax")
      Lmax = std::stoi(next());
    else if (a == "--below")
      below = next();
    else if (a == "--above")
      above = next();
    else if (a == "--primes")
      o.primesOnly = true;
    else if (a == "--engine")
      engine = next();
    else if (a == "--table-max")
      o.tableMax = (u64)parse_num(next()).get_d();
    else if (a == "--k")
      o.forceK = std::stoi(next());
    else if (a == "--t")
      o.forceT = std::stoi(next());
    else if (a == "--threads")
      o.threads = std::stoi(next());
    else if (a == "--out")
      outFile = next();
    else if (a == "--log")
      logFile = next();
    else if (a == "-v")
      o.verbose = true;
    else if (a == "--algo")
      algo = next();
    else if (a == "--khi")
      vo.khi = std::stoi(next());
    else if (a == "--h1")
      vo.h1 = std::stoi(next());
    else if (a == "--r")
      vo.r = std::stoi(next());
    else if (a == "--tp")
      vo.tp = std::stoi(next());
    else if (a == "--batch")
      vo.batchEntries = parse_num(next()).get_d();
    else if (a == "--mem")
      memCapGB = std::stod(next());
    else if (a == "--part") {
      std::string v = next();
      if (sscanf(v.c_str(), "%u/%u", &o.partI, &o.partN) != 2 || o.partN == 0 || o.partI >= o.partN) {
        fprintf(stderr, "bad --part (want i/n with 0 <= i < n)\n");
        return 1;
      }
    }
    else if (a == "-h" || a == "--help") {
      usage();
      return 0;
    } else
      rest.push_back(a);
  }
  if (B1 < 2 || B2 < 2 || B1 == B2 || B1 > 36 || B2 > 36) {
    fprintf(stderr, "bases must be distinct and in 2..36\n");
    return 1;
  }
  // roles: Q = base checked through its top digits (prefer a power of two, else the smaller base)
  auto pow2 = [](u32 b) { return (b & (b - 1)) == 0; };
  if (palBase > 0) {
    o.P = (u32)palBase;
    o.Q = (o.P == B1) ? B2 : B1;
  } else if (pow2(B1) != pow2(B2)) {
    o.Q = pow2(B1) ? B1 : B2;
    o.P = pow2(B1) ? B2 : B1;
  } else {
    o.Q = std::min(B1, B2);
    o.P = std::max(B1, B2);
  }
  if (o.threads <= 0) o.threads = (int)std::thread::hardware_concurrency();

  if (cmd == "check") {
    for (auto& s : rest) {
      mpz_class n = parse_num(s);
      bool p1 = is_pal_base(n, B1), p2 = is_pal_base(n, B2);
      int pr = mpz_probab_prime_p(n.get_mpz_t(), 40);
      printf("%s: base%u %s (%s), base%u %s (%s), %s\n", n.get_str().c_str(), B1, to_base(n, B1).c_str(),
             p1 ? "pal" : "not pal", B2, to_base(n, B2).c_str(), p2 ? "pal" : "not pal",
             pr ? "prime" : "composite");
    }
    return 0;
  }

  // ---- length range ---------------------------------------------------------------------
  if (!below.empty()) {
    mpz_class X = parse_num(below) - 1;
    Lmax = (int)ndigits_base(X, o.P);
  }
  if (!above.empty()) Lmin = (int)ndigits_base(parse_num(above), o.P);
  if (Lmin < 0) Lmin = 1;
  if (Lmax < 0) Lmax = Lmin;

  fprintf(stderr, "# dbpal %s: bases (%u,%u)  enumerate base P=%u, check base Q=%u  lengths L=%d..%d%s\n",
          cmd.c_str(), B1, B2, o.P, o.Q, Lmin, Lmax, o.primesOnly ? "  [primes only]" : "");

  if (cmd == "plan") {
    o.gpu = engine != "cpu";
    if (o.partN > 1) {  // same deterministic planning as a split search
      if (o.tableMax == 0) o.tableMax = 1000000000ull;
      o.memBudget = 12e9;
      o.gpu = true;
    }
    if (vo.k < 0) vo.k = o.forceK;
    for (int L = Lmin; L <= Lmax; L++) {
      bool small = false;
      try {
        Plan pl = make_plan(o, (u32)L);
        print_plan(pl, stdout);
        small = pl.brute || pl.skip;
        if (!small) printf("        table memory %.2f GB\n", (table_offs_bytes(pl) + table_ents_bytes(pl)) / 1e9);
      } catch (std::exception& e) {
        printf("  L=%d v1: %s\n", L, e.what());
      }
      if (!small) {
        SearchOpts o2 = o;
        o2.forceK = -1;
        Plan2 p2 = plan2_choose(o2, vo, (u32)L);
        print_plan2(p2, stdout);
      }
    }
    return 0;
  }

  bool useGpu = false;
  GpuEngine* gpu = nullptr;
  if ((cmd == "search" || cmd == "selftest2") && engine != "cpu") {
    gpu = gpu_create(o.verbose);
    useGpu = gpu != nullptr;
    if (!useGpu && engine == "gpu") {
      fprintf(stderr, "GPU engine unavailable\n");
      return 1;
    }
  }
  o.gpu = useGpu;
  if (memCapGB > 0) {
    o.memBudget = memCapGB * 1e9;
    if (o.tableMax == 0) o.tableMax = (u64)(memCapGB * 1e9 / 9.0);
  }
  if (o.partN > 1) {
    // every machine searching a part of the same length must choose the same plan:
    // fixed budgets and the GPU cost model, independent of this machine's RAM / engine
    if (o.tableMax == 0) o.tableMax = 1000000000ull;
    o.memBudget = 12e9;
    o.gpu = true;
  }
  fprintf(stderr, "# engine: %s, %d CPU threads\n", useGpu ? gpu_name(gpu) : "cpu", o.threads);

  FILE* fout = outFile.empty() ? nullptr : fopen(outFile.c_str(), "a");
  FILE* flog = logFile.empty() ? nullptr : fopen(logFile.c_str(), "a");
  std::vector<Hit> all;
  int failures = 0;

  for (int L = Lmin; L <= Lmax; L++) {
    if (cmd == "selftest") {
      // table method with every feasible k, compared with brute force
      Plan bp = make_plan(o, (u32)L);
      if (bp.skip) {
        printf("L=%d skipped (%s)\n", L, bp.skipReason.c_str());
        continue;
      }
      std::vector<NumQ> ref;
      run_brute(bp, o.threads, ref);
      std::set<std::string> refs;
      for (auto& n : ref) refs.insert(from_numq(n, bp.prm).get_str());
      for (int k = 1; 2 * k < L; k++) {
        SearchOpts o2 = o;
        o2.forceK = k;
        Plan pl;
        try {
          pl = make_plan(o2, (u32)L);
        } catch (std::exception& e) {
          continue;
        }
        if (pl.skip || pl.brute) continue;
        std::vector<u32> offs(pl.prm.Dt + 1);
        std::vector<u64> ents(pl.tableEntries);
        build_table(pl, TableView{offs.data(), ents.data()}, o.threads);
        std::vector<NumQ> got;
        RunStats st = run_cpu(pl, TableView{offs.data(), ents.data()}, o.threads, got);
        std::set<std::string> gs;
        for (auto& n : got) gs.insert(from_numq(n, pl.prm).get_str());
        bool ok = gs == refs && st.overflow == 0;
        if (!ok) failures++;
        printf("L=%2d k=%2d t=%2u table=%-10llu nodes=%-12llu cands=%-10llu found=%zu ref=%zu %s\n", L, k,
               pl.prm.t, (unsigned long long)pl.tableEntries, (unsigned long long)st.nodes,
               (unsigned long long)st.cands, gs.size(), refs.size(), ok ? "OK" : "MISMATCH");
        if (!ok) {
          for (auto& s : refs)
            if (!gs.count(s)) printf("   missing %s\n", s.c_str());
          for (auto& s : gs)
            if (!refs.count(s)) printf("   extra   %s\n", s.c_str());
        }
      }
      continue;
    }

    if (cmd == "selftest2") {
      Plan bp = make_plan(o, (u32)L);
      if (bp.skip) continue;
      std::vector<NumQ> ref;
      run_brute(bp, o.threads, ref);
      std::set<std::string> refs;
      for (auto& n : ref) refs.insert(from_numq(n, bp.prm).get_str());
      int tried = 0;
      for (int k = 2; 2 * k < L; k++)
        for (int khi = 1; khi < k; khi++)
          for (int h1 = 0; h1 <= (L - 2 * k) / 2; h1++)
            for (int r = 1; r <= 12; r += 3) {
              V2Opts v2 = vo;
              v2.k = k;
              v2.khi = khi;
              v2.h1 = h1;
              v2.r = r;
              v2.batchEntries = 50 + 37 * r;
              Plan2 p2 = plan2_choose(o, v2, (u32)L);
              if (!p2.ok) continue;
              plan2_build(p2, o.threads);
              std::vector<NumQ> got;
              RunStats st = useGpu ? run_gpu2(gpu, p2, o.threads, got, false) : run_cpu2(p2, o.threads, got);
              std::set<std::string> gs;
              for (auto& n : got) gs.insert(from_numq(n, p2.base.prm).get_str());
              bool ok = gs == refs && st.overflow == 0;
              tried++;
              if (!ok) {
                failures++;
                printf("L=%d k=%d khi=%d h1=%d r=%d tp=%u w=%u: found %zu ref %zu overflow %llu MISMATCH\n", L, k,
                       khi, h1, r, p2.base.prm.tp, p2.base.prm.w, gs.size(), refs.size(),
                       (unsigned long long)st.overflow);
                for (auto& s : refs)
                  if (!gs.count(s)) printf("   missing %s\n", s.c_str());
                for (auto& s : gs)
                  if (!refs.count(s)) printf("   extra   %s\n", s.c_str());
              }
            }
      printf("L=%d: %d v2 configurations checked against brute force (%zu solutions)\n", L, tried, refs.size());
      continue;
    }

    Plan pl = make_plan(o, (u32)L);
    if (o.verbose || cmd == "search") print_plan(pl, stderr);
    if (pl.skip) {
      if (flog)
        fprintf(flog,
                "{\"bases\":[%u,%u],\"P\":%u,\"Q\":%u,\"L\":%d,\"primesOnly\":%s,\"status\":\"skipped\",\"reason\":\"%s\","
                "\"time\":\"%s\"}\n",
                B1, B2, o.P, o.Q, L, o.primesOnly ? "true" : "false", pl.skipReason.c_str(), now_iso().c_str());
      continue;
    }
    std::vector<NumQ> raw;
    RunStats st;
    double tbuild = 0;
    const char* method = "brute";
    bool useV2 = false;
    Plan2 p2;
    if (!pl.brute && cmd == "search" && algo != "v1" && pl.prm.gmode == 0) {
      p2 = plan2_choose(o, vo, (u32)L);
      useV2 = p2.ok && (algo == "v2" || p2.estSeconds < pl.costEst);
      if (algo == "v2" && !p2.ok) fprintf(stderr, "  v2 not applicable: %s\n", p2.why.c_str());
      if (p2.ok && (o.verbose || useV2)) print_plan2(p2, stderr);
    }
    std::string planId;
    if (useV2) {
      const Params& q = p2.base.prm;
      char b[256];
      snprintf(b, sizeof b, "v2 k=%u khi=%u h1=%u r0=%u r=%u tp=%u w=%u Qr=%llu", q.k, q.khi, q.h1, q.r0, q.r, q.tp,
               q.w, (unsigned long long)q.Qr);
      planId = b;
    } else if (!pl.brute) {
      char b[128];
      snprintf(b, sizeof b, "v1 k=%u t=%u", pl.prm.k, pl.prm.t);
      planId = b;
    }
    if (o.partN > 1 && !pl.brute && !useV2) {
      fprintf(stderr, "dbpal: --part needs the v2 method (L=%d would use v1: its prefix split depends on the engine)\n", L);
      return 4;
    }
    if (pl.brute || cmd == "brute") {
      st = run_brute(pl.brute ? pl : make_plan(o, (u32)L), o.threads, raw);
    } else if (useV2) {
      method = useGpu ? "gpu-v2" : "cpu-v2";
      auto tb0 = std::chrono::steady_clock::now();
      plan2_build(p2, o.threads);
      tbuild = std::chrono::duration<double>(std::chrono::steady_clock::now() - tb0).count();
      st = useGpu ? run_gpu2(gpu, p2, o.threads, raw, o.verbose) : run_cpu2(p2, o.threads, raw);
    } else if (useGpu) {
      method = "gpu";
      st = run_gpu(gpu, pl, o.threads, raw, &tbuild, o.verbose);
    } else {
      method = "cpu";
      std::vector<u32> offs(pl.prm.Dt + 1);
      std::vector<u64> ents(pl.tableEntries);
      tbuild = build_table(pl, TableView{offs.data(), ents.data()}, o.threads);
      st = run_cpu(pl, TableView{offs.data(), ents.data()}, o.threads, raw);
    }
    // verify independently with GMP
    std::set<std::string> seen;
    std::vector<Hit> hits;
    for (auto& r : raw) {
      mpz_class n = from_numq(r, pl.prm);
      std::string key = n.get_str();
      if (seen.count(key)) continue;
      seen.insert(key);
      bool ok = is_pal_base(n, o.P) && is_pal_base(n, o.Q) && ndigits_base(n, o.P) == (u32)L;
      if (ok && o.verbose && useV2) fprintf(stderr, "    hit %s\n", key.c_str());
      if (!ok) {
        fprintf(stderr, "!! engine reported a non-solution: %s\n", key.c_str());
        failures++;
        continue;
      }
      bool pr = mpz_probab_prime_p(n.get_mpz_t(), 40) > 0;
      if (o.primesOnly && !pr) continue;
      hits.push_back({n, pr});
    }
    std::sort(hits.begin(), hits.end(), [](const Hit& a, const Hit& b) { return a.n < b.n; });
    for (auto& h : hits) {
      printf("%s\t%s\tL%u=%u\tL%u=%u\t%s\t%s\n", h.n.get_str().c_str(), h.prime ? "prime" : "-", o.P,
             ndigits_base(h.n, o.P), o.Q, ndigits_base(h.n, o.Q), to_base(h.n, B1).c_str(),
             to_base(h.n, B2).c_str());
      if (fout)
        fprintf(fout, "%s\t%d\t%u\t%u\t%s\t%s\n", h.n.get_str().c_str(), h.prime ? 1 : 0, B1, B2,
                to_base(h.n, B1).c_str(), to_base(h.n, B2).c_str());
      all.push_back(h);
    }
    fflush(stdout);
    if (fout) fflush(fout);
    fprintf(stderr,
            "  L=%d done [%s]: %zu hits  nodes=%.3g probes=%.3g cands=%.3g pruned=%.3g  table %.1fs  search %.1fs%s\n",
            L, method, hits.size(), (double)st.nodes, (double)st.probes, (double)st.cands, (double)st.pruned,
            tbuild, st.seconds, st.overflow ? "  (OVERFLOW!)" : "");
    if (st.overflow) failures++;
    if (flog) {
      bool isPart = o.partN > 1 && !pl.brute;
      char partStr[64] = "";
      if (isPart) snprintf(partStr, sizeof partStr, ",\"part\":\"%u/%u\"", o.partI, o.partN);
      std::string partPlan = isPart ? ",\"plan\":\"" + planId + "\"" : std::string();
      fprintf(flog,
              "{\"bases\":[%u,%u],\"P\":%u,\"Q\":%u,\"L\":%d,\"primesOnly\":%s,\"status\":\"%s\"%s,\"method\":\"%s\","
              "\"k\":%u,\"hits\":%zu,\"nodes\":%llu,\"cands\":%llu,\"tableSeconds\":%.2f,\"seconds\":%.2f,"
              "\"upto\":\"%s\",\"time\":\"%s\"}\n",
              B1, B2, o.P, o.Q, L, o.primesOnly ? "true" : "false", isPart ? "part" : "done",
              (std::string(partStr) + partPlan).c_str(), method, pl.prm.k, hits.size(),
              (unsigned long long)st.nodes, (unsigned long long)st.cands, tbuild, st.seconds,
              pl.highN.get_str().c_str(), now_iso().c_str());
      fflush(flog);
    }
  }
  if (cmd == "selftest" || cmd == "selftest2") printf("%s: %s\n", cmd.c_str(), failures ? "FAILURES" : "all OK");
  if (fout) fclose(fout);
  if (flog) fclose(flog);
  if (gpu) gpu_destroy(gpu);
  return failures ? 2 : 0;
}
