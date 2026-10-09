/* ipn_search.c — exhaustive search for infinitary perfect numbers N <= X.
 *
 * N is infinitary perfect iff, writing N as a product of distinct
 * Fermi-Dirac primes q = p^(2^k) (p prime, k >= 0), prod(q+1) = 2*prod(q).
 *
 * Branch-and-bound over FD factorizations with exact rational bookkeeping:
 * a node holds chosen components (strictly increasing), the needed remaining
 * ratio r = a/b = prod_future (q+1)/q as a signed factor list, and the budget
 * B = X / prod(chosen).  a | 2*prod(chosen) <= 2X so a,b fit in u64.
 *
 * Branch/prune rules (each sound, see README):
 *  - a == b  -> complete solution;  a < b -> dead.
 *  - closure: if (a-b) | b, q* = b/(a-b) is the unique 1-component finish.
 *  - every prime P | b needs a future component P^(2^k) in (f, B]; the
 *    product of those minima over distinct P must be <= B.
 *  - every prime P | a needs a future q ≡ -1 (mod P), q in (f, B].
 *  - qmin = max(f+1, ceil(b/(a-b))); if qmin^3 > B then <= 2 components
 *    remain: solved exactly (closure + 2-component divisor/scan method).
 *  - else scan smallest-next-component q in [b/(a-b)+1 .. mroof*a/(a-b)]
 *    where mroof = floor(log_qmin B)  (from (1+1/q)^m >= r necessity).
 *
 * Build: cc -O3 -I/opt/homebrew/include -L/opt/homebrew/lib \
 *        -o ipn_search ipn_search.c -lprimesieve -lpthread
 * Usage: ./ipn_search -x <bound> [-t threads] [-c ckpt] [-g task_target]
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <inttypes.h>
#include <math.h>
#include <time.h>
#include <signal.h>
#include <pthread.h>
#include <stdatomic.h>
#include <primesieve.h>

typedef uint64_t u64;
typedef uint8_t u8;
typedef unsigned __int128 u128;

#define MAXCOMP 64
#define MAXR    256
#define PLIMIT  16777216ULL

static u128 X;                      /* inclusive search bound */
static u64 *PR; static size_t NPR;  /* primes <= PLIMIT */

static char *sprint128(char *buf, u128 v){   /* buf >= 48 bytes */
    char t[48]; int i=0;
    do { t[i++]='0'+(int)(v%10); v/=10; } while(v);
    char *p=buf;
    while(i) *p++ = t[--i];
    *p=0;
    return buf;
}
static u128 parse128(const char *s){
    u128 v=0;
    while(*s>='0' && *s<='9'){ v = v*10 + (u128)(*s-'0'); s++; }
    return v;
}
static u64 isqrt128(u128 n){
    u64 r=(u64)sqrtl((long double)n);
    while(r>0 && (u128)r*r > n) r--;
    while((u128)(r+1)*(r+1) <= n) r++;
    return r;
}
static volatile sig_atomic_t g_stop = 0;

static _Atomic long long st_nodes, st_closure, st_twocomp, st_divenum,
                         st_fallback, st_warn;

/* --------------------------- arithmetic --------------------------- */
static inline u64 mulmod(u64 a, u64 b, u64 m){ return (u64)((u128)a*b % m); }
static u64 powmod(u64 a, u64 e, u64 m){
    u64 r=1; a%=m;
    while(e){ if(e&1) r=mulmod(r,a,m); a=mulmod(a,a,m); e>>=1; }
    return r;
}
static u64 gcd64(u64 a, u64 b){ while(b){ u64 t=a%b; a=b; b=t; } return a; }

static int is_prime_u64(u64 n){
    static const u64 W[12]={2,3,5,7,11,13,17,19,23,29,31,37};
    if(n<2) return 0;
    for(int i=0;i<12;i++){ if(n==W[i]) return 1; if(n%W[i]==0) return 0; }
    u64 d=n-1; int s=0; while(!(d&1)){ d>>=1; s++; }
    for(int i=0;i<12;i++){
        u64 x=powmod(W[i],d,n);
        if(x==1||x==n-1) continue;
        int ok=0;
        for(int j=1;j<s;j++){ x=mulmod(x,x,n); if(x==n-1){ ok=1; break; } }
        if(!ok) return 0;
    }
    return 1;
}

static u64 brent_rho(u64 n){
    if(!(n&1)) return 2;
    if(n%3==0) return 3;
    for(u64 c=1;;c++){
        u64 x=2, y=2, d=1, q=1, ys=y, r=1;
        const u64 m=128;
        while(d==1){
            x=y;
            for(u64 i=0;i<r;i++) y=(mulmod(y,y,n)+c)%n;
            for(u64 k=0;k<r && d==1;k+=m){
                ys=y;
                u64 lim = (m < r-k) ? m : r-k;
                for(u64 i=0;i<lim;i++){
                    y=(mulmod(y,y,n)+c)%n;
                    u64 t = x>y ? x-y : y-x;
                    if(t) q=mulmod(q,t,n);
                }
                d=gcd64(q,n);
            }
            r<<=1;
        }
        if(d==n){
            d=1;
            while(d==1){
                ys=(mulmod(ys,ys,n)+c)%n;
                u64 t = x>ys ? x-ys : ys-x;
                d = t ? gcd64(t,n) : n;
            }
        }
        if(d!=n) return d;
        /* cycle degenerated: retry with next c */
    }
}

/* factorization of n into (prime, exp) pairs, sorted by prime */
typedef struct { u64 p; int e; } pe_t;

static void fac_insert(pe_t *f, int *nf, u64 p, int e){
    int i=0;
    while(i<*nf && f[i].p<p) i++;
    if(i<*nf && f[i].p==p){ f[i].e+=e; return; }
    memmove(f+i+1, f+i, (size_t)(*nf-i)*sizeof(pe_t));
    f[i].p=p; f[i].e=e; (*nf)++;
}

static void factor_rec(u64 n, pe_t *f, int *nf){
    if(n==1) return;
    if(is_prime_u64(n)){ fac_insert(f,nf,n,1); return; }
    u64 d=brent_rho(n);
    factor_rec(d,f,nf);
    factor_rec(n/d,f,nf);
}

static void factor64(u64 n, pe_t *f, int *nf){
    *nf=0;
    for(size_t i=0;i<NPR && PR[i]<3000;i++){
        u64 p=PR[i];
        if(p*p>n) break;
        if(n%p==0){ int e=0; while(n%p==0){ n/=p; e++; } fac_insert(f,nf,p,e); }
    }
    if(n>1) factor_rec(n,f,nf);
}

static u64 isqrt64(u64 n){
    u64 r=(u64)sqrtl((long double)n);
    while(r>0 && r>n/r) r--;
    while((r+1) <= n/(r+1)) r++;
    return r;
}

/* v is a Fermi-Dirac prime (p^(2^k))?  If yes, return base p and k. */
static int fd_decompose(u64 v, u64 *p, int *k){
    if(v<2) return 0;
    int kk=0;
    for(;;){
        u64 s=isqrt64(v);
        if(s*s==v && v>1){ v=s; kk++; if(kk>6) return 0; }
        else break;
    }
    if(!is_prime_u64(v)) return 0;
    *p=v; *k=kk;
    return 1;
}
static int is_fd_prime(u64 v){ u64 p; int k; return fd_decompose(v,&p,&k); }

/* --------------------- 256-bit helpers (for big X) ----------------- */
typedef struct { u128 hi, lo; } u256;

static u256 mul_128_128(u128 a, u128 b){
    u64 a0=(u64)a, a1=(u64)(a>>64), b0=(u64)b, b1=(u64)(b>>64);
    u128 p00=(u128)a0*b0, p01=(u128)a0*b1, p10=(u128)a1*b0, p11=(u128)a1*b1;
    u128 mid = p01 + p10;
    u128 carry_mid = mid < p01 ? ((u128)1<<64) : 0;
    u256 r;
    r.lo = p00 + (mid<<64);
    u128 carry_lo = r.lo < p00 ? 1 : 0;
    r.hi = p11 + (mid>>64) + carry_mid + carry_lo;
    return r;
}

/* n / d and n % d for n < d*2^128 (i.e. n.hi < d); quotient fits u128 */
static int divmod_256_128(u256 n, u128 d, u128 *q, u128 *r){
    if(n.hi >= d) return 0;      /* quotient >= 2^128: caller treats as invalid */
    u128 rem=n.hi, quot=0;
    for(int i=127;i>=0;i--){
        int top = (rem >> 127) != 0;
        rem = (rem<<1) | ((n.lo>>i)&1);
        quot <<= 1;
        if(top || rem >= d){ rem -= d; quot |= 1; }
    }
    *q=quot; *r=rem;
    return 1;
}

/* 128-bit primality: MR with 40 fixed bases via 256-bit mulmod.
 * (Deterministic below 3.3e24; beyond that error < 4^-40 — reported
 *  solutions get independent verification before publishing.) */
static u128 mulmod128(u128 a, u128 b, u128 m){
    u256 p = mul_128_128(a%m, b%m);
    u128 q,r;
    divmod_256_128(p, m, &q, &r);   /* p.hi < m since both factors < m */
    return r;
}
static u128 powmod128(u128 a, u128 e, u128 m){
    u128 r=1; a%=m;
    while(e){ if(e&1) r=mulmod128(r,a,m); a=mulmod128(a,a,m); e>>=1; }
    return r;
}
static int is_prime_128(u128 n){
    if(n <= UINT64_MAX) return is_prime_u64((u64)n);
    static const u64 W[40]={2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,
        67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,
        163,167,173};
    for(int i=0;i<40;i++) if(n%W[i]==0) return 0;
    u128 d=n-1; int s=0; while(!(d&1)){ d>>=1; s++; }
    for(int i=0;i<40;i++){
        u128 x=powmod128(W[i],d,n);
        if(x==1||x==n-1) continue;
        int ok=0;
        for(int j=1;j<s;j++){ x=mulmod128(x,x,n); if(x==n-1){ ok=1; break; } }
        if(!ok) return 0;
    }
    return 1;
}
static u128 isqrt128w(u128 n){
    u128 r=(u128)(u64)sqrtl((long double)n);
    if(n > ((u128)1<<126)) r=(u128)1<<63;
    while(r>0 && r > n/r) r--;
    while(r+1 <= n/(r+1)) r++;
    return r;
}
static int is_fd_prime_128(u128 v){
    if(v <= UINT64_MAX) return is_fd_prime((u64)v);
    int kk=0;
    for(;;){
        u128 s=isqrt128w(v);
        if(s*s==v && v>1){ v=s; kk++; if(kk>7) return 0; }
        else break;
    }
    return is_prime_128(v);
}

/* ----------------------------- state ------------------------------ */
typedef struct {
    u64  comps[MAXCOMP]; int nc;    /* chosen components, increasing */
    pe_t r[MAXR]; int nr;           /* needed ratio prod p^e (e signed) */
    u128 prod;                      /* product of comps (<= X) */
    u64  f;                         /* floor: all future comps > f */
} state_t;

/* mode: emit callback either recurses (search) or collects tasks (gen) */
struct ctx;
typedef void (*emit_fn)(struct ctx *c, const state_t *child);
typedef struct ctx {
    emit_fn emit;
    int gen_limit;                  /* task depth in gen mode, -1 = search */
    void *gen_sink;
    size_t root;                    /* root task index (search mode) — carried
                                       in ctx so scan-helper threads inherit it */
    void *sb;                       /* per-worker sblock scratch (may be NULL) */
} ctx_t;

/* ------------------------- solution store ------------------------- */
static pthread_mutex_t sol_mu = PTHREAD_MUTEX_INITIALIZER;
static struct { u128 n; u64 comps[MAXCOMP]; int nc; u128 extra; } sols[1024];
static int nsols = 0;

/* extra: optional final component > 2^64 (0 = none); it is provably the
 * largest, so it never needs to live inside the u64 search state. */
static void record_solution2(const u64 *comps, int nc, u128 extra, int quiet){
    u64 sorted[MAXCOMP];
    memcpy(sorted, comps, (size_t)nc*sizeof(u64));
    for(int i=1;i<nc;i++){ u64 v=sorted[i]; int j=i;
        while(j>0 && sorted[j-1]>v){ sorted[j]=sorted[j-1]; j--; } sorted[j]=v; }
    u128 n=1, s=1;
    for(int i=0;i<nc;i++){
        if(i && sorted[i]==sorted[i-1]){ fprintf(stderr,"BUG dup comp\n"); abort(); }
        n *= sorted[i]; s *= sorted[i]+1;
    }
    if(extra){ n *= extra; s *= extra+1; }
    if(s != 2*n || n > X){ fprintf(stderr,"BUG bad solution\n"); abort(); }
    pthread_mutex_lock(&sol_mu);
    if(nsols==1024){ pthread_mutex_unlock(&sol_mu); return; }
    for(int i=0;i<nsols;i++)
        if(sols[i].n==n){ pthread_mutex_unlock(&sol_mu); return; }
    sols[nsols].n=n; sols[nsols].nc=nc; sols[nsols].extra=extra;
    memcpy(sols[nsols].comps, sorted, (size_t)nc*sizeof(u64));
    nsols++;
    if(!quiet){
        char nb[48];
        printf("FOUND %s =", sprint128(nb,n));
        for(int i=0;i<nc;i++) printf(" %" PRIu64, sorted[i]);
        if(extra){ char eb[48]; printf(" %s", sprint128(eb,extra)); }
        printf("\n");
        fflush(stdout);
    }
    pthread_mutex_unlock(&sol_mu);
}
static void record_solution(const u64 *comps, int nc, int quiet){
    record_solution2(comps, nc, 0, quiet);
}

/* a,b from factor list; both guaranteed <= 2X (abort = logic bug) */
static void compute_ab(const state_t *st, u128 *pa, u128 *pb){
    u128 a=1, b=1;
    for(int i=0;i<st->nr;i++){
        u64 p=st->r[i].p; int e=st->r[i].e;
        if(e>0){ for(int j=0;j<e;j++){ a*=p; if(a > 4*X){ fprintf(stderr,"BUG a overflow\n"); abort(); } } }
        else if(e<0){ for(int j=0;j<-e;j++){ b*=p; if(b > 4*X){ fprintf(stderr,"BUG b overflow\n"); abort(); } } }
    }
    *pa=a; *pb=b;
}

/* child = parent + component v (= p0^(2^k)); factors of v+1 go negative.
 * make_child2 takes the precomputed factorization of v+1. */
static void make_child2(const state_t *par, u64 v, u64 p0, int k, state_t *ch,
                        const pe_t *vf, int nvf){
    *ch = *par;
    if(ch->nc>=MAXCOMP){ fprintf(stderr,"BUG comp overflow\n"); abort(); }
    ch->comps[ch->nc++]=v;
    ch->prod = par->prod * (u128)v;
    ch->f = v;
    /* merge: +2^k at p0, minus factors of v+1 */
    pe_t tmp[MAXR]; int nt=0;
    memcpy(tmp, par->r, (size_t)par->nr*sizeof(pe_t)); nt=par->nr;
    if(nt+nvf+1 > MAXR){ fprintf(stderr,"BUG r overflow\n"); abort(); }
    fac_insert(tmp,&nt,p0,1<<k);
    for(int i=0;i<nvf;i++) fac_insert(tmp,&nt,vf[i].p,-vf[i].e);
    int m=0;
    for(int i=0;i<nt;i++) if(tmp[i].e!=0) tmp[m++]=tmp[i];
    memcpy(ch->r, tmp, (size_t)m*sizeof(pe_t)); ch->nr=m;
}
static void make_child(const state_t *par, u64 v, u64 p0, int k, state_t *ch){
    pe_t vf[64]; int nvf; factor64(v+1, vf, &nvf);
    make_child2(par, v, p0, k, ch, vf, nvf);
}

/* cheap pre-check before building a child: the largest prime L of v+1,
 * if not already positive in r, lands in the child's denominator and
 * needs an L-power component in (v, B/v].  Reject v when impossible. */
static int child_viable(const state_t *par, u64 v, u128 B, const pe_t *vf, int nvf){
    u64 L = vf[nvf-1].p;
    int pos=0;
    for(int i=0;i<par->nr;i++) if(par->r[i].p==L && par->r[i].e>0){ pos=1; break; }
    if(pos) return 1;
    u128 budget = B / v;
    u128 w = L;
    while(w <= v){ w = w*w; }
    return w <= budget;
}

/* ---------------- exactly-two-components completion ----------------
 * Need q1<q2, both FD primes > f, q1*q2 <= B, (q1+1)(q2+1)*b = a*q1*q2.
 * Identity: ((a-b)q1 - b)((a-b)q2 - b) = a*b, so d1=(a-b)q1-b runs over
 * divisors of a*b with d1 <= b+(a-b) (since q1 <= 2b/(a-b)+1).           */
static void try_pair(const state_t *st, u64 q1, u128 q2, u128 B){
    if(q1<=st->f || q2<=q1) return;
    if(q2 > B/q1) return;
    if(!is_fd_prime(q1) || !is_fd_prime_128(q2)) return;
    u64 c[MAXCOMP]; memcpy(c, st->comps, (size_t)st->nc*sizeof(u64));
    c[st->nc]=q1;
    if(q2 <= UINT64_MAX){
        c[st->nc+1]=(u64)q2;
        record_solution(c, st->nc+2, 0);
    } else
        record_solution2(c, st->nc+1, q2, 0);
}

/* given q1, solve (q1+1)(q2+1)*b = a*q1*q2 exactly:
 * q2 = b(q1+1) / (dq*q1 - b); returns 1 and *q2 on exact integer. */
static int solve_q2(u64 q1, u128 b, u128 dq, u128 *q2){
    u128 den = dq*q1 - b;            /* fits: q1 <= 3b/dq+2 by caller */
    u256 num = mul_128_128(b, (u128)q1+1);
    u128 q,r;
    if(!divmod_256_128(num, den, &q, &r)) return 0;
    if(r) return 0;
    *q2=q;
    return 1;
}

/* divisor DFS over ab (entries |e|), partial products capped at dmax.
 * Returns 0 if the enumeration cap was hit (caller falls back to scan). */
static int div_dfs(const state_t *st, const pe_t *fl, int nfl, int idx,
                   u128 cur, u128 b, u128 dq, u128 B, long long *cnt){
    if(*cnt > (1<<24)) return 0;
    if(idx==nfl){
        (*cnt)++;
        if((cur + b) % dq == 0){
            u128 q1w=(cur+b)/dq;
            if(q1w > st->f && q1w <= UINT64_MAX && q1w*q1w <= B){
                u128 q2;
                if(solve_q2((u64)q1w, b, dq, &q2) && q2 <= X)
                    try_pair(st, (u64)q1w, q2, B);
            }
        }
        return 1;
    }
    u64 p=fl[idx].p; int e=fl[idx].e<0 ? -fl[idx].e : fl[idx].e;
    u128 v=cur;
    u128 dmax = b + dq;   /* d1 bound */
    for(int j=0;j<=e;j++){
        if(v > dmax) break;
        if(!div_dfs(st, fl, nfl, idx+1, v, b, dq, B, cnt)) return 0;
        v *= p;
    }
    return 1;
}

/* wide two-component scans are embarrassingly parallel: candidates are
 * independent and record_solution is mutex-protected.  A global helper
 * budget (init nthreads in main) lets tail-phase scans use idle cores. */
static _Atomic int g_par_avail;

typedef struct { const state_t *st; u128 b, dq, B; u64 hi, start, stride; } tcs_arg_t;
static void *tcs_worker(void *a){
    tcs_arg_t *t=(tcs_arg_t*)a;
    for(u64 q1=t->start; q1<=t->hi && !g_stop; q1+=t->stride){
        if(!is_fd_prime(q1)) continue;
        u128 q2;
        if(solve_q2(q1, t->b, t->dq, &q2) && q2 <= X)
            try_pair(t->st, q1, q2, t->B);
    }
    return NULL;
}

static void two_comp_scan(const state_t *st, u128 b, u128 dq, u128 B, u64 lo, u64 hi){
    if(hi < lo) return;
    u64 width = hi - lo + 1;
    int extra = 0;
    if(width > 2000000){
        int cur = atomic_load(&g_par_avail);
        while(cur > 0){
            int take = cur < 15 ? cur : 15;
            if(atomic_compare_exchange_weak(&g_par_avail, &cur, cur - take)){ extra = take; break; }
        }
        fprintf(stderr,"WIDE scan width=%" PRIu64 " helpers=%d\n", width, extra);
    }
    if(extra == 0){
        for(u64 q1=lo; q1<=hi && !g_stop; q1++){
            if(!is_fd_prime(q1)) continue;
            u128 q2;
            if(solve_q2(q1, b, dq, &q2) && q2 <= X)
                try_pair(st, q1, q2, B);
        }
        return;
    }
    int n = extra + 1;
    pthread_t th[16]; tcs_arg_t args[16];
    for(int k=0;k<n;k++)
        args[k] = (tcs_arg_t){ st, b, dq, B, hi, lo + (u64)k, (u64)n };
    for(int k=1;k<n;k++) pthread_create(&th[k], NULL, tcs_worker, &args[k]);
    tcs_worker(&args[0]);
    for(int k=1;k<n;k++) pthread_join(th[k], NULL);
    atomic_fetch_add(&g_par_avail, extra);
}

static void two_comp_resolve(const state_t *st, u128 a, u128 b, u128 dq, u128 B){
    atomic_fetch_add_explicit(&st_twocomp,1,memory_order_relaxed);
    u128 lo128 = b/dq;                   /* q1 > lo */
    if(st->f > lo128) lo128 = st->f;
    u128 hi128 = 3*b/dq + 2;
    u64 sB = isqrt128(B);
    if(hi128 > sB) hi128 = sB;
    if(hi128 <= lo128) return;
    if(lo128 >= UINT64_MAX){
        atomic_fetch_add_explicit(&st_warn,1,memory_order_relaxed);
        return;
    }
    u64 lo=(u64)lo128, hi=(u64)hi128;
    u64 width = hi - lo;
    /* estimated divisor count of a*b */
    u128 tau=1;
    for(int i=0;i<st->nr;i++){
        int e=st->r[i].e<0?-st->r[i].e:st->r[i].e;
        tau *= (u64)(e+1);
        if(tau > (u128)1<<40) break;
    }
    (void)a;
    if(width <= 200000 || (u128)width*4 < tau){
        two_comp_scan(st, b, dq, B, lo+1, hi);
    } else {
        atomic_fetch_add_explicit(&st_divenum,1,memory_order_relaxed);
        long long cnt=0;
        if(!div_dfs(st, st->r, st->nr, 0, 1, b, dq, B, &cnt)){
            atomic_fetch_add_explicit(&st_warn,1,memory_order_relaxed);
            two_comp_scan(st, b, dq, B, lo+1, hi);
        }
    }
}

/* ------------------------- main node logic ------------------------ */
static void expand_node(ctx_t *cx, const state_t *st);
static void scan_candidate(ctx_t *cx, const state_t *st, u64 v, u128 B);

static void scan_candidate(ctx_t *cx, const state_t *st, u64 v, u128 B);


/* callers always know the FD decomposition v = p0^(2^k); no re-check */
static void scan_candidate2(ctx_t *cx, const state_t *st, u64 v, u64 p0, int k, u128 B){
    pe_t vf[64]; int nvf; factor64(v+1, vf, &nvf);
    if(!child_viable(st, v, B, vf, nvf)) return;
    state_t ch;
    make_child2(st, v, p0, k, &ch, vf, nvf);
    cx->emit(cx, &ch);
}
static void scan_candidate(ctx_t *cx, const state_t *st, u64 v, u128 B){
    scan_candidate2(cx, st, v, v, 0, B);   /* v is a known prime */
}

/* ------------- batched scan: interval-factor v+1 for a block ----------
 * For prime candidates v in [A, Aend), factor all v+1 together by sieving
 * small primes over the interval, then finishing residuals with MR/rho.
 * 10-50x cheaper than per-candidate rho on the giant near-miss scans.   */
#define SBW 262144            /* block width */
#define SP0 65536             /* sieve small primes below this */
#define MAXVF 24
typedef struct {
    u64  res[SBW];            /* remaining cofactor of v+1 */
    u8   cand[SBW];
    pe_t vf[SBW/8][MAXVF];    /* factor lists, indexed per candidate slot */
    int  nvf[SBW/8];
    int  slot[SBW];           /* position -> candidate slot */
} sblock_t;

/* candidates listed in cands[nc]; emits scan_candidate-equivalent children */
static void process_block(ctx_t *cx, const state_t *st, u128 B,
                          u64 A, u64 Aend, const u64 *cands, int nc, sblock_t *sb){
    if(nc==0) return;
    u64 lo1 = A+1, hi1 = Aend;    /* v+1 range: [A+1, Aend] inclusive */
    memset(sb->cand, 0, (size_t)(hi1-lo1+1));
    for(int c=0;c<nc;c++){
        size_t off = cands[c]+1-lo1;
        sb->cand[off]=1; sb->slot[off]=c; sb->nvf[c]=0; sb->res[off]=cands[c]+1;
    }
    for(size_t i=0;i<NPR && PR[i]<SP0;i++){
        u64 p=PR[i];
        u64 m = ((lo1 + p - 1)/p)*p;
        for(; m<=hi1; m+=p){
            size_t off=m-lo1;
            if(!sb->cand[off]) continue;
            int c=sb->slot[off], e=0;
            while(sb->res[off]%p==0){ sb->res[off]/=p; e++; }
            if(e && sb->nvf[c]<MAXVF){ sb->vf[c][sb->nvf[c]].p=p; sb->vf[c][sb->nvf[c]].e=e; sb->nvf[c]++; }
        }
    }
    for(int c=0;c<nc && !g_stop;c++){
        u64 v=cands[c];
        size_t off=v+1-lo1;
        u64 r=sb->res[off];
        pe_t *vf=sb->vf[c]; int nvf=sb->nvf[c];
        if(r>1){
            if(is_prime_u64(r)){ vf[nvf].p=r; vf[nvf].e=1; nvf++; }
            else {
                pe_t extra[16]; int ne=0;
                factor_rec(r, extra, &ne);
                for(int j=0;j<ne && nvf<MAXVF;j++){ vf[nvf++]=extra[j]; }
            }
        }
        /* sort ascending (small lists; child_viable wants largest last) */
        for(int i=1;i<nvf;i++){ pe_t t=vf[i]; int j=i;
            while(j>0 && vf[j-1].p>t.p){ vf[j]=vf[j-1]; j--; } vf[j]=t; }
        if(!child_viable(st, v, B, vf, nvf)) continue;
        state_t ch;
        make_child2(st, v, v, 0, &ch, vf, nvf);
        cx->emit(cx, &ch);
    }

}

/* per-thread sblock freelist: reentrant (nested scans check out distinct
 * blocks) and allocation-free in steady state */
static __thread sblock_t *t_sbpool[8];
static __thread int t_nsb;
static sblock_t *sb_get(void){
    if(t_nsb) return t_sbpool[--t_nsb];
    sblock_t *s = malloc(sizeof *s);
    if(!s){ fprintf(stderr,"sblock alloc failed\n"); abort(); }
    return s;
}
static void sb_put(sblock_t *s){ if(t_nsb<8) t_sbpool[t_nsb++]=s; else free(s); }
static void sb_drain(void){ while(t_nsb) free(t_sbpool[--t_nsb]); }

/* block-strided range scan: candidates gathered per SBW-value block, then
 * batch interval-factored in process_block. */
static void scan_range_blocks(ctx_t *cx, const state_t *st, u128 B,
                              u64 lo, u64 Q, int idx, int nthr){
    sblock_t *sb = sb_get();
    u64 cands[SBW/8];
    u64 nblocks = (Q - lo)/SBW + 1;
    for(u64 b=(u64)idx; b<nblocks && !g_stop; b+=(u64)nthr){
        u64 A = lo + b*SBW;
        u64 Aend = (Q - A >= SBW) ? A + SBW : Q + 1;    /* [A, Aend) */
        if(A > Q) break;
        int nc=0;
        size_t i0=0, i1=NPR;
        while(i0<i1){ size_t m=(i0+i1)/2; if(PR[m]<A) i0=m+1; else i1=m; }
        for(size_t i=i0; i<NPR && PR[i]<Aend; i++){
            if(nc>=SBW/8){ fprintf(stderr,"BUG cand overflow\n"); abort(); }
            cands[nc++]=PR[i];
        }
        u64 fbs = PR[NPR-1]+1; if(fbs<A) fbs=A;
        if(Aend > fbs){
            u64 v=fbs; if(!(v&1)) v++;
            for(; v<Aend && !g_stop; v+=2)
                if(is_prime_u64(v)){
                    if(nc>=SBW/8){ fprintf(stderr,"BUG cand overflow\n"); abort(); }
                    cands[nc++]=v;
                }
        }
        process_block(cx, st, B, A, Aend, cands, nc, sb);
    }
    sb_put(sb);
}

/* parallel main-scan worker.  release_budget: helper returns its budget
 * slot the moment it finishes, so a straggling teammate does not hold
 * idle capacity hostage. */
typedef struct { ctx_t *cx; const state_t *st; u128 B; u64 lo, Q;
                 int idx, nthr, release_budget; } mscan_arg_t;
static void *mscan_worker(void *a){
    mscan_arg_t *m=(mscan_arg_t*)a;
    scan_range_blocks(m->cx, m->st, m->B, m->lo, m->Q, m->idx, m->nthr);
    sb_drain();
    if(m->release_budget) atomic_fetch_add(&g_par_avail, 1);
    return NULL;
}

static void expand_node(ctx_t *cx, const state_t *st){
    if(g_stop) return;
    atomic_fetch_add_explicit(&st_nodes,1,memory_order_relaxed);
    u128 a,b; compute_ab(st,&a,&b);
    if(a==b){ record_solution(st->comps, st->nc, 0); return; }
    if(a<b) return;
    u128 B = X / st->prod;
    u128 dq = a-b;
    /* closure: single-component completion */
    if(b % dq == 0){
        u128 qs = b/dq;
        if(qs > st->f && qs <= B && is_fd_prime_128(qs)){
            atomic_fetch_add_explicit(&st_closure,1,memory_order_relaxed);
            u64 c[MAXCOMP]; memcpy(c, st->comps, (size_t)st->nc*sizeof(u64));
            if(qs <= UINT64_MAX){
                c[st->nc]=(u64)qs;
                record_solution(c, st->nc+1, 0);
            } else
                record_solution2(c, st->nc, qs, 0);
        }
    }
    /* prune: each denominator prime P needs P^(2^k) in (f,B]; product <= B */
    u128 need=1;
    for(int i=0;i<st->nr;i++){
        if(st->r[i].e >= 0) continue;
        u64 P=st->r[i].p;
        u128 w=P;
        while(w <= st->f) w=w*w;
        if(w > B) return;
        need *= w;
        if(need > B) return;
    }
    /* prune: numerator prime P needs q ≡ -1 (mod P) in (f,B] */
    for(int i=0;i<st->nr;i++){
        if(st->r[i].e <= 0) continue;
        u64 P=st->r[i].p;
        u64 m0 = (st->f+1)/P + 1;
        u128 v0 = (u128)m0*P - 1;
        if(v0 > B) return;
    }
    u128 qmin128 = st->f+1;
    u128 cb = (b + dq - 1)/dq;        /* ceil(b/dq) */
    if(cb > qmin128) qmin128 = cb;
    if(qmin128 >= ((u128)1<<42) || qmin128*qmin128*qmin128 > B){
        two_comp_resolve(st, a, b, dq, B);   /* at most 2 more components */
        return;
    }
    u64 qmin=(u64)qmin128;
    /* scan bound: mroof = floor(log_qmin B); Q = min(B, mroof*a/dq) */
    int mroof=0; { u128 t=1; while(t*qmin <= B){ t*=qmin; mroof++; } }
    u128 Q128 = (u128)mroof*a/dq;
    if(Q128 > B) Q128 = B;
    if(Q128 > UINT64_MAX){
        atomic_fetch_add_explicit(&st_warn,1,memory_order_relaxed);
        Q128 = UINT64_MAX;
    }
    u64 Q = (u64)Q128;
    u128 lo128 = b/dq + 1;
    if(st->f+1 > lo128) lo128 = st->f+1;
    if(lo128 > Q) return;
    u64 lo = (u64)lo128;
    /* primes in [lo, Q]; wide ranges go parallel in search mode */
    {
        size_t i0=0, i1=NPR;
        while(i0<i1){ size_t m=(i0+i1)/2; if(PR[m]<lo) i0=m+1; else i1=m; }
        u64 fb_lo=0;
        if(Q > PR[NPR-1]){
            atomic_fetch_add_explicit(&st_fallback,1,memory_order_relaxed);
            fb_lo = PR[NPR-1]+1; if(fb_lo<lo) fb_lo=lo;
            if(!(fb_lo&1)) fb_lo++;
        }
        size_t iq0=i0, iq1=NPR;   /* upper bound: first index with PR > Q */
        while(iq0<iq1){ size_t m=(iq0+iq1)/2; if(PR[m]<=Q) iq0=m+1; else iq1=m; }
        u64 est = (u64)(iq0-i0) + (fb_lo ? (Q-fb_lo)/2 : 0);
        int extra=0;
        /* helpers whenever budget is free: cheap near the subtree top (each
         * candidate may carry a whole child subtree), higher floor deeper */
        u64 floor_est = st->nc <= 6 ? 1000 : 20000;
        if(cx->gen_limit < 0 && est > floor_est){
            int cur = atomic_load(&g_par_avail);
            while(cur > 0){
                int take = cur < 15 ? cur : 15;
                if(atomic_compare_exchange_weak(&g_par_avail, &cur, cur - take)){ extra=take; break; }
            }
            if((extra && est > 1000000) || (getenv("IPN_TRACE") && est > 100000))
                fprintf(stderr,"WIDEMAIN nc=%d est=%" PRIu64 " helpers=%d avail=%d\n",
                        st->nc, est, extra, atomic_load(&g_par_avail));
        }
        if(extra==0){
            if(est >= 64){
                scan_range_blocks(cx, st, B, lo, Q, 0, 1);
            } else {
                for(size_t i=i0; i<NPR && PR[i]<=Q; i++){
                    if(g_stop) return;
                    scan_candidate(cx, st, PR[i], B);
                }
                if(fb_lo)
                    for(u64 v=fb_lo; v<=Q && !g_stop; v+=2)
                        if(is_prime_u64(v)) scan_candidate(cx, st, v, B);
            }
        } else {
            int n = extra + 1;
            pthread_t th[16]; mscan_arg_t ma[16];
            pthread_attr_t at; pthread_attr_init(&at);
            pthread_attr_setstacksize(&at, 16*1024*1024);
            for(int k=0;k<n;k++)
                ma[k] = (mscan_arg_t){ .cx=cx, .st=st, .B=B, .lo=lo, .Q=Q,
                                       .idx=k, .nthr=n,
                                       .release_budget = (k>0) };
            for(int k=1;k<n;k++) pthread_create(&th[k], &at, mscan_worker, &ma[k]);
            mscan_worker(&ma[0]);   /* main thread's slot is not from budget */
            for(int k=1;k<n;k++) pthread_join(th[k], NULL);
            pthread_attr_destroy(&at);
        }
    }
    /* FD prime powers p^(2^k), k>=1, in [lo, Q] */
    for(int k=1;k<=6;k++){
        u64 pmax = Q;
        for(int j=0;j<k;j++) pmax = isqrt64(pmax);
        if(pmax < 2) break;
        for(size_t i=0;i<NPR && PR[i]<=pmax;i++){
            u128 v=PR[i];
            for(int j=0;j<k;j++) v=v*v;
            if(v>=lo && v<=Q) scan_candidate2(cx, st, (u64)v, PR[i], k, B);
        }
        if(pmax > PR[NPR-1]){ fprintf(stderr,"WARN power base beyond sieve\n");
            atomic_fetch_add_explicit(&st_warn,1,memory_order_relaxed); }
    }
}

/* --------------------------- emit modes --------------------------- */
static void emit_search(ctx_t *cx, const state_t *ch);   /* work-queue, below */

typedef struct { u64 *buf; size_t len, cap; size_t *off; int *cnt; size_t ntasks, tcap; } tasklist_t;

/* optional prefix restriction (-P): gen only descends into matching subtrees */
static u64 (*g_pfx)[16]; static int *g_pfxlen; static size_t g_npfx = 0;

static void emit_gen(ctx_t *cx, const state_t *ch){
    tasklist_t *tl = (tasklist_t*)cx->gen_sink;
    if(g_npfx){
        int m=0;
        for(size_t p=0;p<g_npfx && !m;p++){
            int k = ch->nc < g_pfxlen[p] ? ch->nc : g_pfxlen[p];
            if(!memcmp(ch->comps, g_pfx[p], (size_t)k*sizeof(u64))) m=1;
        }
        if(!m) return;
    }
    if((int)ch->nc >= cx->gen_limit){
        if(tl->ntasks==tl->tcap){
            tl->tcap = tl->tcap? tl->tcap*2 : 1024;
            tl->off = realloc(tl->off, tl->tcap*sizeof(size_t));
            tl->cnt = realloc(tl->cnt, tl->tcap*sizeof(int));
        }
        if(tl->len + MAXCOMP > tl->cap){
            tl->cap = tl->cap? tl->cap*2 : 65536;
            tl->buf = realloc(tl->buf, tl->cap*sizeof(u64));
        }
        tl->off[tl->ntasks]=tl->len; tl->cnt[tl->ntasks]=ch->nc;
        memcpy(tl->buf+tl->len, ch->comps, (size_t)ch->nc*sizeof(u64));
        tl->len += (size_t)ch->nc; tl->ntasks++;
        return;
    }
    expand_node(cx, ch);
}

static void root_state(state_t *st){
    memset(st,0,sizeof *st);
    st->nc=0; st->prod=1; st->f=1;
    st->r[0].p=2; st->r[0].e=1; st->nr=1;   /* r = 2/1 */
}

/* rebuild full state from a component list (must be increasing) */
static void rebuild_state(const u64 *comps, int nc, state_t *st){
    root_state(st);
    for(int i=0;i<nc;i++){
        u64 p0; int k;
        if(!fd_decompose(comps[i],&p0,&k)){ fprintf(stderr,"BUG rebuild\n"); abort(); }
        state_t ch; make_child(st, comps[i], p0, k, &ch);
        *st = ch;
    }
}

/* --------------------------- worker pool -------------------------- */
static tasklist_t g_tasks;
static _Atomic size_t g_next;
static uint8_t *g_done;
static _Atomic size_t g_ndone;

/* work-stealing queue: subtree states parked for any idle worker.
 * outstanding[t] counts live nodes of root task t; when it drops to 0
 * the task is complete.  increment happens BEFORE push, decrement after
 * the node is fully expanded, so the count can never falsely hit 0. */
#define QCAP 8192
static struct { state_t st; size_t root; } *g_q;
static size_t g_qhead, g_qtail;            /* ring indices, guarded */
static _Atomic int g_qlen;
static pthread_mutex_t q_mu = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t  q_cv = PTHREAD_COND_INITIALIZER;
static _Atomic long long *g_outstanding;   /* per task */
static _Atomic int g_active;               /* workers currently expanding */

static void task_node_done(size_t root){
    if(atomic_fetch_sub(&g_outstanding[root],1)==1){
        if(!g_stop){ g_done[root]=1; atomic_fetch_add(&g_ndone,1); }
        pthread_mutex_lock(&q_mu);         /* wake sleepers to re-check exit */
        pthread_cond_broadcast(&q_cv);
        pthread_mutex_unlock(&q_mu);
    }
}

/* try to park a child for another worker; returns 0 if queue full */
static int q_try_push(const state_t *ch, size_t root){
    if(atomic_load_explicit(&g_qlen, memory_order_relaxed) >= QCAP-1) return 0;
    atomic_fetch_add(&g_outstanding[root],1);
    pthread_mutex_lock(&q_mu);
    if(g_qlen >= QCAP-1){
        pthread_mutex_unlock(&q_mu);
        atomic_fetch_sub(&g_outstanding[root],1);
        return 0;
    }
    g_q[g_qtail].st = *ch; g_q[g_qtail].root = root;
    g_qtail = (g_qtail+1) % QCAP;
    atomic_fetch_add(&g_qlen,1);
    pthread_cond_signal(&q_cv);
    pthread_mutex_unlock(&q_mu);
    return 1;
}

static void emit_search(ctx_t *cx, const state_t *ch){
    /* park the child if any worker might be idle; else recurse inline */
    if(atomic_load_explicit(&g_qlen, memory_order_relaxed) < 64 &&
       q_try_push(ch, cx->root))
        return;
    expand_node(cx, ch);
}

static void *worker(void *arg){
    (void)arg;
    ctx_t cx = { emit_search, -1, NULL, 0, NULL };
    for(;;){
        if(g_stop) break;
        /* prefer parked subtrees */
        pthread_mutex_lock(&q_mu);
        while(g_qlen == 0 && !g_stop){
            size_t nx = atomic_load(&g_next);
            if(nx < g_tasks.ntasks) break;              /* fresh tasks remain */
            if(atomic_load(&g_active)==0){              /* all drained: exit */
                pthread_mutex_unlock(&q_mu);
                pthread_cond_broadcast(&q_cv);
                sb_drain();
                return NULL;
            }
            struct timespec ts; clock_gettime(CLOCK_REALTIME,&ts);
            ts.tv_sec += 1;
            pthread_cond_timedwait(&q_cv,&q_mu,&ts);
        }
        if(g_stop){ pthread_mutex_unlock(&q_mu); break; }
        if(g_qlen > 0){
            state_t st = g_q[g_qhead].st;
            size_t root = g_q[g_qhead].root;
            g_qhead = (g_qhead+1) % QCAP;
            atomic_fetch_sub(&g_qlen,1);
            atomic_fetch_add(&g_active,1);
            pthread_mutex_unlock(&q_mu);
            cx.root = root;
            expand_node(&cx, &st);
            task_node_done(root);
            atomic_fetch_sub(&g_active,1);
            continue;
        }
        pthread_mutex_unlock(&q_mu);
        /* pull a fresh root task (active held across the pull so exiting
         * workers cannot miss us) */
        atomic_fetch_add(&g_active,1);
        size_t i = atomic_fetch_add(&g_next, 1);
        if(i >= g_tasks.ntasks){ atomic_fetch_sub(&g_active,1); continue; }
        if(g_done[i]){ atomic_fetch_add(&g_ndone,1); atomic_fetch_sub(&g_active,1); continue; }
        cx.root = i;
        atomic_fetch_add(&g_outstanding[i],1);
        state_t st;
        rebuild_state(g_tasks.buf+g_tasks.off[i], g_tasks.cnt[i], &st);
        expand_node(&cx, &st);
        task_node_done(i);
        atomic_fetch_sub(&g_active,1);
    }
    return NULL;
}

/* --------------------------- checkpoint --------------------------- */
static u64 tasks_hash(void){
    u64 h=1469598103934665603ULL;
    for(size_t i=0;i<g_tasks.ntasks;i++)
        for(int j=0;j<g_tasks.cnt[i];j++){
            u64 v=g_tasks.buf[g_tasks.off[i]+j];
            for(int s=0;s<64;s+=8){ h ^= (v>>s)&0xff; h *= 1099511628211ULL; }
        }
    return h;
}

static void ckpt_save(const char *path, int gen_limit){
    char tmp[1024]; snprintf(tmp,sizeof tmp,"%s.tmp",path);
    FILE *f=fopen(tmp,"w");
    if(!f){ perror("ckpt"); return; }
    char xb[48];
    fprintf(f,"IPN2 X=%s L=%d tasks=%zu hash=%016" PRIx64 "\n",
            sprint128(xb,X), gen_limit, g_tasks.ntasks, tasks_hash());
    fprintf(f,"B ");
    for(size_t i=0;i<g_tasks.ntasks;i+=4){
        int nib=0;
        for(size_t j=0;j<4 && i+j<g_tasks.ntasks;j++) if(g_done[i+j]) nib |= 1<<j;
        fputc("0123456789abcdef"[nib], f);
    }
    fputc('\n', f);
    pthread_mutex_lock(&sol_mu);
    for(int i=0;i<nsols;i++){
        fprintf(f,"S");
        for(int j=0;j<sols[i].nc;j++) fprintf(f," %" PRIu64, sols[i].comps[j]);
        if(sols[i].extra){ char eb[48]; fprintf(f," %s", sprint128(eb,sols[i].extra)); }
        fputc('\n', f);
    }
    pthread_mutex_unlock(&sol_mu);
    fclose(f);
    rename(tmp,path);
}

/* returns gen_limit from ckpt, or -1 if none/invalid */
static int ckpt_read_header(const char *path, u128 *cx, int *cl, size_t *ct, u64 *ch){
    FILE *f=fopen(path,"r");
    if(!f) return 0;
    char xs[48];
    int ok = fscanf(f,"IPN2 X=%47s L=%d tasks=%zu hash=%" SCNx64, xs,cl,ct,ch)==4;
    fclose(f);
    if(ok) *cx = parse128(xs);
    return ok;
}

static void ckpt_apply(const char *path){
    FILE *f=fopen(path,"r");
    if(!f) return;
    /* the bitmap line is ntasks/4 chars; buffer must hold it whole */
    size_t cap = g_tasks.ntasks/4 + 65536;
    char *line = malloc(cap);
    if(!line){ fclose(f); return; }
    if(!fgets(line,(int)cap,f)){ free(line); fclose(f); return; }
    while(fgets(line,(int)cap,f)){
        if(line[0]=='B'){
            char *s=line+2;
            for(size_t i=0; *s && *s!='\n'; s++, i+=4){
                int nib = *s<='9' ? *s-'0' : *s-'a'+10;
                for(size_t j=0;j<4 && i+j<g_tasks.ntasks;j++)
                    if(nib & (1<<j)) g_done[i+j]=1;
            }
        } else if(line[0]=='S'){
            u64 c[MAXCOMP]; int nc=0; u128 extra=0;
            char *s=line+1;
            for(;;){
                while(*s==' ') s++;
                if(*s<'0'||*s>'9') break;
                char *e=s; while(*e>='0'&&*e<='9') e++;
                char sv=*e; *e=0;
                u128 v=parse128(s);
                *e=sv; s=e;
                if(v > UINT64_MAX) extra=v;
                else if(nc<MAXCOMP) c[nc++]=(u64)v;
            }
            if(nc) record_solution2(c,nc,extra,1);
        }
    }
    free(line);
    fclose(f);
}

/* ------------------------------ main ------------------------------ */
static void on_sig(int s){ (void)s; g_stop=1; }

static void gen_tasks(int gen_limit){
    free(g_tasks.buf); free(g_tasks.off); free(g_tasks.cnt);
    memset(&g_tasks,0,sizeof g_tasks);
    ctx_t cx = { emit_gen, gen_limit, &g_tasks };
    state_t root; root_state(&root);
    expand_node(&cx, &root);
}

int main(int argc, char **argv){
    X = 1342989770695372800ULL;      /* smallest known IPN above a(17) */
    int nthreads=12, target=4096;
    const char *ckpt=NULL;
    int force_L=-1, shard_m=0, shard_k=0;
    const char *pfx_file=NULL, *dump_file=NULL;
    for(int i=1;i<argc;i++){
        if(!strcmp(argv[i],"-x") && i+1<argc){
            ++i;
            if(strchr(argv[i],'e')||strchr(argv[i],'E')) X=(u128)strtold(argv[i],NULL);
            else X=parse128(argv[i]);
        }
        else if(!strcmp(argv[i],"-t") && i+1<argc) nthreads=atoi(argv[++i]);
        else if(!strcmp(argv[i],"-c") && i+1<argc) ckpt=argv[++i];
        else if(!strcmp(argv[i],"-g") && i+1<argc) target=atoi(argv[++i]);
        else if(!strcmp(argv[i],"-L") && i+1<argc) force_L=atoi(argv[++i]);
        else if(!strcmp(argv[i],"-P") && i+1<argc) pfx_file=argv[++i];
        else if(!strcmp(argv[i],"-M") && i+1<argc) shard_m=atoi(argv[++i]);
        else if(!strcmp(argv[i],"-K") && i+1<argc) shard_k=atoi(argv[++i]);
        else if(!strcmp(argv[i],"-D") && i+1<argc) dump_file=argv[++i];
        else { fprintf(stderr,"usage: %s [-x bound] [-t threads] [-c ckpt] [-g tasks]\n"
                       "  [-L depth] [-P prefixfile] [-M shards -K shard] [-D dumpfile]\n",argv[0]); return 1; }
    }
    char xb[48]; sprint128(xb,X);
    char defck[256];
    if(!ckpt){ snprintf(defck,sizeof defck,"ipn_%s.ckpt",xb); ckpt=defck; }
    /* audited overflow-safe to 1e32 (256-bit two-component arithmetic;
       den = dq*q1 <= 3b+2dq <= 10X; q1*q1 checks bounded by u64max^2) */
    if(X > (u128)100000000000000000ULL*1000000000000000ULL){
        fprintf(stderr,"bound too large (max 1e32)\n"); return 1; }

    if(pfx_file){
        FILE *pf=fopen(pfx_file,"r");
        if(!pf){ perror(pfx_file); return 1; }
        size_t pcap=1024;
        g_pfx=malloc(pcap*sizeof *g_pfx); g_pfxlen=malloc(pcap*sizeof *g_pfxlen);
        char pl[1024];
        while(fgets(pl,sizeof pl,pf)){
            int n=0; char *s=pl;
            while(n<16){ char *e; u64 v=strtoull(s,&e,10); if(e==s) break; g_pfx[g_npfx][n++]=v; s=e; }
            if(n){ g_pfxlen[g_npfx]=n; g_npfx++;
                if(g_npfx==pcap){ pcap*=2; g_pfx=realloc(g_pfx,pcap*sizeof *g_pfx);
                    g_pfxlen=realloc(g_pfxlen,pcap*sizeof *g_pfxlen); } }
        }
        fclose(pf);
    }
    PR = (u64*)primesieve_generate_primes(0, PLIMIT, &NPR, UINT64_PRIMES);
    fprintf(stderr,"X=%s  threads=%d  ckpt=%s  (%zu primes sieved)\n",
            xb, nthreads, ckpt, NPR);
    signal(SIGINT,on_sig); signal(SIGTERM,on_sig);

    /* resume header (fixes gen_limit) or adaptive gen */
    u128 hx; u64 hh; int hl; size_t ht;
    int gen_limit=-1;
    int have_ck = ckpt_read_header(ckpt,&hx,&hl,&ht,&hh);
    if(have_ck){
        if(hx!=X){ fprintf(stderr,"ckpt is for a different X\n"); return 1; }
        gen_limit=hl;
        gen_tasks(gen_limit);
        if(g_tasks.ntasks!=ht || tasks_hash()!=hh){
            fprintf(stderr,"ckpt task list mismatch (rebuilt %zu, ckpt %zu)\n",g_tasks.ntasks,ht);
            return 1;
        }
        fprintf(stderr,"resuming from %s\n",ckpt);
        if(force_L>0 && force_L!=gen_limit){
            fprintf(stderr,"-L %d conflicts with ckpt depth %d\n",force_L,gen_limit); return 1; }
    } else if(force_L>0){
        gen_limit=force_L;
        gen_tasks(gen_limit);
    } else {
        for(gen_limit=2; gen_limit<=10; gen_limit++){
            gen_tasks(gen_limit);
            if((int)g_tasks.ntasks >= target) break;
        }
    }
    g_done = calloc(g_tasks.ntasks?g_tasks.ntasks:1, 1);
    /* -P: only tasks whose leading components match a listed prefix stay live */
    if(pfx_file){
        size_t kept=0;
        for(size_t i=0;i<g_tasks.ntasks;i++){
            int match=0;
            for(size_t p=0;p<g_npfx && !match;p++){
                if(g_tasks.cnt[i] >= g_pfxlen[p] &&
                   !memcmp(g_tasks.buf+g_tasks.off[i], g_pfx[p], (size_t)g_pfxlen[p]*sizeof(u64)))
                    match=1;
            }
            if(!match) g_done[i]=1; else kept++;
        }
        fprintf(stderr,"prefix filter: %zu prefixes, %zu tasks live\n",g_npfx,kept);
    }
    /* -M/-K: keep only every M-th live task (rank % M == K) */
    if(shard_m>1){
        size_t r=0, kept=0;
        for(size_t i=0;i<g_tasks.ntasks;i++){
            if(g_done[i]) continue;
            if((int)(r % (size_t)shard_m) != shard_k) g_done[i]=1; else kept++;
            r++;
        }
        fprintf(stderr,"shard %d/%d: %zu tasks live\n",shard_k,shard_m,kept);
    }
    if(have_ck) ckpt_apply(ckpt);
    /* -D: write pending task component lists and exit */
    if(dump_file){
        FILE *df=fopen(dump_file,"w");
        if(!df){ perror(dump_file); return 1; }
        size_t nd=0;
        for(size_t i=0;i<g_tasks.ntasks;i++){
            if(g_done[i]) continue;
            for(int j=0;j<g_tasks.cnt[i];j++)
                fprintf(df,"%" PRIu64 "%c", g_tasks.buf[g_tasks.off[i]+j],
                        j+1==g_tasks.cnt[i]?'\n':' ');
            nd++;
        }
        fclose(df);
        fprintf(stderr,"dumped %zu pending tasks to %s\n",nd,dump_file);
        return 0;
    }
    size_t already=0;
    for(size_t i=0;i<g_tasks.ntasks;i++) if(g_done[i]) already++;
    fprintf(stderr,"tasks: %zu at depth %d (%zu already done)\n",
            g_tasks.ntasks, gen_limit, already);

    atomic_store(&g_next,0); atomic_store(&g_ndone,0);
    atomic_store(&g_par_avail, nthreads);
    g_q = calloc(QCAP, sizeof *g_q);
    g_outstanding = calloc(g_tasks.ntasks?g_tasks.ntasks:1, sizeof *g_outstanding);
    if(!g_q || !g_outstanding){ fprintf(stderr,"alloc failed\n"); return 1; }
    pthread_t tid[256]; if(nthreads>256) nthreads=256;
    pthread_attr_t at; pthread_attr_init(&at);
    pthread_attr_setstacksize(&at, 16*1024*1024);
    for(int i=0;i<nthreads;i++) pthread_create(&tid[i],&at,worker,NULL);

    time_t t0=time(NULL), last_ck=t0;
    long long lastn=0; time_t lastt=t0;
    while(atomic_load(&g_ndone) < g_tasks.ntasks && !g_stop){
        struct timespec ts={1,0}; nanosleep(&ts,NULL);
        time_t now=time(NULL);
        if(now-lastt >= 5){
            long long n=atomic_load(&st_nodes);
            fprintf(stderr,"[%llds] tasks %zu/%zu  nodes %.3e (%.1fM/s)  sols %d  cl %lld 2c %lld dv %lld fb %lld wn %lld\n",
                (long long)(now-t0), atomic_load(&g_ndone), g_tasks.ntasks,
                (double)n, (double)(n-lastn)/1e6/(double)(now-lastt),
                nsols, (long long)st_closure,(long long)st_twocomp,
                (long long)st_divenum,(long long)st_fallback,(long long)st_warn);
            lastn=n; lastt=now;
        }
        if(now-last_ck >= 30){ ckpt_save(ckpt,gen_limit); last_ck=now; }
    }
    for(int i=0;i<nthreads;i++) pthread_join(tid[i],NULL);
    ckpt_save(ckpt,gen_limit);

    if(g_stop){
        fprintf(stderr,"interrupted at %zu/%zu tasks — checkpoint saved, rerun to resume\n",
                atomic_load(&g_ndone), g_tasks.ntasks);
        return 2;
    }
    /* final report, sorted */
    pthread_mutex_lock(&sol_mu);
    for(int i=1;i<nsols;i++){
        int j=i; while(j>0 && sols[j-1].n>sols[j].n){
            __typeof__(sols[0]) t=sols[j-1]; sols[j-1]=sols[j]; sols[j]=t; j--; }
    }
    printf("\n== COMPLETE: all infinitary perfect numbers <= %s ==\n", xb);
    for(int i=0;i<nsols;i++){
        char nb[48];
        printf("a(%d) = %s  =", i+1, sprint128(nb,sols[i].n));
        for(int j=0;j<sols[i].nc;j++) printf(" %" PRIu64, sols[i].comps[j]);
        if(sols[i].extra){ char eb[48]; printf(" %s", sprint128(eb,sols[i].extra)); }
        printf("\n");
    }
    pthread_mutex_unlock(&sol_mu);
    fprintf(stderr,"done in %llds, %lld nodes, warnings %lld\n",
            (long long)(time(NULL)-t0), (long long)st_nodes, (long long)st_warn);
    return 0;
}
