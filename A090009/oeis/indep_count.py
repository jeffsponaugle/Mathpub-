# Independent count of the numbers <= X whose digits admit a chain of L primes
# (all L terms odd and prime to 3 and 5), written from scratch: its own walk of
# the low four digits, its own digit-sum counting.  No code shared with dspattern.
import sys
from functools import lru_cache
L, X = int(sys.argv[1]), int(sys.argv[2])
B = 10**4
def ds(n): return sum(map(int, str(n)))

def admissible(S, t, x0):
    """walk p = H*B + x0 with ds(H) = S and exactly t trailing 9s in H;
       returns (ok, carries)"""
    if x0 % 2 == 0 or (S + ds(x0)) % 3 == 0: return (False, False)
    S1, u, carried = S + 1 - 9 * t, x0, False
    for k in range(L):
        low, hs = (u, S) if u < B else (u - B, S1)
        if u >= B: carried = True
        if low % 10 in (0, 5): return (False, carried)
        if k == L - 1: break
        s = hs + ds(low)
        if s % 2: return (False, carried)
        u += s
    return (True, carried)

@lru_cache(maxsize=None)
def strings(n, s):              # n-digit strings (leading zeros) with digit sum s
    if s < 0: return 0
    if n == 0: return 1 if s == 0 else 0
    return sum(strings(n - 1, s - d) for d in range(10))

def count_upto(Y, m):           # 0 <= y <= Y with digit sum m
    if Y < 0 or m < 0: return 0
    digs, total, rem = str(Y), 0, m
    for i, ch in enumerate(digs):
        left = len(digs) - i - 1
        for d in range(int(ch)): total += strings(left, rem - d)
        rem -= int(ch)
        if rem < 0: return total
    return total + (rem == 0)

def count_exact_t(Hlo, Hhi, S, t):   # H in [Hlo,Hhi], ds(H)=S, exactly t trailing 9s
    # H = Q*10^t + (10^t - 1), Q not ending in 9, ds(Q) = S - 9t
    p, c, m = 10**t, 10**t - 1, S - 9 * t
    if m < 0 or Hhi < c: return 0
    qlo = 0 if Hlo <= c else -(-(Hlo - c) // p)
    qhi = (Hhi - c) // p
    def cnt(Y):                 # Q <= Y, ds = m, last digit != 9
        if Y < 0: return 0
        return count_upto(Y, m) - (count_upto((Y - 9) // 10, m - 9) if Y >= 9 else 0)
    return cnt(qhi) - cnt(qlo - 1)

Hmax = X // B + 1
Smax, tmax = 9 * len(str(Hmax)), len(str(Hmax))
total, npat = 0, 0
for S in range(Smax + 1):
    for x0 in range(1, B, 2):
        # non-carrying chains do not depend on t
        ok, carried = admissible(S, 0, x0)
        hlo, hhi = 1, (X - x0) // B if X >= x0 else -1
        if ok and not carried:
            npat += 1
            total += count_upto(hhi, S) - count_upto(hlo - 1, S)
        for t in range(0, min(tmax, S // 9) + 1):
            ok, carried = admissible(S, t, x0)
            if ok and carried:
                npat += 1
                total += count_exact_t(hlo, hhi, S, t)
print(f"L={L}, X={X}: {npat} patterns, {total:,} candidates")
