#!/usr/bin/env python3
"""
piverify — verify a file of pi digits using the BBP digit-extraction formula.

The BBP formula (Bailey-Borwein-Plouffe, 1996) computes hexadecimal digits
of pi at an arbitrary position without computing any earlier digits:

    pi = sum_{k>=0} 16^-k * ( 4/(8k+1) - 2/(8k+4) - 1/(8k+5) - 1/(8k+6) )

frac(16^n * pi) is computed with modular exponentiation, giving the hex
digits starting at position n+1 in O(n) small operations and O(1) memory.

Verification modes (auto-detected by `verify`):

  hex file  -> true random spot checks: sample positions, compare the file's
               digits against BBP at each position.

  dec file  -> "anchor" verification: there is no practical base-10 digit
               extraction formula, so instead we compute
               frac(16^m * x_file) from ALL N decimal digits in the file
               (one big modular exponentiation) and compare it to BBP's
               frac(16^m * pi). Base conversion mixes every decimal digit
               into the result, so with m chosen near its maximum usable
               value (~0.8305*N), a single wrong digit ANYWHERE in the file
               (except the last ~few dozen digits) scrambles the anchor
               with probability ~1 - 16^-t. This reads the whole file but
               is a genuine whole-file integrity proof, which random
               decimal spot checks cannot be.

Commands:
    verify FILE     verify a digit file (auto-detects dec/hex)
    digit POS       print hex digits of pi at 1-based position POS
    gen N           generate a reference digit file (for testing)

Requires only the standard library; uses gmpy2 automatically if installed
(large speedup for big decimal files).
"""

import argparse
import math
import os
import random
import sys
import time
from multiprocessing import Pool

try:
    sys.set_int_max_str_digits(0)  # allow huge int<->str conversions
except AttributeError:
    pass

try:
    import gmpy2
    from gmpy2 import mpz
    _powmod = gmpy2.powmod
    _isqrt = gmpy2.isqrt
    HAVE_GMPY2 = True
except ImportError:
    mpz = int
    _powmod = pow
    _isqrt = math.isqrt
    HAVE_GMPY2 = False


class PiVerifyError(Exception):
    pass


# First fractional digits of pi, used to sanity-check file format/offset.
KNOWN_DEC = (b"1415926535897932384626433832795028841971693993751"
             b"058209749445923078164062862089986280348253421170679")
KNOWN_HEX = b"243f6a8885a308d313198a2e03707344a4093822299f31d0082efa98ec4e6c89"

WS = b" \t\r\n\x0b\x0c"
STRIP = WS + b".,_"          # characters ignored when reading digit files
ALLOWED = {"dec": b"0123456789", "hex": b"0123456789abcdefABCDEF"}
KNOWN = {"dec": KNOWN_DEC, "hex": KNOWN_HEX}

# log16(10) = 1.2041199826559248, as rational bounds (upper, lower)
L16_10_HI = 12041199827
L16_10_LO = 12041199826
L16_10_DEN = 10**10


def log(msg):
    print(msg, flush=True)


# ---------------------------------------------------------------------------
# BBP hex digit extraction (integer fixed-point, optionally parallel)
# ---------------------------------------------------------------------------

def _bbp_chunk(args):
    """Partial BBP sums for k in [a, b), as B-bit fixed-point integers."""
    n, a, b, B = args
    s1 = s4 = s5 = s6 = 0
    pm = _powmod
    for k in range(a, b):
        e = n - k
        d = 8 * k + 1
        s1 += (int(pm(16, e, d)) << B) // d
        d += 3
        s4 += (int(pm(16, e, d)) << B) // d
        d += 1
        s5 += (int(pm(16, e, d)) << B) // d
        d += 1
        s6 += (int(pm(16, e, d)) << B) // d
    m = (1 << B) - 1
    return s1 & m, s4 & m, s5 & m, s6 & m


_BBP_BIN = ["unset"]


def _bbp_binary():
    """Path to the compiled C engine (bbp.c), if present."""
    if _BBP_BIN[0] == "unset":
        cand = os.environ.get("PIVERIFY_BBP")
        if not (cand and os.access(cand, os.X_OK)):
            cand = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "bbp")
            if not os.access(cand, os.X_OK):
                cand = None
        _BBP_BIN[0] = cand
    return _BBP_BIN[0]


def bbp_hex(pos0, t=8, procs=1, verbose=False):
    """Return t hex digits of pi starting at 0-based fractional position pos0.

    bbp_hex(0, 8) == '243f6a88' (pi = 3.243f6a88... in hex).
    Uses the compiled C engine (make bbp) when available and t <= 16;
    otherwise exact Python integer fixed-point with 96 guard bits.
    """
    if t <= 16 and _bbp_binary():
        import subprocess
        cmd = [_bbp_binary(), "-d", str(t), "-T", str(max(procs, 1))]
        if verbose:
            cmd.append("-v")
        r = subprocess.run(cmd + [str(pos0 + 1)],
                           stdout=subprocess.PIPE, text=True)
        if r.returncode == 0:
            out = r.stdout.split()
            if len(out) == 2 and out[0] == str(pos0 + 1):
                return out[1]
        log("      (C engine failed; falling back to Python BBP)")
    B = 4 * t + 96
    K = pos0 + 1
    start = time.time()
    if procs > 1 and K > 200_000:
        njobs = procs * 8
        step = -(-K // njobs)
        jobs = [(pos0, a, min(a + step, K), B) for a in range(0, K, step)]
        sums = [0, 0, 0, 0]
        with Pool(procs) as pool:
            done = 0
            for part in pool.imap_unordered(_bbp_chunk, jobs):
                for i in range(4):
                    sums[i] += part[i]
                done += 1
                if verbose and done % max(1, len(jobs) // 8) == 0:
                    log("      ... BBP %d/%d chunks (%.1fs)"
                        % (done, len(jobs), time.time() - start))
        s1, s4, s5, s6 = sums
    else:
        s1, s4, s5, s6 = _bbp_chunk((pos0, 0, K, B))

    # tail terms k > pos0 (16^(pos0-k) is a right shift)
    k = pos0 + 1
    while True:
        shift = 4 * (k - pos0)
        if shift > B:
            break
        s1 += (1 << B) // ((8 * k + 1) << shift)
        s4 += (1 << B) // ((8 * k + 4) << shift)
        s5 += (1 << B) // ((8 * k + 5) << shift)
        s6 += (1 << B) // ((8 * k + 6) << shift)
        k += 1

    total = (4 * s1 - 2 * s4 - s5 - s6) % (1 << B)
    return format(total >> (B - 4 * t), "0{}x".format(t))


# ---------------------------------------------------------------------------
# Chudnovsky (binary splitting) — reference digits for `gen`
# ---------------------------------------------------------------------------

_C3_24 = 640320**3 // 24


def _bs(a, b):
    if b - a == 1:
        if a == 0:
            P = Q = mpz(1)
        else:
            P = mpz(6 * a - 5) * (2 * a - 1) * (6 * a - 1)
            Q = mpz(a) ** 3 * _C3_24
        T = P * (13591409 + 545140134 * a)
        if a & 1:
            T = -T
        return P, Q, T
    m = (a + b) // 2
    p1, q1, t1 = _bs(a, m)
    p2, q2, t2 = _bs(m, b)
    return p1 * p2, q1 * q2, t1 * q2 + p1 * t2


def pi_scaled(F):
    """floor(pi * F) for a large integer scale F."""
    F = mpz(F)
    prec_digits = F.bit_length() * 0.30103 + 10
    terms = int(prec_digits / 14.181647462725477) + 2
    _, Q, T = _bs(0, max(terms, 1))
    G = 64
    sq = _isqrt(10005 * (F << G) ** 2)
    result = (426880 * sq * Q) // (T << G)
    if result <= 0:
        raise PiVerifyError("internal error: Chudnovsky sign")
    return result


# ---------------------------------------------------------------------------
# Digit file access
# ---------------------------------------------------------------------------

def sniff_base(path):
    with open(path, "rb") as f:
        head = f.read(8192)
    s = head.translate(None, STRIP).lower()
    if s.startswith(b"3" + KNOWN_HEX[:8]) or s.startswith(KNOWN_HEX[:8]):
        return "hex"
    if s.startswith(b"3" + KNOWN_DEC[:7]) or s.startswith(KNOWN_DEC[:7]):
        return "dec"
    raise PiVerifyError(
        "file does not begin with the digits of pi in either base 10 or 16\n"
        "  (expected 3.1415926535... or 3.243f6a8885...; got %r...)"
        % head[:32])


def _bad_char_error(chunk, base):
    bad = chunk.translate(None, ALLOWED[base] + STRIP)
    return PiVerifyError(
        "unexpected character %r in digit file (base %s)" % (bad[:1], base))


class DigitFile:
    """Random access to the fractional digits of a pi digit file.

    Handles an optional '3.' / '3' prefix and ignores whitespace.  Files
    that are one contiguous run of digits are read with O(1) seeks; files
    with line breaks or other formatting fall back to streaming scans.
    """

    def __init__(self, path, base):
        self.path = path
        self.base = base
        self.allowed = ALLOWED[base]
        size = os.path.getsize(path)
        with open(path, "rb") as f:
            head = f.read(8192)
        if not head:
            raise PiVerifyError("file is empty")

        known = KNOWN[base]
        stripped = head.translate(None, STRIP).lower()
        if stripped.startswith(b"3" + known[:8]):
            self.skip = 1  # stream of digits includes the leading '3'
        elif stripped.startswith(known[:8]):
            self.skip = 0
        else:
            raise PiVerifyError(
                "file does not start with the digits of pi (base %s): %r..."
                % (base, head[:24]))

        # byte-exact prefix for contiguous ("pure") files
        i = 0
        while i < len(head) and head[i] in WS:
            i += 1
        if head[i:i + 2] == b"3.":
            data0 = i + 2
        elif self.skip and head[i:i + 1] == b"3":
            data0 = i + 1
        else:
            data0 = i
        body = head[data0:]
        self.pure = len(body.translate(None, self.allowed)) == 0

        if self.pure:
            self.data0 = data0
            with open(path, "rb") as f:
                f.seek(max(0, size - 4096))
                tail = f.read()
            trail = len(tail) - len(tail.rstrip(WS))
            core = tail.rstrip(WS)
            if core and core[-16:].translate(None, self.allowed):
                raise PiVerifyError("file tail contains non-digit data; "
                                    "cannot use fast seek mode")
            self.n = size - data0 - trail
        else:
            self.data0 = None
            self.n = self._count()

    def _count(self):
        log("  (formatted file: counting digits with a full pass...)")
        n = 0
        with open(self.path, "rb") as f:
            while True:
                chunk = f.read(1 << 22)
                if not chunk:
                    break
                d = chunk.translate(None, STRIP)
                if d.translate(None, self.allowed):
                    raise _bad_char_error(chunk, self.base)
                n += len(d)
        return n - self.skip

    def fetch(self, ranges):
        """ranges: iterable of (start0, length) in fractional-digit indices.
        Returns {(start0, length): lowercase digit string}."""
        ranges = sorted(set(ranges))
        out = {}
        if self.pure:
            with open(self.path, "rb") as f:
                for s, ln in ranges:
                    f.seek(self.data0 + s)
                    b = f.read(ln)
                    if len(b) != ln or b.translate(None, self.allowed):
                        raise PiVerifyError(
                            "could not read %d digits at position %d" % (ln, s))
                    out[(s, ln)] = b.lower().decode("ascii")
            return out

        # streaming scan for formatted files
        results = {r: bytearray() for r in ranges}
        idx = 0
        ri = 0
        skip_left = self.skip
        with open(self.path, "rb") as f:
            while ri < len(ranges):
                chunk = f.read(1 << 22)
                if not chunk:
                    raise PiVerifyError("EOF before requested digit position")
                d = chunk.translate(None, STRIP)
                if d.translate(None, self.allowed):
                    raise _bad_char_error(chunk, self.base)
                if skip_left:
                    drop = min(skip_left, len(d))
                    d = d[drop:]
                    skip_left -= drop
                clen = len(d)
                for j in range(ri, len(ranges)):
                    s, ln = ranges[j]
                    if s >= idx + clen:
                        break
                    e = s + ln
                    if e <= idx:
                        continue
                    lo = max(s, idx) - idx
                    hi = min(e, idx + clen) - idx
                    results[(s, ln)] += d[lo:hi]
                while (ri < len(ranges)
                       and len(results[ranges[ri]]) == ranges[ri][1]):
                    ri += 1
                idx += clen
        return {r: bytes(v).lower().decode("ascii")
                for r, v in results.items()}


def load_dec_digits(path):
    """Load ALL fractional decimal digits as one bytes object."""
    parts = []
    with open(path, "rb") as f:
        while True:
            chunk = f.read(1 << 24)
            if not chunk:
                break
            d = chunk.translate(None, STRIP)
            if d.translate(None, ALLOWED["dec"]):
                raise _bad_char_error(chunk, "dec")
            parts.append(d)
    s = b"".join(parts)
    if s.startswith(b"3" + KNOWN_DEC[:8]):
        s = s[1:]
    if not s.startswith(KNOWN_DEC[:8]):
        raise PiVerifyError("file does not start with 3.14159265...")
    return s


# ---------------------------------------------------------------------------
# verify: decimal files (whole-file anchor verification)
# ---------------------------------------------------------------------------

def verify_dec(path, args, rng):
    t = args.digits if args.digits else 12
    log("mode: DECIMAL file -> whole-file anchor verification against BBP")
    log("  (base 10 has no practical digit-extraction formula; instead the")
    log("   file's full value is base-converted at random hex anchor points")
    log("   and compared with BBP, which catches an error at ANY position)")
    log("")

    t0 = time.time()
    log("[1/4] reading file...")
    digits = load_dec_digits(path)
    n = len(digits)
    log("      %s fractional decimal digits (%.1fs)" % (f"{n:,}", time.time() - t0))
    if n < 2000:
        raise PiVerifyError("need at least 2000 digits for anchor verification")

    log("[2/4] checking first %d digits against known value..." % min(n, 100))
    if digits[:100] != KNOWN_DEC[:min(n, 100)]:
        log("      FAIL: leading digits are wrong")
        return False
    log("      ok")

    log("[3/4] converting file to a %s-digit integer..." % f"{n:,}")
    t0 = time.time()
    if HAVE_GMPY2:
        d_int = mpz(digits.decode("ascii"))
    else:
        d_int = int(digits)
    mod = mpz(10) ** n
    log("      done (%.1fs)%s" % (time.time() - t0,
        "" if HAVE_GMPY2 else "  [tip: pip install gmpy2 for large files]"))

    # Highest anchor m with 16^(m+t+20) <= 10^n, so file truncation error
    # (< 10^-n) stays far below the t compared hex digits.
    m_cap = n * L16_10_DEN // L16_10_HI - t - 20
    anchors = [m_cap]
    while len(anchors) < max(1, args.anchors):
        m = m_cap - rng.randrange(1, 512)
        if m not in anchors:
            anchors.append(m)

    log("[4/4] checking %d anchor(s) near hex position %s..."
        % (len(anchors), f"{m_cap:,}"))
    ok = True
    for m in anchors:
        t0 = time.time()
        y = (_powmod(2, 4 * m, mod) * d_int) % mod
        file_hex = format(int((y << (4 * t)) // mod), "0{}x".format(t))
        t1 = time.time()
        bbp = bbp_hex(m, t, procs=args.processes, verbose=args.verbose)
        status = "PASS" if file_hex == bbp else "FAIL"
        ok = ok and (file_hex == bbp)
        log("      anchor hex pos %s: file->%s  bbp->%s  %s"
            "  (convert %.1fs, bbp %.1fs)"
            % (f"{m:,}", file_hex, bbp, status, t1 - t0, time.time() - t1))

    m_min = min(anchors)
    covered = (m_min + t) * L16_10_LO // L16_10_DEN
    covered = min(covered, n)
    log("")
    if ok:
        log("RESULT: PASS")
        log("  A single wrong digit at any decimal position 1..%s would have"
            % f"{covered:,}")
        log("  scrambled every anchor (miss probability ~16^-%d per anchor)."
            % t)
        if covered < n:
            log("  Note: the final %s digits are beyond the anchors' resolution"
                % f"{n - covered:,}")
            log("  and are not covered by this check.")
    else:
        log("RESULT: FAIL — the file's digits do not match pi.")
        log("  (An anchor mismatch means at least one digit somewhere in the")
        log("   covered range is wrong, or digits are missing/inserted.)")
    return ok


# ---------------------------------------------------------------------------
# verify: hex files (true random spot checks)
# ---------------------------------------------------------------------------

def verify_hex(path, args, rng):
    t = args.digits if args.digits else 8
    log("mode: HEX file -> random spot checks with the BBP formula")
    log("")

    log("[1/3] scanning file...")
    df = DigitFile(path, "hex")
    n = df.n
    log("      %s fractional hex digits (%s access)"
        % (f"{n:,}", "seek" if df.pure else "streaming"))
    if n < 2 * t:
        raise PiVerifyError("file too small")

    log("[2/3] checking first %d digits against known value..." % min(n, 64))
    lead = df.fetch([(0, min(n, 64))])[(0, min(n, 64))]
    if lead != KNOWN_HEX[:min(n, 64)].decode():
        log("      FAIL: leading digits are wrong")
        return False
    log("      ok")

    positions = {n - t}                      # always check the very end:
    for _ in range(args.checks):             # computation errors propagate,
        positions.add(rng.randrange(0, n - t + 1))   # so the tail is decisive
    ranges = sorted((p, t) for p in positions)

    log("[3/3] running %d spot checks (positions chosen at random%s)..."
        % (len(ranges), "" if args.seed is None else ", seed=%d" % args.seed))
    file_digits = df.fetch(ranges)
    ok = True
    for pos, ln in ranges:
        t0 = time.time()
        expect = bbp_hex(pos, ln, procs=args.processes, verbose=args.verbose)
        got = file_digits[(pos, ln)]
        status = "PASS" if got == expect else "FAIL"
        ok = ok and (got == expect)
        log("      hex digits %s..%s: file->%s  bbp->%s  %s  (%.1fs)"
            % (f"{pos + 1:,}", f"{pos + ln:,}", got, expect, status,
               time.time() - t0))

    log("")
    if ok:
        log("RESULT: PASS — all sampled positions match BBP.")
        log("  Note: spot checks verify the sampled digits. The end-of-file")
        log("  check is the strong one: any error in a pi *computation*")
        log("  propagates to all later digits, so a matching tail plus")
        log("  random interior samples is the standard verification method.")
    else:
        log("RESULT: FAIL — at least one sampled position does not match pi.")
    return ok


# ---------------------------------------------------------------------------
# commands
# ---------------------------------------------------------------------------

def cmd_verify(args):
    base = args.base
    if base == "auto":
        base = sniff_base(args.file)
    rng = random.Random(args.seed)
    log("piverify: %s" % args.file)
    log("gmpy2: %s | processes: %d | BBP engine: %s" %
        ("yes" if HAVE_GMPY2 else "no (stdlib big-int fallback)",
         args.processes,
         "C (%s)" % _bbp_binary() if _bbp_binary() else "Python"))
    log("")
    if base == "dec":
        ok = verify_dec(args.file, args, rng)
    else:
        ok = verify_hex(args.file, args, rng)
    return 0 if ok else 1


def cmd_digit(args):
    pos1 = args.position
    if pos1 < 1:
        raise PiVerifyError("position is 1-based; must be >= 1")
    t0 = time.time()
    d = bbp_hex(pos1 - 1, args.count, procs=args.processes,
                verbose=args.verbose)
    log("hex digits %s..%s of pi (after the point): %s   [%.1fs]"
        % (f"{pos1:,}", f"{pos1 + args.count - 1:,}", d, time.time() - t0))
    return 0


def cmd_gen(args):
    n = args.digits_count
    if n < 100:
        raise PiVerifyError("gen needs at least 100 digits")
    if n > 5_000_000:
        log("warning: generating %s digits in pure Python may take a while"
            % f"{n:,}")
    t0 = time.time()
    if args.base == "dec":
        s = str(pi_scaled(mpz(10) ** (n + 8)) // mpz(10) ** 8)
        digits = s[1:]
        check = KNOWN_DEC[:min(n, 100)].decode()
    else:
        v = pi_scaled(mpz(2) ** (4 * n + 32)) >> 32
        s = format(int(v), "x")
        digits = s[1:]
        check = KNOWN_HEX[:min(n, 64)].decode()
    if len(digits) != n or not digits.startswith(check):
        raise PiVerifyError("internal self-check failed generating digits")
    with open(args.out, "w") as f:
        f.write("3.")
        w = args.line_width
        if w <= 0:
            f.write(digits)
        else:
            f.write("\n")
            for i in range(0, n, w):
                f.write(digits[i:i + w])
                f.write("\n")
    log("wrote %s %s digits of pi to %s (%.1fs)"
        % (f"{n:,}", args.base, args.out, time.time() - t0))
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="piverify",
        description="Verify files of pi digits with the BBP formula.")
    ap.add_argument("--processes", type=int,
                    default=max(1, os.cpu_count() or 1),
                    help="worker processes for BBP (default: all cores)")
    ap.add_argument("--verbose", action="store_true",
                    help="show BBP progress for long computations")
    sub = ap.add_subparsers(dest="cmd", required=True)

    v = sub.add_parser("verify", help="verify a pi digit file")
    v.add_argument("file")
    v.add_argument("--base", choices=["auto", "dec", "hex"], default="auto")
    v.add_argument("--checks", type=int, default=8,
                   help="random spot checks for hex files (default 8)")
    v.add_argument("--anchors", type=int, default=3,
                   help="anchor points for decimal files (default 3)")
    v.add_argument("--digits", type=int, default=0,
                   help="hex digits compared per check (default 8 hex / 12 dec)")
    v.add_argument("--seed", type=int, default=None,
                   help="seed the random position generator (reproducible)")
    v.set_defaults(fn=cmd_verify)

    d = sub.add_parser("digit", help="hex digits of pi at a position")
    d.add_argument("position", type=int, help="1-based hex digit position")
    d.add_argument("--count", type=int, default=8)
    d.set_defaults(fn=cmd_digit)

    g = sub.add_parser("gen", help="generate a reference digit file (testing)")
    g.add_argument("digits_count", type=int)
    g.add_argument("--base", choices=["dec", "hex"], default="dec")
    g.add_argument("-o", "--out", required=True)
    g.add_argument("--line-width", type=int, default=0,
                   help="wrap digits at this width (default: one line)")
    g.set_defaults(fn=cmd_gen)

    args = ap.parse_args(argv)
    try:
        return args.fn(args)
    except PiVerifyError as e:
        log("error: %s" % e)
        return 2


if __name__ == "__main__":
    sys.exit(main())
