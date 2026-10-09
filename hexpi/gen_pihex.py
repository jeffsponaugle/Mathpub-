#!/usr/bin/env python3
"""Generate the first N fractional hexadecimal digits of Pi (for testing
hexcount). Machin's formula with integer fixed-point arithmetic.

Usage: gen_pihex.py N outfile [--plain]

Default output starts with "3." and wraps lines at 64 chars, to exercise
hexcount's prefix-skipping and whitespace handling; --plain writes the raw
digit stream only.
"""
import sys


def arctan_inv(x, one):
    """arctan(1/x) scaled by `one`, by Taylor series."""
    total, n = 0, 0
    xpow, xsq = one // x, x * x
    while xpow:
        term = xpow // (2 * n + 1)
        total += -term if n & 1 else term
        xpow //= xsq
        n += 1
    return total


def pi_hex_digits(n):
    prec = 4 * n + 64  # guard bits
    one = 1 << prec
    pi = 16 * arctan_inv(5, one) - 4 * arctan_inv(239, one)
    frac = pi - 3 * one
    return format((frac << (4 * n)) >> prec, "0%dx" % n)


def main():
    args = [a for a in sys.argv[1:] if a != "--plain"]
    plain = "--plain" in sys.argv
    if len(args) != 2:
        sys.exit("usage: gen_pihex.py N outfile [--plain]")
    n, out = int(args[0]), args[1]
    digits = pi_hex_digits(n)
    with open(out, "w") as f:
        if plain:
            f.write(digits)
        else:
            f.write("3.\n")
            for i in range(0, n, 64):
                f.write(digits[i:i + 64] + "\n")
    print("wrote %d hex digits to %s (first 32: %s)" % (n, out, digits[:32]))


if __name__ == "__main__":
    main()
