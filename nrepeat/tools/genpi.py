#!/usr/bin/env python3
"""Generate pi digits with gmpy2/MPFR for testing nrepeat.

usage: genpi.py NDIGITS > pi.txt        (writes "3." + NDIGITS fractional digits)
requires: pip install gmpy2
40M digits ~45s, 600M digits ~14min on an M-series Mac.
"""
import sys
import gmpy2

nd = int(sys.argv[1])
gmpy2.get_context().precision = int(nd * 3.3220) + 64
digits, exp, _ = gmpy2.const_pi().digits(10)
assert digits.startswith("3") and exp == 1
sys.stdout.write("3." + digits[1:1 + nd])
