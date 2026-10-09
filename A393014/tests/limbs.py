#!/usr/bin/env python3
"""Multi-limb arithmetic test: shrink the limbs (DBPAL_LIMB_DIGITS) so that small numbers use every
limb of ./dbpal (3 limbs) and ./dbpal4 (4 limbs), and compare the hits with ./dbpal at full-size
limbs.  Covers v1, v2 and the non-coprime (gcd) mode on both engines; this is what certifies the
searches beyond ~10^54, which only ./dbpal4 can do."""
import os
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# bases, base-P length range, check base Q, algorithm
CASES = [("10,15", 17, 10, "auto"), ("6,9", 26, 6, "auto"), ("2,10", 24, 2, "auto"), ("10,12", 22, 10, "auto"),
         ("10,11", 23, 10, "v1"), ("10,11", 23, 10, "v2"), ("2,9", 25, 2, "v2"), ("2,3", 50, 2, "v2"), ("7,10", 24, 7, "v2")]


def hits(binary, b, Lmax, algo, engine, limb_digits=None):
    with tempfile.NamedTemporaryFile(suffix=".tsv") as f:
        env = dict(os.environ)
        if limb_digits:
            env["DBPAL_LIMB_DIGITS"] = str(limb_digits)
        r = subprocess.run([os.path.join(ROOT, binary), "search", "-b", b, "--Lmin", "1", "--Lmax", str(Lmax),
                            "--algo", algo, "--engine", engine, "--out", f.name], capture_output=True, text=True, env=env)
        if r.returncode != 0:
            return None, r.stderr.strip().splitlines()[-1]
        return {l.split("\t")[0] for l in open(f.name) if l[0].isdigit()}, ""


def main():
    fail = 0
    probe = subprocess.run([os.path.join(ROOT, "dbpal"), "search", "-b", "2,9", "--Lmin", "1", "--Lmax", "1",
                            "--engine", "gpu"], capture_output=True, text=True)
    engines = ["cpu", "gpu"] if probe.returncode == 0 else ["cpu"]  # GPU: Metal or CUDA, when present
    for b, Lmax, Q, algo in CASES:
        ref, err = hits("dbpal", b, Lmax, algo, "cpu")
        P = int(b.split(",")[0]) if int(b.split(",")[1]) == Q else int(b.split(",")[1])
        for binary, nl in (("dbpal", 3), ("dbpal4", 4)):
            if not os.path.exists(os.path.join(ROOT, binary)):
                continue
            c = 1  # smallest limb with Q^(c*nl) > P^Lmax, so that the top limb is in use
            while Q ** (c * nl) <= P ** Lmax:
                c += 1
            for eng in engines:
                got, err = hits(binary, b, Lmax, algo, eng, c)
                ok = got is not None and got == ref
                fail += not ok
                print(f"   {binary:6} ({b}) {algo:4} {eng}: limbs of {c} base-{Q} digits, L<={Lmax}: "
                      + (f"{len(got)} hits, {'OK' if ok else 'MISMATCH'}" if got is not None else f"MISMATCH ({err})"))
    return 1 if fail else 0


if __name__ == "__main__":
    sys.exit(main())
