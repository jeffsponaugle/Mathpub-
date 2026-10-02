#!/bin/bash
# wait for the missing11 queue, then validate the <353> code path on known
# size-10 shapes and re-verify Deleuran's size-11 values that are cheap for us
cd /Users/Jeff.Sponaugle/src/math/A112389
while pgrep -f "results/missing11.txt" >/dev/null || pgrep -f "lego -r" >/dev/null; do sleep 30; done
./lego -r 352 -a 02 >> results/validate10.txt 2>/dev/null
./lego -r 343 -a 02 >> results/validate10.txt 2>/dev/null
python3 verify_lasse.py results/deleuran-sum-for-size.py 11 1800 \
  32222 23231 23222 222221 132221 22322 22421 42221 4241 434 443 3242 3341 4331 4232 4223 14231 \
  123221 13331 22331 23321 33221 3332 3323 32231 > results/verify11.txt 2>&1
echo FOLLOWUP_DONE >> results/verify11.txt
