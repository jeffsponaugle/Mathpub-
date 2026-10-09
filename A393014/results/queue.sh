#!/bin/bash
# Prioritised GPU campaign queue (sequential, resumable: completed lengths are skipped).
# Order: A393014's cross-referenced pairs and A393014 itself first, then deeper / family searches.
cd "$(dirname "$0")/.."
C="python3 tools/campaign.py run"
$C --bases 3,10 --below 10^40     # A046473 / A007633
$C --bases 8,10 --below 10^44     # A046477 / A029804
$C --bases 2,10 --below 10^40     # A046472 / A007632 (cross-check of the known b-file)
$C --bases 2,9  --below 9^49      # A393014 / A259385
$C --bases 9,10 --below 10^44     # A046478 / A029965
$C --bases 3,10 --below 10^44
$C --bases 2,3  --below 3^101     # A060792 (known to 3^93)
$C --bases 7,10 --below 10^42     # A046476 / A029964
$C --bases 10,11 --below 10^42    # A046479 / A029966
$C --bases 10,13 --below 10^42    # A046481 / A029968
$C --bases 4,10 --below 10^40     # A046474 / A029961
$C --bases 6,10 --below 10^40     # A046475 / A029963
$C --bases 10,12 --below 10^40    # A046480 / A029967
$C --bases 10,14 --below 10^40    # A046482 / A029969
$C --bases 10,15 --below 10^40    # A046483 / A029970
$C --bases 10,16 --below 10^40    # A046484 / A029731
echo QUEUE DONE
