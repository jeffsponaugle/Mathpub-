#!/bin/bash
# Builds perft and naive for this machine, stamping perft with the source hash, compiler
# and architecture (printed at the top of every run log).
cd "$(dirname "$0")" || exit 1
if [ "$(uname -m)" = arm64 ] || [ "$(uname -m)" = aarch64 ]; then ARCH=-mcpu=native; else ARCH=-march=native; fi
CXX=${CXX:-$(command -v clang++ || command -v g++)}
CC=${CC:-$(command -v clang || command -v gcc)}
SRC=$( (shasum perft.cpp 2>/dev/null || sha1sum perft.cpp) | cut -c1-12)
ID="$SRC $(basename "$CXX") $("$CXX" -dumpversion) $(uname -m)"
"$CXX" -O3 -std=c++20 $ARCH -pthread -Wall -Wextra -DBUILD_ID="\"$ID\"" -o perft perft.cpp || exit 1
"$CC" -O3 $ARCH -pthread -Wall -o naive naive.c || exit 1
echo "built: $ID"
