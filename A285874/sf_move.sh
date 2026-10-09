#!/bin/bash
# Prints "<move>: <perft>" from Stockfish for the no-pawn start position after one first move.
FEN="rnbqkbnr/8/8/8/8/8/8/RNBQKBNR w KQkq - 0 1"
n=$(printf "position fen %s moves %s\ngo perft %s\nquit\n" "$FEN" "$1" "$2" | /opt/homebrew/bin/stockfish | awk "/Nodes searched/ {print \$3}")
echo "$1: $n"
