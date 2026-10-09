/* bbp.h — multithreaded BBP hex digit extraction for pi. */
#ifndef PIVERIFY_BBP_H
#define PIVERIFY_BBP_H

#include <stdint.h>

typedef unsigned __int128 u128;
typedef uint64_t u64;

/* frac(16^n * pi) as a 128-bit fixed-point value (see bbp.c for the
 * error analysis: reliable to well below the top 64 bits for n < 2^59). */
u128 bbp_frac(u64 n, int threads, int verbose);

/* t (1..16) hex digits of pi starting at 0-based fractional position
 * pos0, written to out as a NUL-terminated lowercase string. */
void bbp_hex_digits(u64 pos0, int t, int threads, int verbose, char out[17]);

#endif
