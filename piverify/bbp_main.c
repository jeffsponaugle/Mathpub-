/* bbp_main.c — standalone CLI for the BBP engine.
 *
 * Usage: bbp [-d digits<=16] [-t threads] [-v] position [position ...]
 *        (positions are 1-based; prints "position hexdigits" per line)
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "bbp.h"

int main(int argc, char **argv) {
    int digits = 8;
    long threads = sysconf(_SC_NPROCESSORS_ONLN);
    int verbose = 0;
    int opt;
    while ((opt = getopt(argc, argv, "d:t:T:v")) != -1) {
        switch (opt) {
        case 'd': digits = atoi(optarg); break;
        case 't':
        case 'T': threads = atol(optarg); break;
        case 'v': verbose = 1; break;
        default:
            fprintf(stderr,
                    "usage: bbp [-d digits<=16] [-t threads] [-v] pos...\n");
            return 2;
        }
    }
    if (digits < 1 || digits > 16) {
        fprintf(stderr, "error: -d must be 1..16\n");
        return 2;
    }
    if (optind >= argc) {
        fprintf(stderr, "usage: bbp [-d digits<=16] [-t threads] [-v] pos...\n");
        return 2;
    }
    for (int i = optind; i < argc; i++) {
        char *end;
        unsigned long long p = strtoull(argv[i], &end, 10);
        if (*end || p < 1 || p > (1ULL << 59)) {
            fprintf(stderr, "error: bad position %s\n", argv[i]);
            return 2;
        }
        char buf[17];
        bbp_hex_digits((u64)p - 1, digits, (int)threads, verbose, buf);
        printf("%llu %s\n", p, buf);
        fflush(stdout);
    }
    return 0;
}
