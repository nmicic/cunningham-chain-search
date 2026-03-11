/*
 * Copyright (c) 2026 Nenad Mićić <nenad@micic.be>
 * SPDX-License-Identifier: Apache-2.0
 *
 * validate_chains.c — Validate Cunningham Chain (first kind) files
 * Build: gcc -O2 -o validate_chains validate_chains.c -lgmp
 * Usage: ./validate_chains
 *
 * For each all_CCN.txt file, verifies every hex entry starts a
 * first-kind chain of length >= N using GMP's mpz_probab_prime_p.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <gmp.h>

static int cc_length(const mpz_t p) {
    mpz_t q;
    mpz_init_set(q, p);
    int len = 0;
    /* mpz_probab_prime_p with 25 rounds is essentially certain */
    while (mpz_probab_prime_p(q, 25)) {
        len++;
        mpz_mul_2exp(q, q, 1);
        mpz_add_ui(q, q, 1);
    }
    mpz_clear(q);
    return len;
}

int main(void) {
    const char *dir = ".";  /* directory containing CC_8.txt .. CC_16.txt */
    int targets[] = {16, 15, 14, 13, 12, 11, 10, 9, 8};
    int ntargets = sizeof(targets) / sizeof(targets[0]);

    int total_ok = 0, total_fail = 0, total_longer = 0;
    mpz_t p;
    mpz_init(p);
    char line[4096];

    for (int ti = 0; ti < ntargets; ti++) {
        int expected = targets[ti];
        char fname[512];
        snprintf(fname, sizeof(fname), "%s/all_CC%d.txt", dir, expected);

        FILE *fp = fopen(fname, "r");
        if (!fp) {
            printf("Cannot open %s\n", fname);
            continue;
        }

        printf("=== Validating CC%d ===\n", expected);
        fflush(stdout);
        int ok = 0, fail = 0, longer = 0, count = 0;
        clock_t t0 = clock();

        while (fgets(line, sizeof(line), fp)) {
            /* strip whitespace */
            char *s = line;
            while (*s && isspace((unsigned char)*s)) s++;
            size_t len = strlen(s);
            while (len > 0 && isspace((unsigned char)s[len-1])) s[--len] = '\0';
            if (len == 0) continue;

            /* strip 0x prefix */
            if (len >= 2 && s[0] == '0' && (s[1] == 'x' || s[1] == 'X')) {
                s += 2;
                len -= 2;
            }
            if (len == 0) continue;

            /* parse hex */
            if (mpz_set_str(p, s, 16) != 0) continue;
            if (mpz_cmp_ui(p, 1) <= 0) continue;
            count++;

            int clen = cc_length(p);

            if (clen >= expected) {
                ok++;
                if (clen > expected) {
                    longer++;
                    printf("  #%d: chain=%d (>%d!) p=0x%s\n", count, clen, expected, s);
                }
            } else {
                fail++;
                printf("  FAIL #%d: chain=%d (expected >=%d) p=0x%s\n", count, clen, expected, s);
            }

            if (count % 50000 == 0) {
                double elapsed = (double)(clock() - t0) / CLOCKS_PER_SEC;
                printf("  ... %d checked (%.1fs)\n", count, elapsed);
                fflush(stdout);
            }
        }
        fclose(fp);

        double elapsed = (double)(clock() - t0) / CLOCKS_PER_SEC;
        printf("  Results: %d checked, %d OK, %d LONGER, %d FAIL (%.1fs)\n\n",
               count, ok, longer, fail, elapsed);
        fflush(stdout);

        total_ok += ok;
        total_fail += fail;
        total_longer += longer;
    }

    printf("=== TOTAL: %d OK, %d longer-than-expected, %d FAIL ===\n",
           total_ok, total_longer, total_fail);
    if (total_fail == 0)
        printf("ALL CHAINS VALID!\n");
    else
        printf("*** %d FAILURES DETECTED ***\n", total_fail);

    mpz_clear(p);
    return total_fail > 0 ? 1 : 0;
}
