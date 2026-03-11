/*
 * Copyright (c) 2026 Nenad Mićić <nenad@micic.be>
 * SPDX-License-Identifier: Apache-2.0
 *
 * cc_immunization_heatmap.c — Immunization heatmap per CC-length × bitsize
 *
 * Reads combined CC file (format: "CC10 0xHEX 25") and for each root p
 * checks immunization against small primes q < 37:
 *   Immune to q  ⟺  (p + 1) ≡ 0 (mod q)
 *   Kill position ⟺  smallest j where 2^j·(p+1) ≡ 1 (mod q)
 *
 * Output: CSV with immunization rates per (cc_length, bits, prime).
 *
 * No GMP needed — computes p mod q from hex digits directly.
 *
 * Build:
 *   gcc -O2 -o cc_immunization_heatmap cc_immunization_heatmap.c -lm
 *
 * Usage:
 *   ./cc_immunization_heatmap /path/to/CC_x.txt [--kill-csv kill_positions.csv]
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

/* Small primes < 37 */
static const int PRIMES[] = {3, 5, 7, 11, 13, 17, 19, 23, 29, 31};
#define NPRIMES 10

/* Max CC length and bit size we track */
#define MAX_CC  25
#define MAX_BITS 300

/* Bucket: per (cc, bits) pair */
typedef struct {
    int count;                     /* total roots in this bucket */
    int immune[NPRIMES];           /* count immune to each prime: (p+1) ≡ 0 mod q */
    int safe[NPRIMES];             /* count safe: immune OR kill position >= cc length */
    int kill_hist[NPRIMES][MAX_CC]; /* kill_hist[p][j] = how many roots killed at position j */
} Bucket;

static Bucket buckets[MAX_CC + 1][MAX_BITS + 1];

/* Parse a hex digit */
static inline int hexval(char c) {
    if (c >= '0' && c <= '9') return c - '0';
    if (c >= 'a' && c <= 'f') return c - 'a' + 10;
    if (c >= 'A' && c <= 'F') return c - 'A' + 10;
    return -1;
}

/* Compute hex_string mod q without bignum — process digit by digit */
static int hex_mod(const char *hex, int len, int q) {
    int r = 0;
    for (int i = 0; i < len; i++) {
        int d = hexval(hex[i]);
        if (d < 0) return -1;
        r = (r * 16 + d) % q;
    }
    return r;
}

/* Find kill position: smallest j >= 0 where 2^j * (p+1) ≡ 1 (mod q)
 * Returns j, or -1 if immune (p+1 ≡ 0 mod q), or -2 if not killed within MAX_CC */
static int kill_position(int p_mod_q, int q) {
    int pp1 = (p_mod_q + 1) % q;
    if (pp1 == 0) return -1;  /* immune */

    /* Find j where 2^j * pp1 ≡ 1 (mod q) */
    int val = pp1;
    for (int j = 0; j < MAX_CC; j++) {
        if (val == 1) return j;       /* 2^j * (p+1) ≡ 1 mod q → q | chain member j */
        val = (val * 2) % q;
    }
    return -2;  /* not killed within range */
}

int main(int argc, char **argv) {
    const char *input = NULL;
    const char *kill_csv = NULL;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--kill-csv") == 0 && i + 1 < argc) {
            kill_csv = argv[++i];
        } else if (!input) {
            input = argv[i];
        }
    }

    if (!input) {
        fprintf(stderr, "Usage: %s <CC_x.txt> [--kill-csv kill_positions.csv]\n", argv[0]);
        return 1;
    }

    FILE *fp = fopen(input, "r");
    if (!fp) {
        fprintf(stderr, "Cannot open %s\n", input);
        return 1;
    }

    memset(buckets, 0, sizeof(buckets));

    char line[4096];
    int total = 0, skipped = 0;

    while (fgets(line, sizeof(line), fp)) {
        /* Trim */
        char *s = line;
        while (*s && isspace((unsigned char)*s)) s++;
        size_t len = strlen(s);
        while (len > 0 && isspace((unsigned char)s[len - 1])) s[--len] = '\0';
        if (len == 0) continue;

        /* Parse: CC<len> 0x<hex> <bits> */
        char *p1 = s;
        if (p1[0] != 'C' || p1[1] != 'C') { skipped++; continue; }
        p1 += 2;
        char *end;
        int cc = (int)strtol(p1, &end, 10);
        if (end == p1 || cc < 1 || cc > MAX_CC) { skipped++; continue; }

        /* Skip whitespace to hex */
        char *p2 = end;
        while (*p2 && isspace((unsigned char)*p2)) p2++;
        if (*p2 == '0' && (p2[1] == 'x' || p2[1] == 'X')) p2 += 2;

        /* Find end of hex */
        char *hex_start = p2;
        int hex_len = 0;
        while (hexval(p2[hex_len]) >= 0) hex_len++;
        if (hex_len == 0) { skipped++; continue; }

        /* Compute actual bit length from hex value.
         * Third column is hex digit count, NOT bit size. */
        int lead = hexval(hex_start[0]);
        int lead_bits = 0;
        if (lead > 0) { int v = lead; while (v) { lead_bits++; v >>= 1; } }
        int bits = (hex_len - 1) * 4 + lead_bits;
        if (bits > MAX_BITS) bits = MAX_BITS;

        /* Process immunization for each small prime */
        Bucket *b = &buckets[cc][bits];
        b->count++;

        for (int pi = 0; pi < NPRIMES; pi++) {
            int q = PRIMES[pi];
            int p_mod = hex_mod(hex_start, hex_len, q);
            if (p_mod < 0) continue;

            int kp = kill_position(p_mod, q);
            if (kp == -1) {
                b->immune[pi]++;
                b->safe[pi]++;
            } else if (kp == -2) {
                /* Never killed (residue not in <2> mod q) */
                b->safe[pi]++;
            } else if (kp >= 0) {
                if (kp >= cc) {
                    /* Kill position beyond chain length — safe for this chain */
                    b->safe[pi]++;
                }
                if (kp < MAX_CC) {
                    b->kill_hist[pi][kp]++;
                }
            }
        }

        total++;
        if (total % 200000 == 0) {
            fprintf(stderr, "  [%d roots processed]\n", total);
        }
    }
    fclose(fp);

    fprintf(stderr, "Processed %d roots (%d skipped)\n", total, skipped);

    /* ---- Output 1: Immunization rate heatmap CSV ---- */
    /* Header: immune = (p+1)≡0 mod q; safe = immune OR kill_pos >= cc OR not in <2> */
    printf("cc,bits,count");
    for (int pi = 0; pi < NPRIMES; pi++)
        printf(",imm_%d,imm_%d_pct,safe_%d,safe_%d_pct",
               PRIMES[pi], PRIMES[pi], PRIMES[pi], PRIMES[pi]);
    printf(",safe_all,safe_all_pct\n");

    for (int cc = 1; cc <= MAX_CC; cc++) {
        for (int bits = 1; bits <= MAX_BITS; bits++) {
            Bucket *b = &buckets[cc][bits];
            if (b->count == 0) continue;

            printf("%d,%d,%d", cc, bits, b->count);

            int min_safe = b->count;
            for (int pi = 0; pi < NPRIMES; pi++) {
                double imm_pct = 100.0 * b->immune[pi] / b->count;
                double safe_pct = 100.0 * b->safe[pi] / b->count;
                printf(",%d,%.2f,%d,%.2f",
                       b->immune[pi], imm_pct, b->safe[pi], safe_pct);
                if (b->safe[pi] < min_safe) min_safe = b->safe[pi];
            }
            /* Lower bound on all-safe (min of individual safe counts) */
            printf(",%d,%.2f", min_safe, 100.0 * min_safe / b->count);
            printf("\n");
        }
    }

    /* ---- Output 2: Kill position CSV (optional) ---- */
    if (kill_csv) {
        FILE *kf = fopen(kill_csv, "w");
        if (!kf) {
            fprintf(stderr, "Cannot write %s\n", kill_csv);
            return 1;
        }

        fprintf(kf, "cc,bits,prime,position,count,pct_of_bucket\n");
        for (int cc = 1; cc <= MAX_CC; cc++) {
            for (int bits = 1; bits <= MAX_BITS; bits++) {
                Bucket *b = &buckets[cc][bits];
                if (b->count == 0) continue;
                for (int pi = 0; pi < NPRIMES; pi++) {
                    for (int j = 0; j < MAX_CC; j++) {
                        if (b->kill_hist[pi][j] > 0) {
                            fprintf(kf, "%d,%d,%d,%d,%d,%.2f\n",
                                    cc, bits, PRIMES[pi], j,
                                    b->kill_hist[pi][j],
                                    100.0 * b->kill_hist[pi][j] / b->count);
                        }
                    }
                }
            }
        }
        fclose(kf);
        fprintf(stderr, "Kill positions written to %s\n", kill_csv);
    }

    /* ---- Summary to stderr ---- */
    fprintf(stderr, "\n=== SUMMARY ===\n");
    fprintf(stderr, "  immune = (p+1) ≡ 0 mod q\n");
    fprintf(stderr, "  safe   = immune OR residue not in <2> mod q OR kill beyond chain\n");
    fprintf(stderr, "  random = expected 1/q for random primes\n\n");
    fprintf(stderr, "  %5s  %8s %8s %8s  %6s %6s\n",
            "q", "immune%", "safe%", "random%", "ord2q", "|safe|/q");
    fprintf(stderr, "  %5s  %8s %8s %8s  %6s %6s\n",
            "-----", "--------", "--------", "--------", "------", "------");
    for (int pi = 0; pi < NPRIMES; pi++) {
        int q = PRIMES[pi];
        int tot_imm = 0, tot_safe = 0;
        for (int cc = 1; cc <= MAX_CC; cc++)
            for (int bits = 1; bits <= MAX_BITS; bits++) {
                tot_imm += buckets[cc][bits].immune[pi];
                tot_safe += buckets[cc][bits].safe[pi];
            }
        /* Compute ord(2, q) */
        int ord = 1;
        int v = 2;
        while (v != 1) { v = (v * 2) % q; ord++; }
        /* Safe residues: q - |<2>| = q - ord (excluding 0 which is immune) */
        int safe_residues = q - ord;  /* includes 0 */
        fprintf(stderr, "  q=%2d: %7.2f%% %7.2f%% %7.2f%%  %5d  %d/%d\n",
                q,
                100.0 * tot_imm / total,
                100.0 * tot_safe / total,
                100.0 / q,
                ord, safe_residues, q);
    }

    return 0;
}
