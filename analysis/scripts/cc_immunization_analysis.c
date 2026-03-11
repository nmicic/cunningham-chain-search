/*
 * Copyright (c) 2026 Nenad Mićić <nenad@micic.be>
 * SPDX-License-Identifier: Apache-2.0
 *
 * cc_immunization_analysis.c — Full immunization analysis → JSON
 *
 * Reads combined CC file (format: "CC10 0xHEX 25") and produces a single
 * JSON file with all immunization data for the dashboard.
 *
 * Analyses:
 *   1. Per (cc, bits) bucket: immune & safe counts per prime
 *   2. Kill position histogram
 *   3. Immune combination fingerprints (which primes each root is immune to)
 *   4. Residue class distribution per prime
 *   5. Group theory constants (ord(2,q), subgroup, safe set)
 *
 * No GMP needed — computes p mod q from hex digits directly.
 *
 * Build:
 *   gcc -O2 -o cc_immunization_analysis cc_immunization_analysis.c -lm
 *
 * Usage:
 *   ./cc_immunization_analysis CC_x.txt -o immunization.json
 *   ./cc_immunization_analysis CC_x.txt -o immunization.json --min-bits 80
 *   ./cc_immunization_analysis CC_x.txt -o immunization.json --min-bits 80 --min-cc 10
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>

/* Small primes < 37 */
static const int PRIMES[] = {3, 5, 7, 11, 13, 17, 19, 23, 29, 31};
#define NPRIMES 10

#define MAX_CC   25
#define MAX_BITS 512

/* 2^NPRIMES = 1024 possible immune combinations */
#define NCOMBOS (1 << NPRIMES)

/* ---- Per (cc, bits) bucket ---- */
typedef struct {
    int count;
    int immune[NPRIMES];
    int safe[NPRIMES];
    int kill_hist[NPRIMES][MAX_CC];
    int residue[NPRIMES][32];  /* residue[(p+1) mod q] counts, max q=31 */
} Bucket;

/* ---- Per CC length: combo tracking ---- */
typedef struct {
    int count;
    int combo_count[NCOMBOS];  /* combo_count[bitmask] = how many roots */
} CCCombo;

static Bucket  buckets[MAX_CC + 1][MAX_BITS + 1];
static CCCombo combos[MAX_CC + 1];

/* ---- Helpers ---- */
static inline int hexval(char c) {
    if (c >= '0' && c <= '9') return c - '0';
    if (c >= 'a' && c <= 'f') return c - 'a' + 10;
    if (c >= 'A' && c <= 'F') return c - 'A' + 10;
    return -1;
}

static int hex_mod(const char *hex, int len, int q) {
    int r = 0;
    for (int i = 0; i < len; i++) {
        int d = hexval(hex[i]);
        if (d < 0) return -1;
        r = (r * 16 + d) % q;
    }
    return r;
}

static int kill_position(int p_mod_q, int q) {
    int pp1 = (p_mod_q + 1) % q;
    if (pp1 == 0) return -1;   /* immune */
    int val = pp1;
    for (int j = 0; j < MAX_CC; j++) {
        if (val == 1) return j;
        val = (val * 2) % q;
    }
    return -2;  /* never killed (not in <2> mod q) */
}

static int ord2(int q) {
    int v = 2, o = 1;
    while (v != 1) { v = (v * 2) % q; o++; }
    return o;
}

/* ---- JSON output helpers ---- */
static FILE *jf;
static int j_first;  /* for comma tracking */

static void j_obj_start(void)  { fprintf(jf, "{"); j_first = 1; }
static void j_obj_end(void)    { fprintf(jf, "}"); }
static void j_arr_start(void)  { fprintf(jf, "["); j_first = 1; }
static void j_arr_end(void)    { fprintf(jf, "]"); }
static void j_comma(void)      { if (!j_first) fprintf(jf, ","); j_first = 0; }

static void j_key(const char *k)           { j_comma(); fprintf(jf, "\"%s\":", k); }
static void j_key_int(const char *k, int v) { j_comma(); fprintf(jf, "\"%s\":%d", k, v); }
static void j_key_dbl(const char *k, double v) { j_comma(); fprintf(jf, "\"%s\":%.4f", k, v); }
static void j_key_str(const char *k, const char *v) { j_comma(); fprintf(jf, "\"%s\":\"%s\"", k, v); }
static void j_int(int v)     { j_comma(); fprintf(jf, "%d", v); }
static void j_dbl(double v)  { j_comma(); fprintf(jf, "%.4f", v); }

/* ---- Main ---- */
int main(int argc, char **argv) {
    const char *input = NULL;
    const char *output = "immunization.json";
    int min_bits = 0;
    int min_cc = 0;

    for (int i = 1; i < argc; i++) {
        if ((strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "--output") == 0) && i + 1 < argc) {
            output = argv[++i];
        } else if (strcmp(argv[i], "--min-bits") == 0 && i + 1 < argc) {
            min_bits = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--min-cc") == 0 && i + 1 < argc) {
            min_cc = atoi(argv[++i]);
        } else if (argv[i][0] != '-' && !input) {
            input = argv[i];
        }
    }

    if (!input) {
        fprintf(stderr, "Usage: %s <CC_x.txt> [-o output.json] [--min-bits N] [--min-cc N]\n", argv[0]);
        return 1;
    }

    FILE *fp = fopen(input, "r");
    if (!fp) {
        fprintf(stderr, "Cannot open %s\n", input);
        return 1;
    }

    fprintf(stderr, "Reading %s (min-bits=%d, min-cc=%d)...\n", input, min_bits, min_cc);

    memset(buckets, 0, sizeof(buckets));
    memset(combos, 0, sizeof(combos));

    char line[4096];
    int total = 0, skipped = 0, filtered = 0;
    clock_t t0 = clock();

    while (fgets(line, sizeof(line), fp)) {
        char *s = line;
        while (*s && isspace((unsigned char)*s)) s++;
        size_t len = strlen(s);
        while (len > 0 && isspace((unsigned char)s[len - 1])) s[--len] = '\0';
        if (len == 0) continue;

        /* Parse CC label */
        if (s[0] != 'C' || s[1] != 'C') { skipped++; continue; }
        char *end;
        int cc = (int)strtol(s + 2, &end, 10);
        if (end == s + 2 || cc < 1 || cc > MAX_CC) { skipped++; continue; }

        /* Skip whitespace to hex */
        char *p2 = end;
        while (*p2 && isspace((unsigned char)*p2)) p2++;
        if (*p2 == '0' && (p2[1] == 'x' || p2[1] == 'X')) p2 += 2;

        char *hex_start = p2;
        int hex_len = 0;
        while (hexval(p2[hex_len]) >= 0) hex_len++;
        if (hex_len == 0) { skipped++; continue; }

        /* Compute actual bit length from hex value.
         * The third column in CC_x.txt is the hex digit count, NOT bit size.
         * Bit length = (hex_len - 1) * 4 + bits_in_leading_digit */
        int lead = hexval(hex_start[0]);
        int lead_bits = 0;
        if (lead > 0) { int v = lead; while (v) { lead_bits++; v >>= 1; } }
        int bits = (hex_len - 1) * 4 + lead_bits;
        /* Skip the third column (hex digit count) — we don't use it */
        char *p3 = p2 + hex_len;
        while (*p3 && isspace((unsigned char)*p3)) p3++;
        if (*p3) {
            /* just consume it, we computed bits from hex above */
            (void)strtol(p3, &end, 10);
        } else {
            bits = hex_len * 4;
        }
        if (bits > MAX_BITS) bits = MAX_BITS;

        /* Apply filters */
        if (bits < min_bits || cc < min_cc) { filtered++; continue; }

        /* Process */
        Bucket *b = &buckets[cc][bits];
        b->count++;

        int immune_mask = 0;

        for (int pi = 0; pi < NPRIMES; pi++) {
            int q = PRIMES[pi];
            int p_mod = hex_mod(hex_start, hex_len, q);
            if (p_mod < 0) continue;

            int pp1_mod = (p_mod + 1) % q;
            b->residue[pi][pp1_mod]++;

            int kp = kill_position(p_mod, q);
            if (kp == -1) {
                b->immune[pi]++;
                b->safe[pi]++;
                immune_mask |= (1 << pi);
            } else if (kp == -2) {
                b->safe[pi]++;
            } else if (kp >= 0) {
                if (kp >= cc) b->safe[pi]++;
                if (kp < MAX_CC) b->kill_hist[pi][kp]++;
            }
        }

        combos[cc].count++;
        combos[cc].combo_count[immune_mask]++;

        total++;
        if (total % 200000 == 0) {
            fprintf(stderr, "  [%d roots processed]\n", total);
        }
    }
    fclose(fp);

    double elapsed = (double)(clock() - t0) / CLOCKS_PER_SEC;
    fprintf(stderr, "Processed %d roots (%d skipped, %d filtered) in %.2fs\n",
            total, skipped, filtered, elapsed);

    /* ================================================================
     * Write JSON
     * ================================================================ */
    jf = fopen(output, "w");
    if (!jf) {
        fprintf(stderr, "Cannot write %s\n", output);
        return 1;
    }

    j_obj_start();

    /* ---- meta ---- */
    j_key("meta"); j_obj_start();
    j_key_str("source_file", input);
    j_key_int("total_roots", total);
    j_key_int("skipped_lines", skipped);
    j_key_int("filtered_out", filtered);
    j_key_int("min_bits", min_bits);
    j_key_int("min_cc", min_cc);
    j_key("primes"); j_arr_start();
    for (int pi = 0; pi < NPRIMES; pi++) j_int(PRIMES[pi]);
    j_arr_end();
    j_obj_end();

    /* ---- group_theory ---- */
    j_key("group_theory"); j_arr_start();
    for (int pi = 0; pi < NPRIMES; pi++) {
        int q = PRIMES[pi];
        int o = ord2(q);
        /* Compute subgroup <2> mod q */
        int sg[32], sg_len = 0;
        int v = 1;
        for (int i = 0; i < q; i++) {
            sg[sg_len++] = v;
            v = (v * 2) % q;
            if (v == 1) break;
        }
        /* Safe set = {0} union residues NOT in <2> */
        int safe_set[32], safe_len = 0;
        safe_set[safe_len++] = 0;
        for (int r = 1; r < q; r++) {
            int in_sg = 0;
            for (int k = 0; k < sg_len; k++)
                if (sg[k] == r) { in_sg = 1; break; }
            if (!in_sg) safe_set[safe_len++] = r;
        }

        j_comma(); j_obj_start(); j_first = 1;
        j_key_int("q", q);
        j_key_int("ord2q", o);
        j_key_int("subgroup_size", sg_len);
        j_key("subgroup"); j_arr_start();
        for (int k = 0; k < sg_len; k++) j_int(sg[k]);
        j_arr_end();
        j_key_int("safe_count", safe_len);
        j_key("safe_set"); j_arr_start();
        for (int k = 0; k < safe_len; k++) j_int(safe_set[k]);
        j_arr_end();
        j_key_dbl("expected_immune_among_safe", safe_len > 0 ? 1.0 / safe_len : 0);
        j_obj_end();
    }
    j_arr_end();

    /* ---- per_cc: aggregated per chain length ---- */
    j_key("per_cc"); j_arr_start();
    for (int cc = 1; cc <= MAX_CC; cc++) {
        int cc_count = 0;
        int cc_imm[NPRIMES] = {0}, cc_safe[NPRIMES] = {0};
        for (int bits = 0; bits <= MAX_BITS; bits++) {
            Bucket *b = &buckets[cc][bits];
            if (b->count == 0) continue;
            cc_count += b->count;
            for (int pi = 0; pi < NPRIMES; pi++) {
                cc_imm[pi] += b->immune[pi];
                cc_safe[pi] += b->safe[pi];
            }
        }
        if (cc_count == 0) continue;

        j_comma(); j_obj_start(); j_first = 1;
        j_key_int("cc", cc);
        j_key_int("count", cc_count);
        j_key("immune_pct"); j_arr_start();
        for (int pi = 0; pi < NPRIMES; pi++)
            j_dbl(100.0 * cc_imm[pi] / cc_count);
        j_arr_end();
        j_key("safe_pct"); j_arr_start();
        for (int pi = 0; pi < NPRIMES; pi++)
            j_dbl(100.0 * cc_safe[pi] / cc_count);
        j_arr_end();
        j_obj_end();
    }
    j_arr_end();

    /* ---- per_bucket: full detail ---- */
    j_key("per_bucket"); j_arr_start();
    for (int cc = 1; cc <= MAX_CC; cc++) {
        for (int bits = 0; bits <= MAX_BITS; bits++) {
            Bucket *b = &buckets[cc][bits];
            if (b->count == 0) continue;

            j_comma(); j_obj_start(); j_first = 1;
            j_key_int("cc", cc);
            j_key_int("bits", bits);
            j_key_int("count", b->count);
            j_key("immune"); j_arr_start();
            for (int pi = 0; pi < NPRIMES; pi++) j_int(b->immune[pi]);
            j_arr_end();
            j_key("safe"); j_arr_start();
            for (int pi = 0; pi < NPRIMES; pi++) j_int(b->safe[pi]);
            j_arr_end();
            j_obj_end();
        }
    }
    j_arr_end();

    /* ---- kill_positions ---- */
    j_key("kill_positions"); j_arr_start();
    for (int cc = 1; cc <= MAX_CC; cc++) {
        for (int bits = 0; bits <= MAX_BITS; bits++) {
            Bucket *b = &buckets[cc][bits];
            if (b->count == 0) continue;
            for (int pi = 0; pi < NPRIMES; pi++) {
                for (int j = 0; j < MAX_CC; j++) {
                    if (b->kill_hist[pi][j] == 0) continue;
                    j_comma(); j_obj_start(); j_first = 1;
                    j_key_int("cc", cc);
                    j_key_int("bits", bits);
                    j_key_int("q", PRIMES[pi]);
                    j_key_int("pos", j);
                    j_key_int("count", b->kill_hist[pi][j]);
                    j_obj_end();
                }
            }
        }
    }
    j_arr_end();

    /* ---- combinations: immune fingerprints per CC ---- */
    j_key("combinations"); j_arr_start();
    for (int cc = 1; cc <= MAX_CC; cc++) {
        if (combos[cc].count == 0) continue;

        /* Collect non-zero combos, sort by count descending */
        typedef struct { int mask; int count; } MC;
        MC mc[NCOMBOS];
        int nmc = 0;
        for (int m = 0; m < NCOMBOS; m++) {
            if (combos[cc].combo_count[m] > 0) {
                mc[nmc].mask = m;
                mc[nmc].count = combos[cc].combo_count[m];
                nmc++;
            }
        }
        /* Simple insertion sort (max 1024 entries) */
        for (int i = 1; i < nmc; i++) {
            MC tmp = mc[i];
            int j = i - 1;
            while (j >= 0 && mc[j].count < tmp.count) {
                mc[j + 1] = mc[j];
                j--;
            }
            mc[j + 1] = tmp;
        }

        j_comma(); j_obj_start(); j_first = 1;
        j_key_int("cc", cc);
        j_key_int("total", combos[cc].count);
        j_key_int("unique_combos", nmc);

        /* Output all combos (or top 50 if too many) */
        int out_n = nmc < 50 ? nmc : 50;
        j_key("top"); j_arr_start();
        for (int i = 0; i < out_n; i++) {
            j_comma(); j_obj_start(); j_first = 1;

            /* Build fingerprint string: list of immune primes */
            char fp_str[128];
            int fp_len = 0;
            fp_str[0] = '\0';
            for (int pi = 0; pi < NPRIMES; pi++) {
                if (mc[i].mask & (1 << pi)) {
                    if (fp_len > 0)
                        fp_len += snprintf(fp_str + fp_len, sizeof(fp_str) - fp_len, ",");
                    fp_len += snprintf(fp_str + fp_len, sizeof(fp_str) - fp_len, "%d", PRIMES[pi]);
                }
            }
            j_key_str("primes", fp_str);
            j_key_int("mask", mc[i].mask);
            j_key_int("count", mc[i].count);
            j_key_dbl("pct", 100.0 * mc[i].count / combos[cc].count);

            /* Count how many primes in this combo */
            int np = 0;
            for (int pi = 0; pi < NPRIMES; pi++)
                if (mc[i].mask & (1 << pi)) np++;
            j_key_int("num_immune", np);

            j_obj_end();
        }
        j_arr_end();
        j_obj_end();
    }
    j_arr_end();

    /* ---- residue_distribution: per prime, aggregated (p+1) mod q distribution ---- */
    j_key("residue_distribution"); j_arr_start();
    for (int pi = 0; pi < NPRIMES; pi++) {
        int q = PRIMES[pi];
        int res_total[32] = {0};
        for (int cc = 1; cc <= MAX_CC; cc++)
            for (int bits = 0; bits <= MAX_BITS; bits++)
                for (int r = 0; r < q; r++)
                    res_total[r] += buckets[cc][bits].residue[pi][r];

        j_comma(); j_obj_start(); j_first = 1;
        j_key_int("q", q);
        j_key("counts"); j_arr_start();
        for (int r = 0; r < q; r++) j_int(res_total[r]);
        j_arr_end();
        j_obj_end();
    }
    j_arr_end();

    j_obj_end();
    fprintf(jf, "\n");
    fclose(jf);

    long fsize = 0;
    FILE *sz = fopen(output, "r");
    if (sz) { fseek(sz, 0, SEEK_END); fsize = ftell(sz); fclose(sz); }
    fprintf(stderr, "Wrote %s (%ld KB)\n", output, fsize / 1024);

    /* ---- Console summary ---- */
    fprintf(stderr, "\n=== SUMMARY ===\n");
    fprintf(stderr, "  %5s  %8s %8s  %6s  %s\n", "q", "immune%", "safe%", "ord2q", "top combo contribution");
    fprintf(stderr, "  %5s  %8s %8s  %6s  %s\n", "-----", "--------", "--------", "------", "---------------------");
    for (int pi = 0; pi < NPRIMES; pi++) {
        int q = PRIMES[pi];
        int tot_imm = 0, tot_safe = 0;
        for (int cc = 1; cc <= MAX_CC; cc++)
            for (int bits = 0; bits <= MAX_BITS; bits++) {
                tot_imm += buckets[cc][bits].immune[pi];
                tot_safe += buckets[cc][bits].safe[pi];
            }
        fprintf(stderr, "  q=%2d: %7.2f%% %7.2f%%  %5d\n",
                q, 100.0 * tot_imm / total, 100.0 * tot_safe / total, ord2(q));
    }

    /* Top 5 global combos */
    fprintf(stderr, "\n  Top immune combinations (all CCs):\n");
    int g_combo[NCOMBOS] = {0};
    for (int cc = 1; cc <= MAX_CC; cc++)
        for (int m = 0; m < NCOMBOS; m++)
            g_combo[m] += combos[cc].combo_count[m];

    for (int rank = 0; rank < 10; rank++) {
        int best_m = -1, best_c = 0;
        for (int m = 0; m < NCOMBOS; m++) {
            if (g_combo[m] > best_c) { best_c = g_combo[m]; best_m = m; }
        }
        if (best_m < 0) break;

        char fp_str[128];
        int fp_len = 0;
        fp_str[0] = '\0';
        int np = 0;
        for (int pi = 0; pi < NPRIMES; pi++) {
            if (best_m & (1 << pi)) {
                if (fp_len > 0) fp_len += snprintf(fp_str + fp_len, sizeof(fp_str) - fp_len, ",");
                fp_len += snprintf(fp_str + fp_len, sizeof(fp_str) - fp_len, "%d", PRIMES[pi]);
                np++;
            }
        }
        fprintf(stderr, "  #%d: {%s} (%d primes) — %d roots (%.2f%%)\n",
                rank + 1, fp_str, np, best_c, 100.0 * best_c / total);
        g_combo[best_m] = 0;
    }

    return 0;
}
