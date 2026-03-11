/*
 * Copyright (c) 2026 Nenad Mićić <nenad@micic.be>
 * SPDX-License-Identifier: Apache-2.0
 *
 * =============================================================================
 * Cunningham Chain Constructor - GMP v33-03 (CC18a Optimized)
 *
 * Based on v33-02. Final hardening: dynamic storage, CLI validation.
 * Target hardware: AMD Ryzen 9 7950X3D (128 MB V-Cache), RTX 5090 (future GPU).
 *
 * OPTIMIZATIONS (inherited from v31):
 *   - OPT-19: Precomputed kM_mod LUT tables (5.2 MB, fits V-Cache)
 *   - OPT-20: Per-tile-per-base residue precomputation
 *   - OPT-21: Tile*M hoist outside base loop
 *   - CC18a-ONLY: First-kind chains only
 *   - Default line-depth=17
 *
 * NEW in v32:
 *   - OPT-22: Software prefetch for kM_mod line-sieve row after L2 pass.
 *   - Simplified is_prime_pocklington_native to single witness a=2.
 *
 * NEW in v33-01:
 *   - Removed 2^127 hard ceiling on chain-follow arithmetic.
 *   - Dual-path proving: native u128 (<=127-bit chains) / GMP (>127-bit).
 *   - g_wide_mode flag computed once at startup; zero regression for <=127-bit.
 *   - L2/ext-L2/line-sieve screening unchanged (u32 residues) in all modes.
 *   - Wide-mode uses mpz_probab_prime_p (GMP BPSW) for chain links.
 *
 * FIXES in v33-02 (bulletproofing):
 *   - FIX-10: Dynamic hex string allocation (mpz_get_str(NULL,...)) —
 *             no buffer overflow at any bit size (was 256-byte stack buffers).
 *   - FIX-11: Guard mpz_urandomm against wide_num_tiles==0 (SIGFPE crash).
 *   - FIX-12: Guard prefix_bits > target_bits (was negative shift → OOM).
 *   - FIX-13: Checkpoint hex buffer enlarged to 4096 (supports ~16000-bit).
 *   - FIX-14: max_val generous: 2^(bits + 2*target) for longer chain discovery.
 *   - FIX-15: Enlarged result_strings/best_start (was fixed 256-char).
 *   - FIX-16: queue_found_result uses dynamic allocation.
 *   - FIX-17: Removed dead code (result == -1 check).
 *   - FIX-18: Removed undocumented --no-ext-l2/--ext-l2 from help.
 *
 * FIXES in v33-03:
 *   - FIX-19: CLI validation for --bits (>=2) and --target (>=1). Prevents
 *             memory exhaustion from mpz_ui_pow_ui with negative/zero exponent.
 *   - FIX-20: Arbitrary-length checkpoint tile loading via mpz_inp_str
 *             (was fscanf %4095s — truncation at >16380-bit tiles).
 *   - FIX-21: Fully dynamic result storage (char* with strdup). No fixed-size
 *             limit on result hex strings (was 1024-byte truncation at ~4092 bits).
 *   - FIX-22: Defensive guard for rc_num_tiles==0 in non-wide random mode
 *             (xoshiro_bounded already safe, but explicit guard for clarity).
 *
 * Build:
 *   gcc -O3 -march=native -flto -o cc_v33_03 cc_gmp_v33_03.c \
 *       -lgmp -lpthread -lm
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <signal.h>
#include <pthread.h>
#include <gmp.h>
#include <unistd.h>
#include <sys/time.h>
#ifdef __linux__
#include <sched.h>
#endif

typedef uint64_t u64;
typedef uint32_t u32;
typedef uint16_t u16;
typedef uint8_t  u8;

extern int g_quiet_mode;
#define QPRINTF(...) do { if (!g_quiet_mode) printf(__VA_ARGS__); } while(0)

#define MAX_CHAIN 64
#define MAX_RESULTS 100000
#define MAX_SIEVE_CHAIN_LEN 32
#define LIKELY(x)   __builtin_expect(!!(x), 1)
#define UNLIKELY(x) __builtin_expect(!!(x), 0)

/* Chain kinds — CC18a ONLY (first-kind) */
#define KIND_FIRST  1   /* 2p+1 chains, n ≡ 5 (mod 6) */

/* =============================================================================
 * LATTICE CONFIGURATION — Search profile only (CC18a)
 * ========================================================================== */
#define PROFILE_NAME "cc18a"
#define PROFILE_BANNER "v33-03-CC18a"
#define PROFILE_CRT_DESC "CRT primes: 5, 7, 11, 13, 17, 19"
#define PROFILE_L2_DESC "Bitmask: 37, 41, 43, 47, 53, 59, 61"
#define CRT_PRIMES_COUNT 6
static const u32 CRT_PRIMES[CRT_PRIMES_COUNT] = {5, 7, 11, 13, 17, 19};
static const u64 LATTICE_M = 5ULL * 7 * 11 * 13 * 17 * 19;
#define FILTER_PRIMES_COUNT 7
static const u32 FILTER_PRIMES[FILTER_PRIMES_COUNT] = {37, 41, 43, 47, 53, 59, 61};

#define WHEEL_PRIMES_COUNT 3
static const u32 WHEEL_PRIMES[WHEEL_PRIMES_COUNT] = {23, 29, 31};
#define WHEEL_PERIOD_CONST (23u * 29u * 31u)
static const u64 WHEEL_PERIOD = WHEEL_PERIOD_CONST;

/* OPT-7: Extended L2 filter primes (byte lookup, too large for bitmask) */
#define EXT_FILTER_PRIMES_COUNT 7
static const u32 EXT_FILTER_PRIMES[EXT_FILTER_PRIMES_COUNT] = {67, 71, 73, 79, 83, 89, 97};

/* OPT-6: Pre-filter primes for chain ceiling prediction */
#define PREFILTER_PRIMES_COUNT 15
static const u32 PREFILTER_PRIMES[PREFILTER_PRIMES_COUNT] = {
    67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137
};
#define PREFILTER_MAX_POS 32

#define TRIAL_PRIMES_COUNT 62
static const u32 TRIAL_PRIMES[TRIAL_PRIMES_COUNT] = {
    67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139,
    149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229,
    233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317,
    331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409
};

/* =============================================================================
 * OPT-13: Extended chain screening primes (67-863, 132 primes)
 * v26 used 83 primes (67-547) catching ~33% of chain links.
 * v27 extends to 132 primes (67-863) catching ~37% (+6.8% more survivors caught).
 * Extra cost: ~49ns/step. Savings per catch: ~15us PRP test. Clear net win.
 * ========================================================================== */
#define CHAIN_SCREEN_COUNT 132
#define CHAIN_SCREEN_MAX_PRIME 863
static const u32 CHAIN_SCREEN_PRIMES[CHAIN_SCREEN_COUNT] = {
    67,  71,  73,  79,  83,  89,  97,  101, 103, 107, 109, 113, 127, 131, 137, 139,
    149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229,
    233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317,
    331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421,
    431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521,
    523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619,
    631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733,
    739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839,
    853, 857, 859, 863
};

/* Runtime-configurable screen depth (--screen-primes N, default: all 132) */
static int g_active_screen_count = CHAIN_SCREEN_COUNT;

/* =============================================================================
 * OPT-17: MULTI-POSITION LINE-SIEVE
 * Checks positions 1..line_depth for primes 101-863 BEFORE BPSW.
 * Position 0 is handled by batched trial division (faster).
 * Kill table: line_kill[prime_idx][residue] = 1 if ANY position 1..depth is killed.
 *
 * For first-kind position i: n_i = 2^i*(n_0+1) - 1
 *   n_i ≡ 0 (mod q) when n_0 ≡ 2^(-i) - 1 (mod q)
 * For second-kind position i: n_i = 2^i*(n_0-1) + 1
 *   n_i ≡ 0 (mod q) when n_0 ≡ 1 - 2^(-i) (mod q)
 *
 * Uses CHAIN_SCREEN_PRIMES starting from index 7 (skip 67-97 already in ext-L2).
 * ========================================================================== */
#define LINE_SIEVE_START_IDX 7   /* Skip primes 67-97 (handled by ext-L2) */
#define LINE_SIEVE_COUNT (CHAIN_SCREEN_COUNT - LINE_SIEVE_START_IDX)  /* 125 primes: 101-863 */
#define LINE_SIEVE_MAX_PRIME 864 /* Array dimension: max prime + 1 */

static u8 g_line_kill_first[LINE_SIEVE_COUNT][LINE_SIEVE_MAX_PRIME];
static u32 g_M_mod_line[LINE_SIEVE_COUNT];
static u32 g_wheel_period_mod_line[LINE_SIEVE_COUNT];
static int g_line_depth = 17;  /* Default: check positions 1..17 for CC18 */

/* =============================================================================
 * OPT-19: PRECOMPUTED kM_mod LUT — eliminates Barrett/division from hot paths
 *
 * For line-sieve: g_kM_mod_line[k_offset][ls] = (k_offset % q) * (M % q) % q
 * Transposed layout: row = k_offset (WHEEL_PERIOD entries), col = prime index.
 * Each row is 125 u16 = 250 bytes = 4 cache lines. Sequential access per candidate.
 * Total: 20677 × 125 × 2 = 5.17 MB — fits in 7950X3D V-Cache (128 MB).
 *
 * For L2/ext-L2 filters: g_kmod_filter/ext stores (k%p)*(M%p)%p (repurposed).
 * Same table size as before, just different precomputed values.
 * ========================================================================== */
static u16 g_kM_mod_line[WHEEL_PERIOD_CONST][LINE_SIEVE_COUNT];

/* =============================================================================
 * OPT-18: BARRETT REDUCTION — division-free modular arithmetic
 *
 * For divisor q (u32), precompute magic multiplier M and shift s such that:
 *   n mod q = n - q * ((n * M) >> (32 + s))    (with at most 1 correction)
 *
 * Eliminates hardware divl (~25 cycles) in line-sieve hot path.
 * Replaces with mul+shift (~8 cycles). Two divisions per prime → 2× savings.
 * ========================================================================== */
typedef struct {
    u32 q;         /* the prime */
    u32 shift;     /* extra shift beyond 32 bits */
    u64 magic;     /* ceil(2^(32+shift) / q) */
} BarrettRecip;

static BarrettRecip g_barrett_line[LINE_SIEVE_COUNT];

/* Compute Barrett reciprocal for divisor q.
 * For any n < q*q (guaranteed when q ≤ 863, n < 863*20677+863 < 2^24):
 *   quotient = (u32)(((u64)n * magic) >> (32 + shift))
 *   remainder = n - quotient * q
 *   if (remainder >= q) remainder -= q;
 */
static inline void barrett_init(BarrettRecip* br, u32 q) {
    br->q = q;
    /* Choose shift = ceil(log2(q)) to ensure magic fits in u32 range for small n,
     * but we store magic as u64 to handle the multiply correctly */
    u32 s = 0;
    u32 tmp = q - 1;
    while (tmp > 0) { s++; tmp >>= 1; }
    br->shift = s;
    /* magic = ceil(2^(32+s) / q) */
    br->magic = (((u64)1 << (32 + s)) + q - 1) / q;
}

/* Fast modulo: n mod q using precomputed Barrett reciprocal.
 * n must be < 2^32 (i.e. fits in u32). */
static inline u32 barrett_mod32(u32 n, const BarrettRecip* br) {
    u32 quotient = (u32)(((u64)n * br->magic) >> (32 + br->shift));
    u32 r = n - quotient * br->q;
    if (r >= br->q) r -= br->q;
    return r;
}
#define MAX_TRIAL_BATCHES 12
typedef struct {
    unsigned long product;
    int start;
    int count;
} TrialBatch;

static TrialBatch g_trial_batches_full[MAX_TRIAL_BATCHES];
static int g_num_batches_full = 0;
static TrialBatch g_trial_batches_extl2[MAX_TRIAL_BATCHES];
static int g_num_batches_extl2 = 0;

static void init_trial_batches(void) {
    int i;
    g_num_batches_full = 0;
    i = 0;
    while (i < TRIAL_PRIMES_COUNT) {
        unsigned long prod = 1;
        int start = i;
        while (i < TRIAL_PRIMES_COUNT) {
            __uint128_t test = (__uint128_t)prod * TRIAL_PRIMES[i];
            if (test > 0xFFFFFFFFFFFFFFFFULL) break;
            prod = (unsigned long)test;
            i++;
        }
        g_trial_batches_full[g_num_batches_full].product = prod;
        g_trial_batches_full[g_num_batches_full].start = start;
        g_trial_batches_full[g_num_batches_full].count = i - start;
        g_num_batches_full++;
    }

    g_num_batches_extl2 = 0;
    i = EXT_FILTER_PRIMES_COUNT;
    while (i < TRIAL_PRIMES_COUNT) {
        unsigned long prod = 1;
        int start = i;
        while (i < TRIAL_PRIMES_COUNT) {
            __uint128_t test = (__uint128_t)prod * TRIAL_PRIMES[i];
            if (test > 0xFFFFFFFFFFFFFFFFULL) break;
            prod = (unsigned long)test;
            i++;
        }
        g_trial_batches_extl2[g_num_batches_extl2].product = prod;
        g_trial_batches_extl2[g_num_batches_extl2].start = start;
        g_trial_batches_extl2[g_num_batches_extl2].count = i - start;
        g_num_batches_extl2++;
    }
}

static u32 g_M_mod_filter[FILTER_PRIMES_COUNT];
static u32 g_M_mod_ext_filter[EXT_FILTER_PRIMES_COUNT];

static u32 g_wheel_period_mod_filter[FILTER_PRIMES_COUNT];
static u32 g_wheel_period_mod_ext_filter[EXT_FILTER_PRIMES_COUNT];
static u32 g_wheel_period_mod6;

/* k mod p lookup tables for wheel offsets [0, WHEEL_PERIOD). */
static u16 g_kmod_filter[FILTER_PRIMES_COUNT][WHEEL_PERIOD_CONST];
static u16 g_kmod_ext_filter[EXT_FILTER_PRIMES_COUNT][WHEEL_PERIOD_CONST];

/* =============================================================================
 * GLOBAL STATE
 * ========================================================================== */

volatile sig_atomic_t shutdown_requested = 0;

/* OPT-4+5: BaseWheel — first-kind only, pre-partitioned mod6 buckets */
typedef struct {
    u64 base;
    u32* wheel_by_mod6[6];
    int  wheel_size_by_mod6[6];
    int  wheel_total;
    u32* residues;
    u32* ext_residues;
    u32* line_residues;   /* OPT-17: base mod line-sieve primes */
    u32  target_k_mod6;
} BaseWheel;

BaseWheel* g_base_wheels = NULL;
int g_num_bases = 0;
int g_total_wheel_positions = 0;
double g_avg_wheel_density = 0.0;

/* L2 filter bitmasks — first-kind only */
u64 g_filter_mask_first[FILTER_PRIMES_COUNT][MAX_SIEVE_CHAIN_LEN + 1];

/* OPT-7: Extended L2 filter (byte arrays, primes 67-97) — first-kind only */
u8 g_ext_filter_first[EXT_FILTER_PRIMES_COUNT][MAX_SIEVE_CHAIN_LEN + 1][128];

/* OPT-6: Pre-filter kill-position tables — first-kind only */
u8 g_kill_pos_first[PREFILTER_PRIMES_COUNT][140];

mpz_t g_lattice_m;
mpz_t g_k_range_min;
mpz_t g_k_range_size;
mpz_t g_k_range_max;

int g_sequential_mode = 0;
int g_random_chunk_mode = 0;  /* OPT-16: random-sequential chunks (default for non-sequential) */
u64 g_chunk_tiles = 0;        /* tiles per random chunk (0 = auto-compute) */
int g_quiet_mode = 0;
int g_full_quiet_mode = 0;
int g_use_prefix = 0;
int g_use_ext_l2 = 1;  /* Always ON for CC18a */
int g_sieve_only_mode = 0;  /* Sieve-only: skip primality proving (for benchmarking/GPU portability) */
int g_wide_mode = 0;        /* Set at startup if bits+target > 127; uses GMP for proving path */
static int g_target_bits_global = 89;
static int g_target_length_global = 18;

mpz_t g_n_range_min;
mpz_t g_n_range_max;

pthread_mutex_t g_seq_lock = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t g_print_lock = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t g_log_file_lock = PTHREAD_MUTEX_INITIALIZER;
mpz_t g_current_tile;
int g_search_complete = 0;
int g_pin_threads = 0;
int g_pin_base_cpu = 0;
int g_report_interval_sec = 1;
FILE* g_log_fp = NULL;

/* Checkpoint state */
const char* g_checkpoint_file = NULL;
int g_checkpoint_interval_sec = 60;
double g_last_checkpoint_time = 0.0;

/* =============================================================================
 * RESULTS
 * ========================================================================== */

typedef struct {
    pthread_mutex_t lock;
    u64 tiles_processed;
    u64 wheel_candidates;
    u64 passed_l2_filter;
    u64 passed_ext_filter;
    u64 primes_found;
    u64 chains_found;
    u64 chains_by_length[MAX_CHAIN + 1];
    u64 non_roots_skipped;
    u64 prefilter_skipped;
    u64 chain_screen_saved;
    u64 chain_links_tested;
    u64 line_sieve_passed;    /* OPT-17 */
    u64 line_sieve_rejected;  /* OPT-17 */
    char *result_strings[MAX_RESULTS]; /* FIX-21: dynamic (strdup), no size limit */
    int result_lengths[MAX_RESULTS];
    int num_results;

    int best_length;
    char *best_start;  /* FIX-21: dynamic (strdup), no size limit */
} SearchResults;

SearchResults g_results;

/* =============================================================================
 * THREAD CONFIG WITH MPZ POOL
 * ========================================================================== */

typedef struct {
    int thread_id;
    int target_length;
    int sieve_len;
    int target_bits;
    int log_threshold;
    u64 tiles_per_batch;
    int random_chunk_mode;     /* OPT-16: thread operates in random-chunk mode */
    u64 rng_state[4];          /* OPT-16: xoshiro256** per-thread PRNG */
    gmp_randstate_t rand_state;

    mpz_t n, k, temp, max_val, root;
    mpz_t tile_start, full_k, wheel_period;
    mpz_t pred, current;
    mpz_t temp2;
    mpz_t prev_prime;

    u32 base_residues[FILTER_PRIMES_COUNT];

    char found_buf[8192];
    size_t found_buf_len;
} ThreadConfig;

/* =============================================================================
 * SIGNAL HANDLER & THREAD UTILS
 * ========================================================================== */

void signal_handler(int sig) {
    (void)sig;
    shutdown_requested = 1;
}

static void pin_thread_if_requested(int thread_id) {
#ifdef __linux__
    if (!g_pin_threads) return;
    long ncpu = sysconf(_SC_NPROCESSORS_ONLN);
    if (ncpu <= 0) return;
    int cpu = (g_pin_base_cpu + thread_id) % (int)ncpu;
    cpu_set_t set;
    CPU_ZERO(&set);
    CPU_SET(cpu, &set);
    (void)pthread_setaffinity_np(pthread_self(), sizeof(set), &set);
#else
    (void)thread_id;
#endif
}

static inline void flush_thread_found_buffer(ThreadConfig* cfg) {
    if (!g_full_quiet_mode || !g_log_fp || cfg->found_buf_len == 0) return;
    pthread_mutex_lock(&g_log_file_lock);
    fwrite(cfg->found_buf, 1, cfg->found_buf_len, g_log_fp);
    pthread_mutex_unlock(&g_log_file_lock);
    cfg->found_buf_len = 0;
}

/* FIX-16: Compute line size from hex string length (no fixed truncation). */
static inline void queue_found_result(ThreadConfig* cfg, int length, int kind_unused, const char* n_hex) {
    (void)kind_unused;
    if (!g_full_quiet_mode || !g_log_fp) return;

    size_t hex_len = strlen(n_hex);
    size_t line_sz = hex_len + 32;  /* "CC99a 0x" + hex + "\n" + NUL */
    char line_stack[320];
    char *line = (line_sz <= sizeof(line_stack)) ? line_stack : (char*)malloc(line_sz);
    if (!line) return;

    int n = snprintf(line, line_sz, "CC%da 0x%s\n", length, n_hex);
    if (n <= 0) { if (line != line_stack) free(line); return; }

    if ((size_t)n >= sizeof(cfg->found_buf)) {
        pthread_mutex_lock(&g_log_file_lock);
        fwrite(line, 1, (size_t)n, g_log_fp);
        pthread_mutex_unlock(&g_log_file_lock);
        if (line != line_stack) free(line);
        return;
    }

    if (cfg->found_buf_len + (size_t)n > sizeof(cfg->found_buf)) {
        flush_thread_found_buffer(cfg);
    }

    memcpy(cfg->found_buf + cfg->found_buf_len, line, (size_t)n);
    cfg->found_buf_len += (size_t)n;
    if (line != line_stack) free(line);
}

/* FIX-4: Immediate persistence for critical chains (CC >= target).
 * Writes directly to log file with fflush, bypassing batch buffer.
 * Also writes to stderr as crash-safe backup (stderr is unbuffered). */
static void persist_critical_chain(int length, const char* kind_str, const char* hex_root,
                                    int thread_id) {
    /* Always announce on stderr (unbuffered = crash-safe) */
    fprintf(stderr, "\n*** CRITICAL FIND: CC%d%s root=0x%s [T%d] ***\n",
            length, kind_str, hex_root, thread_id);

    /* Write to log file immediately if available */
    if (g_log_fp) {
        pthread_mutex_lock(&g_log_file_lock);
        fprintf(g_log_fp, "*** CC%d%s 0x%s (immediate)\n", length, kind_str, hex_root);
        fflush(g_log_fp);
        pthread_mutex_unlock(&g_log_file_lock);
    }
}

static int g_max_results_warned = 0;  /* FIX-5: warn once on overflow */

/* OPT-19: Tables store (k%p)*(M%p)%p instead of just k%p.
 * This eliminates the multiplication and modular reduction from the hot path.
 * Hot path becomes: n_mod_p = base_tile_r + table[k_offset]; if >=p, -=p. */
static void init_kmod_tables(void) {
    for (int i = 0; i < FILTER_PRIMES_COUNT; i++) {
        u32 p = FILTER_PRIMES[i];
        u32 M_mod_p = g_M_mod_filter[i];
        for (u32 k = 0; k < WHEEL_PERIOD; k++) {
            g_kmod_filter[i][k] = (u16)((u64)(k % p) * M_mod_p % p);
        }
    }
    for (int i = 0; i < EXT_FILTER_PRIMES_COUNT; i++) {
        u32 p = EXT_FILTER_PRIMES[i];
        u32 M_mod_p = g_M_mod_ext_filter[i];
        for (u32 k = 0; k < WHEEL_PERIOD; k++) {
            g_kmod_ext_filter[i][k] = (u16)((u64)(k % p) * M_mod_p % p);
        }
    }
}

/* =============================================================================
 * CHECKPOINT SUPPORT
 * ========================================================================== */

static void write_checkpoint(void) {
    if (!g_checkpoint_file || !g_sequential_mode) return;
    char tmp_path[512];
    snprintf(tmp_path, sizeof(tmp_path), "%s.tmp", g_checkpoint_file);

    FILE* fp = fopen(tmp_path, "w");
    if (!fp) return;

    pthread_mutex_lock(&g_seq_lock);
    char* tile_hex = mpz_get_str(NULL, 16, g_current_tile);
    pthread_mutex_unlock(&g_seq_lock);

    fprintf(fp, "v33\n%s\n", tile_hex);
    free(tile_hex);
    fclose(fp);
    rename(tmp_path, g_checkpoint_file);
}

static int load_checkpoint(void) {
    if (!g_checkpoint_file) return 0;
    FILE* fp = fopen(g_checkpoint_file, "r");
    if (!fp) return 0;

    /* FIX-20: Read version with fscanf, then use mpz_inp_str for tile.
     * mpz_inp_str handles arbitrary-length hex (no fixed buffer limit). */
    char version[32];
    if (fscanf(fp, "%31s", version) != 1) {
        fclose(fp);
        return 0;
    }

    if (strcmp(version, "v33") != 0 && strcmp(version, "v32") != 0 && strcmp(version, "v31") != 0 && strcmp(version, "v30") != 0 && strcmp(version, "v29") != 0 && strcmp(version, "v28") != 0) {
        fprintf(stderr, "WARNING: Checkpoint version mismatch (%s), ignoring\n", version);
        fclose(fp);
        return 0;
    }

    if (mpz_inp_str(g_current_tile, fp, 16) == 0) {
        fprintf(stderr, "WARNING: Invalid checkpoint tile value, ignoring\n");
        fclose(fp);
        return 0;
    }
    fclose(fp);

    char *tile_hex = mpz_get_str(NULL, 16, g_current_tile);
    QPRINTF("RESUMED from checkpoint: tile 0x%s\n", tile_hex);
    free(tile_hex);
    return 1;
}

/* =============================================================================
 * MATH HELPERS
 * ========================================================================== */

static inline u64 mod_inverse(u64 a, u64 m) {
    if (m == 1) return 0;
    int64_t m0 = (int64_t)m, a0 = (int64_t)(a % m), x0 = 0, x1 = 1;
    if (a0 == 0) return 0;
    while (a0 > 1) {
        if (m0 == 0) return 0;
        int64_t q = a0 / m0, t = m0;
        m0 = a0 % m0; a0 = t;
        t = x0; x0 = x1 - q * x0; x1 = t;
    }
    if (x1 < 0) x1 += (int64_t)m;
    return (u64)x1;
}

/* OPT-5: Compute forbidden residues for given chain kind */
static void compute_forbidden_residues(u32 q, int chain_length, int kind,
                                        u32* forbidden, int* num_forbidden) {
    *num_forbidden = 0;
    int max_i = chain_length;
    if (max_i > (int)q - 1) max_i = q - 1;

    for (int i = 0; i < max_i; i++) {
        u64 pow2i = 1;
        for (int j = 0; j < i; j++) pow2i = (pow2i * 2) % q;
        u64 inv_pow2i = mod_inverse(pow2i, q);

        u32 r;
        if (kind == KIND_FIRST) {
            r = (u32)((inv_pow2i + q - 1) % q);
        } else {
            r = (u32)((1 + q - inv_pow2i) % q);
        }

        int found = 0;
        for (int j = 0; j < *num_forbidden; j++) {
            if (forbidden[j] == r) { found = 1; break; }
        }
        if (!found) forbidden[(*num_forbidden)++] = r;
    }
}

static void compute_valid_residues(u32 q, int chain_length, int kind,
                                    u32* valid, int* num_valid) {
    u32 forbidden[256];
    int num_forbidden;
    compute_forbidden_residues(q, chain_length, kind, forbidden, &num_forbidden);
    *num_valid = 0;
    for (u32 r = 0; r < q; r++) {
        int is_forbidden = 0;
        for (int j = 0; j < num_forbidden; j++) {
            if (forbidden[j] == r) { is_forbidden = 1; break; }
        }
        if (!is_forbidden) valid[(*num_valid)++] = r;
    }
}

static u64 crt_combine(u64 r1, u64 m1, u64 r2, u64 m2) {
    u64 inv = mod_inverse(m1 % m2, m2);
    u64 diff = (r2 >= r1 % m2) ? (r2 - r1 % m2) : (m2 - (r1 % m2 - r2));
    return r1 + m1 * ((diff * inv) % m2);
}

/* =============================================================================
 * LEVEL 0: CRT BASE GENERATION (OPT-5: kind-aware)
 * ========================================================================== */

static int generate_crt_bases(int chain_length, int kind, u64* bases, int max_bases) {
    u32 valid[CRT_PRIMES_COUNT][64];
    int num_valid[CRT_PRIMES_COUNT];
    for (int i = 0; i < CRT_PRIMES_COUNT; i++) {
        compute_valid_residues(CRT_PRIMES[i], chain_length, kind, valid[i], &num_valid[i]);
    }
    int total_combinations = 1;
    for (int i = 0; i < CRT_PRIMES_COUNT; i++) total_combinations *= num_valid[i];
    if (total_combinations > max_bases) return 0;

    int count = 0;
    typedef struct { u64 r, m; int prime_idx; } CRTState;
    CRTState stack[CRT_PRIMES_COUNT * 64];
    int stack_top = 0;

    for (int i = 0; i < num_valid[0] && count < max_bases; i++) {
        stack[stack_top++] = (CRTState){valid[0][i], CRT_PRIMES[0], 1};
    }
    while (stack_top > 0 && count < max_bases) {
        CRTState state = stack[--stack_top];
        if (state.prime_idx >= CRT_PRIMES_COUNT) {
            bases[count++] = state.r;
        } else {
            int p_idx = state.prime_idx;
            for (int i = 0; i < num_valid[p_idx]; i++) {
                u64 combined = crt_combine(state.r, state.m, valid[p_idx][i], CRT_PRIMES[p_idx]);
                stack[stack_top++] = (CRTState){combined, state.m * CRT_PRIMES[p_idx], p_idx + 1};
            }
        }
    }
    /* Sort and deduplicate */
    for (int i = 0; i < count - 1; i++) {
        for (int j = i + 1; j < count; j++) {
            if (bases[j] < bases[i]) { u64 t = bases[i]; bases[i] = bases[j]; bases[j] = t; }
        }
    }
    int unique = 0;
    for (int i = 0; i < count; i++) {
        if (i == 0 || bases[i] != bases[unique-1]) bases[unique++] = bases[i];
    }
    return unique;
}

/* =============================================================================
 * LEVEL 1: WHEEL GENERATION (OPT-4+5: both kinds)
 * ========================================================================== */

static void generate_optimized_wheels(int chain_length) {
    g_base_wheels = NULL;
    g_num_bases = 0;
    g_total_wheel_positions = 0;

    for (int i = 0; i < FILTER_PRIMES_COUNT; i++) {
        g_M_mod_filter[i] = (u32)(LATTICE_M % FILTER_PRIMES[i]);
    }
    for (int i = 0; i < EXT_FILTER_PRIMES_COUNT; i++) {
        g_M_mod_ext_filter[i] = (u32)(LATTICE_M % EXT_FILTER_PRIMES[i]);
    }
    for (int i = 0; i < FILTER_PRIMES_COUNT; i++) {
        g_wheel_period_mod_filter[i] = (u32)(WHEEL_PERIOD % FILTER_PRIMES[i]);
    }
    for (int i = 0; i < EXT_FILTER_PRIMES_COUNT; i++) {
        g_wheel_period_mod_ext_filter[i] = (u32)(WHEEL_PERIOD % EXT_FILTER_PRIMES[i]);
    }
    g_wheel_period_mod6 = (u32)(WHEEL_PERIOD % 6);
    init_kmod_tables();

    u64 temp_bases_first[16384];
    int num_first = generate_crt_bases(chain_length, KIND_FIRST, temp_bases_first, 16384);

    g_base_wheels = (BaseWheel*)malloc(num_first * sizeof(BaseWheel));

    {
        int kind = KIND_FIRST;
        u64* temp_bases = temp_bases_first;
        int num_bases = num_first;

        for (int b = 0; b < num_bases; b++) {
            u64 base = temp_bases[b];

            u32 forbidden[WHEEL_PRIMES_COUNT][64];
            int num_forbidden[WHEEL_PRIMES_COUNT];
            u64 step[WHEEL_PRIMES_COUNT];
            u64 base_r[WHEEL_PRIMES_COUNT];

            for (int i = 0; i < WHEEL_PRIMES_COUNT; i++) {
                u32 q = WHEEL_PRIMES[i];
                compute_forbidden_residues(q, chain_length, kind, forbidden[i], &num_forbidden[i]);
                step[i] = LATTICE_M % q;
                base_r[i] = base % q;
            }

            u32 base_mod6 = (u32)(base % 6);
            u32 target_n_mod6 = (kind == KIND_FIRST) ? 5 : 1;
            u32 inv_m6 = (LATTICE_M % 6 == 1) ? 1 : 5;
            u32 target_k_mod6 = (u32)(((target_n_mod6 - base_mod6 + 6) % 6) * inv_m6 % 6);

            u32* temp_wheel = (u32*)malloc(WHEEL_PERIOD * sizeof(u32));
            int total_valid = 0;

            for (u32 k = 0; k < WHEEL_PERIOD; k++) {
                int valid = 1;
                for (int i = 0; i < WHEEL_PRIMES_COUNT && valid; i++) {
                    u32 q = WHEEL_PRIMES[i];
                    u32 r = (u32)((base_r[i] + (u64)k * step[i]) % q);
                    for (int j = 0; j < num_forbidden[i]; j++) {
                        if (r == forbidden[i][j]) { valid = 0; break; }
                    }
                }
                if (valid) temp_wheel[total_valid++] = k;
            }

            if (total_valid == 0) { free(temp_wheel); continue; }

            BaseWheel* bw = &g_base_wheels[g_num_bases];
            bw->base = base;
            bw->wheel_total = total_valid;
            bw->target_k_mod6 = target_k_mod6;

            bw->residues = (u32*)malloc(FILTER_PRIMES_COUNT * sizeof(u32));
            for (int i = 0; i < FILTER_PRIMES_COUNT; i++) {
                bw->residues[i] = (u32)(base % FILTER_PRIMES[i]);
            }

            bw->ext_residues = (u32*)malloc(EXT_FILTER_PRIMES_COUNT * sizeof(u32));
            for (int i = 0; i < EXT_FILTER_PRIMES_COUNT; i++) {
                bw->ext_residues[i] = (u32)(base % EXT_FILTER_PRIMES[i]);
            }

            /* OPT-17: Line-sieve base residues */
            bw->line_residues = (u32*)malloc(LINE_SIEVE_COUNT * sizeof(u32));
            for (int i = 0; i < LINE_SIEVE_COUNT; i++) {
                bw->line_residues[i] = (u32)(base % CHAIN_SCREEN_PRIMES[LINE_SIEVE_START_IDX + i]);
            }

            for (int m = 0; m < 6; m++) bw->wheel_size_by_mod6[m] = 0;
            for (int w = 0; w < total_valid; w++) bw->wheel_size_by_mod6[temp_wheel[w] % 6]++;
            for (int m = 0; m < 6; m++) {
                if (bw->wheel_size_by_mod6[m] > 0)
                    bw->wheel_by_mod6[m] = (u32*)malloc(bw->wheel_size_by_mod6[m] * sizeof(u32));
                else
                    bw->wheel_by_mod6[m] = NULL;
                bw->wheel_size_by_mod6[m] = 0;
            }
            for (int w = 0; w < total_valid; w++) {
                int m = temp_wheel[w] % 6;
                bw->wheel_by_mod6[m][bw->wheel_size_by_mod6[m]++] = temp_wheel[w];
            }

            free(temp_wheel);
            g_total_wheel_positions += total_valid;
            g_num_bases++;
        }
    }

    g_avg_wheel_density = (g_num_bases > 0) ?
        100.0 * g_total_wheel_positions / (g_num_bases * WHEEL_PERIOD) : 0.0;
}

/* =============================================================================
 * LEVEL 2: FILTER BITMASKS (OPT-5: both kinds)
 * ========================================================================== */

static void compute_filter_bitmasks(int max_chain_len) {
    memset(g_filter_mask_first, 0, sizeof(g_filter_mask_first));
    memset(g_ext_filter_first, 0, sizeof(g_ext_filter_first));

    for (int p_idx = 0; p_idx < FILTER_PRIMES_COUNT; p_idx++) {
        u32 q = FILTER_PRIMES[p_idx];
        for (int L = 1; L <= max_chain_len && L <= MAX_SIEVE_CHAIN_LEN; L++) {
            u64 mask = 0;
            int max_i = L;
            if (max_i > (int)q - 1) max_i = q - 1;
            for (int i = 0; i < max_i; i++) {
                u64 pow2i = 1;
                for (int j = 0; j < i; j++) pow2i = (pow2i * 2) % q;
                u64 inv = mod_inverse(pow2i, q);
                u32 forbidden_r = (u32)((inv + q - 1) % q);
                if (forbidden_r < 64) mask |= (1ULL << forbidden_r);
            }
            g_filter_mask_first[p_idx][L] = mask;
        }
    }

    for (int p_idx = 0; p_idx < EXT_FILTER_PRIMES_COUNT; p_idx++) {
        u32 q = EXT_FILTER_PRIMES[p_idx];
        for (int L = 1; L <= max_chain_len && L <= MAX_SIEVE_CHAIN_LEN; L++) {
            int max_i = L;
            if (max_i > (int)q - 1) max_i = q - 1;
            for (int i = 0; i < max_i; i++) {
                u64 pow2i = 1;
                for (int j = 0; j < i; j++) pow2i = (pow2i * 2) % q;
                u64 inv = mod_inverse(pow2i, q);
                u32 forbidden_r = (u32)((inv + q - 1) % q);
                g_ext_filter_first[p_idx][L][forbidden_r] = 1;
            }
        }
    }
}

/* =============================================================================
 * OPT-6: CHAIN PRE-FILTER INITIALIZATION
 * ========================================================================== */

static void init_prefilter_tables(void) {
    memset(g_kill_pos_first, PREFILTER_MAX_POS, sizeof(g_kill_pos_first));

    for (int p_idx = 0; p_idx < PREFILTER_PRIMES_COUNT; p_idx++) {
        u32 q = PREFILTER_PRIMES[p_idx];
        for (int i = 0; i < PREFILTER_MAX_POS && i < (int)q; i++) {
            u64 pow2i = 1;
            for (int j = 0; j < i; j++) pow2i = (pow2i * 2) % q;
            u64 inv = mod_inverse(pow2i, q);
            u32 kill_r = (u32)((inv + q - 1) % q);
            if (g_kill_pos_first[p_idx][kill_r] == PREFILTER_MAX_POS) {
                g_kill_pos_first[p_idx][kill_r] = (u8)i;
            }
        }
    }
}

/* OPT-6: Predict maximum chain length using pre-filter tables.
 * OPT-8: If out_residues != NULL, stores residues for reuse by chain screening.
 *
 * NOTE: When ext-L2 filter clears positions 0..sieve_len-1 for primes 67-97,
 * and line-sieve clears positions 1..line_depth for primes 101-863, any
 * candidate reaching this function has already survived those positions.
 * If log_threshold <= min(sieve_len, line_depth+1), the prefilter can never
 * produce ceiling < log_threshold, so prefilter_skipped will be 0.
 * This is EXPECTED BEHAVIOR, not a bug. */
static inline int prefilter_predict_ceiling(mpz_t n, int cutoff, u32* out_residues) {
    int ceiling = PREFILTER_MAX_POS;

    u8 (*table)[140] = g_kill_pos_first;

    for (int p_idx = 0; p_idx < PREFILTER_PRIMES_COUNT; p_idx++) {
        u32 q = PREFILTER_PRIMES[p_idx];
        u32 r = (u32)mpz_fdiv_ui(n, q);
        if (out_residues) out_residues[p_idx] = r;
        int kill = table[p_idx][r];
        if (kill < ceiling) {
            ceiling = kill;
            if (ceiling <= 1 || (cutoff > 0 && ceiling < cutoff)) {
                /* FIX-1: Complete remaining residues on early-exit.
                 * Was: `if (out_residues && cutoff <= 0)` — the cutoff <= 0
                 * condition was dead code in production (cutoff = log_threshold > 0),
                 * leaving garbage in out_residues[p_idx+1..14]. */
                if (out_residues) {
                    for (int j = p_idx + 1; j < PREFILTER_PRIMES_COUNT; j++)
                        out_residues[j] = (u32)mpz_fdiv_ui(n, PREFILTER_PRIMES[j]);
                }
                return ceiling;
            }
        }
    }
    return ceiling;
}

/* =============================================================================
 * OPT-17: LINE-SIEVE INITIALIZATION
 * Build kill tables for positions 1..line_depth.
 * Position 0 is handled by batched trial division (more efficient).
 *
 * line_kill[idx][r] = 1 if residue r is killed at ANY position 1..line_depth
 * ========================================================================== */

static void init_line_sieve(int line_depth) {
    if (line_depth <= 0) return;

    memset(g_line_kill_first, 0, sizeof(g_line_kill_first));

    int total_kills_first = 0;

    for (int idx = 0; idx < LINE_SIEVE_COUNT; idx++) {
        u32 q = CHAIN_SCREEN_PRIMES[LINE_SIEVE_START_IDX + idx];

        for (int pos = 1; pos <= line_depth && pos < (int)q; pos++) {
            u64 pow2 = 1;
            for (int j = 0; j < pos; j++) pow2 = (pow2 * 2) % q;
            u64 inv2pos = mod_inverse(pow2, q);

            u32 kill_r_first = (u32)((inv2pos + q - 1) % q);
            if (!g_line_kill_first[idx][kill_r_first]) {
                g_line_kill_first[idx][kill_r_first] = 1;
                total_kills_first++;
            }
        }
    }

    /* Precompute modular constants */
    for (int i = 0; i < LINE_SIEVE_COUNT; i++) {
        u32 q = CHAIN_SCREEN_PRIMES[LINE_SIEVE_START_IDX + i];
        g_M_mod_line[i] = (u32)(LATTICE_M % q);
        g_wheel_period_mod_line[i] = (u32)(WHEEL_PERIOD % q);
        barrett_init(&g_barrett_line[i], q);
    }

    /* OPT-19: Build kM_mod_line LUT — (k%q)*(M%q)%q for all wheel positions.
     * Transposed layout: g_kM_mod_line[k_offset][prime_idx].
     * Row access for a given k_offset is 125 u16 = 250 bytes = 4 cache lines. */
    for (u32 k = 0; k < WHEEL_PERIOD; k++) {
        for (int ls = 0; ls < LINE_SIEVE_COUNT; ls++) {
            u32 q = CHAIN_SCREEN_PRIMES[LINE_SIEVE_START_IDX + ls];
            u32 k_mod_q = k % q;
            g_kM_mod_line[k][ls] = (u16)((u64)k_mod_q * g_M_mod_line[ls] % q);
        }
    }

    QPRINTF("\n=== OPT-17+19: LINE-SIEVE (positions 1..%d, kM_mod LUT) ===\n", line_depth);
    QPRINTF("Primes: %d (101-%d)\n", LINE_SIEVE_COUNT,
            CHAIN_SCREEN_PRIMES[CHAIN_SCREEN_COUNT - 1]);
    QPRINTF("Kill entries (first-kind):  %d\n", total_kills_first);
    QPRINTF("Kill tables: %zu KB, kM_mod LUT: %zu KB\n",
            sizeof(g_line_kill_first) / 1024,
            sizeof(g_kM_mod_line) / 1024);

    if (!g_quiet_mode) {
        for (int i = 0; i < 5 && i < LINE_SIEVE_COUNT; i++) {
            u32 q = CHAIN_SCREEN_PRIMES[LINE_SIEVE_START_IDX + i];
            u32 kills = 0;
            for (u32 r = 0; r < q; r++)
                if (g_line_kill_first[i][r]) kills++;
            printf("  q=%u: %u/%u residues killed (%.1f%%)\n",
                   q, kills, q, 100.0 * kills / q);
        }
    }
}

/* =============================================================================
 * OPT-15: NATIVE 128-BIT MONTGOMERY ARITHMETIC
 *
 * For 89-106 bit numbers (chain links), replaces GMP mpz_powm with a hand-rolled
 * Montgomery multiplication using __uint128_t. Eliminates mpz_t allocation,
 * normalization, and function-call overhead.
 *
 * Montgomery representation: ā = a·R mod n, where R = 2^128.
 * Montgomery product: MonPro(ā, b̄) = ā·b̄·R⁻¹ mod n.
 * Uses REDC algorithm with precomputed n' = -n⁻¹ mod 2^64.
 *
 * Cost per mulmod: 8 mul64 + ~20 add/cmp ≈ 30-40 cycles
 * Cost per powm (89-bit exp): ~133 mulmods ≈ 4000-5000 cycles ≈ 1.5-2 µs
 * ========================================================================== */

typedef unsigned __int128 u128;

/* Convert mpz_t to native u128 (for numbers ≤ 128 bits) */
static inline u128 mpz_to_u128(const mpz_t x) {
    u64 lo = (u64)mpz_getlimbn(x, 0);
    u64 hi = (mpz_size(x) > 1) ? (u64)mpz_getlimbn(x, 1) : 0;
    return ((u128)hi << 64) | lo;
}

/* Convert u128 back to mpz_t */
static inline void u128_to_mpz(mpz_t x, u128 val) {
    u64 lo = (u64)val;
    u64 hi = (u64)(val >> 64);
    if (hi) {
        mpz_set_ui(x, hi);
        mpz_mul_2exp(x, x, 64);
        mpz_add_ui(x, x, lo);
    } else {
        mpz_set_ui(x, lo);
    }
}

/* ==========================================================================
 * OPT-16: Fast thread-local PRNG (xoshiro256** by Blackman & Vigna)
 * Used for random-chunk mode tile position selection.
 * ========================================================================== */

static inline u64 xoshiro_rotl(u64 x, int k) { return (x << k) | (x >> (64 - k)); }

static inline u64 xoshiro_next(u64 s[4]) {
    u64 result = xoshiro_rotl(s[1] * 5, 7) * 9;
    u64 t = s[1] << 17;
    s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3];
    s[2] ^= t; s[3] = xoshiro_rotl(s[3], 45);
    return result;
}

static inline void xoshiro_seed(u64 s[4], u64 seed) {
    /* splitmix64 seeding */
    for (int i = 0; i < 4; i++) {
        seed += 0x9E3779B97F4A7C15ULL;
        u64 z = seed;
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        s[i] = z ^ (z >> 31);
    }
}

/* Generate random u64 in [0, range) without bias (rejection sampling) */
static inline u64 xoshiro_bounded(u64 s[4], u64 range) {
    if (range <= 1) return 0;
    u64 threshold = (-range) % range;  /* = (2^64 - range) % range */
    for (;;) {
        u64 r = xoshiro_next(s);
        if (r >= threshold) return r % range;
    }
}

/* Compute -n⁻¹ mod 2^64 using Hensel lifting (Newton's method) */
static inline u64 compute_ninv(u64 n0) {
    /* n0 must be odd. Start with n⁻¹ ≡ 1 (mod 2), lift to mod 2^64. */
    u64 inv = 1;
    inv *= 2 - n0 * inv;  /* correct to 2 bits  */
    inv *= 2 - n0 * inv;  /* correct to 4 bits  */
    inv *= 2 - n0 * inv;  /* correct to 8 bits  */
    inv *= 2 - n0 * inv;  /* correct to 16 bits */
    inv *= 2 - n0 * inv;  /* correct to 32 bits */
    inv *= 2 - n0 * inv;  /* correct to 64 bits */
    return -inv;           /* -n⁻¹ mod 2^64     */
}

/* Montgomery context for a given modulus */
typedef struct {
    u64 n[2];      /* modulus limbs (lo, hi) */
    u64 ninv;      /* -n⁻¹ mod 2^64         */
    u64 r1[2];     /* R mod n = 1 in Montgomery form */
    u64 r2[2];     /* R² mod n (for converting to Montgomery form) */
} MontCtx;

/* Compute R mod n = 2^128 mod n, AND R² mod n = 2^256 mod n,
 * using repeated doubling. Avoids expensive __udivti3 (software 128-bit div). */
static inline void compute_r_and_r2(u128 n128, u128 *r_out, u128 *r2_out) {
    u128 x = 1;
    for (int i = 0; i < 128; i++) {
        x <<= 1;
        if (x >= n128) x -= n128;
    }
    *r_out = x;
    for (int i = 0; i < 128; i++) {
        x <<= 1;
        if (x >= n128) x -= n128;
    }
    *r2_out = x;
}

/* Initialize Montgomery context for modulus n */
static inline void mont_ctx_init(MontCtx* ctx, u128 n128) {
    ctx->n[0] = (u64)n128;
    ctx->n[1] = (u64)(n128 >> 64);
    ctx->ninv = compute_ninv(ctx->n[0]);

    u128 r, r2;
    compute_r_and_r2(n128, &r, &r2);
    ctx->r1[0] = (u64)r;
    ctx->r1[1] = (u64)(r >> 64);
    ctx->r2[0] = (u64)r2;
    ctx->r2[1] = (u64)(r2 >> 64);
}

/* Montgomery multiplication: r = a * b * R⁻¹ mod n
 * a, b must be in Montgomery form (i.e., a = ā·R mod n).
 * Uses schoolbook 2×2 multiply + 2-step REDC. */
static inline void __attribute__((hot))
mont_mul(u64 r[2], const u64 a[2], const u64 b[2], const MontCtx* ctx) {
    const u64 *n = ctx->n;
    const u64 ninv = ctx->ninv;

    /* Schoolbook 2×2 → 4 limbs */
    u128 p00 = (u128)a[0] * b[0];
    u128 p01 = (u128)a[0] * b[1];
    u128 p10 = (u128)a[1] * b[0];
    u128 p11 = (u128)a[1] * b[1];

    u64 T0 = (u64)p00;
    u128 mid = (p00 >> 64) + (u64)p01 + (u64)p10;
    u64 T1 = (u64)mid;
    u128 hi = (mid >> 64) + (p01 >> 64) + (p10 >> 64) + (u64)p11;
    u64 T2 = (u64)hi;
    u64 T3 = (u64)(hi >> 64) + (u64)(p11 >> 64);

    /* REDC step 0 */
    u64 m0 = T0 * ninv;
    u128 c = (u128)m0 * n[0] + T0;
    c = (c >> 64) + (u128)m0 * n[1] + T1;
    T1 = (u64)c;
    c = (c >> 64) + T2;
    T2 = (u64)c;
    c = (c >> 64) + T3;
    T3 = (u64)c;
    u64 T4 = (u64)(c >> 64);

    /* REDC step 1 */
    u64 m1 = T1 * ninv;
    c = (u128)m1 * n[0] + T1;
    c = (c >> 64) + (u128)m1 * n[1] + T2;
    r[0] = (u64)c;
    c = (c >> 64) + T3;
    r[1] = (u64)c;
    u64 carry = (u64)(c >> 64) + T4;

    /* Conditional subtract: if result >= n */
    if (carry || r[1] > n[1] || (r[1] == n[1] && r[0] >= n[0])) {
        u128 val = ((u128)r[1] << 64) | r[0];
        val -= ((u128)n[1] << 64) | n[0];
        r[0] = (u64)val;
        r[1] = (u64)(val >> 64);
    }
}

/* Montgomery squaring: r = a² * R⁻¹ mod n (3 muls instead of 4) */
static inline void __attribute__((hot))
mont_sqr(u64 r[2], const u64 a[2], const MontCtx* ctx) {
    const u64 *n = ctx->n;
    const u64 ninv = ctx->ninv;

    u128 p00 = (u128)a[0] * a[0];
    u128 p01 = (u128)a[0] * a[1];  /* cross term, doubled */
    u128 p11 = (u128)a[1] * a[1];
    u128 dp01 = p01 << 1;          /* safe: p01 <= 2^127, dp01 <= 2^128-2, fits u128 */

    u64 T0 = (u64)p00;
    u128 mid = (p00 >> 64) + (u64)dp01;
    u64 T1 = (u64)mid;
    u128 hi = (mid >> 64) + (dp01 >> 64) + (u64)p11;
    u64 T2 = (u64)hi;
    u64 T3 = (u64)(hi >> 64) + (u64)(p11 >> 64);

    /* REDC step 0 */
    u64 m0 = T0 * ninv;
    u128 c = (u128)m0 * n[0] + T0;
    c = (c >> 64) + (u128)m0 * n[1] + T1;
    T1 = (u64)c;
    c = (c >> 64) + T2;
    T2 = (u64)c;
    c = (c >> 64) + T3;
    T3 = (u64)c;
    u64 T4 = (u64)(c >> 64);

    /* REDC step 1 */
    u64 m1 = T1 * ninv;
    c = (u128)m1 * n[0] + T1;
    c = (c >> 64) + (u128)m1 * n[1] + T2;
    r[0] = (u64)c;
    c = (c >> 64) + T3;
    r[1] = (u64)c;
    u64 carry = (u64)(c >> 64) + T4;

    if (carry || r[1] > n[1] || (r[1] == n[1] && r[0] >= n[0])) {
        u128 val = ((u128)r[1] << 64) | r[0];
        val -= ((u128)n[1] << 64) | n[0];
        r[0] = (u64)val;
        r[1] = (u64)(val >> 64);
    }
}

/* Convert value to Montgomery form: ā = a * R mod n = MonPro(a, R²) */
static inline void mont_to(u64 r[2], const u64 a[2], const MontCtx* ctx) {
    mont_mul(r, a, ctx->r2, ctx);
}

/* Convert from Montgomery form: a = ā * R⁻¹ mod n = MonPro(ā, 1) */
static inline void mont_from(u64 r[2], const u64 a[2], const MontCtx* ctx) {
    u64 one[2] = {1, 0};
    mont_mul(r, a, one, ctx);
}

/* Modular exponentiation: base^exp mod n, using Montgomery.
 * Returns result as plain (non-Montgomery) value. */
static inline u128 __attribute__((hot))
mont_powm(u128 base, u128 exp, const MontCtx* ctx) {
    if (exp == 0) return 1;

    /* Convert base to Montgomery form */
    u64 b[2] = {(u64)base, (u64)(base >> 64)};
    u64 bm[2];
    mont_to(bm, b, ctx);

    /* Find top bit of exponent to skip leading zeros */
    int top_bit;
    u64 exp_hi = (u64)(exp >> 64);
    u64 exp_lo = (u64)exp;
    if (exp_hi) {
        top_bit = 127 - __builtin_clzll(exp_hi);
    } else {
        top_bit = 63 - __builtin_clzll(exp_lo);
    }

    /* Result starts at base (Montgomery form) — first 1-bit is the top bit */
    u64 rm[2] = {bm[0], bm[1]};
    u64 tmp[2];

    /* Binary exponentiation from second-highest bit down */
    for (int bit = top_bit - 1; bit >= 0; bit--) {
        mont_sqr(tmp, rm, ctx);
        rm[0] = tmp[0]; rm[1] = tmp[1];

        if ((exp >> bit) & 1) {
            mont_mul(tmp, rm, bm, ctx);
            rm[0] = tmp[0]; rm[1] = tmp[1];
        }
    }

    /* Convert back from Montgomery form */
    u64 result[2];
    mont_from(result, rm, ctx);
    return ((u128)result[1] << 64) | result[0];
}

/* Optimized: keep x in Montgomery form throughout the squaring loop */
static inline int __attribute__((hot))
miller_rabin_mont(u128 n, u64 a_val, const MontCtx* ctx) {
    u128 nm1 = n - 1;
    int s = 0;
    {
        u128 tmp = nm1;
        while (!(tmp & 1)) { tmp >>= 1; s++; }
    }
    u128 d = nm1 >> s;

    /* Compute a^d mod n using Montgomery powm, get result in plain form */
    u128 x = mont_powm(a_val, d, ctx);

    if (x == 1 || x == nm1) return 1;

    /* Keep squaring in Montgomery form for efficiency */
    u64 xm[2], tmp2[2];
    u64 xv[2] = {(u64)x, (u64)(x >> 64)};
    mont_to(xm, xv, ctx);

    /* Precompute nm1 in Montgomery form for comparison */
    u64 nm1_v[2] = {(u64)nm1, (u64)(nm1 >> 64)};
    u64 nm1_m[2];
    mont_to(nm1_m, nm1_v, ctx);

    for (int r = 1; r < s; r++) {
        mont_sqr(tmp2, xm, ctx);
        xm[0] = tmp2[0]; xm[1] = tmp2[1];

        /* Check x == n-1 (in Montgomery form) */
        if (xm[0] == nm1_m[0] && xm[1] == nm1_m[1]) return 1;
        /* Check x == 1 (in Montgomery form: r1) */
        if (xm[0] == ctx->r1[0] && xm[1] == ctx->r1[1]) return 0;
    }

    return 0;
}

/* =============================================================================
 * OPT-15: NATIVE PRIMALITY TESTS (replace GMP in chain hot path)
 * ========================================================================== */

/* Native Pocklington proof for first-kind chain links.
 * N = 2p + 1, so N-1 = 2p. Known factor: p.
 * Pocklington: if a^(N-1) ≡ 1 (mod N) AND gcd(a^((N-1)/p) - 1, N) = 1,
 * then N is prime. Since (N-1)/p = 2, the GCD check is gcd(a²-1, N).
 *
 * v32: Single witness a=2. For CC18a first-kind candidates:
 *   - n ≡ 5 (mod 6) guaranteed by CRT lattice construction
 *   - Therefore gcd(3, n) = 1 (n is not divisible by 3)
 *   - a=2: gcd(a²-1, n) = gcd(3, n) = 1 — always succeeds for valid primes
 *   - Eliminates the 3-witness loop and all GCD branching */
static inline int __attribute__((hot))
is_prime_pocklington_native(u128 n, const MontCtx* ctx) {
    if (n <= 5) return (n == 2 || n == 3 || n == 5);

    /* Fermat test: 2^(n-1) ≡ 1 (mod n) */
    u128 result = mont_powm(2, n - 1, ctx);
    if (result != 1) return 0;  /* definitely composite */

    /* GCD check: gcd(2²-1, n) = gcd(3, n).
     * For CC18a: n ≡ 5 (mod 6), so n is not divisible by 3.
     * gcd(3, n) = 1. Pocklington proof succeeds. */
    return 1;
}

/* Native Miller-Rabin for second-kind chain links (no Pocklington).
 * Uses bases 2, 3, 5: for candidates that passed heavy screening,
 * the probability of a 3-base MR pseudoprime is astronomically low. */
static inline int __attribute__((hot))
is_prime_mr_native(u128 n, const MontCtx* ctx) {
    /* Handle small n: MR witnesses must be < n */
    if (n <= 5) {
        return (n == 2 || n == 3 || n == 5);
    }
    if (!miller_rabin_mont(n, 2, ctx)) return 0;
    if (n > 3 && !miller_rabin_mont(n, 3, ctx)) return 0;
    if (n > 5 && !miller_rabin_mont(n, 5, ctx)) return 0;
    return 1;
}

/* OPT-15: Native root primality test for candidates that already passed
 * L2 filter + trial division. Replaces mpz_probab_prime_p(n, 1) which
 * calls GMP BPSW. Uses 3-base native Montgomery MR.
 *
 * For search mode: false positives just produce a chain that ends at len=1.
 * No correctness impact, and probability is negligible (<10^-10 at 89 bits).
 * For construct mode: chains are independently verified, so this is safe. */
static inline int __attribute__((hot))
is_prime_root_native(u128 n) {
    if (n < 2) return 0;
    if (n == 2 || n == 3 || n == 5) return 1;
    if (!(n & 1)) return 0;

    MontCtx mctx;
    mont_ctx_init(&mctx, n);

    if (!miller_rabin_mont(n, 2, &mctx)) return 0;
    if (!miller_rabin_mont(n, 3, &mctx)) return 0;
    if (!miller_rabin_mont(n, 5, &mctx)) return 0;
    return 1;
}

/* =============================================================================
 * PRIMALITY TESTING (OPT-2: reduced MR rounds)
 * ========================================================================== */

static inline int is_prime_fast(mpz_t n, mpz_t temp) {
    if (UNLIKELY(mpz_cmp_ui(n, 2) < 0)) return 0;
    if (mpz_cmp_ui(n, 2) == 0) return 1;
    if (mpz_even_p(n)) return 0;

    for (int i = 0; i < FILTER_PRIMES_COUNT; i++) {
        if (mpz_divisible_ui_p(n, FILTER_PRIMES[i])) {
            return mpz_cmp_ui(n, FILTER_PRIMES[i]) == 0;
        }
    }
    for (int i = 0; i < TRIAL_PRIMES_COUNT; i++) {
        if (mpz_divisible_ui_p(n, TRIAL_PRIMES[i])) {
            return mpz_cmp_ui(n, TRIAL_PRIMES[i]) == 0;
        }
    }

    return mpz_probab_prime_p(n, 1) > 0;
}

/* Root primality for candidates that already passed L2 (and optionally ext-L2) */
static inline int __attribute__((hot))
is_prime_root_filtered(mpz_t n, mpz_t temp, int ext_l2_applied) {
    if (UNLIKELY(mpz_cmp_ui(n, 2) < 0)) return 0;
    if (mpz_cmp_ui(n, 2) == 0) return 1;
    if (mpz_even_p(n)) return 0;

    /* OPT-15: Convert to u128 once, use for both trial div and MR */
    u128 n128 = mpz_to_u128(n);

    /* OPT-12+15: Batched trial division using native u128 % u64 */
    const TrialBatch *batches = ext_l2_applied ? g_trial_batches_extl2 : g_trial_batches_full;
    int num_batches = ext_l2_applied ? g_num_batches_extl2 : g_num_batches_full;
    for (int b = 0; b < num_batches; b++) {
        unsigned long r = (unsigned long)(n128 % batches[b].product);
        for (int j = 0; j < batches[b].count; j++) {
            if (r % TRIAL_PRIMES[batches[b].start + j] == 0) {
                return 0;
            }
        }
    }

    /* OPT-15: Native Montgomery MR instead of GMP BPSW */
    return is_prime_root_native(n128);
}

/* =============================================================================
 * OPT-1: POCKLINGTON PRIMALITY PROVING FOR CHAIN ELEMENTS
 *
 * For first-kind chain: N = 2p + 1, so N-1 = 2p, known prime factor p.
 * (N-1)/p = 2, so GCD check uses a^2 - 1 (small constant, no modexp needed).
 * ========================================================================== */

static inline int is_prime_pocklington(mpz_t n, mpz_t p_factor, mpz_t temp, mpz_t temp2) {
    if (mpz_cmp_ui(n, 2) < 0) return 0;
    if (mpz_even_p(n)) return 0;

    for (int i = 0; i < FILTER_PRIMES_COUNT; i++) {
        if (mpz_divisible_ui_p(n, FILTER_PRIMES[i])) {
            return mpz_cmp_ui(n, FILTER_PRIMES[i]) == 0;
        }
    }
    for (int i = 0; i < TRIAL_PRIMES_COUNT; i++) {
        if (mpz_divisible_ui_p(n, TRIAL_PRIMES[i])) {
            return mpz_cmp_ui(n, TRIAL_PRIMES[i]) == 0;
        }
    }

    static const unsigned long witnesses[] = {2, 3, 5};

    for (int w = 0; w < 3; w++) {
        unsigned long a = witnesses[w];

        mpz_sub_ui(temp, n, 1);
        mpz_set_ui(temp2, a);
        mpz_powm(temp2, temp2, temp, n);

        if (mpz_cmp_ui(temp2, 1) != 0) {
            return 0;
        }

        mpz_set_ui(temp, (unsigned long)(a * a - 1));
        mpz_gcd(temp, temp, n);

        if (mpz_cmp_ui(temp, 1) == 0) {
            return 1;
        }
    }

    return mpz_probab_prime_p(n, 1) > 0;
}

/* OPT-8: Pocklington without trial division -- residue screening already done */
static inline int __attribute__((hot))
is_prime_pocklington_screened(mpz_t n, mpz_t p_factor, mpz_t temp, mpz_t temp2) {
    static const unsigned long witnesses[] = {2, 3, 5};

    for (int w = 0; w < 3; w++) {
        unsigned long a = witnesses[w];

        mpz_sub_ui(temp, n, 1);
        mpz_set_ui(temp2, a);
        mpz_powm(temp2, temp2, temp, n);

        if (mpz_cmp_ui(temp2, 1) != 0) {
            return 0;
        }

        mpz_set_ui(temp, (unsigned long)(a * a - 1));
        mpz_gcd(temp, temp, n);

        if (mpz_cmp_ui(temp, 1) == 0) {
            return 1;
        }
    }

    return mpz_probab_prime_p(n, 1) > 0;
}

/* OPT-8: is_prime_fast without trial division (for screened second-kind elements) */
static inline int __attribute__((hot))
is_prime_fast_screened(mpz_t n, mpz_t temp) {
    if (UNLIKELY(mpz_cmp_ui(n, 2) < 0)) return 0;
    if (mpz_cmp_ui(n, 2) == 0) return 1;
    if (mpz_even_p(n)) return 0;
    return mpz_probab_prime_p(n, 1) > 0;
}

/* =============================================================================
 * CHAIN FOLLOWING - FIRST KIND (2p+1)
 * OPT-8:  Residue Screening (67-547 in v26, 67-863 in v27)
 * OPT-13: Extended to 132 primes
 * OPT-14: Division-free recurrence (conditional subtract replaces % q)
 *
 * DOC-1: Why screen_caught can be 0 in production:
 * Ext-L2 clears positions 0..sieve_len-1 for primes 67-97.
 * Line-sieve clears positions 1..line_depth for primes 101-863.
 * Together, ALL screening primes (67-863) are cleared for positions
 * 0..min(sieve_len-1, line_depth). If no chain exceeds that length,
 * screening can never detect a composite (the small factors are guaranteed
 * absent). Example: with sieve_len=18 and line_depth=12, positions 0-12
 * are cleared for all primes. Only chains reaching position 13+ would
 * benefit from screening. This is correct behavior.
 * ========================================================================== */

static inline int is_chain_root_first(mpz_t n, mpz_t pred, mpz_t temp) {
    if (mpz_cmp_ui(n, 3) < 0) return 0;
    mpz_sub_ui(pred, n, 1);
    mpz_fdiv_q_2exp(pred, pred, 1);
    if (mpz_cmp_ui(pred, 2) < 0) return 1;
    return !is_prime_fast(pred, temp);
}

/*
 * OPT-14: Division-free recurrence explanation
 *
 * First-kind: r_next = (2*r + 1) mod q
 *   Since r < q, we have 2*r + 1 < 2q + 1, and for q >= 67, 2*r + 1 <= 2*(q-1)+1 = 2q-1 < 2q.
 *   So: 2r+1 < 2q, meaning at most ONE subtraction of q suffices.
 *   Result: nr = 2*r + 1; if (nr >= q) nr -= q;
 *   This replaces a ~26-cycle hardware div with a ~2-cycle cmp+cmov.
 *
 * Second-kind: r_next = (2*r - 1) mod q
 *   We compute nr = 2*r + q - 1 (add q to avoid unsigned underflow).
 *   Range: nr in [q-1, 3q-3]. Need at most TWO subtractions of q.
 */

static int __attribute__((hot))
follow_chain_first_screened(mpz_t n, mpz_t max_val, mpz_t current,
                            mpz_t temp, mpz_t temp2, mpz_t prev_prime,
                            u32* prefilter_residues, int root_known_prime,
                            u64* screen_saved, u64* links_tested) {
    int len = 0;
    const int scr_count = g_active_screen_count;

    /* OPT-15: Convert to native u128 for entire chain follow */
    u128 n128 = mpz_to_u128(n);
    u128 max128 = mpz_to_u128(max_val);
    u128 cur = n128;

    /* Initialize screening residues */
    u32 res[CHAIN_SCREEN_COUNT];
    if (prefilter_residues) {
        for (int i = 0; i < PREFILTER_PRIMES_COUNT && i < scr_count; i++)
            res[i] = prefilter_residues[i];
        for (int i = PREFILTER_PRIMES_COUNT; i < scr_count; i++)
            res[i] = (u32)(n128 % CHAIN_SCREEN_PRIMES[i]);
    } else {
        for (int i = 0; i < scr_count; i++)
            res[i] = (u32)(n128 % CHAIN_SCREEN_PRIMES[i]);
    }

    if (root_known_prime) {
        len = 1;
    } else {
        if (!is_prime_fast(n, temp)) return 0;
        len = 1;
    }

    while (len < MAX_CHAIN) {
        u128 next = (cur << 1) + 1;  /* 2p + 1 */

        if (next >= max128) break;

        if (links_tested) (*links_tested)++;

        /* OPT-13+14: Update residues with division-free recurrence and screen */
        int composite = 0;
        for (int i = 0; i < scr_count; i++) {
            u32 q = CHAIN_SCREEN_PRIMES[i];
            u32 nr = 2 * res[i] + 1;
            if (nr >= q) nr -= q;
            res[i] = nr;
            if (UNLIKELY(nr == 0)) { composite = 1; break; }
        }
        if (composite) {
            if (screen_saved) (*screen_saved)++;
            break;
        }

        /* OPT-15: Native Pocklington primality test.
         * FIX-17: Removed dead code (result == -1 path);
         * is_prime_pocklington_native only returns 0 or 1. */
        MontCtx mctx;
        mont_ctx_init(&mctx, next);
        if (!is_prime_pocklington_native(next, &mctx)) break;
        cur = next;
        len++;
    }

    /* Set mpz current to final chain value (for potential logging) */
    u128_to_mpz(current, cur);
    return len;
}

/* Legacy wrapper */
static int follow_chain_first(mpz_t n, mpz_t max_val, mpz_t current, mpz_t temp,
                               mpz_t temp2, mpz_t prev_prime) {
    return follow_chain_first_screened(n, max_val, current, temp, temp2, prev_prime,
                                       NULL, 0, NULL, NULL);
}

/* =============================================================================
 * GMP WIDE-MODE CHAIN FOLLOW (>127-bit chain members)
 *
 * Clone of follow_chain_first_screened() using mpz_t instead of u128 for
 * chain arithmetic. Screening residues remain u32 (division-free recurrence).
 * Used only when g_wide_mode == 1 (bits+target > 127).
 * ========================================================================== */

static int __attribute__((hot))
follow_chain_first_screened_gmp(mpz_t n, mpz_t max_val, mpz_t current,
                                mpz_t temp, mpz_t temp2, mpz_t prev_prime,
                                u32* prefilter_residues, int root_known_prime,
                                u64* screen_saved, u64* links_tested) {
    int len = 0;
    const int scr_count = g_active_screen_count;

    /* Initialize screening residues via GMP (no u128) */
    u32 res[CHAIN_SCREEN_COUNT];
    if (prefilter_residues) {
        for (int i = 0; i < PREFILTER_PRIMES_COUNT && i < scr_count; i++)
            res[i] = prefilter_residues[i];
        for (int i = PREFILTER_PRIMES_COUNT; i < scr_count; i++)
            res[i] = (u32)mpz_fdiv_ui(n, CHAIN_SCREEN_PRIMES[i]);
    } else {
        for (int i = 0; i < scr_count; i++)
            res[i] = (u32)mpz_fdiv_ui(n, CHAIN_SCREEN_PRIMES[i]);
    }

    if (root_known_prime) {
        len = 1;
    } else {
        if (!is_prime_fast(n, temp)) return 0;
        len = 1;
    }

    /* Work with mpz_t current for chain follow */
    mpz_set(current, n);

    while (len < MAX_CHAIN) {
        /* next = 2*cur + 1 */
        mpz_mul_2exp(temp, current, 1);
        mpz_add_ui(temp, temp, 1);

        if (mpz_cmp(temp, max_val) >= 0) break;

        if (links_tested) (*links_tested)++;

        /* OPT-13+14: Update residues with division-free recurrence and screen */
        int composite = 0;
        for (int i = 0; i < scr_count; i++) {
            u32 q = CHAIN_SCREEN_PRIMES[i];
            u32 nr = 2 * res[i] + 1;
            if (nr >= q) nr -= q;
            res[i] = nr;
            if (UNLIKELY(nr == 0)) { composite = 1; break; }
        }
        if (composite) {
            if (screen_saved) (*screen_saved)++;
            break;
        }

        /* Primality via GMP BPSW */
        if (mpz_probab_prime_p(temp, 1) <= 0) break;

        mpz_set(current, temp);
        len++;
    }

    return len;
}

/* GMP wide-mode wrapper (no screening) */
static int follow_chain_first_gmp(mpz_t n, mpz_t max_val, mpz_t current,
                                   mpz_t temp, mpz_t temp2, mpz_t prev_prime) {
    return follow_chain_first_screened_gmp(n, max_val, current, temp, temp2,
                                            prev_prime, NULL, 0, NULL, NULL);
}

/* Backward chain walkers — first-kind only */
static void follow_chain_backward_first(mpz_t n, mpz_t root_out, mpz_t pred, mpz_t temp) {
    mpz_set(root_out, n);
    int steps = 0;
    while (steps < MAX_CHAIN) {
        mpz_sub_ui(pred, root_out, 1);
        mpz_fdiv_q_2exp(pred, pred, 1);
        if (mpz_cmp_ui(pred, 2) < 0) break;
        if (!is_prime_fast(pred, temp)) break;
        mpz_set(root_out, pred);
        steps++;
    }
}

static int get_full_chain_length(mpz_t n, mpz_t max_val, mpz_t root, mpz_t pred,
                                  mpz_t current, mpz_t temp, mpz_t temp2, mpz_t prev_prime) {
    follow_chain_backward_first(n, root, pred, temp);
    if (g_wide_mode)
        return follow_chain_first_gmp(root, max_val, current, temp, temp2, prev_prime);
    return follow_chain_first(root, max_val, current, temp, temp2, prev_prime);
}

/* =============================================================================
 * WORKER: SEQUENTIAL TILE-WALK
 * ========================================================================== */

void* __attribute__((hot)) worker_sequential(void* arg) {
    ThreadConfig* cfg = (ThreadConfig*)arg;
    pin_thread_if_requested(cfg->thread_id);

    /* v33: Compute max_val dynamically based on target bits + chain length.
     * In wide mode, this exceeds 2^127; in native mode, capped at 2^127. */
    /* FIX-14: More generous max_val to avoid premature chain truncation.
     * Chain members grow as 2p+1 per step — worst case ~doubles each step.
     * Use bits + 2*target to provide ample headroom. */
    if (g_wide_mode)
        mpz_ui_pow_ui(cfg->max_val, 2, g_target_bits_global + 2 * g_target_length_global);
    else
        mpz_ui_pow_ui(cfg->max_val, 2, 127);
    mpz_set_ui(cfg->wheel_period, WHEEL_PERIOD);

    u64 local_tiles = 0, local_wheel = 0, local_filter = 0;
    u64 local_ext_filter = 0, local_primes = 0, local_chains = 0;
    u64 local_chain_counts[MAX_CHAIN + 1] = {0};
    u64 local_non_roots = 0, local_prefilter_skip = 0;
    u64 local_screen_saved = 0, local_links_tested = 0;
    u64 local_line_passed = 0, local_line_rejected = 0;
    int local_best_length = 0;  /* FIX-2: track best chain locally */

    mpz_t my_tile_start, my_tile_end;
    mpz_init(my_tile_start);
    mpz_init(my_tile_end);

    /* Precompute k-range as u128 for random-chunk tile generation AND inner loop.
     * In wide mode, u128 cannot represent these values; use mpz_t instead. */
    u128 rc_k_min = 0, rc_k_max = 0, rc_k_range = 0;
    u64 rc_num_tiles = 0;
    u128 k_min128 = 0, k_max128 = 0, n_min128 = 0, n_max128 = 0;
    /* Wide mode: full-precision tile count (may exceed u64) */
    mpz_t wide_num_tiles;
    mpz_init(wide_num_tiles);
    if (!g_wide_mode) {
        rc_k_min = mpz_to_u128(g_k_range_min);
        rc_k_max = mpz_to_u128(g_k_range_max);
        rc_k_range = rc_k_max - rc_k_min;
        rc_num_tiles = (u64)(rc_k_range / WHEEL_PERIOD);
        k_min128 = rc_k_min;
        k_max128 = rc_k_max;
        n_min128 = mpz_to_u128(g_n_range_min);
        n_max128 = mpz_to_u128(g_n_range_max);
    } else {
        /* Wide mode: tile count can exceed u64 — keep as mpz_t.
         * For 135-bit + prefix 0b1: ~2^99 tiles, far beyond u64. */
        mpz_sub(wide_num_tiles, g_k_range_max, g_k_range_min);
        mpz_fdiv_q_ui(wide_num_tiles, wide_num_tiles, WHEEL_PERIOD);
    }

    while (!shutdown_requested && !g_search_complete) {
        if (cfg->random_chunk_mode) {
            /* OPT-16: Pick random tile position, walk chunk_size tiles */
            if (!g_wide_mode) {
                /* FIX-22: Explicit guard for rc_num_tiles==0 (degenerate narrow range).
                 * xoshiro_bounded(state,0) already returns 0 safely, but this makes
                 * the intent explicit and avoids relying on that subtle behavior. */
                u64 rand_tile = (rc_num_tiles > 0)
                    ? xoshiro_bounded(cfg->rng_state, rc_num_tiles) : 0;
                u128 chunk_start = rc_k_min + (u128)rand_tile * WHEEL_PERIOD;
                u128 chunk_end = chunk_start + (u128)cfg->tiles_per_batch * WHEEL_PERIOD;
                if (chunk_end > rc_k_max) chunk_end = rc_k_max;
                u128_to_mpz(my_tile_start, chunk_start);
                u128_to_mpz(my_tile_end, chunk_end);
            } else {
                /* FIX-11: Guard against wide_num_tiles == 0 (degenerate range).
                 * mpz_urandomm is undefined for n==0; fall back to k_range_min. */
                if (mpz_sgn(wide_num_tiles) == 0) {
                    mpz_set(my_tile_start, g_k_range_min);
                    mpz_set(my_tile_end, my_tile_start);
                    mpz_addmul_ui(my_tile_end, cfg->wheel_period, cfg->tiles_per_batch);
                    if (mpz_cmp(my_tile_end, g_k_range_max) > 0)
                        mpz_set(my_tile_end, g_k_range_max);
                } else {
                /* Wide mode: tile count can exceed u64 (e.g. 2^99 for 135-bit).
                 * Use GMP random to pick a tile index with full precision. */
                mpz_t rand_tile_mpz;
                mpz_init(rand_tile_mpz);
                mpz_urandomm(rand_tile_mpz, cfg->rand_state, wide_num_tiles);
                /* chunk_start = k_range_min + rand_tile * WHEEL_PERIOD */
                mpz_mul_ui(rand_tile_mpz, rand_tile_mpz, WHEEL_PERIOD);
                mpz_add(my_tile_start, g_k_range_min, rand_tile_mpz);
                mpz_set(my_tile_end, my_tile_start);
                mpz_addmul_ui(my_tile_end, cfg->wheel_period, cfg->tiles_per_batch);
                if (mpz_cmp(my_tile_end, g_k_range_max) > 0)
                    mpz_set(my_tile_end, g_k_range_max);
                mpz_clear(rand_tile_mpz);
            } /* end else wide_num_tiles > 0 */
            }
        } else {
            /* Sequential: grab from shared cursor */
            pthread_mutex_lock(&g_seq_lock);
            if (mpz_cmp(g_current_tile, g_k_range_max) >= 0) {
                g_search_complete = 1;
                pthread_mutex_unlock(&g_seq_lock);
                break;
            }
            mpz_set(my_tile_start, g_current_tile);
            mpz_add_ui(my_tile_end, g_current_tile, cfg->tiles_per_batch * WHEEL_PERIOD);
            if (mpz_cmp(my_tile_end, g_k_range_max) > 0) mpz_set(my_tile_end, g_k_range_max);
            mpz_set(g_current_tile, my_tile_end);
            pthread_mutex_unlock(&g_seq_lock);
        }

        int first_tile = 1;
        u32 tile_mod6 = 0;
        u32 tile_residues[FILTER_PRIMES_COUNT];
        u32 tile_ext_residues[EXT_FILTER_PRIMES_COUNT];
        u32 tile_line_residues[LINE_SIEVE_COUNT];

        for (mpz_set(cfg->tile_start, my_tile_start);
             mpz_cmp(cfg->tile_start, my_tile_end) < 0 && !shutdown_requested;
             mpz_add(cfg->tile_start, cfg->tile_start, cfg->wheel_period)) {

            local_tiles++;

            if (UNLIKELY(first_tile)) {
                tile_mod6 = (u32)mpz_fdiv_ui(cfg->tile_start, 6);
                for (int i = 0; i < FILTER_PRIMES_COUNT; i++)
                    tile_residues[i] = (u32)mpz_fdiv_ui(cfg->tile_start, FILTER_PRIMES[i]);
                for (int i = 0; i < EXT_FILTER_PRIMES_COUNT; i++)
                    tile_ext_residues[i] = (u32)mpz_fdiv_ui(cfg->tile_start, EXT_FILTER_PRIMES[i]);
                if (g_line_depth > 0) {
                    for (int i = 0; i < LINE_SIEVE_COUNT; i++)
                        tile_line_residues[i] = (u32)mpz_fdiv_ui(cfg->tile_start,
                            CHAIN_SCREEN_PRIMES[LINE_SIEVE_START_IDX + i]);
                }
                first_tile = 0;
            } else {
                tile_mod6 = (tile_mod6 + g_wheel_period_mod6) % 6;
                for (int i = 0; i < FILTER_PRIMES_COUNT; i++) {
                    tile_residues[i] += g_wheel_period_mod_filter[i];
                    if (tile_residues[i] >= FILTER_PRIMES[i])
                        tile_residues[i] -= FILTER_PRIMES[i];
                }
                for (int i = 0; i < EXT_FILTER_PRIMES_COUNT; i++) {
                    tile_ext_residues[i] += g_wheel_period_mod_ext_filter[i];
                    if (tile_ext_residues[i] >= EXT_FILTER_PRIMES[i])
                        tile_ext_residues[i] -= EXT_FILTER_PRIMES[i];
                }
                if (g_line_depth > 0) {
                    for (int i = 0; i < LINE_SIEVE_COUNT; i++) {
                        tile_line_residues[i] += g_wheel_period_mod_line[i];
                        u32 q = CHAIN_SCREEN_PRIMES[LINE_SIEVE_START_IDX + i];
                        if (tile_line_residues[i] >= q)
                            tile_line_residues[i] -= q;
                    }
                }
            }

            /* OPT-16: Compute tile128 once per tile (constant across all bases).
             * In wide mode, tile128 is unused; proving path uses cfg->tile_start directly. */
            u128 tile128 = g_wide_mode ? 0 : mpz_to_u128(cfg->tile_start);

            /* OPT-21: Hoist tile*M mod p outside base loop.
             * tile_residues[i] * g_M_mod_filter[i] % p is the same for all bases
             * within a tile. Precompute once, then base loop is just add+branch.
             * Saves g_num_bases × 139 multiplies per tile (~10× fewer muls). */
            u32 tile_M_l2[FILTER_PRIMES_COUNT];
            for (int i = 0; i < FILTER_PRIMES_COUNT; i++)
                tile_M_l2[i] = (u32)((u64)tile_residues[i] * g_M_mod_filter[i] % FILTER_PRIMES[i]);

            u32 tile_M_ext[EXT_FILTER_PRIMES_COUNT];
            for (int i = 0; i < EXT_FILTER_PRIMES_COUNT; i++)
                tile_M_ext[i] = (u32)((u64)tile_ext_residues[i] * g_M_mod_ext_filter[i] % EXT_FILTER_PRIMES[i]);

            u16 tile_M_line[LINE_SIEVE_COUNT];
            if (g_line_depth > 0) {
                for (int ls = 0; ls < LINE_SIEVE_COUNT; ls++) {
                    u32 q = CHAIN_SCREEN_PRIMES[LINE_SIEVE_START_IDX + ls];
                    tile_M_line[ls] = (u16)((u64)tile_line_residues[ls] * g_M_mod_line[ls] % q);
                }
            }

            for (int b = 0; b < g_num_bases; b++) {
                BaseWheel* bw = &g_base_wheels[b];

                u32 valid_k_mod6 = (bw->target_k_mod6 + 6 - tile_mod6) % 6;
                int bucket_size = bw->wheel_size_by_mod6[valid_k_mod6];
                u32* bucket_wheel = bw->wheel_by_mod6[valid_k_mod6];

                /* OPT-19+20+21: base+tile residues. Multiply already hoisted. */
                u32 base_tile_l2[FILTER_PRIMES_COUNT];
                for (int i = 0; i < FILTER_PRIMES_COUNT; i++) {
                    u32 t = bw->residues[i] + tile_M_l2[i];
                    if (t >= FILTER_PRIMES[i]) t -= FILTER_PRIMES[i];
                    base_tile_l2[i] = t;
                }

                u32 base_tile_ext[EXT_FILTER_PRIMES_COUNT];
                for (int i = 0; i < EXT_FILTER_PRIMES_COUNT; i++) {
                    u32 t = bw->ext_residues[i] + tile_M_ext[i];
                    if (t >= EXT_FILTER_PRIMES[i]) t -= EXT_FILTER_PRIMES[i];
                    base_tile_ext[i] = t;
                }

                /* OPT-19+20+21: base+tile for line-sieve. No multiplies in base loop. */
                u16 base_tile_line[LINE_SIEVE_COUNT];
                if (g_line_depth > 0) {
                    for (int ls = 0; ls < LINE_SIEVE_COUNT; ls++) {
                        u32 t = bw->line_residues[ls] + tile_M_line[ls];
                        u32 q = CHAIN_SCREEN_PRIMES[LINE_SIEVE_START_IDX + ls];
                        if (t >= q) t -= q;
                        base_tile_line[ls] = (u16)t;
                    }
                }

                for (int w = 0; w < bucket_size; w++) {
                    u32 k_offset = bucket_wheel[w];
                    local_wheel++;

                    /* OPT-19: L2 bitmask filter — no multiply, no division.
                     * n_mod_p = base_tile_l2[i] + kmod_filter[i][k_offset]
                     * kmod_filter now stores (k%p)*(M%p)%p precomputed. */
                    int passed_l2 = 1;
                    for (int i = 0; i < FILTER_PRIMES_COUNT && passed_l2; i++) {
                        u32 n_mod_p = base_tile_l2[i] + g_kmod_filter[i][k_offset];
                        if (n_mod_p >= FILTER_PRIMES[i]) n_mod_p -= FILTER_PRIMES[i];
                        if ((g_filter_mask_first[i][cfg->sieve_len] >> n_mod_p) & 1ULL) {
                            passed_l2 = 0;
                        }
                    }
                    if (!passed_l2) continue;
                    local_filter++;

                    /* OPT-22: Prefetch kM_mod line-sieve row for this k_offset.
                     * Fires at 3.88% (L2 pass rate). Hides ~10-12 cycle L3 latency
                     * for the 0.8% of candidates that reach line-sieve.
                     * Row = 125 u16 = 250 bytes. First cache line (64 bytes = 32 primes)
                     * covers most early-exit cases (avg ~15 primes before kill). */
                    if (g_line_depth > 0) {
                        __builtin_prefetch(&g_kM_mod_line[k_offset], 0, 1);
                    }

                    /* OPT-19: Extended L2 filter — no multiply, no division. */
                    {
                        int passed_ext = 1;
                        for (int i = 0; i < EXT_FILTER_PRIMES_COUNT && passed_ext; i++) {
                            u32 n_mod_p = base_tile_ext[i] + g_kmod_ext_filter[i][k_offset];
                            if (n_mod_p >= EXT_FILTER_PRIMES[i]) n_mod_p -= EXT_FILTER_PRIMES[i];
                            if (g_ext_filter_first[i][cfg->sieve_len][n_mod_p]) {
                                passed_ext = 0;
                            }
                        }
                        if (!passed_ext) continue;
                    }
                    local_ext_filter++;

                    /* OPT-19: Line-sieve — NO Barrett, NO division in hot path.
                     * n_mod_q = base_tile_line[ls] + g_kM_mod_line[k_offset][ls]
                     * One add + one conditional subtract per prime (~6 cycles).
                     * Was: 2× Barrett + multiply (~19 cycles). ~3× faster. */
                    if (g_line_depth > 0) {
                        int line_fail = 0;
                        const u16 *kM_row = g_kM_mod_line[k_offset];
                        for (int ls = 0; ls < LINE_SIEVE_COUNT; ls++) {
                            u32 n_mod_q = base_tile_line[ls] + kM_row[ls];
                            u32 q = CHAIN_SCREEN_PRIMES[LINE_SIEVE_START_IDX + ls];
                            if (n_mod_q >= q) n_mod_q -= q;
                            if (UNLIKELY(g_line_kill_first[ls][n_mod_q])) {
                                line_fail = 1;
                                break;
                            }
                        }
                        if (LIKELY(line_fail)) {
                            local_line_rejected++;
                            continue;
                        }
                        local_line_passed++;
                    }

                    if (UNLIKELY(g_sieve_only_mode)) continue;

                    if (UNLIKELY(g_wide_mode)) {
                    /* ========================================================
                     * GMP WIDE-MODE PROVING PATH (>127-bit chains)
                     * Candidate construction, trial division, primality, and
                     * chain follow all use mpz_t / GMP. Screening residues
                     * remain u32 (unchanged).
                     * ======================================================== */

                    /* Candidate construction via GMP */
                    mpz_set(cfg->temp, cfg->tile_start);
                    mpz_add_ui(cfg->temp, cfg->temp, k_offset);
                    mpz_mul_ui(cfg->n, cfg->temp, LATTICE_M);
                    mpz_add_ui(cfg->n, cfg->n, bw->base);

                    /* Range checks */
                    if (mpz_cmp(cfg->n, g_n_range_min) < 0) continue;
                    if (mpz_cmp(cfg->n, g_n_range_max) > 0) continue;

                    /* Trial division via GMP */
                    {
                        int trial_fail = 0;
                        for (int tb = 0; tb < g_num_batches_extl2; tb++) {
                            unsigned long r = mpz_fdiv_ui(cfg->n, g_trial_batches_extl2[tb].product);
                            for (int j = 0; j < g_trial_batches_extl2[tb].count; j++) {
                                if (r % TRIAL_PRIMES[g_trial_batches_extl2[tb].start + j] == 0) {
                                    trial_fail = 1; break;
                                }
                            }
                            if (trial_fail) break;
                        }
                        if (trial_fail) continue;
                    }

                    /* Root primality via GMP BPSW */
                    if (mpz_probab_prime_p(cfg->n, 1) <= 0) continue;
                    local_primes++;

                    int is_root = is_chain_root_first(cfg->n, cfg->pred, cfg->temp);

                    if (!is_root) {
                        local_non_roots++;
                        int full_len = get_full_chain_length(cfg->n, cfg->max_val,
                            cfg->root, cfg->pred, cfg->current, cfg->temp,
                            cfg->temp2, cfg->prev_prime);
                        if (full_len > 0 && full_len <= MAX_CHAIN)
                            local_chain_counts[full_len]++;
                        if (full_len >= cfg->target_length) local_chains++;
                        if (full_len >= cfg->log_threshold) {
                            char *root_str = mpz_get_str(NULL, 16, cfg->root);
                            if (full_len >= cfg->target_length)
                                persist_critical_chain(full_len, "a", root_str, cfg->thread_id);
                            if (g_full_quiet_mode) {
                                queue_found_result(cfg, full_len, KIND_FIRST, root_str);
                                if (full_len >= cfg->target_length)
                                    flush_thread_found_buffer(cfg);
                            } else {
                                if (!g_quiet_mode) {
                                    pthread_mutex_lock(&g_print_lock);
                                    printf("  [T%d] [NON-ROOT] CC%da elem, root: 0x%s\n",
                                           cfg->thread_id, full_len, root_str);
                                    fflush(stdout);
                                    pthread_mutex_unlock(&g_print_lock);
                                }
                                pthread_mutex_lock(&g_results.lock);
                                if (g_results.num_results < MAX_RESULTS) {
                                    g_results.result_strings[g_results.num_results] = strdup(root_str);
                                    g_results.result_lengths[g_results.num_results] = full_len;
                                    g_results.num_results++;
                                } else if (!g_max_results_warned) {
                                    fprintf(stderr, "WARNING: result buffer full (%d)\n", MAX_RESULTS);
                                    g_max_results_warned = 1;
                                }
                                if (full_len > g_results.best_length) {
                                    g_results.best_length = full_len;
                                    free(g_results.best_start); g_results.best_start = strdup(root_str);
                                }
                                pthread_mutex_unlock(&g_results.lock);
                            }
                            free(root_str);
                        }
                        continue;
                    }

                    /* Root chain follow via GMP */
                    u32 pf_residues_w[PREFILTER_PRIMES_COUNT];
                    int ceiling_w = prefilter_predict_ceiling(cfg->n, cfg->log_threshold, pf_residues_w);
                    if (ceiling_w < cfg->log_threshold) {
                        local_prefilter_skip++;
                        continue;
                    }

                    int length_w = follow_chain_first_screened_gmp(cfg->n, cfg->max_val, cfg->current,
                                    cfg->temp, cfg->temp2, cfg->prev_prime,
                                    pf_residues_w, 1, &local_screen_saved, &local_links_tested);

                    if (length_w > 0 && length_w <= MAX_CHAIN)
                        local_chain_counts[length_w]++;
                    if (length_w > local_best_length)
                        local_best_length = length_w;

                    if (length_w >= cfg->log_threshold) {
                        char *n_str = mpz_get_str(NULL, 16, cfg->n);
                        if (length_w >= cfg->target_length)
                            persist_critical_chain(length_w, "a", n_str, cfg->thread_id);
                        if (g_full_quiet_mode) {
                            queue_found_result(cfg, length_w, KIND_FIRST, n_str);
                            if (length_w >= cfg->target_length)
                                flush_thread_found_buffer(cfg);
                        } else {
                            pthread_mutex_lock(&g_print_lock);
                            if (g_quiet_mode)
                                printf("CC%da 0x%s\n", length_w, n_str);
                            else
                                printf("  [T%d] *** CC%da FOUND! 0x%s\n",
                                       cfg->thread_id, length_w, n_str);
                            fflush(stdout);
                            pthread_mutex_unlock(&g_print_lock);
                            pthread_mutex_lock(&g_results.lock);
                            if (g_results.num_results < MAX_RESULTS) {
                                g_results.result_strings[g_results.num_results] = strdup(n_str);
                                g_results.result_lengths[g_results.num_results] = length_w;
                                g_results.num_results++;
                            } else if (!g_max_results_warned) {
                                fprintf(stderr, "WARNING: result buffer full (%d), chain stats still tracked\n", MAX_RESULTS);
                                g_max_results_warned = 1;
                            }
                            if (length_w > g_results.best_length) {
                                g_results.best_length = length_w;
                                free(g_results.best_start); g_results.best_start = strdup(n_str);
                            }
                            pthread_mutex_unlock(&g_results.lock);
                        }
                        free(n_str);
                    }

                    if (length_w >= cfg->target_length) local_chains++;

                    } else {
                    /* ========================================================
                     * NATIVE u128 PROVING PATH (<=127-bit chains, unchanged)
                     * ======================================================== */

                    /* Build candidate n = base + full_k * M */
                    u128 full_k128 = tile128 + k_offset;
                    if (full_k128 < k_min128) continue;
                    if (full_k128 >= k_max128) continue;

                    u128 n128 = full_k128 * (u128)LATTICE_M + bw->base;
                    if (n128 < n_min128) continue;
                    if (n128 > n_max128) continue;

                    /* OPT-15: Batched trial division + native MR */
                    {
                        int trial_fail = 0;
                        for (int tb = 0; tb < g_num_batches_extl2; tb++) {
                            unsigned long r = (unsigned long)(n128 % g_trial_batches_extl2[tb].product);
                            for (int j = 0; j < g_trial_batches_extl2[tb].count; j++) {
                                if (r % TRIAL_PRIMES[g_trial_batches_extl2[tb].start + j] == 0) {
                                    trial_fail = 1; break;
                                }
                            }
                            if (trial_fail) break;
                        }
                        if (trial_fail) continue;
                    }

                    if (!is_prime_root_native(n128)) continue;
                    local_primes++;

                    u128_to_mpz(cfg->n, n128);

                    int is_root = is_chain_root_first(cfg->n, cfg->pred, cfg->temp);

                    if (!is_root) {
                        local_non_roots++;
                        int full_len = get_full_chain_length(cfg->n, cfg->max_val,
                            cfg->root, cfg->pred, cfg->current, cfg->temp,
                            cfg->temp2, cfg->prev_prime);
                        /* FIX-3: Count ALL non-root chains, not just >= log_threshold */
                        if (full_len > 0 && full_len <= MAX_CHAIN)
                            local_chain_counts[full_len]++;
                        if (full_len >= cfg->target_length) local_chains++;
                        if (full_len >= cfg->log_threshold) {
                            char *root_str = mpz_get_str(NULL, 16, cfg->root);

                            /* FIX-4: Immediate persistence for target-length chains */
                            if (full_len >= cfg->target_length) {
                                persist_critical_chain(full_len, "a", root_str, cfg->thread_id);
                            }

                            if (g_full_quiet_mode) {
                                queue_found_result(cfg, full_len, KIND_FIRST, root_str);
                                if (full_len >= cfg->target_length)
                                    flush_thread_found_buffer(cfg);
                            } else {
                                if (!g_quiet_mode) {
                                    pthread_mutex_lock(&g_print_lock);
                                    printf("  [T%d] [NON-ROOT] CC%da elem, root: 0x%s\n",
                                           cfg->thread_id, full_len, root_str);
                                    fflush(stdout);  /* FIX-8 */
                                    pthread_mutex_unlock(&g_print_lock);
                                }
                                pthread_mutex_lock(&g_results.lock);
                                if (g_results.num_results < MAX_RESULTS) {
                                    g_results.result_strings[g_results.num_results] = strdup(root_str);
                                    g_results.result_lengths[g_results.num_results] = full_len;
                                    g_results.num_results++;
                                } else if (!g_max_results_warned) {
                                    fprintf(stderr, "WARNING: result buffer full (%d)\n", MAX_RESULTS);
                                    g_max_results_warned = 1;
                                }
                                /* FIX-6: Track best chain root */
                                if (full_len > g_results.best_length) {
                                    g_results.best_length = full_len;
                                    free(g_results.best_start); g_results.best_start = strdup(root_str);
                                }
                                pthread_mutex_unlock(&g_results.lock);
                            }
                            free(root_str);
                        }
                        continue;
                    }

                    u32 pf_residues[PREFILTER_PRIMES_COUNT];
                    int ceiling = prefilter_predict_ceiling(cfg->n, cfg->log_threshold, pf_residues);
                    if (ceiling < cfg->log_threshold) {
                        local_prefilter_skip++;
                        continue;
                    }

                    int length = follow_chain_first_screened(cfg->n, cfg->max_val, cfg->current,
                                    cfg->temp, cfg->temp2, cfg->prev_prime,
                                    pf_residues, 1, &local_screen_saved, &local_links_tested);

                    if (length > 0 && length <= MAX_CHAIN)
                        local_chain_counts[length]++;
                    /* FIX-2: Track best chain length for ALL chains */
                    if (length > local_best_length)
                        local_best_length = length;

                    if (length >= cfg->log_threshold) {
                        char *n_str = mpz_get_str(NULL, 16, cfg->n);

                        /* FIX-4: Immediate persistence for target-length chains */
                        if (length >= cfg->target_length) {
                            persist_critical_chain(length, "a", n_str, cfg->thread_id);
                        }

                        if (g_full_quiet_mode) {
                            queue_found_result(cfg, length, KIND_FIRST, n_str);
                            /* FIX-4: Immediate flush for target chains */
                            if (length >= cfg->target_length)
                                flush_thread_found_buffer(cfg);
                        } else {
                            pthread_mutex_lock(&g_print_lock);
                            if (g_quiet_mode)
                                printf("CC%da 0x%s\n", length, n_str);
                            else
                                printf("  [T%d] *** CC%da FOUND! 0x%s\n",
                                       cfg->thread_id, length, n_str);
                            /* FIX-8: Flush stdout immediately for qualifying chains */
                            fflush(stdout);
                            pthread_mutex_unlock(&g_print_lock);

                            pthread_mutex_lock(&g_results.lock);
                            if (g_results.num_results < MAX_RESULTS) {
                                g_results.result_strings[g_results.num_results] = strdup(n_str);
                                g_results.result_lengths[g_results.num_results] = length;
                                g_results.num_results++;
                            } else if (!g_max_results_warned) {
                                /* FIX-5: Warn once when results buffer is full */
                                fprintf(stderr, "WARNING: result buffer full (%d), chain stats still tracked\n", MAX_RESULTS);
                                g_max_results_warned = 1;
                            }
                            /* FIX-6: Always update best_start for best chain */
                            if (length > g_results.best_length) {
                                g_results.best_length = length;
                                free(g_results.best_start); g_results.best_start = strdup(n_str);
                            }
                            pthread_mutex_unlock(&g_results.lock);
                        }
                        free(n_str);  /* FIX-10: dynamic alloc from mpz_get_str */
                    }

                    if (length >= cfg->target_length) local_chains++;

                    } /* end if/else g_wide_mode */
                }
            }
        }

        /* Flush local counters to global */
        pthread_mutex_lock(&g_results.lock);
        g_results.tiles_processed += local_tiles;
        g_results.wheel_candidates += local_wheel;
        g_results.passed_l2_filter += local_filter;
        g_results.passed_ext_filter += local_ext_filter;
        g_results.primes_found += local_primes;
        g_results.chains_found += local_chains;
        g_results.non_roots_skipped += local_non_roots;
        g_results.prefilter_skipped += local_prefilter_skip;
        g_results.chain_screen_saved += local_screen_saved;
        g_results.chain_links_tested += local_links_tested;
        g_results.line_sieve_passed += local_line_passed;
        g_results.line_sieve_rejected += local_line_rejected;
        /* FIX-2: Flush local best chain length */
        if (local_best_length > g_results.best_length)
            g_results.best_length = local_best_length;
        for (int i = 0; i <= MAX_CHAIN; i++) {
            g_results.chains_by_length[i] += local_chain_counts[i];
            local_chain_counts[i] = 0;
        }
        pthread_mutex_unlock(&g_results.lock);

        flush_thread_found_buffer(cfg);

        local_tiles = 0; local_wheel = 0; local_filter = 0;
        local_ext_filter = 0; local_primes = 0; local_chains = 0;
        local_non_roots = 0; local_prefilter_skip = 0;
        local_screen_saved = 0; local_links_tested = 0;
        local_line_passed = 0; local_line_rejected = 0;
        local_best_length = 0;
    }

    mpz_clear(my_tile_start);
    mpz_clear(my_tile_end);
    mpz_clear(wide_num_tiles);
    flush_thread_found_buffer(cfg);
    return NULL;
}

/* CC18a: worker_random removed. Only sequential/random-chunk mode. */

/* =============================================================================
 * UTILITIES
 * ========================================================================== */

static int parse_prefix(const char* str, mpz_t value, int* bits) {
    mpz_set_ui(value, 0);
    *bits = 0;
    if (str[0] != '0' || (str[1] != 'b' && str[1] != 'B')) {
        fprintf(stderr, "ERROR: Prefix must be binary format (0b...)\n");
        return -1;
    }
    const char* p = str + 2;
    if (*p != '1') { fprintf(stderr, "ERROR: Binary prefix must start with 1\n"); return -1; }
    while (*p == '0' || *p == '1') {
        mpz_mul_2exp(value, value, 1);
        if (*p == '1') mpz_add_ui(value, value, 1);
        (*bits)++;
        p++;
    }
    return 0;
}

static void init_thread_config(ThreadConfig* cfg, int thread_id, int target_len,
                                int target_bits, int log_thresh, u64 tiles_batch) {
    cfg->thread_id = thread_id;
    cfg->target_length = target_len;
    cfg->sieve_len = (target_len <= MAX_SIEVE_CHAIN_LEN) ? target_len : MAX_SIEVE_CHAIN_LEN;
    cfg->target_bits = target_bits;
    cfg->log_threshold = log_thresh;
    cfg->tiles_per_batch = tiles_batch;
    cfg->random_chunk_mode = g_random_chunk_mode;
    cfg->found_buf_len = 0;

    /* Seed fast PRNG for random-chunk mode */
    xoshiro_seed(cfg->rng_state,
        time(NULL) ^ ((u64)thread_id * 0x9E3779B97F4A7C15ULL) ^ clock() ^ (u64)getpid());

    gmp_randinit_mt(cfg->rand_state);
    gmp_randseed_ui(cfg->rand_state, time(NULL) ^ (thread_id * 0x9E3779B97F4A7C15ULL) ^ clock());

    mpz_init(cfg->n); mpz_init(cfg->k); mpz_init(cfg->temp);
    mpz_init(cfg->max_val); mpz_init(cfg->root);
    mpz_init(cfg->tile_start); mpz_init(cfg->full_k); mpz_init(cfg->wheel_period);
    mpz_init(cfg->pred); mpz_init(cfg->current);
    mpz_init(cfg->temp2); mpz_init(cfg->prev_prime);
}

static void cleanup_thread_config(ThreadConfig* cfg) {
    gmp_randclear(cfg->rand_state);
    mpz_clear(cfg->n); mpz_clear(cfg->k); mpz_clear(cfg->temp);
    mpz_clear(cfg->max_val); mpz_clear(cfg->root);
    mpz_clear(cfg->tile_start); mpz_clear(cfg->full_k); mpz_clear(cfg->wheel_period);
    mpz_clear(cfg->pred); mpz_clear(cfg->current);
    mpz_clear(cfg->temp2); mpz_clear(cfg->prev_prime);
}

static void print_usage(const char* prog) {
    printf("Usage: %s [options]\n\n", prog);
    printf("Build profile: %s\n\n", PROFILE_NAME);
    printf("Options:\n");
    printf("  --target N         Target chain length (default: 18)\n");
    printf("  --bits N           Bit size for candidates (default: 89)\n");
    printf("  --threads N        Number of threads (default: auto)\n");
    printf("  --log N            Log chains >= N (default: 6)\n");
    printf("  --report N         Progress update interval in seconds (default: 1)\n");
    printf("  --continuous       Run until interrupted\n");
    printf("  --prefix PAT       Binary prefix (e.g., 0b101)\n");
    printf("  --sequential       Sequential tile-walk (requires --prefix)\n");
    printf("  --output FILE      Append results to file\n");
    printf("  --quiet            Minimal output (CC chains only)\n");
    printf("  --full-quiet       No console output; batch-write found CC lines to --output\n");
    printf("  --pin              Pin worker threads to CPU cores (Linux only)\n");
    printf("  --pin-base N       Base CPU index for pinning (default: 0)\n");
    /* FIX-18: --no-ext-l2/--ext-l2 removed from help — not parsed in CC18a build */
    printf("  --screen-primes N  Number of chain screening primes, 1-%d (default: %d)\n",
           CHAIN_SCREEN_COUNT, CHAIN_SCREEN_COUNT);
    printf("  --checkpoint FILE  Write checkpoint to FILE every 60s (sequential mode)\n");
    printf("  --resume           Resume from checkpoint file\n");
    printf("  --ckpt-interval N  Checkpoint interval in seconds (default: 60)\n");
    printf("  --chunk-tiles N    Tiles per random chunk (default: 500)\n");
    printf("  --line-depth N     Line-sieve: check positions 1..N for primes 101-863 (0=off, default: 17)\n");
    printf("                     Recommended: 5 for CC14+, 8 for CC16+, 17 for CC18\n");
    printf("  --sieve-only       Run only sieve filters (no primality proving); for benchmarking\n");
    printf("  --test             Run unit tests\n");
    printf("\nv33: +GMP wide-mode (>127-bit chains) +OPT-22 (prefetch) +OPT-19 (kM_mod LUT)\n");
    printf("Line-sieve checks chain positions 1..N for primes 101-863 before BPSW\n");
    printf("Searches >128-bit candidates work without truncation (GMP wide-mode)\n");
    printf("Output format: CC<len><kind> 0x<hex>  (e.g., CC7a, CC8b)\n");
}

/* =============================================================================
 * UNIT TESTS
 * ========================================================================== */

int run_unit_tests(void) {
    printf("=== %s UNIT TESTS ===\n\n", PROFILE_BANNER);
    int passed = 0, failed = 0;
    init_trial_batches();
    mpz_t n, temp, temp2, pred, prev_prime, current, max_val;
    mpz_init(n); mpz_init(temp); mpz_init(temp2); mpz_init(pred);
    mpz_init(prev_prime); mpz_init(current); mpz_init(max_val);
    mpz_ui_pow_ui(max_val, 2, 127);  /* Existing tests use native u128 path */

    /* Test 1: Known CC13 first-kind root */
    mpz_set_str(n, "3d3aa382c918303b5470bad", 16);
    if (is_prime_fast(n, temp)) { printf("  OK CC13 root is prime\n"); passed++; }
    else { printf("  FAIL CC13 root is_prime\n"); failed++; }

    if (is_chain_root_first(n, pred, temp)) { printf("  OK CC13 root is first-kind root\n"); passed++; }
    else { printf("  FAIL CC13 root check\n"); failed++; }

    int len1 = follow_chain_first(n, max_val, current, temp, temp2, prev_prime);
    if (len1 == 13) { printf("  OK CC13 chain length = %d\n", len1); passed++; }
    else { printf("  FAIL CC13 chain length = %d (expected 13)\n", len1); failed++; }

    /* Test 2: First-kind chain from small prime 2 (2 -> 5 -> 11 -> 23 -> 47, CC5) */
    mpz_set_ui(n, 2);
    int len2 = follow_chain_first(n, max_val, current, temp, temp2, prev_prime);
    if (len2 == 5) { printf("  OK First-kind from 2: CC%d (2->5->11->23->47)\n", len2); passed++; }
    else { printf("  FAIL First-kind from 2: CC%d (expected 5)\n", len2); failed++; }

    /* Test 3: Pre-filter tables */
    init_prefilter_tables();
    int prefilter_ok = 1;
    for (int i = 0; i < PREFILTER_PRIMES_COUNT; i++) {
        u32 q = PREFILTER_PRIMES[i];
        if (g_kill_pos_first[i][q-1] != PREFILTER_MAX_POS) {
            printf("  FAIL First-kind immunity: prime %u, residue %u kill_pos=%d (expected %d)\n",
                   q, q-1, g_kill_pos_first[i][q-1], PREFILTER_MAX_POS);
            prefilter_ok = 0;
            break;
        }
    }
    if (prefilter_ok) { printf("  OK First-kind immunizing residues (q-1) all have kill_pos=MAX\n"); passed++; }
    else failed++;

    /* Test 4: First-kind immunity at residue 0 (r=0 means n%q=0, killed at pos 0) */
    prefilter_ok = 1;
    for (int i = 0; i < PREFILTER_PRIMES_COUNT; i++) {
        u32 q = PREFILTER_PRIMES[i];
        /* residue q-1 is the immunizing residue for first-kind (already tested in Test 3) */
        /* Check residue 0: should have kill_pos = 0 (divisible at position 0) */
        if (g_kill_pos_first[i][0] != 0) {
            printf("  FAIL First-kind r=0: prime %u, kill_pos=%d (expected 0)\n",
                   q, g_kill_pos_first[i][0]);
            prefilter_ok = 0;
            break;
        }
    }
    if (prefilter_ok) { printf("  OK First-kind residue 0 correctly has kill_pos=0 for all primes\n"); passed++; }
    else failed++;

    /* Test 5: Pre-filter predicts CC13 ceiling >= 13 */
    mpz_set_str(n, "3d3aa382c918303b5470bad", 16);
    int ceiling = prefilter_predict_ceiling(n, 0, NULL);
    if (ceiling >= 13) { printf("  OK CC13 pre-filter ceiling = %d (>= 13)\n", ceiling); passed++; }
    else { printf("  FAIL CC13 pre-filter ceiling = %d (should be >= 13)\n", ceiling); failed++; }

    /* Test 6: Filter bitmasks -- immunizing residues not forbidden */
    compute_filter_bitmasks(18);
    int mask_ok = 1;
    for (int i = 0; i < FILTER_PRIMES_COUNT; i++) {
        u32 q = FILTER_PRIMES[i];
        if (q > 64) continue;
        if ((g_filter_mask_first[i][18] >> (q-1)) & 1ULL) {
            printf("  FAIL First-kind mask: prime %u forbids immunizing residue %u\n", q, q-1);
            mask_ok = 0;
        }
    }
    if (mask_ok) { printf("  OK Filter bitmasks: immunizing residues not forbidden\n"); passed++; }
    else failed++;

    /* Test 7: OPT-8+13 chain screening gives same result as unscreened */
    mpz_set_str(n, "3d3aa382c918303b5470bad", 16);
    u32 test_pf_res[PREFILTER_PRIMES_COUNT];
    prefilter_predict_ceiling(n, 0, test_pf_res);
    u64 test_saved = 0, test_links = 0;
    int screened_len = follow_chain_first_screened(n, max_val, current, temp, temp2, prev_prime,
                                                    test_pf_res, 1, &test_saved, &test_links);
    int plain_len = follow_chain_first(n, max_val, current, temp, temp2, prev_prime);
    if (screened_len == plain_len && screened_len == 13) {
        printf("  OK OPT-13 screened CC13 = %d (matches unscreened %d, %d primes)\n",
               screened_len, plain_len, g_active_screen_count);
        passed++;
    } else {
        printf("  FAIL OPT-13 screened=%d, unscreened=%d (expected 13)\n", screened_len, plain_len);
        failed++;
    }

    /* Test 8: OPT-14 division-free recurrence correctness */
    {
        mpz_set_str(n, "3d3aa382c918303b5470bad", 16);
        u32 r = (u32)mpz_fdiv_ui(n, 67);
        mpz_mul_2exp(temp, n, 1);
        mpz_add_ui(temp, temp, 1);
        u32 r_next = (u32)mpz_fdiv_ui(temp, 67);
        /* Division-free: nr = 2*r + 1; if (nr >= 67) nr -= 67; */
        u32 nr = 2 * r + 1;
        if (nr >= 67) nr -= 67;
        if (r_next == nr) {
            printf("  OK OPT-14 div-free recurrence: (2*%u+1) mod 67 = %u\n", r, nr);
            passed++;
        } else {
            printf("  FAIL OPT-14 div-free recurrence: actual=%u, calc=%u (%% gives %u)\n",
                   r_next, nr, (2*r+1) % 67);
            failed++;
        }
    }

    /* Test 9: OPT-14 second-kind division-free recurrence */
    {
        mpz_set_ui(n, 89);  /* small prime for easy verification */
        u32 r = (u32)mpz_fdiv_ui(n, 67);  /* 89 mod 67 = 22 */
        mpz_mul_2exp(temp, n, 1);
        mpz_sub_ui(temp, temp, 1);  /* 177 */
        u32 r_next = (u32)mpz_fdiv_ui(temp, 67);  /* 177 mod 67 = 43 */
        /* Division-free second kind: nr = 2*r + q - 1 = 44 + 66 = 110; 110-67=43 */
        u32 nr = 2 * r + 67 - 1;
        if (nr >= 67) nr -= 67;
        if (nr >= 67) nr -= 67;
        if (r_next == nr) {
            printf("  OK OPT-14 second-kind div-free: (2*%u-1) mod 67 = %u\n", r, nr);
            passed++;
        } else {
            printf("  FAIL OPT-14 second-kind: actual=%u, calc=%u\n", r_next, nr);
            failed++;
        }
    }

    /* Test 10: Division-free edge case: r=0 for second kind */
    {
        u32 q = 67;
        u32 r = 0;
        /* (2*0 - 1) mod 67 = 66 = q-1 */
        u32 nr = 2 * r + q - 1;  /* = 66 */
        if (nr >= q) nr -= q;
        if (nr >= q) nr -= q;
        if (nr == q - 1) {
            printf("  OK OPT-14 edge case r=0: (2*0-1) mod %u = %u\n", q, nr);
            passed++;
        } else {
            printf("  FAIL OPT-14 edge case r=0: got %u, expected %u\n", nr, q-1);
            failed++;
        }
    }

    /* Test 11: Exhaustive div-free verification for all residues of prime 67 */
    {
        u32 q = 67;
        int ok = 1;
        for (u32 r = 0; r < q; r++) {
            /* First kind */
            u32 expected_f = (2 * r + 1) % q;
            u32 nr_f = 2 * r + 1;
            if (nr_f >= q) nr_f -= q;
            if (nr_f != expected_f) { ok = 0; printf("  FAIL first r=%u\n", r); break; }
            /* Second kind */
            u32 expected_s = (2 * r + q - 1) % q;
            u32 nr_s = 2 * r + q - 1;
            if (nr_s >= q) nr_s -= q;
            if (nr_s >= q) nr_s -= q;
            if (nr_s != expected_s) { ok = 0; printf("  FAIL second r=%u\n", r); break; }
        }
        if (ok) { printf("  OK OPT-14 exhaustive verification for q=67: all %u residues match\n", q); passed++; }
        else failed++;
    }

    /* Test 12: Verify extended screening primes are correct */
    {
        int ok = 1;
        for (int i = 0; i < CHAIN_SCREEN_COUNT; i++) {
            u32 p = CHAIN_SCREEN_PRIMES[i];
            /* Check primality of each screening prime */
            mpz_set_ui(temp, p);
            if (mpz_probab_prime_p(temp, 25) == 0) {
                printf("  FAIL Screen prime %u at index %d is not prime!\n", p, i);
                ok = 0; break;
            }
            /* Check ordering */
            if (i > 0 && CHAIN_SCREEN_PRIMES[i] <= CHAIN_SCREEN_PRIMES[i-1]) {
                printf("  FAIL Screen primes not in order at index %d\n", i);
                ok = 0; break;
            }
        }
        if (ok) { printf("  OK All %d screening primes verified prime and ordered\n", CHAIN_SCREEN_COUNT); passed++; }
        else failed++;
    }

    /* Test 15: OPT-18 Barrett reduction correctness */
    {
        int ok = 1;
        /* Test all 125 line-sieve primes with representative inputs */
        u32 test_vals[] = {0, 1, 2, 100, 999, 5000, 10000, 20676, 50000, 100000, 500000, 745832};
        int nvals = sizeof(test_vals) / sizeof(test_vals[0]);
        for (int i = 0; i < LINE_SIEVE_COUNT && ok; i++) {
            u32 q = CHAIN_SCREEN_PRIMES[LINE_SIEVE_START_IDX + i];
            BarrettRecip br;
            barrett_init(&br, q);
            for (int j = 0; j < nvals && ok; j++) {
                u32 n = test_vals[j];
                u32 expected = n % q;
                u32 got = barrett_mod32(n, &br);
                if (got != expected) {
                    printf("  FAIL OPT-18: barrett_mod32(%u, %u) = %u (expected %u)\n",
                           n, q, got, expected);
                    ok = 0;
                }
            }
            /* Also exhaustive test for small range (0..q-1 and q..2*q) */
            for (u32 n = 0; n < 2 * q && ok; n++) {
                u32 expected = n % q;
                u32 got = barrett_mod32(n, &br);
                if (got != expected) {
                    printf("  FAIL OPT-18: barrett_mod32(%u, %u) = %u (expected %u)\n",
                           n, q, got, expected);
                    ok = 0;
                }
            }
        }
        if (ok) { printf("  OK OPT-18 Barrett reduction matches hardware div for all 125 primes\n"); passed++; }
        else failed++;
    }

    /* Test 16: OPT-15 Montgomery powm correctness vs GMP */
    {
        int ok = 1;
        /* Test case: 2^13 mod 97 = 44 */
        u128 mod97 = 97;
        MontCtx mctx;
        mont_ctx_init(&mctx, mod97);
        u128 result = mont_powm(2, 13, &mctx);
        if (result != 44) { printf("  FAIL OPT-15: 2^13 mod 97 = %llu (expected 44)\n", (unsigned long long)result); ok = 0; }

        /* Test: 3^50 mod 1000000007 */
        u128 big_mod = 1000000007ULL;
        mont_ctx_init(&mctx, big_mod);
        result = mont_powm(3, 50, &mctx);
        mpz_t gmp_r, gmp_b, gmp_e, gmp_m;
        mpz_inits(gmp_r, gmp_b, gmp_e, gmp_m, NULL);
        mpz_set_ui(gmp_b, 3); mpz_set_ui(gmp_e, 50); mpz_set_ui(gmp_m, 1000000007ULL);
        mpz_powm(gmp_r, gmp_b, gmp_e, gmp_m);
        u64 gmp_val = mpz_get_ui(gmp_r);
        if ((u64)result != gmp_val) {
            printf("  FAIL OPT-15: 3^50 mod 1e9+7: native=%llu gmp=%llu\n",
                   (unsigned long long)result, (unsigned long long)gmp_val);
            ok = 0;
        }

        /* Test with 90-bit number: CC13 root (hex 3d3aa382c918303b5470bad) */
        mpz_set_str(n, "3d3aa382c918303b5470bad", 16);
        u128 n89 = mpz_to_u128(n);
        mont_ctx_init(&mctx, n89);
        /* 2^(n-1) mod n should be 1 if n is prime */
        u128 nm1 = n89 - 1;
        result = mont_powm(2, nm1, &mctx);
        if (result != 1) {
            printf("  FAIL OPT-15: Fermat test on CC13 root gave %llu (expected 1)\n",
                   (unsigned long long)result);
            ok = 0;
        }

        /* Test Pocklington on first link: 2*root+1 */
        u128 link1 = (n89 << 1) + 1;
        mont_ctx_init(&mctx, link1);
        int pock = is_prime_pocklington_native(link1, &mctx);
        if (pock != 1) {
            printf("  FAIL OPT-15: Pocklington on CC13 link 1 returned %d (expected 1)\n", pock);
            ok = 0;
        }

        /* Test composite detection */
        u128 composite = 15;
        mont_ctx_init(&mctx, composite);
        result = mont_powm(2, 14, &mctx);
        if (result == 1) {
            printf("  FAIL OPT-15: Fermat test should fail for n=15\n");
            ok = 0;
        }

        if (ok) { printf("  OK OPT-15 Montgomery powm matches GMP (small, 1e9+7, 90-bit, composite)\n"); passed++; }
        else failed++;
        mpz_clears(gmp_r, gmp_b, gmp_e, gmp_m, NULL);
    }

    /* Test 16: OPT-15 Miller-Rabin native matches GMP for known primes/composites */
    {
        int ok = 1;
        /* Test known 64-bit primes */
        u64 primes64[] = {1000000007ULL, 998244353ULL, 999999999989ULL, 18446744073709551557ULL};
        for (int i = 0; i < 4; i++) {
            u128 p = primes64[i];
            MontCtx mctx;
            mont_ctx_init(&mctx, p);
            if (!is_prime_mr_native(p, &mctx)) {
                printf("  FAIL OPT-15 MR: %llu reported composite\n", (unsigned long long)primes64[i]);
                ok = 0;
            }
        }
        /* Test known composites: Carmichael numbers */
        u64 composites[] = {561, 1105, 1729, 2821};
        for (int i = 0; i < 4; i++) {
            u128 c = composites[i];
            MontCtx mctx;
            mont_ctx_init(&mctx, c);
            if (is_prime_mr_native(c, &mctx)) {
                printf("  FAIL OPT-15 MR: Carmichael %llu reported prime\n", (unsigned long long)composites[i]);
                ok = 0;
            }
        }
        if (ok) { printf("  OK OPT-15 native MR correct on 4 primes + 4 Carmichael numbers\n"); passed++; }
        else failed++;
    }

    /* Test 18: v33 GMP wide-mode chain follow matches u128 path */
    {
        mpz_set_str(n, "3d3aa382c918303b5470bad", 16);
        /* Both paths use 2^127 max_val (fits u128) for fair comparison */
        int gmp_len = follow_chain_first_gmp(n, max_val, current, temp, temp2, prev_prime);
        int native_len = follow_chain_first(n, max_val, current, temp, temp2, prev_prime);
        if (gmp_len == 13 && gmp_len == native_len) {
            printf("  OK v33 GMP chain follow CC13 = %d (matches native %d)\n", gmp_len, native_len);
            passed++;
        } else {
            printf("  FAIL v33 GMP chain follow = %d, native = %d (expected 13)\n", gmp_len, native_len);
            failed++;
        }
    }

    /* Test 19: v33 GMP wide-mode — known 135-bit CC15a root
     * Root: 27353790674175627273118204975428644651729
     * Hex:  0x5062b474874cd337fe0f061aaef7c442d1
     * Chain spans 135-149 bits, well beyond the u128 ceiling. */
    {
        int ok = 1;
        mpz_t cc15_root, big_max;
        mpz_init(cc15_root);
        mpz_init(big_max);

        mpz_set_str(cc15_root, "5062b474874cd337fe0f061aaef7c442d1", 16);
        mpz_ui_pow_ui(big_max, 2, 200);

        /* Verify it's a 135-bit prime */
        if (mpz_sizeinbase(cc15_root, 2) != 135) {
            printf("  FAIL v33 CC15: root is %zu bits (expected 135)\n",
                   mpz_sizeinbase(cc15_root, 2));
            ok = 0;
        }
        if (mpz_probab_prime_p(cc15_root, 25) <= 0) {
            printf("  FAIL v33 CC15: root is not prime\n");
            ok = 0;
        }

        /* Verify it's a first-kind root ((n-1)/2 is NOT prime) */
        if (ok) {
            int is_root = is_chain_root_first(cc15_root, pred, temp);
            if (!is_root) {
                printf("  FAIL v33 CC15: root is not a first-kind chain root\n");
                ok = 0;
            }
        }

        /* Follow chain with GMP wide-mode, verify length == 15 */
        if (ok) {
            int cc15_len = follow_chain_first_gmp(cc15_root, big_max, current,
                                                   temp, temp2, prev_prime);
            if (cc15_len != 15) {
                printf("  FAIL v33 CC15: GMP chain follow returned %d (expected 15)\n", cc15_len);
                ok = 0;
            }
        }

        /* Also verify screened version matches */
        if (ok) {
            u32 pf_res15[PREFILTER_PRIMES_COUNT];
            for (int i = 0; i < PREFILTER_PRIMES_COUNT; i++)
                pf_res15[i] = (u32)mpz_fdiv_ui(cc15_root, PREFILTER_PRIMES[i]);
            u64 ss15 = 0, lt15 = 0;
            int scr_len = follow_chain_first_screened_gmp(cc15_root, big_max, current,
                            temp, temp2, prev_prime, pf_res15, 1, &ss15, &lt15);
            if (scr_len != 15) {
                printf("  FAIL v33 CC15: screened GMP chain = %d (expected 15)\n", scr_len);
                ok = 0;
            }
        }

        if (ok) {
            printf("  OK v33 CC15a 135-bit root: GMP chain follow = 15 (screened matches)\n");
            passed++;
        } else {
            failed++;
        }

        mpz_clear(cc15_root);
        mpz_clear(big_max);
    }

    printf("\n--- %d passed, %d failed ---\n", passed, failed);

    mpz_clear(n); mpz_clear(temp); mpz_clear(temp2); mpz_clear(pred);
    mpz_clear(prev_prime); mpz_clear(current); mpz_clear(max_val);

    return (failed == 0) ? 0 : 1;
}

/* =============================================================================
 * MAIN
 * ========================================================================== */

int main(int argc, char** argv) {
    int target_length = 18;
    int target_bits = 89;
    int num_threads = 0;
    int log_threshold = 6;
    int continuous_mode = 0;
    const char* prefix_str = NULL;
    const char* output_file = NULL;
    int resume_mode = 0;

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--help") || !strcmp(argv[i], "-h")) { print_usage(argv[0]); return 0; }
        else if (!strcmp(argv[i], "--target") && i+1 < argc) target_length = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--bits") && i+1 < argc) target_bits = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--threads") && i+1 < argc) num_threads = atoi(argv[++i]);
        else if ((!strcmp(argv[i], "--log") || !strcmp(argv[i], "-l")) && i+1 < argc)
            log_threshold = atoi(argv[++i]);
        else if ((!strcmp(argv[i], "--report") || !strcmp(argv[i], "-r")) && i+1 < argc)
            g_report_interval_sec = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--continuous") || !strcmp(argv[i], "-c")) continuous_mode = 1;
        else if ((!strcmp(argv[i], "--prefix") || !strcmp(argv[i], "-p")) && i+1 < argc)
            prefix_str = argv[++i];
        else if ((!strcmp(argv[i], "--sequential") || !strcmp(argv[i], "-s")))
            g_sequential_mode = 1;
        else if ((!strcmp(argv[i], "--output") || !strcmp(argv[i], "-o")) && i+1 < argc)
            output_file = argv[++i];
        else if (!strcmp(argv[i], "--quiet") || !strcmp(argv[i], "-q")) g_quiet_mode = 1;
        else if (!strcmp(argv[i], "--full-quiet")) { g_full_quiet_mode = 1; g_quiet_mode = 1; }
        else if (!strcmp(argv[i], "--pin")) g_pin_threads = 1;
        else if (!strcmp(argv[i], "--pin-base") && i+1 < argc) g_pin_base_cpu = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--screen-primes") && i+1 < argc) {
            g_active_screen_count = atoi(argv[++i]);
            if (g_active_screen_count < 1) g_active_screen_count = 1;
            if (g_active_screen_count > CHAIN_SCREEN_COUNT) g_active_screen_count = CHAIN_SCREEN_COUNT;
        }
        else if (!strcmp(argv[i], "--checkpoint") && i+1 < argc) g_checkpoint_file = argv[++i];
        else if (!strcmp(argv[i], "--resume")) resume_mode = 1;
        else if (!strcmp(argv[i], "--ckpt-interval") && i+1 < argc)
            g_checkpoint_interval_sec = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--chunk-tiles") && i+1 < argc)
            g_chunk_tiles = (u64)atoll(argv[++i]);
        else if (!strcmp(argv[i], "--line-depth") && i+1 < argc)
            g_line_depth = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--sieve-only")) g_sieve_only_mode = 1;
        else if (!strcmp(argv[i], "--test")) { return run_unit_tests(); }
        else if (argv[i][0] == '-') {
            if (!strcmp(argv[i], "-t") && i+1 < argc) num_threads = atoi(argv[++i]);
            else if (!strcmp(argv[i], "-b") && i+1 < argc) target_bits = atoi(argv[++i]);
            else {
                fprintf(stderr, "Unknown option: %s (try --help)\n", argv[i]);
                return 1;
            }
        }
    }

    /* FIX-19: Validate critical CLI parameters to prevent undefined behavior.
     * --bits 0 or negative → mpz_ui_pow_ui with huge unsigned exponent → OOM.
     * --target 0 or negative → invalid loop bounds and wheel generation. */
    if (target_bits < 2) {
        fprintf(stderr, "ERROR: --bits must be >= 2 (got %d)\n", target_bits);
        return 1;
    }
    if (target_length < 1) {
        fprintf(stderr, "ERROR: --target must be >= 1 (got %d)\n", target_length);
        return 1;
    }

    if (g_report_interval_sec <= 0) g_report_interval_sec = 1;
    if (g_pin_base_cpu < 0) g_pin_base_cpu = 0;

    if (num_threads <= 0) {
        num_threads = sysconf(_SC_NPROCESSORS_ONLN);
        if (num_threads <= 0) num_threads = 4;
    }

    if (g_sequential_mode && !prefix_str) {
        printf("ERROR: --sequential requires --prefix\n");
        return 1;
    }

    /* OPT-16: Auto-select random-chunk mode when not sequential */
    if (!g_sequential_mode) {
        g_random_chunk_mode = 1;
        if (g_chunk_tiles == 0) g_chunk_tiles = 500;
    }

    signal(SIGINT, signal_handler);
    signal(SIGTERM, signal_handler);

    g_log_fp = NULL;
    if (output_file) {
        g_log_fp = fopen(output_file, "a");
        if (g_log_fp) {
            setvbuf(g_log_fp, NULL, _IOFBF, 1 << 20);
        }
    }
    if (g_full_quiet_mode && !g_log_fp) {
        fprintf(stderr, "WARNING: --full-quiet without --output; results won't be persisted\n");
    }

    if (!g_quiet_mode) {
        printf("╔═══════════════════════════════════════════════════════════════════╗\n");
        printf("║  CUNNINGHAM CHAIN CONSTRUCTOR - GMP %-29s ║\n", PROFILE_BANNER);
        printf("║  CC18a-ONLY | OPT-19..22 | GMP wide-mode for >127-bit chains    ║\n");
        printf("╚═══════════════════════════════════════════════════════════════════╝\n\n");
    }

    mpz_init(g_lattice_m); mpz_init(g_k_range_min); mpz_init(g_k_range_size);
    mpz_init(g_k_range_max); mpz_init(g_current_tile);
    mpz_init(g_n_range_min); mpz_init(g_n_range_max);
    mpz_set_ui(g_lattice_m, LATTICE_M);

    mpz_t prefix;
    int prefix_bits = 0;
    mpz_init(prefix);
    if (prefix_str) {
        if (parse_prefix(prefix_str, prefix, &prefix_bits) != 0) {
            mpz_clear(prefix); return 1;
        }
        g_use_prefix = 1;
        QPRINTF("PREFIX: %s (%d bits)\n", prefix_str, prefix_bits);
    }

    const char* kind_str = "FIRST-KIND ONLY (CC18a)";

    QPRINTF("\n=== LEVEL 0+1: CRT+WHEEL (%s) ===\n", kind_str);
    QPRINTF("%s (M = %llu)\n", PROFILE_CRT_DESC, (unsigned long long)LATTICE_M);
    QPRINTF("Wheel primes: 23, 29, 31 (W = %llu)\n", (unsigned long long)WHEEL_PERIOD);

    generate_optimized_wheels(target_length);

    QPRINTF("Valid bases: %d (first-kind)\n", g_num_bases);

    if (g_num_bases == 0) {
        fprintf(stderr, "ERROR: No valid CRT bases found for CC%d\n", target_length);
        return 1;
    }

    QPRINTF("Total wheel positions: %d\n", g_total_wheel_positions);
    QPRINTF("Average density: %.2f%%\n", g_avg_wheel_density);

    QPRINTF("\n=== LEVEL 2: FILTER (%s) ===\n", kind_str);
    QPRINTF("%s\n", PROFILE_L2_DESC);
    if (g_use_ext_l2)
        QPRINTF("Extended (OPT-7): 67, 71, 73, 79, 83, 89, 97\n");
    else
        QPRINTF("Extended (OPT-7): DISABLED\n");
    compute_filter_bitmasks(target_length > MAX_SIEVE_CHAIN_LEN ? MAX_SIEVE_CHAIN_LEN : target_length);

    QPRINTF("\n=== OPT-6: CHAIN PRE-FILTER ===\n");
    QPRINTF("Primes: 67-137 (%d primes), kill-position tables\n", PREFILTER_PRIMES_COUNT);
    init_prefilter_tables();
    init_trial_batches();
    QPRINTF("Pre-filter tables computed\n");

    /* OPT-17: Initialize line-sieve if enabled */
    if (g_line_depth > 0) {
        if (g_line_depth > target_length - 1)
            g_line_depth = target_length - 1;
        init_line_sieve(g_line_depth);
    }

    QPRINTF("\n=== LEVEL 3: PRIMALITY ===\n");
    QPRINTF("Trial division: 67-409 (%d primes, %d/%d batches full/extl2)\n",
            TRIAL_PRIMES_COUNT, g_num_batches_full, g_num_batches_extl2);
    QPRINTF("Root test: BPSW + 1 MR (OPT-2)\n");
    QPRINTF("Chain links: Pocklington proving (OPT-1, first-kind only)\n");
    QPRINTF("Chain screening (OPT-13): %d/%d primes (67-%u), div-free recurrence (OPT-14)\n",
            g_active_screen_count, CHAIN_SCREEN_COUNT,
            CHAIN_SCREEN_PRIMES[g_active_screen_count - 1]);

    /* DOC-1: Inform about expected prefilter/screening behavior */
    {
        int sieve_len_used = (target_length <= MAX_SIEVE_CHAIN_LEN) ? target_length : MAX_SIEVE_CHAIN_LEN;
        int min_cleared = g_line_depth;  /* line-sieve clears pos 1..line_depth for primes 101-863 */
        /* ext-L2 clears pos 0..sieve_len-1 for primes 67-97 */
        int ext_cleared = g_use_ext_l2 ? (sieve_len_used - 1) : 0;
        QPRINTF("NOTE: Ext-L2 clears pos 0..%d (primes 67-97), line-sieve clears pos 1..%d (primes 101-863)\n",
                ext_cleared, g_line_depth);
        if (log_threshold <= min_cleared + 1 && log_threshold <= ext_cleared + 1) {
            QPRINTF("  -> Prefilter (cutoff=%d) guaranteed pass: PF=0 is expected (positions already cleared)\n",
                    log_threshold);
        }
        if (min_cleared >= 1) {
            QPRINTF("  -> Screen catches only at position %d+ (chains CC%d+ only)\n",
                    min_cleared + 1, min_cleared + 1);
        }
    }

    /* Compute k ranges */
    mpz_t two_pow_bits, two_pow_min;
    mpz_init(two_pow_bits); mpz_init(two_pow_min);
    mpz_ui_pow_ui(two_pow_bits, 2, target_bits);
    mpz_fdiv_q(g_k_range_max, two_pow_bits, g_lattice_m);
    mpz_ui_pow_ui(two_pow_min, 2, target_bits - 1);
    mpz_fdiv_q(g_k_range_min, two_pow_min, g_lattice_m);
    mpz_sub(g_k_range_size, g_k_range_max, g_k_range_min);
    mpz_add_ui(g_k_range_size, g_k_range_size, 1);

    mpz_set(g_n_range_min, two_pow_min);
    mpz_sub_ui(g_n_range_max, two_pow_bits, 1);

    if (g_use_prefix) {
        mpz_t cand_min, cand_max, min_valid, max_valid;
        mpz_init(cand_min); mpz_init(cand_max); mpz_init(min_valid); mpz_init(max_valid);

        /* FIX-12: Guard against prefix_bits > target_bits (negative shift → OOM) */
        if (prefix_bits > target_bits) {
            fprintf(stderr, "\n*** ERROR: Prefix has %d bits but target is only %d bits ***\n",
                    prefix_bits, target_bits);
            fprintf(stderr, "    Prefix must be shorter than --bits value.\n");
            mpz_clear(cand_min); mpz_clear(cand_max); mpz_clear(min_valid); mpz_clear(max_valid);
            return 1;
        }
        int shift = target_bits - prefix_bits;
        mpz_mul_2exp(cand_min, prefix, shift);
        mpz_add_ui(cand_max, prefix, 1);
        mpz_mul_2exp(cand_max, cand_max, shift);
        mpz_sub_ui(cand_max, cand_max, 1);

        mpz_ui_pow_ui(min_valid, 2, target_bits - 1);
        mpz_ui_pow_ui(max_valid, 2, target_bits);
        mpz_sub_ui(max_valid, max_valid, 1);

        if (mpz_cmp(cand_max, min_valid) < 0 || mpz_cmp(cand_min, max_valid) > 0) {
            fprintf(stderr, "\n*** ERROR: Prefix does not overlap %d-bit range ***\n", target_bits);
            mpz_clear(cand_min); mpz_clear(cand_max); mpz_clear(min_valid); mpz_clear(max_valid);
            return 1;
        }

        if (mpz_cmp(cand_min, min_valid) < 0) mpz_set(cand_min, min_valid);
        if (mpz_cmp(cand_max, max_valid) > 0) mpz_set(cand_max, max_valid);

        mpz_set(g_n_range_min, cand_min);
        mpz_set(g_n_range_max, cand_max);

        u64 largest_base = g_base_wheels[g_num_bases - 1].base;
        u64 smallest_base = g_base_wheels[0].base;

        if (mpz_cmp_ui(cand_min, largest_base) > 0) {
            mpz_sub_ui(g_k_range_min, cand_min, largest_base);
            mpz_cdiv_q(g_k_range_min, g_k_range_min, g_lattice_m);
        } else {
            mpz_set_ui(g_k_range_min, 0);
        }

        mpz_sub_ui(g_k_range_max, cand_max, smallest_base);
        mpz_fdiv_q(g_k_range_max, g_k_range_max, g_lattice_m);

        mpz_sub(g_k_range_size, g_k_range_max, g_k_range_min);
        mpz_add_ui(g_k_range_size, g_k_range_size, 1);

        mpz_clear(cand_min); mpz_clear(cand_max); mpz_clear(min_valid); mpz_clear(max_valid);
    }

    mpz_fdiv_q_ui(g_current_tile, g_k_range_min, WHEEL_PERIOD);
    mpz_mul_ui(g_current_tile, g_current_tile, WHEEL_PERIOD);

    /* Resume from checkpoint if requested */
    if (resume_mode && g_checkpoint_file) {
        load_checkpoint();
    }

    QPRINTF("\n=== SEARCH CONFIGURATION ===\n");
    QPRINTF("Target:         CC%d+\n", target_length);
    QPRINTF("Target bits:    %d\n", target_bits);
    QPRINTF("Threads:        %d\n", num_threads);
    QPRINTF("Mode:           %s\n", g_sequential_mode ? "SEQUENTIAL" : "RANDOM");
    QPRINTF("Pinning:        %s\n", g_pin_threads ? "ON" : "OFF");
    if (g_pin_threads) QPRINTF("Pin base CPU:   %d\n", g_pin_base_cpu);
    QPRINTF("Report every:   %d s\n", g_report_interval_sec);
    QPRINTF("Chain kind:     FIRST-KIND ONLY (CC18a)\n");
    QPRINTF("Ext L2:         ON (67-97)\n");
    QPRINTF("Screen primes:  %d (67-%u)\n", g_active_screen_count,
            CHAIN_SCREEN_PRIMES[g_active_screen_count - 1]);
    if (g_line_depth > 0)
        QPRINTF("Line-sieve:     ON (pos 1..%d, %d primes 101-863)\n", g_line_depth, LINE_SIEVE_COUNT);
    else
        QPRINTF("Line-sieve:     OFF (enable with --line-depth N)\n");
    if (g_sieve_only_mode)
        QPRINTF("Sieve-only:     ON (no primality proving)\n");
    if (g_checkpoint_file)
        QPRINTF("Checkpoint:     %s (every %ds)\n", g_checkpoint_file, g_checkpoint_interval_sec);
    QPRINTF("Bases:          %d (first-kind)\n", g_num_bases);

    int k_bits = mpz_sizeinbase(g_k_range_size, 2);
    QPRINTF("K range bits:   %d\n", k_bits);

    /* v33: Wide-mode selection. If chain members exceed 127 bits, use GMP
     * for the proving path instead of native u128 Montgomery arithmetic. */
    g_target_bits_global = target_bits;
    g_target_length_global = target_length;
    {
        int max_chain_bits = target_bits + target_length;
        if (max_chain_bits > 127) {
            g_wide_mode = 1;
            QPRINTF("Wide mode:      ENABLED (GMP path for %d-bit chain members)\n", max_chain_bits);
            QPRINTF("                L2/ext-L2/line-sieve screening unchanged (u32 residues)\n");
            QPRINTF("                Primality: GMP BPSW | Chain follow: GMP arithmetic\n");
        } else {
            QPRINTF("Chain safety:   %d-bit root + CC%d = max %d bits (< 127, native u128 OK)\n",
                    target_bits, target_length, max_chain_bits);
        }
    }

    memset(&g_results, 0, sizeof(g_results));
    pthread_mutex_init(&g_results.lock, NULL);

    pthread_t* threads = (pthread_t*)malloc(num_threads * sizeof(pthread_t));
    ThreadConfig* configs = (ThreadConfig*)malloc(num_threads * sizeof(ThreadConfig));

    const char* mode_str = g_sequential_mode ? "SEQUENTIAL" : "RANDOM-CHUNK";
    QPRINTF("\n=== STARTING SEARCH (%s: +OPT-13..18, mode=%s, line=%d) ===\n", PROFILE_BANNER, mode_str, g_line_depth);
    if (g_random_chunk_mode)
        QPRINTF("  Chunk size: %llu tiles (~%.1f candidates/chunk)\n",
                (unsigned long long)g_chunk_tiles,
                (double)g_chunk_tiles * g_num_bases * 5.0);
    QPRINTF("\n");

    struct timeval start_time, now;
    gettimeofday(&start_time, NULL);
    g_last_checkpoint_time = start_time.tv_sec + start_time.tv_usec / 1e6;

    for (int i = 0; i < num_threads; i++) {
        u64 batch = g_random_chunk_mode ? g_chunk_tiles :
                    g_sequential_mode ? 500 : 500;
        init_thread_config(&configs[i], i, target_length, target_bits, log_threshold, batch);
        pthread_create(&threads[i], NULL, worker_sequential, &configs[i]);
    }

    while (!shutdown_requested && !g_search_complete) {
        sleep((unsigned int)g_report_interval_sec);
        gettimeofday(&now, NULL);
        double elapsed = (now.tv_sec - start_time.tv_sec) + (now.tv_usec - start_time.tv_usec) / 1e6;
        double now_d = now.tv_sec + now.tv_usec / 1e6;

        /* Periodic checkpoint */
        if (g_checkpoint_file && g_sequential_mode &&
            (now_d - g_last_checkpoint_time) >= g_checkpoint_interval_sec) {
            write_checkpoint();
            g_last_checkpoint_time = now_d;
        }

        pthread_mutex_lock(&g_results.lock);
        u64 wheel = g_results.wheel_candidates;
        u64 primes = g_results.primes_found;
        u64 chains = g_results.chains_found;
        u64 prefilter_skip = g_results.prefilter_skipped;
        u64 screen_saved = g_results.chain_screen_saved;
        u64 links_tested = g_results.chain_links_tested;
        int best = g_results.best_length;
        pthread_mutex_unlock(&g_results.lock);

        double wheel_rate = wheel / (elapsed > 0 ? elapsed : 1) / 1e6;
        double primes_per_sec = primes / (elapsed > 0 ? elapsed : 1);
        double screen_pct = (links_tested > 0) ? 100.0 * screen_saved / links_tested : 0;

        if (!g_quiet_mode) {
            if (g_sequential_mode && !g_random_chunk_mode) {
                pthread_mutex_lock(&g_seq_lock);
                mpz_t progress; mpz_init(progress);
                mpz_sub(progress, g_current_tile, g_k_range_min);
                double pct = 100.0 * mpz_get_d(progress) / mpz_get_d(g_k_range_size);
                mpz_clear(progress);
                pthread_mutex_unlock(&g_seq_lock);

                /* ETA calculation */
                double eta_sec = (pct > 0.0001) ? (elapsed / pct * (100.0 - pct)) : 0;
                char eta_buf[64] = "";
                if (pct > 0.001 && eta_sec > 0 && eta_sec < 1e9) {
                    if (eta_sec < 3600)
                        snprintf(eta_buf, sizeof(eta_buf), " ETA:%.0fm", eta_sec / 60);
                    else if (eta_sec < 86400)
                        snprintf(eta_buf, sizeof(eta_buf), " ETA:%.1fh", eta_sec / 3600);
                    else
                        snprintf(eta_buf, sizeof(eta_buf), " ETA:%.1fd", eta_sec / 86400);
                }

                printf("\r%.4f%% | W: %llu (%.2fM/s) | P: %llu (%.0f/s) | PF: %llu | Scr: %llu/%.0f%% | CC%d+: %llu | Best: CC%d | %.0fs%s   ",
                       pct,
                       (unsigned long long)wheel, wheel_rate,
                       (unsigned long long)primes, primes_per_sec,
                       (unsigned long long)prefilter_skip,
                       (unsigned long long)screen_saved, screen_pct,
                       target_length, (unsigned long long)chains,
                       best, elapsed, eta_buf);
            } else if (g_random_chunk_mode) {
                u64 tiles = 0;
                pthread_mutex_lock(&g_results.lock);
                tiles = g_results.tiles_processed;
                pthread_mutex_unlock(&g_results.lock);
                printf("\r[RChunk] W: %llu (%.2fM/s) | P: %llu (%.0f/s) | T: %llu | PF: %llu | Scr: %llu/%.0f%% | CC%d+: %llu | Best: CC%d | %.0fs   ",
                       (unsigned long long)wheel, wheel_rate,
                       (unsigned long long)primes, primes_per_sec,
                       (unsigned long long)tiles,
                       (unsigned long long)prefilter_skip,
                       (unsigned long long)screen_saved, screen_pct,
                       target_length, (unsigned long long)chains,
                       best, elapsed);
            } else {
                printf("\rW: %llu (%.2fM/s) | P: %llu (%.0f/s) | PF: %llu | Scr: %llu/%.0f%% | CC%d+: %llu | Best: CC%d | %.0fs   ",
                       (unsigned long long)wheel, wheel_rate,
                       (unsigned long long)primes, primes_per_sec,
                       (unsigned long long)prefilter_skip,
                       (unsigned long long)screen_saved, screen_pct,
                       target_length, (unsigned long long)chains,
                       best, elapsed);
            }
            fflush(stdout);
        }

        if (!continuous_mode && !g_sequential_mode && elapsed > 10) break;
    }

    /* Final checkpoint before exit */
    if (g_checkpoint_file && g_sequential_mode) write_checkpoint();

    if (!g_search_complete) shutdown_requested = 1;

    for (int i = 0; i < num_threads; i++) {
        pthread_join(threads[i], NULL);
        cleanup_thread_config(&configs[i]);
    }

    QPRINTF("\n\n");
    if (g_search_complete) QPRINTF("*** SEARCH COMPLETE ***\n\n");
    else QPRINTF("*** STOPPED ***\n\n");

    gettimeofday(&now, NULL);
    double total_time = (now.tv_sec - start_time.tv_sec) + (now.tv_usec - start_time.tv_usec) / 1e6;

    QPRINTF("=== FINAL RESULTS ===\n");
    QPRINTF("Wheel candidates (L1): %llu\n", (unsigned long long)g_results.wheel_candidates);
    QPRINTF("Passed L2 bitmask:     %llu (%.2f%%)\n",
           (unsigned long long)g_results.passed_l2_filter,
           100.0 * g_results.passed_l2_filter / (g_results.wheel_candidates + 1));
    QPRINTF("Passed ext L2 (OPT-7): %llu (%.2f%% of L2)\n",
           (unsigned long long)g_results.passed_ext_filter,
           100.0 * g_results.passed_ext_filter / (g_results.passed_l2_filter + 1));
    if (g_line_depth > 0) {
        u64 line_total = g_results.line_sieve_passed + g_results.line_sieve_rejected;
        QPRINTF("Line-sieve (OPT-17):   %llu passed, %llu rejected (%.2f%% kill, pos 1..%d, %d primes)\n",
               (unsigned long long)g_results.line_sieve_passed,
               (unsigned long long)g_results.line_sieve_rejected,
               100.0 * g_results.line_sieve_rejected / (line_total + 1),
               g_line_depth, LINE_SIEVE_COUNT);
    }
    QPRINTF("Primes found:          %llu\n", (unsigned long long)g_results.primes_found);
    QPRINTF("Prefilter skipped:     %llu (%.1f%% of roots)\n",
           (unsigned long long)g_results.prefilter_skipped,
           100.0 * g_results.prefilter_skipped /
           (g_results.primes_found - g_results.non_roots_skipped + 1));
    QPRINTF("Non-roots skipped:     %llu\n", (unsigned long long)g_results.non_roots_skipped);
    QPRINTF("Chain links tested:    %llu\n", (unsigned long long)g_results.chain_links_tested);
    QPRINTF("Screen caught (OPT-13):%llu (%.1f%% of links, %d primes 67-%u, div-free)\n",
           (unsigned long long)g_results.chain_screen_saved,
           100.0 * g_results.chain_screen_saved / (g_results.chain_links_tested + 1),
           g_active_screen_count, CHAIN_SCREEN_PRIMES[g_active_screen_count - 1]);
    QPRINTF("Chains >= CC%d:        %llu\n", target_length, (unsigned long long)g_results.chains_found);

    QPRINTF("\nChain distribution (first-kind only):\n");
    for (int i = 1; i <= MAX_CHAIN; i++) {
        if (g_results.chains_by_length[i] > 0) {
            QPRINTF("  CC%2da: %llu%s\n", i,
                   (unsigned long long)g_results.chains_by_length[i],
                   (i >= target_length) ? " *** TARGET" : (i >= log_threshold ? " *" : ""));
        }
    }

    if (g_results.num_results > 0) {
        QPRINTF("\n=== TOP RESULTS ===\n");
        for (int i = 0; i < g_results.num_results && i < 20; i++) {
            QPRINTF("  CC%2da: 0x%s\n", g_results.result_lengths[i],
                   g_results.result_strings[i]);
        }
    }

    QPRINTF("\nTime: %.1f seconds\n", total_time);
    QPRINTF("Wheel rate: %.2f M/sec\n", g_results.wheel_candidates / total_time / 1e6);
    QPRINTF("Primes/sec: %.0f\n", g_results.primes_found / total_time);

    if (g_log_fp) {
        time_t now_t = time(NULL);
        fprintf(g_log_fp, "\n=== CC%d SEARCH (%s) %s", target_length, PROFILE_BANNER, ctime(&now_t));
        fprintf(g_log_fp, "Bits: %d, Threads: %d, Time: %.1fs, Screen: %d primes\n",
                target_bits, num_threads, total_time, g_active_screen_count);
        fprintf(g_log_fp, "Wheel: %llu, L2: %llu, ExtL2: %llu, Primes: %llu, PF_skip: %llu, Screen: %llu/%llu\n",
                (unsigned long long)g_results.wheel_candidates,
                (unsigned long long)g_results.passed_l2_filter,
                (unsigned long long)g_results.passed_ext_filter,
                (unsigned long long)g_results.primes_found,
                (unsigned long long)g_results.prefilter_skipped,
                (unsigned long long)g_results.chain_screen_saved,
                (unsigned long long)g_results.chain_links_tested);
        for (int i = 0; i < g_results.num_results; i++) {
            fprintf(g_log_fp, "CC%da 0x%s\n", g_results.result_lengths[i],
                   g_results.result_strings[i]);
        }
        fflush(g_log_fp);  /* FIX-9: flush before close to ensure all data persisted */
        fclose(g_log_fp);
        g_log_fp = NULL;
        QPRINTF("Results saved to: %s\n", output_file);
    }

    /* Cleanup */
    pthread_mutex_destroy(&g_results.lock);
    pthread_mutex_destroy(&g_print_lock);
    pthread_mutex_destroy(&g_log_file_lock);
    free(threads); free(configs);
    for (int b = 0; b < g_num_bases; b++) {
        for (int m = 0; m < 6; m++) free(g_base_wheels[b].wheel_by_mod6[m]);
        free(g_base_wheels[b].residues);
        free(g_base_wheels[b].ext_residues);
        free(g_base_wheels[b].line_residues);
    }
    free(g_base_wheels);
    /* FIX-21: Free dynamic result strings */
    for (int i = 0; i < g_results.num_results; i++)
        free(g_results.result_strings[i]);
    free(g_results.best_start);
    mpz_clear(g_lattice_m); mpz_clear(g_k_range_min); mpz_clear(g_k_range_size);
    mpz_clear(g_k_range_max); mpz_clear(g_current_tile);
    mpz_clear(g_n_range_min); mpz_clear(g_n_range_max);
    mpz_clear(prefix); mpz_clear(two_pow_bits); mpz_clear(two_pow_min);

    return 0;
}
