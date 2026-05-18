/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2026 Nenad Micic <nenad@micic.be>
 *
 * cc20_first_kind_immune_v16.cu
 *
 * First-kind Cunningham chain search over multiplicative immune families.
 *
 *   p     = A * D * 2^E - 1
 *   p_j   = A * D * 2^(E + j) - 1     for j = 0..target_len - 1
 *
 * D is the immune-factor product (every q | D forces p_j ≡ -1 (mod q) for
 * all j). The search variable is A; we build a CRT-join coefficient wheel
 * over small primes coprime to D, then a packed line sieve over a wider
 * band of primes also coprime to D.
 *
 * Forbidden-residue math (for any sieve prime q with D % q != 0):
 *   A * D * 2^(E+j) - 1 ≡ 0 (mod q)
 *   A ≡ inv(D mod q) * inv(2^(E+j) mod q)  (mod q)   for some j ∈ [0..target)
 * Family is inadmissible at q iff |F_q| == q (dedup count).
 *
 * Primes q | D are immune: no sieve, no forbidden set.
 *
 * Reverse-depth proving: test top-position chain term first
 *   top = p * 2^(target-1) + (2^(target-1) - 1)
 * then forward chain-follow only on candidates whose top survives MR.
 * (Same idea as v14/v15.)
 *
 * Build: see Makefile (sm_120 default for RTX 5090).
 */

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <signal.h>
#include <unistd.h>
#include <getopt.h>
#include <sys/time.h>
#include <fcntl.h>
#include <errno.h>
#include <gmp.h>
#include <pthread.h>

typedef uint64_t u64;
typedef uint32_t u32;
typedef uint16_t u16;
typedef uint8_t  u8;
typedef unsigned __int128 u128;

#define CUDA_CHECK(call) do {                                              \
    cudaError_t err = (call);                                              \
    if (err != cudaSuccess) {                                              \
        fprintf(stderr, "CUDA error at %s:%d: %s\n",                       \
                __FILE__, __LINE__, cudaGetErrorString(err));              \
        exit(1);                                                           \
    }                                                                      \
} while (0)

#define V16_VERSION       "v16.0.0-first-kind-immune"
#define V16_VERSION_Q_ITER "v16.0.3-drained-cursor"
#define V16_MAX_TARGET    24
#define V16_MAX_IMMUNE    64
#define V16_MAX_WHEEL     16     /* wheel-prime count (primorial must fit u64) */
#define V16_MAX_SIEVE     256
#define V16_MAX_WHEEL_OFFSETS (1u << 22)   /* 4M admissible offsets max */
#define V16_THREADS       256
#define V16_MAX_RESULTS   (1u << 16)
#define V16_MAX_SURVIVORS (1u << 25)        /* 32M survivors/batch; 256MB at 8B each */
#define V16_Q_MAX_SURVIVORS (1u << 22)      /* 4M survivors/batch for Q-iter (16 B each, 64 MiB/slot) */
#define V16_MAX_STREAMS   16                /* hard cap on N-deep dispatch stream pool */
#define V16_PACKED_WORDS_MAX 14  /* ceil(max sieve prime / 64). 14*64 = 896, matches v15 line band 101..863 */
#define V16_TRIAL_PRIMES_MAX 64
#define V16_Q_FORBID_SENTINEL 0xFFFFFFFFu  /* q | D_V → no forbidden residue */
/* Per-(V,q) forbidden-residue bitmask.
 * Max sieve prime is 503 -> ceil(503/64) = 8 u64 words covers everything.
 * Compile-time guard against future --sieve-max bumps below. */
#define V16_Q_MASK_WORDS  8
#define V16_HIST_MAX      64                /* chain-length histogram bins */

/* Survivor record: identifies a (tile_in_launch, variant, wheel_offset_idx)
 * triple that passed the K1 line sieve and is queued for K2 proving.
 * 8 bytes — tile_in_launch fits in u16 since main() clamps launch_tiles to
 * 65535 (CUDA gridDim.y limit); woff_idx stays u32 since n_admissible can
 * reach V16_MAX_WHEEL_OFFSETS = 1<<22. variant_idx is 0 when
 * --sieve-variants is OFF. */
typedef struct __align__(8) {
    u16 tile_in_launch;
    u16 variant_idx;
    u32 woff_idx;
} V16Survivor;

/* Q-iter survivor: one (Q, V) pair that passed the per-variant sieve.
 * 8 B record: split Q into u128 host base + u32 thread offset, matching v15
 * survivor layout. The full Q is reconstructed CPU-side as
 *   Q_full = Q_base_for_batch + (u128)thread_offset
 * where Q_base_for_batch is the per-launch u128 base recorded alongside the
 * survivor batch (see ProveWorker.q_base_lo/q_base_hi). Halves device memory
 * vs. the old 16 B u64-Q record and unlocks Q_max > 2^63. */
typedef struct __align__(8) {
    u32 thread_offset;  /* Q = Q_base + thread_offset (per-launch base passed
                         * separately via the prove-batch handoff) */
    u16 seed_idx;       /* index into cfg->pool_seeds[] */
    u16 _pad;
} QSurvivor;

/* K1.5 survivor: identical wire layout to QSurvivor so the
 * CPU prove worker can consume either type via a const-cast (we tag the mode
 * through ProveWorker.k15_passed). Once the host already discriminates via
 * the worker flag, we save the byte and keep the record 8 B = QSurvivor footprint.
 * Aliased through CPU dispatch. */
typedef struct __align__(8) {
    u32 thread_offset;
    u16 seed_idx;
    u16 _pad;
} K15Survivor;

/* Chain hit reported by the CPU prove pool (mirror of the device-side
 * d_chain_* arrays, used only when --cpu-prove is on). */
typedef struct {
    u64 A_lo;
    u64 A_hi;
    u32 chain_len;
    u32 top_pass;
} HostChainResult;

/* Pool seed: a named D variant we attribute hits to at prove time.
 * The sieve runs with the primary --immune-* config; pool members are
 * tested only via "does D_pool | (p+1)" on prime roots. */
/* Raised 1024 → 8192 for Q-iter: auto:71#:1-3 emits 2516 V, auto:71#:1-4 emits
 * 2516+1820=2516+C(16,4)=2516+1820. Stays in BSS (g_cfg is global); the
 * v16_cfg.pool_seeds array adds ~3.5 MiB to the binary, negligible. */
#define V16_MAX_POOL_SEEDS 8192
#define V16_MAX_AUTO_POOL  16
typedef struct {
    char name[64];
    u128 D;
    char D_dec[80];
    u32  primes[V16_MAX_IMMUNE];
    int  n_primes;
    u64  hits;             /* atomic accumulator (set via __atomic_add_fetch) */
} v16_pool_seed;

/* ============================================================================
 * Small-prime table (used to enumerate D candidates and wheel/sieve primes).
 * ============================================================================ */
static const u32 V16_SMALL_PRIMES[] = {
      2,   3,   5,   7,  11,  13,  17,  19,  23,  29,
     31,  37,  41,  43,  47,  53,  59,  61,  67,  71,
     73,  79,  83,  89,  97, 101, 103, 107, 109, 113,
    127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
    179, 181, 191, 193, 197, 199, 211, 223, 227, 229,
    233, 239, 241, 251, 257, 263, 269, 271, 277, 281,
    283, 293, 307, 311, 313, 317, 331, 337, 347, 349,
    353, 359, 367, 373, 379, 383, 389, 397, 401, 409,
    419, 421, 431, 433, 439, 443, 449, 457, 461, 463,
    467, 479, 487, 491, 499, 503
};
#define V16_SMALL_PRIME_COUNT (sizeof(V16_SMALL_PRIMES) / sizeof(V16_SMALL_PRIMES[0]))

/* 2026-05-14: extended-sieve primes 509..1009 for the CPU-side depth filter.
 * The GPU forbid_mask is capped at V16_Q_MASK_WORDS*64 = 512 (prime ≤ 503),
 * so primes > 503 can't go into the GPU kernel without doubling the mask
 * memory. Instead we precompute a per-(V, q) forbid bitmask on the host
 * (allocated via calloc — see build_ext_forbid_table; the "pinned" wording
 * in earlier revisions was aspirational, not literal) and check it at the
 * *top* of cpu_prove_worker before the expensive Fermat/MR. Math is
 * identical to the GPU sieve — Q must not match the forbid set for any
 * active extended prime. Mirrors v34 OPT-C's packed u64 kill table
 * (cc_gmp_v34_bit-vector_10.c:224). */
static const u32 V16_EXT_PRIMES[] = {
     509,  521,  523,  541,  547,  557,  563,  569,
     571,  577,  587,  593,  599,  601,  607,  613,
     617,  619,  631,  641,  643,  647,  653,  659,
     661,  673,  677,  683,  691,  701,  709,  719,
     727,  733,  739,  743,  751,  757,  761,  769,
     773,  787,  797,  809,  811,  821,  823,  827,
     829,  839,  853,  857,  859,  863,  877,  881,
     883,  887,  907,  911,  919,  929,  937,  941,
     947,  953,  967,  971,  977,  983,  991,  997,
    1009
};
#define V16_EXT_PRIME_COUNT (sizeof(V16_EXT_PRIMES) / sizeof(V16_EXT_PRIMES[0]))
#define V16_EXT_MASK_WORDS  16   /* 16 * 64 = 1024 ≥ max ext prime 1009 */

/* ============================================================================
 * Host configuration.
 * ============================================================================ */
typedef struct {
    int       target_len;
    int       sieve_depth;       /* forbidden-residue positions per prime;
                                  * defaults to target_len. Setting depth >
                                  * target (v15 style) gives a much stronger
                                  * sieve at the cost of losing chains that
                                  * would die at position > target. */
    int       target_bits;       /* if > 0, restrict to this exact bit-size */
    int       bits_min;          /* q-iter bit-band lower edge (default 90) */
    int       bits_max;          /* q-iter bit-band upper edge (default 127) */
    int       bits_rotate_step;  /* width of bit-band sub-stripe for rotation
                                  *      across [bits_min, bits_max]. 0 = off (one
                                  *      static band, current behavior). >0 = rotate
                                  *      through sub-bands of this width every
                                  *      bits_rotate_period_sec, recomputing per-V
                                  *      Q windows on each rotation. Forces equal
                                  *      GPU time per bit-stratum instead of the
                                  *      volume-bias toward bits_max. */
    int       bits_rotate_period_sec; /* rotation period in seconds (default 30). */
    u32       immune_prime_cap;  /* D primorial cap; 0 = explicit set used */
    u32       omit_primes[V16_MAX_IMMUNE];
    int       n_omit;
    u32       immune_set[V16_MAX_IMMUNE];
    int       n_immune_set;
    u32       extra_primes[V16_MAX_IMMUNE];
    int       n_extra;
    int       exp_start;
    char      seed_name[128];
    int       validate_seed;
    u64       validate_A;        /* if non-zero, validate this specific A */
    char      validate_A_str[80];/* parsed from --validate-A or seed default */
    u32       wheel_prime_max;
    u32       sieve_prime_max;
    u64       tiles_per_launch;
    int       tiles_explicit;     /* 1 iff user passed --tiles on the CLI */
    int       sieve_variants;     /* 1 iff --sieve-variants (multi-D K1 path) */
    u64       max_tiles;
    u64       start_tile;        /* explicit A_tile cursor start (overrides --bits when > 0) */
    u64       end_tile;          /* explicit A_tile end (exclusive); 0 = derived from --bits */
    int       min_report_len;
    int       prove_forward;     /* 0=reverse top-down (default), 1=forward */
    int       cpu_prove;         /* 1 -> CPU GMP prove pool, skip K2 */
    int       prove_threads;     /* pthread workers when cpu_prove==1 */
    /* Bounded ProveJob queue depth (default 8).
     * Each slot owns a pre-allocated pinned host buffer (~32 MiB / slot
     * q-iter, ~256 MiB / slot legacy). CLI: --prove-queue-depth. */
    int       prove_queue_depth;
    u32       pool_primes[V16_MAX_IMMUNE];
    int       n_pool_primes;     /* fingerprint primes to test against p+1 */
    v16_pool_seed pool_seeds[V16_MAX_POOL_SEEDS];
    int       n_pool_seeds;
    /* Deferred 'auto:Z#:K' specs; expanded after build_D in main(). */
    char      pool_auto_specs[V16_MAX_AUTO_POOL][64];
    int       n_pool_auto_specs;
    int       gpu_id;
    int       time_limit_sec;
    int       report_sec;
    char      checkpoint_file[256];
    char      resume_file[256];
    char      log_file[256];
    /* Q-iter mode (opt-in). When 0 the binary behaves byte-identical to the
     * commit-5cd449f baseline (multi-D-in-kernel via --sieve-variants). */
    int       q_iter;            /* 1 iff --q-iter */
    int       q_order_random;    /* 0=sequential (default), 1=random */
    u64       q_seed;             /* PRNG seed when q_order_random; default 0 */
    u32       q_chunk;            /* Q values per kernel launch (default 2^24) */
    /* --q-band-mode exhaustive auto-picks sequential per-V when
     * Q-range bits <= exhaustive_max_q_bits, random otherwise. The global
     * q_order_random is ignored in exhaustive mode (per-V flag wins). */
    int       q_band_mode_exhaustive;  /* 0=fixed (default), 1=exhaustive */
    int       exhaustive_max_q_bits;   /* default 40 */
    /* Emit-dedup on p_0. Two causes for duplicate emits in a primorial soak:
     * (1) random Q-order samples with replacement, (2) structural multiplicity
     * — same p_0 reached via different (Q,V) factorizations of p_0+1. Both
     * produce identical p_0; dedup keeps the first verbose HIT, collapses
     * subsequent emits to compact HIT-DUP lines. Default ON. */
    int       dedup_p0;                /* 0=off, 1=on (default) */
    /* K1.5 GPU top-MR pre-filter.
     * k15_enabled defaults OFF for the first deployment so smoke runs can A/B
     * against the baseline; flip to default ON after the regression smoke is
     * green. k15_mr_reps in {2, 3}: witnesses {2,3} at reps=2 mirrors CPU's
     * `mpz_probab_prime_p(topv, 2)` exactly (no-stronger-than-CPU); reps=3
     * adds {5} for a tighter drop ratio. */
    int       k15_enabled;        /* 0 = off (legacy K1→CPU), 1 = K1→K1.5→CPU */
    int       k15_mr_reps;        /* 2 or 3; default 2 (matches CPU reps=2) */
    /* Per-V active-prime list skips sieve primes
     * whose forbid mask is all-zero (q | D_V). Default ON: ~4-5/92 primes
     * pruned for the 3433V exp_d pool; ~10% K1 throughput. */
    int       active_primes;      /* 1 = enable per-V active list (default) */
    /* N-deep stream pool: N GPU streams in a ring; drain happens
     * when we wrap around to reuse a slot. With N=4 (default), up to 4 K1
     * launches are in flight simultaneously, hiding per-launch host overhead
     * that capped gpu= at ~85% on the 2-stream pipeline. Clamped to
     * [2, V16_MAX_STREAMS]. */
    int       n_streams;          /* default 4 */
    /* Startup engine-integrity self-test. Default ON: runs a CPU GMP suite of
     * known Cunningham-chain roots (CC5..CC17) plus q-iter-flavor Q*D_V=p_0+1
     * round-trip vectors before the main dispatch. --no-test skips. Modeled
     * on v15 run_self_test() (cc18_filter_cuda_CpC_v15.cu lines 2345-2470).
     * Added 2026-05-14 after the 7-immune-only pool structural bug: every
     * campaign asserts the chain-follow machinery before burning GPU time. */
    int       run_self_test;      /* 1 = run startup self-test (default), 0 = skip */
    /* Extended CPU sieve (primes 509..1009, depth-aware, pre-MR filter).
     * Default ON. Builds a per-V forbid bitmask in pinned host RAM,
     * checked at the top of cpu_prove_worker before any Fermat/MR. */
    int       use_ext_sieve;      /* 1 = enable extended CPU sieve (default), 0 = skip */
    /* Derived. */
    u32       D_primes[V16_MAX_IMMUNE];
    int       n_D_primes;
    u128      D;                 /* u128 product */
    char      D_dec[80];
    /* Computed wheel + sieve primes. */
    u32       wheel_primes[V16_MAX_WHEEL];
    int       n_wheel_primes;
    u32       sieve_primes[V16_MAX_SIEVE];
    int       n_sieve_primes;
} v16_cfg;

static v16_cfg g_cfg;
static volatile sig_atomic_t g_stop = 0;
static void on_sigint(int s) { (void)s; g_stop = 1; }

/* Extended-sieve forbid table (CPU side, depth-aware).
 * Layout: g_ext_forbid_table[v * V16_EXT_PRIME_COUNT * V16_EXT_MASK_WORDS
 *                            + p * V16_EXT_MASK_WORDS + w]
 *         = u64 bitmask word `w` of variant `v`'s forbid set for ext
 *         prime index `p`. Bit `(r & 63)` of word `r >> 6` is 1 iff Q ≡ r
 *         (mod V16_EXT_PRIMES[p]) would make some p_k (k ∈ [0, depth)) be
 *         divisible by that prime. NULL when --no-ext-sieve. */
static u64 *g_ext_forbid_table = NULL;
static int  g_ext_enabled = 0;        /* set after table is built */
static int  g_ext_total_bytes = 0;    /* for banner reporting */
static int  g_ext_hits = 0;           /* count of post-K1 survivors that hit ext sieve (atomic-incremented from prove-pool workers; read by the report path) */

/* ============================================================================
 * xorshift64* PRNG for random Q sampling (mirrors v15:788-825).
 * In random mode the dispatch loop picks Q_base = Q_min + rand128 % range
 * per-launch instead of walking the sequential cursor. This is the structural
 * fix for the "everything clusters at the bottom bit" phenomenon — random
 * samples uniformly across [bits_min, bits_max] per variant.
 * ============================================================================ */
static u64 rng_state = 1;
static void rng_seed(u64 seed) { rng_state = seed ? seed : 1; }
static u64 rng_next(void) {
    u64 x = rng_state;
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    rng_state = x;
    return x;
}
/* Best-effort high-entropy seed: /dev/urandom first, fallback to
 * CLOCK_MONOTONIC nanoseconds + getpid() + gpu_id. Returns the seed used
 * so the host can print it in the banner for reproducibility logging.
 * Never returns 0 (xorshift64* on 0 produces 0). */
static u64 rng_seed_entropy(int gpu_id) {
    u64 seed = 0;
    int fd = open("/dev/urandom", O_RDONLY);
    if (fd >= 0) {
        ssize_t r = read(fd, &seed, sizeof(seed));
        close(fd);
        if (r == (ssize_t)sizeof(seed) && seed != 0) {
            rng_seed(seed);
            return seed;
        }
    }
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    seed = (u64)ts.tv_sec * 1000000000ULL + (u64)ts.tv_nsec;
    seed ^= (u64)getpid() * 0x517CC1B727220A95ULL;
    seed ^= (u64)gpu_id * 0x9E3779B97F4A7C15ULL;
    if (seed == 0) seed = 1;
    rng_seed(seed);
    return seed;
}

/* ============================================================================
 * Host math helpers.
 * ============================================================================ */
static u32 mod_inv_u32(u32 a, u32 q) {
    long long m0 = (long long)q, x0 = 0, x1 = 1;
    long long aa = (long long)(a % q);
    if (q == 1) return 0;
    while (aa > 1) {
        long long qq = aa / m0;
        long long t = m0;
        m0 = aa % m0;
        aa = t;
        t = x0;
        x0 = x1 - qq * x0;
        x1 = t;
    }
    if (x1 < 0) x1 += (long long)q;
    return (u32)x1;
}

/* Barrett-style modulus for u32 thread_offset % q.
 *
 *   For prime q >= 2 with q < 2^32, precompute magic = floor(2^32 / q).
 *   Then est = (x * magic) >> 32 is an under-estimate of x/q (off by at most
 *   one), so r = x - est*q lies in [0, 2q). A single correction step
 *   (`if (r >= q) r -= q`) returns the true remainder.
 *
 *   This replaces the slow HW divide in the K1 inner loop. Per-prime cost
 *   becomes one 32x32->64 mul + one shift + one mul-sub + one cmp-sub.
 *
 *   Magic always fits in u32: floor(2^32 / q) <= 2^32 / 2 = 2^31 for q>=2.
 *   (For q == 1 the formula is undefined but we never sieve q == 1.)
 *
 *   Layout matches kernel shared-mem footprint (8 B per prime):
 *     - q     : u32 (4 B)
 *     - magic : u32 (4 B)
 */
typedef struct {
    u32 q;
    u32 magic;
} BarrettQ;

static inline u32 barrett_magic_u32(u32 q) {
    /* floor(2^32 / q). Use 64-bit divide. q must be >= 2. */
    return (u32)(((u64)1 << 32) / (u64)q);
}

static inline u32 barrett_mod_u32(u32 x, u32 q, u32 magic) {
    u32 hi = (u32)(((u64)x * (u64)magic) >> 32);
    u32 r  = x - hi * q;
    if (r >= q) r -= q;
    return r;
}

static u32 pow2_mod_u32(int e, u32 q) {
    u64 r = 1 % q;
    u64 base = 2 % q;
    while (e > 0) {
        if (e & 1) r = (r * base) % q;
        base = (base * base) % q;
        e >>= 1;
    }
    return (u32)r;
}

static int u128_bits(u128 v) {
    if (v == 0) return 0;
    int b = 0;
    u128 t = v;
    while (t) { b++; t >>= 1; }
    return b;
}
static void u128_to_dec(u128 v, char *out, size_t cap) {
    char buf[64];
    int n = 0;
    if (v == 0) { snprintf(out, cap, "0"); return; }
    while (v > 0 && n < (int)sizeof(buf)) {
        buf[n++] = (char)('0' + (int)(v % 10));
        v /= 10;
    }
    if (n >= (int)cap) n = (int)cap - 1;
    for (int i = 0; i < n; i++) out[i] = buf[n - 1 - i];
    out[n] = '\0';
}

/* ============================================================================
 * Seed / D construction.
 * ============================================================================ */
static int parse_csv_primes(const char *s, u32 *out, int max) {
    int n = 0;
    const char *p = s;
    while (*p) {
        while (*p == ' ' || *p == ',') p++;
        if (!*p) break;
        char *end;
        unsigned long v = strtoul(p, &end, 10);
        if (end == p) {
            fprintf(stderr, "ERROR: malformed prime in list near '%s'\n", p);
            return -1;
        }
        if (n >= max) {
            fprintf(stderr, "ERROR: too many primes in list (max %d)\n", max);
            return -1;
        }
        out[n++] = (u32)v;
        p = end;
    }
    return n;
}

static int is_prime_u32(u32 n) {
    if (n < 2) return 0;
    if (n < 4) return 1;
    if ((n & 1) == 0) return 0;
    for (u32 i = 3; (u64)i * i <= n; i += 2) {
        if (n % i == 0) return 0;
    }
    return 1;
}

static int build_D(v16_cfg *c) {
    /* Sources of D primes: explicit immune_set, or all small primes <= cap
     * minus omit set, plus extra primes. */
    u32 tmp[V16_MAX_IMMUNE];
    int n = 0;

    if (c->n_immune_set > 0) {
        for (int i = 0; i < c->n_immune_set; i++) {
            if (!is_prime_u32(c->immune_set[i])) {
                fprintf(stderr, "ERROR: immune-factor-set entry %u not prime\n",
                        c->immune_set[i]);
                return -1;
            }
            if (n < V16_MAX_IMMUNE) tmp[n++] = c->immune_set[i];
        }
    } else if (c->immune_prime_cap > 0) {
        for (size_t i = 0; i < V16_SMALL_PRIME_COUNT; i++) {
            u32 q = V16_SMALL_PRIMES[i];
            if (q > c->immune_prime_cap) break;
            int omit = 0;
            for (int j = 0; j < c->n_omit; j++) {
                if (c->omit_primes[j] == q) { omit = 1; break; }
            }
            if (!omit && n < V16_MAX_IMMUNE) tmp[n++] = q;
        }
    } else {
        fprintf(stderr, "ERROR: must specify --immune-prime or --immune-factor-set\n");
        return -1;
    }

    for (int i = 0; i < c->n_extra; i++) {
        u32 q = c->extra_primes[i];
        if (!is_prime_u32(q)) {
            fprintf(stderr, "ERROR: --extra-prime %u not prime\n", q);
            return -1;
        }
        int dup = 0;
        for (int j = 0; j < n; j++) if (tmp[j] == q) { dup = 1; break; }
        if (!dup && n < V16_MAX_IMMUNE) tmp[n++] = q;
    }

    /* Sort ascending for determinism. */
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (tmp[j] < tmp[i]) { u32 t = tmp[i]; tmp[i] = tmp[j]; tmp[j] = t; }

    c->n_D_primes = n;
    memcpy(c->D_primes, tmp, (size_t)n * sizeof(u32));

    c->D = 1;
    for (int i = 0; i < n; i++) {
        u128 nd = c->D * (u128)c->D_primes[i];
        if (nd / c->D_primes[i] != c->D) {
            fprintf(stderr, "ERROR: D overflows u128 at prime %u\n", c->D_primes[i]);
            return -1;
        }
        c->D = nd;
    }
    u128_to_dec(c->D, c->D_dec, sizeof(c->D_dec));
    return 0;
}

/* ============================================================================
 * Pool-spec parser: "Z#", "Z#/X", "Z#/X/Y", "Z#/X/Y/W", ...
 * Builds a list of immune primes (= primorial-Z primes minus omits).
 * Also accepts a bare comma-separated prime list (treated as explicit immune
 * set without primorial parsing). Returns 0 on success, -1 on error.
 * ============================================================================ */
static int parse_pool_spec(const char *spec, v16_pool_seed *out) {
    memset(out, 0, sizeof(*out));
    strncpy(out->name, spec, sizeof(out->name) - 1);
    out->name[sizeof(out->name) - 1] = '\0';

    /* Find Z# in the prefix. Spec like "53#" or "53#/X/Y" — find '#'. */
    const char *hash = strchr(spec, '#');
    if (hash) {
        /* primorial-style spec. */
        long Z = strtol(spec, NULL, 10);
        if (Z < 2 || Z > 500) {
            fprintf(stderr, "ERROR: bad primorial '%s'\n", spec);
            return -1;
        }
        /* Collect primes <= Z. */
        u32 tmp[V16_MAX_IMMUNE];
        int n = 0;
        for (size_t i = 0; i < V16_SMALL_PRIME_COUNT; i++) {
            if ((long)V16_SMALL_PRIMES[i] > Z) break;
            if (n < V16_MAX_IMMUNE) tmp[n++] = V16_SMALL_PRIMES[i];
        }
        /* Parse omits after each '/'. */
        const char *p = hash + 1;
        while (*p) {
            while (*p == '/' || *p == ' ') p++;
            if (!*p) break;
            char *end;
            u32 omit = (u32)strtoul(p, &end, 10);
            if (end == p) {
                fprintf(stderr, "ERROR: malformed omit in '%s'\n", spec);
                return -1;
            }
            int found = -1;
            for (int i = 0; i < n; i++) if (tmp[i] == omit) { found = i; break; }
            if (found < 0) {
                fprintf(stderr, "ERROR: prime %u not in %ld# (spec '%s')\n", omit, Z, spec);
                return -1;
            }
            for (int i = found; i < n - 1; i++) tmp[i] = tmp[i + 1];
            n--;
            p = end;
        }
        for (int i = 0; i < n; i++) out->primes[i] = tmp[i];
        out->n_primes = n;
    } else {
        /* '*'-separated list. Each token is "P" or "P^k" (k>=1).
         * P^k expands to P repeated k times in the primes array (so D picks up
         * P^k as a factor).
         *
         * Note: prime POWERS in D are a primorial-construction tweak — they
         * shrink the A search space by forcing extra small-prime factors into
         * (p+1) — they do NOT add immunity. Immunity at q only requires
         * q | (p+1), independent of higher powers. Use sparingly. */
        int n = 0;
        const char *p = spec;
        while (*p) {
            while (*p == ' ' || *p == '*') p++;
            if (!*p) break;
            char *end;
            unsigned long pr = strtoul(p, &end, 10);
            if (end == p) {
                fprintf(stderr, "ERROR: malformed prime in '%s'\n", spec);
                return -1;
            }
            int exp = 1;
            if (*end == '^') {
                p = end + 1;
                long e = strtol(p, &end, 10);
                /* Observed max exponent in the CC10+ pool data is 16.
                 * Defensive cap at 16. */
                if (end == p || e < 1 || e > 16) {
                    fprintf(stderr, "ERROR: bad exponent in '%s'\n", spec);
                    return -1;
                }
                exp = (int)e;
            }
            for (int i = 0; i < exp; i++) {
                if (n >= V16_MAX_IMMUNE) {
                    fprintf(stderr, "ERROR: too many primes in '%s'\n", spec);
                    return -1;
                }
                out->primes[n++] = (u32)pr;
            }
            p = end;
        }
        out->n_primes = n;
    }

    /* Compute D and check primality. */
    out->D = 1;
    for (int i = 0; i < out->n_primes; i++) {
        if (!is_prime_u32(out->primes[i])) {
            fprintf(stderr, "ERROR: %u in '%s' not prime\n", out->primes[i], spec);
            return -1;
        }
        u128 nd = out->D * (u128)out->primes[i];
        if (out->D != 0 && nd / out->primes[i] != out->D) {
            fprintf(stderr, "ERROR: D overflows u128 in '%s'\n", spec);
            return -1;
        }
        out->D = nd;
    }
    u128_to_dec(out->D, out->D_dec, sizeof(out->D_dec));
    return 0;
}

/* Parse a comma-separated list of pool specs into cfg->pool_seeds.
 * 'auto:...' specs are saved for deferred expansion after build_D. */
static int parse_seed_pool(const char *list, v16_cfg *c) {
    size_t L = strlen(list);
    char *buf = (char *)malloc(L + 1);
    if (!buf) return -1;
    memcpy(buf, list, L + 1);
    char *save = NULL;
    for (char *tok = strtok_r(buf, ",", &save); tok;
              tok = strtok_r(NULL, ",", &save)) {
        while (*tok == ' ') tok++;
        size_t tl = strlen(tok);
        while (tl > 0 && (tok[tl - 1] == ' ' || tok[tl - 1] == '\n')) tok[--tl] = '\0';
        if (!*tok) continue;
        if (strncmp(tok, "auto:", 5) == 0) {
            if (c->n_pool_auto_specs >= V16_MAX_AUTO_POOL) {
                fprintf(stderr, "ERROR: too many auto pool specs (max %d)\n",
                        V16_MAX_AUTO_POOL);
                free(buf); return -1;
            }
            strncpy(c->pool_auto_specs[c->n_pool_auto_specs], tok,
                    sizeof(c->pool_auto_specs[0]) - 1);
            c->pool_auto_specs[c->n_pool_auto_specs]
                              [sizeof(c->pool_auto_specs[0]) - 1] = '\0';
            c->n_pool_auto_specs++;
            continue;
        }
        if (c->n_pool_seeds >= V16_MAX_POOL_SEEDS) {
            fprintf(stderr, "ERROR: too many pool seeds (max %d)\n", V16_MAX_POOL_SEEDS);
            free(buf); return -1;
        }
        if (parse_pool_spec(tok, &c->pool_seeds[c->n_pool_seeds]) != 0) {
            free(buf); return -1;
        }
        c->n_pool_seeds++;
    }
    free(buf);
    return 0;
}

/* Load pool specs from a text file: one spec per line. Lines beginning with
 * '#' and blank lines are ignored; inline trailing '#' comments are
 * stripped, as is leading/trailing whitespace. Each surviving token is
 * routed through parse_seed_pool() so it supports the same syntaxes as the
 * --seed-pool CLI flag (primorial Z#/X/Y, P*Q*R lists, and auto:Z#:K
 * expansion deferred to expand_auto_pools). Returns 0 on success, -1 on
 * the first parse error. */
static int parse_seed_pool_file(const char *path, v16_cfg *c) {
    FILE *fp = fopen(path, "r");
    if (!fp) {
        fprintf(stderr, "ERROR: cannot open --seed-pool-file '%s': %s\n",
                path, strerror(errno));
        return -1;
    }
    char line[1024];
    int lineno = 0;
    int n_added = 0;
    while (fgets(line, sizeof(line), fp)) {
        lineno++;
        /* strip inline comment + trailing whitespace */
        char *hash = strchr(line, '#');
        if (hash) *hash = '\0';
        size_t L = strlen(line);
        while (L > 0 && (line[L - 1] == '\n' || line[L - 1] == '\r' ||
                          line[L - 1] == ' '  || line[L - 1] == '\t')) {
            line[--L] = '\0';
        }
        /* strip leading whitespace */
        char *tok = line;
        while (*tok == ' ' || *tok == '\t') tok++;
        if (!*tok) continue;
        if (parse_seed_pool(tok, c) != 0) {
            fprintf(stderr, "ERROR: --seed-pool-file '%s' line %d: failed "
                    "to parse '%s'\n", path, lineno, tok);
            fclose(fp);
            return -1;
        }
        n_added++;
    }
    fclose(fp);
    fprintf(stderr,
            "INFO: --seed-pool-file '%s' loaded %d pool entries "
            "(non-comment lines)\n", path, n_added);
    return 0;
}

/* Recursive helper: emit pool seeds for every k-subset of drop[] of size k.
 * `chosen` accumulates the current selection. */
static int emit_auto_combo(v16_cfg *c, long Z,
                           const u32 *drop, int n_drop,
                           int start, int k_remain,
                           u32 *chosen, int n_chosen) {
    if (k_remain == 0) {
        if (c->n_pool_seeds >= V16_MAX_POOL_SEEDS) {
            fprintf(stderr, "ERROR: pool overflow at auto-enum (max %d)\n",
                    V16_MAX_POOL_SEEDS);
            return -1;
        }
        v16_pool_seed *ps = &c->pool_seeds[c->n_pool_seeds];
        memset(ps, 0, sizeof(*ps));
        /* Build name "Z#" or "Z#/p1/p2/...". */
        int wpos = snprintf(ps->name, sizeof(ps->name), "%ld#", Z);
        for (int i = 0; i < n_chosen && wpos < (int)sizeof(ps->name) - 1; i++) {
            wpos += snprintf(ps->name + wpos, sizeof(ps->name) - (size_t)wpos,
                             "/%u", chosen[i]);
        }
        /* Build prime list = primes ≤ Z minus chosen. */
        int n_primes = 0;
        for (size_t i = 0; i < V16_SMALL_PRIME_COUNT; i++) {
            u32 q = V16_SMALL_PRIMES[i];
            if ((long)q > Z) break;
            int dropped = 0;
            for (int j = 0; j < n_chosen; j++)
                if (chosen[j] == q) { dropped = 1; break; }
            if (!dropped && n_primes < V16_MAX_IMMUNE) ps->primes[n_primes++] = q;
        }
        ps->n_primes = n_primes;
        ps->D = 1;
        for (int i = 0; i < n_primes; i++) {
            u128 nd = ps->D * (u128)ps->primes[i];
            if (ps->D != 0 && nd / ps->primes[i] != ps->D) {
                fprintf(stderr, "ERROR: D overflows u128 for %s\n", ps->name);
                return -1;
            }
            ps->D = nd;
        }
        u128_to_dec(ps->D, ps->D_dec, sizeof(ps->D_dec));
        c->n_pool_seeds++;
        return 0;
    }
    for (int i = start; i <= n_drop - k_remain; i++) {
        chosen[n_chosen] = drop[i];
        if (emit_auto_combo(c, Z, drop, n_drop, i + 1, k_remain - 1,
                            chosen, n_chosen + 1) != 0) return -1;
    }
    return 0;
}

/* Expand 'auto:Z#:K' or 'auto:Z#:Kmin-Kmax' into combinatorial pool seeds.
 * Must-keep set = c->D_primes (the active sieve immune set); these primes
 * are never dropped, so every emitted variant satisfies D_V ⊇ D_sieve and
 * is reachable from the sieve hits. */
static int expand_auto_pool_one(v16_cfg *c, const char *spec) {
    /* spec format: "auto:Z#:K" or "auto:Z#:Kmin-Kmax" */
    if (strncmp(spec, "auto:", 5) != 0) return -1;
    const char *p = spec + 5;
    char *end;
    long Z = strtol(p, &end, 10);
    if (Z < 2 || Z > 500 || end == p) {
        fprintf(stderr, "ERROR: bad auto spec '%s' (expected auto:Z#:K)\n", spec);
        return -1;
    }
    if (*end != '#') {
        fprintf(stderr, "ERROR: bad auto spec '%s' (expected # after Z)\n", spec);
        return -1;
    }
    end++;
    if (*end != ':') {
        fprintf(stderr, "ERROR: bad auto spec '%s' (expected : after Z#)\n", spec);
        return -1;
    }
    end++;
    long Kmin = strtol(end, &end, 10);
    long Kmax = Kmin;
    if (*end == '-') {
        end++;
        Kmax = strtol(end, &end, 10);
    }
    if (Kmin < 0 || Kmax < Kmin) {
        fprintf(stderr, "ERROR: bad K range in '%s'\n", spec);
        return -1;
    }

    /* Drop-eligible = primes ≤ Z minus must-keep (c->D_primes). */
    u32 drop[V16_MAX_IMMUNE];
    int n_drop = 0;
    for (size_t i = 0; i < V16_SMALL_PRIME_COUNT; i++) {
        u32 q = V16_SMALL_PRIMES[i];
        if ((long)q > Z) break;
        int must_keep = 0;
        for (int j = 0; j < c->n_D_primes; j++) {
            if (c->D_primes[j] == q) { must_keep = 1; break; }
        }
        if (!must_keep && n_drop < V16_MAX_IMMUNE) drop[n_drop++] = q;
    }
    if (Kmax > n_drop) Kmax = n_drop;

    u32 chosen[V16_MAX_IMMUNE];
    for (long k = Kmin; k <= Kmax; k++) {
        if (emit_auto_combo(c, Z, drop, n_drop, 0, (int)k, chosen, 0) != 0)
            return -1;
    }
    return 0;
}

static int expand_auto_pools(v16_cfg *c) {
    for (int i = 0; i < c->n_pool_auto_specs; i++) {
        if (expand_auto_pool_one(c, c->pool_auto_specs[i]) != 0) return -1;
    }
    return 0;
}

/* ============================================================================
 * Parse a '*'-separated prime-power spec ("2*3*5*11^2*13*19") into an mpz_t.
 * Lifts the inner-loop of parse_pool_spec()'s bare-list branch but accumulates
 * into GMP so it can express D_V values for the self-test (some > 64 bits).
 *
 *   parse_v_spec_mpz("2*3*5*11^2*13*17*19*37", out)
 *     -> mpz_set(out, 2 * 3 * 5 * 11^2 * 13 * 17 * 19 * 37)
 *
 * Returns 0 on success, -1 on parse error. Caller must mpz_init() out first.
 * ============================================================================ */
static int parse_v_spec_mpz(const char *spec, mpz_t out) {
    mpz_set_ui(out, 1);
    const char *p = spec;
    int n_factors = 0;
    while (*p) {
        while (*p == ' ' || *p == '*') p++;
        if (!*p) break;
        char *end;
        unsigned long pr = strtoul(p, &end, 10);
        if (end == p || pr < 2) {
            fprintf(stderr, "ERROR: parse_v_spec_mpz: bad prime in '%s'\n", spec);
            return -1;
        }
        int exp = 1;
        if (*end == '^') {
            p = end + 1;
            long e = strtol(p, &end, 10);
            if (end == p || e < 1 || e > 64) {
                fprintf(stderr, "ERROR: parse_v_spec_mpz: bad exponent in '%s'\n", spec);
                return -1;
            }
            exp = (int)e;
        }
        if (!is_prime_u32((u32)pr)) {
            fprintf(stderr, "ERROR: parse_v_spec_mpz: %lu in '%s' is not prime\n", pr, spec);
            return -1;
        }
        for (int i = 0; i < exp; i++) mpz_mul_ui(out, out, pr);
        n_factors++;
        p = end;
    }
    if (n_factors == 0) {
        fprintf(stderr, "ERROR: parse_v_spec_mpz: empty spec '%s'\n", spec);
        return -1;
    }
    return 0;
}

/* ============================================================================
 * Startup engine-integrity self-test (modeled on v15 run_self_test(),
 * cc18_filter_cuda_CpC_v15.cu:2345-2470). Default-ON in v16; --no-test skips.
 *
 *   Suite 1: classic first-kind CC roots from CC5..CC17 — root prime, root
 *            check ((root-1)/2 not prime), chain-follow length assert.
 *   Suite 2: q-iter Q*D_V=p_0+1 round-trip vectors. Exercises the per-V
 *            decomposition the 7-immune-only-pool structural bug touched.
 *            Includes the 2 known CC18 world-record chains.
 *
 * CPU/GMP only. Returns 0 on full pass, 2 on any failure (main returns 2
 * cleanly — no mid-init exit). cfg is read for optional pool-presence check.
 * ============================================================================ */
static int run_v16_self_test(const v16_cfg *cfg) {
    fprintf(stderr, "\n=== v16 startup self-test ===\n");
    int passed = 0, failed = 0;

    struct {
        const char *root_dec;
        int expected_len;
        const char *desc;
    } classics[] = {
        { "2",                       5,  "CC5  classic (2->5->11->23->47)" },
        { "89",                      6,  "CC6  classic" },
        { "1122659",                 7,  "CC7  1st-kind" },
        { "19099919",                8,  "CC8  1st-kind" },
        { "85864769",                9,  "CC9  1st-kind" },
        { "26089808579",             10, "CC10 1st-kind" },
        { "665043081119",            11, "CC11 1st-kind" },
        { "554688278429",            12, "CC12 1st-kind" },
        { "4090932431513069",        13, "CC13 1st-kind" },
        { "95405042230542329",       14, "CC14 1st-kind" },
        { "2759832934171386593519",  17, "CC17 1st-kind (Wroblewski 2008)" },
    };
    int n_classics = (int)(sizeof(classics) / sizeof(classics[0]));

    mpz_t root, current, next, pred;
    mpz_inits(root, current, next, pred, NULL);

    for (int t = 0; t < n_classics; t++) {
        if (mpz_set_str(root, classics[t].root_dec, 10) != 0) {
            fprintf(stderr, "  FAIL [%s]: failed to parse root\n", classics[t].desc);
            failed++;
            continue;
        }
        if (mpz_probab_prime_p(root, 25) == 0) {
            fprintf(stderr, "  FAIL [%s]: root is not prime\n", classics[t].desc);
            failed++;
            continue;
        }
        if (mpz_cmp_ui(root, 2) > 0) {
            mpz_sub_ui(pred, root, 1);
            mpz_fdiv_q_2exp(pred, pred, 1);
            if (mpz_probab_prime_p(pred, 25) > 0) {
                fprintf(stderr, "  FAIL [%s]: (root-1)/2 is prime — not a chain root\n",
                        classics[t].desc);
                failed++;
                continue;
            }
        }
        int chain_len = 1;
        mpz_set(current, root);
        while (chain_len < 100) {
            mpz_mul_2exp(next, current, 1);
            mpz_add_ui(next, next, 1);
            if (mpz_probab_prime_p(next, 25) == 0) break;
            chain_len++;
            mpz_set(current, next);
        }
        if (chain_len == classics[t].expected_len) {
            fprintf(stderr, "  PASS [%s]: CC%d verified (%zu bits)\n",
                    classics[t].desc, chain_len, mpz_sizeinbase(root, 2));
            passed++;
        } else {
            fprintf(stderr, "  FAIL [%s]: expected CC%d, got CC%d\n",
                    classics[t].desc, classics[t].expected_len, chain_len);
            failed++;
        }
    }
    mpz_clears(root, current, next, pred, NULL);

    struct {
        const char *p0_hex;
        int         expected_len;
        const char *Q_hex;
        const char *DV_spec;
    } qvecs[] = {
        { "0xb149165114f43bf4f72249", 18, "0x54623df9bab82c1",
          "2*3*5*11^2*13*17*19*37" },
        { "0x57c4647809417bf1819285", 18, "0x469120b85a014e6db1",
          "2*3*5*11*13*19" },
        { "0xfe5f018b29d68951358f47", 17, "0xd3f90029b7b33d94",
          "2*3*5*11*13^2*19^2" },
        { "0x476a71909e048bca2dd4f5", 17, "0x396b8637a8edad2799",
          "2*3*5*11*13*19" },
        { "0x57654a4fdf5a8c81cbc835", 17, "0x4644a9d315efbec779",
          "2*3*5*11*13*19" },
    };
    int n_qvecs = (int)(sizeof(qvecs) / sizeof(qvecs[0]));

    mpz_t Q, DV, p0_synth, p0_decl, qcur, qnext;
    mpz_inits(Q, DV, p0_synth, p0_decl, qcur, qnext, NULL);

    for (int t = 0; t < n_qvecs; t++) {
        if (parse_v_spec_mpz(qvecs[t].DV_spec, DV) != 0) {
            fprintf(stderr, "  FAIL [q-iter v%d]: bad D_V spec\n", t);
            failed++;
            continue;
        }
        if (mpz_set_str(Q, qvecs[t].Q_hex, 0) != 0) {
            fprintf(stderr, "  FAIL [q-iter v%d]: bad Q hex\n", t);
            failed++;
            continue;
        }
        mpz_mul(p0_synth, Q, DV);
        mpz_sub_ui(p0_synth, p0_synth, 1);
        if (mpz_set_str(p0_decl, qvecs[t].p0_hex, 0) != 0) {
            fprintf(stderr, "  FAIL [q-iter v%d]: bad p_0 hex\n", t);
            failed++;
            continue;
        }
        if (mpz_cmp(p0_synth, p0_decl) != 0) {
            char synth_hex[256], decl_hex[256];
            gmp_snprintf(synth_hex, sizeof(synth_hex), "0x%ZX", p0_synth);
            gmp_snprintf(decl_hex,  sizeof(decl_hex),  "0x%ZX", p0_decl);
            fprintf(stderr,
                "  FAIL [q-iter v%d]: Q*D_V-1 = %s != declared p_0 %s\n",
                t, synth_hex, decl_hex);
            failed++;
            continue;
        }
        if (mpz_probab_prime_p(p0_synth, 25) == 0) {
            fprintf(stderr, "  FAIL [q-iter v%d]: p_0 is not prime\n", t);
            failed++;
            continue;
        }
        int chain_len = 1;
        mpz_set(qcur, p0_synth);
        while (chain_len < 100) {
            mpz_mul_2exp(qnext, qcur, 1);
            mpz_add_ui(qnext, qnext, 1);
            if (mpz_probab_prime_p(qnext, 25) == 0) break;
            chain_len++;
            mpz_set(qcur, qnext);
        }
        if (chain_len != qvecs[t].expected_len) {
            fprintf(stderr, "  FAIL [q-iter v%d]: expected CC%d, got CC%d\n",
                    t, qvecs[t].expected_len, chain_len);
            failed++;
            continue;
        }
        int pool_status = -1;
        if (cfg && cfg->n_pool_seeds > 0) {
            pool_status = 1;
            mpz_t pool_D;
            mpz_init(pool_D);
            for (int s = 0; s < cfg->n_pool_seeds; s++) {
                if (mpz_set_str(pool_D, cfg->pool_seeds[s].D_dec, 10) != 0) continue;
                if (mpz_cmp(pool_D, DV) == 0) { pool_status = 0; break; }
            }
            mpz_clear(pool_D);
        }
        size_t p0_bits = mpz_sizeinbase(p0_synth, 2);
        if (pool_status == 0) {
            fprintf(stderr,
                "  PASS [q-iter v%d]: CC%d round-trip (%zu bits, D_V in pool)\n",
                t, chain_len, p0_bits);
        } else if (pool_status == 1) {
            fprintf(stderr,
                "  PASS [q-iter v%d]: CC%d round-trip (%zu bits, D_V NOT in pool)\n",
                t, chain_len, p0_bits);
        } else {
            fprintf(stderr,
                "  PASS [q-iter v%d]: CC%d round-trip (%zu bits)\n",
                t, chain_len, p0_bits);
        }
        passed++;
    }
    mpz_clears(Q, DV, p0_synth, p0_decl, qcur, qnext, NULL);

    int total = n_classics + n_qvecs;
    if (failed == 0) {
        fprintf(stderr, "INFO: self-test PASSED (%d/%d vectors)\n", passed, total);
        return 0;
    } else {
        fprintf(stderr, "ERROR: self-test FAILED (%d failed of %d)\n",
                failed, total);
        return 2;
    }
}

/* Forward decl — definition appears later in file (the forbidden-residue
 * math used both by GPU forbid_mask and the CPU ext-sieve table). */
static int v16_forbidden_residues(u32 q, u128 D, int exp_start, int target,
                                  u32 *out, int max_out);

/* ============================================================================
 * Extended-sieve forbid table builder.
 *
 * Allocates g_ext_forbid_table = (n_pool_seeds × n_ext_primes × mask_words)
 * u64 entries and populates each (V, q_ext) mask with the depth-aware
 * forbidden Q residue set. Bit r of mask[V][q] is 1 iff Q ≡ r (mod q) would
 * make any of p_0..p_{depth-1} (with p_k = 2^k · Q · D_V − 1) divisible by q.
 *
 * Math is identical to the GPU forbid_mask, just stored CPU-side because
 * V16_Q_MASK_WORDS=8 caps the GPU mask at q ≤ 512. Called once at startup
 * after select_wheel_sieve() and pool load.
 *
 * Returns 0 on success, -1 on alloc failure. Sets g_ext_enabled = 1 on
 * success. Caller respects c->use_ext_sieve (no-op when 0).
 * ============================================================================ */
static int build_ext_forbid_table(const v16_cfg *c) {
    if (!c->use_ext_sieve) {
        fprintf(stderr, "INFO: ext-sieve disabled (--no-ext-sieve)\n");
        return 0;
    }
    int n_v = c->n_pool_seeds;
    int n_p = (int)V16_EXT_PRIME_COUNT;
    if (n_v == 0) {
        /* q-iter mode requires pool. Skip silently when no pool loaded. */
        return 0;
    }
    size_t entries = (size_t)n_v * (size_t)n_p * (size_t)V16_EXT_MASK_WORDS;
    size_t bytes = entries * sizeof(u64);
    g_ext_forbid_table = (u64*)calloc(entries, sizeof(u64));
    if (g_ext_forbid_table == NULL) {
        fprintf(stderr, "ERROR: ext-sieve: calloc(%zu bytes) failed\n", bytes);
        return -1;
    }
    g_ext_total_bytes = (int)bytes;
    int total_forbid_bits = 0;
    for (int v = 0; v < n_v; v++) {
        u128 Dv = c->pool_seeds[v].D;
        for (int p = 0; p < n_p; p++) {
            u32 q = V16_EXT_PRIMES[p];
            u32 fbd[V16_MAX_TARGET + 1];
            int nf = v16_forbidden_residues(q, Dv, c->exp_start,
                                            c->sieve_depth, fbd,
                                            (int)(sizeof(fbd)/sizeof(fbd[0])));
            if (nf < 0) continue;  /* q | D_V → no forbids (immune) */
            u64 *mask = &g_ext_forbid_table[((size_t)v * n_p + p)
                                            * V16_EXT_MASK_WORDS];
            for (int i = 0; i < nf; i++) {
                u32 r = fbd[i];
                mask[r >> 6] |= (1ULL << (r & 63));
                total_forbid_bits++;
            }
        }
    }
    g_ext_enabled = 1;
    fprintf(stderr,
        "#   ext-sieve  : ON  primes=%d range=[%u,%u]  table=%.2f MiB  "
        "avg forbid bits/V=%.2f\n",
        n_p, V16_EXT_PRIMES[0], V16_EXT_PRIMES[n_p - 1],
        (double)bytes / (1024.0 * 1024.0),
        (double)total_forbid_bits / (double)n_v);
    return 0;
}

/* Validate every pool seed satisfies the --sieve-variants preconditions:
 *   (1) D_V is divisible by cfg->D (V's immune set is a superset of sieve D)
 *   (2) V adds no primes that collide with active wheel_primes or sieve_primes
 * Must be called AFTER select_wheel_sieve() so c->wheel_primes / sieve_primes
 * are populated. Returns 0 on success, -1 on the first violation. */
static int validate_sieve_variants(const v16_cfg *c) {
    if (c->n_pool_seeds == 0) {
        fprintf(stderr,
            "ERROR: --sieve-variants requires at least one expanded pool seed\n");
        return -1;
    }
    for (int v = 0; v < c->n_pool_seeds; v++) {
        const v16_pool_seed *ps = &c->pool_seeds[v];
        /* (1) D_V % D_sieve == 0 (u128). */
        if (c->D == 0 || (ps->D % c->D) != 0) {
            fprintf(stderr,
                "ERROR: --sieve-variants: pool seed %s D=%s not divisible by "
                "sieve D=%s (V's primes must include the sieve's D primes)\n",
                ps->name, ps->D_dec, c->D_dec);
            return -1;
        }
        /* (2) V's added primes (those in V.primes not in cfg->D_primes)
         * colliding with active wheel/sieve primes used to be a hard
         * ERROR. It's now a WARN: the wheel/sieve was tuned for cfg.D, so
         * at each colliding prime q where q | D_V, V is fully immune but
         * the wheel only enumerates (q - nf_cfg)/q of the q residue
         * classes — variant V is under-counted by factor nf_cfg/q at that
         * prime (multiplicative across collisions). Hits found are real
         * CC chains (CPU prove trial-divides at every prime); coverage is
         * partial. Acceptable for CC-hunting; not for exhaustive
         * enumeration. See multi_d_design.md "Coverage bias". */
        for (int j = 0; j < ps->n_primes; j++) {
            u32 q = ps->primes[j];
            int in_sieve_D = 0;
            for (int k = 0; k < c->n_D_primes; k++) {
                if (c->D_primes[k] == q) { in_sieve_D = 1; break; }
            }
            if (in_sieve_D) continue;
            for (int k = 0; k < c->n_wheel_primes; k++) {
                if (c->wheel_primes[k] == q) {
                    fprintf(stderr,
                        "WARN: --sieve-variants: pool seed %s adds prime "
                        "%u (active wheel prime); variant V will be "
                        "under-counted at q=%u by factor (q-nf)/q. Hits "
                        "found are real CC chains; coverage is partial. "
                        "Continuing.\n",
                        ps->name, q, q);
                    continue;
                }
            }
            for (int k = 0; k < c->n_sieve_primes; k++) {
                if (c->sieve_primes[k] == q) {
                    fprintf(stderr,
                        "WARN: --sieve-variants: pool seed %s adds prime "
                        "%u (active sieve prime); variant V will be "
                        "under-counted at q=%u by factor (q-nf)/q. Hits "
                        "found are real CC chains; coverage is partial. "
                        "Continuing.\n",
                        ps->name, q, q);
                    continue;
                }
            }
        }
    }
    return 0;
}

/* ============================================================================
 * Forbidden-residue math for v16.
 *   For sieve prime q with D % q != 0, and for j in [0..target-1):
 *     F_q ⊇ { inv(D mod q) * inv(2^(E+j) mod q) mod q }
 *   Dedup. If |F_q| == q, the family is inadmissible at q.
 * ============================================================================ */
static int v16_forbidden_residues(u32 q, u128 D, int exp_start, int target,
                                  u32 *out, int max_out) {
    u32 D_mod_q = (u32)(D % q);
    if (D_mod_q == 0) return -1;   /* immune at q — no forbidden residues */
    u32 inv_D = mod_inv_u32(D_mod_q, q);
    int n = 0;
    for (int j = 0; j < target; j++) {
        u32 p2 = pow2_mod_u32(exp_start + j, q);
        u32 inv_p2 = mod_inv_u32(p2, q);
        u32 r = (u32)(((u64)inv_D * (u64)inv_p2) % q);
        int dup = 0;
        for (int k = 0; k < n; k++) if (out[k] == r) { dup = 1; break; }
        if (!dup) {
            if (n >= max_out) return -2;
            out[n++] = r;
        }
    }
    return n;
}

/* ============================================================================
 * Coefficient wheel (CRT-join, adapted from kt_wheel.c).
 *
 *   For each wheel prime q (coprime to D, ordered ascending):
 *     fbd = v16_forbidden_residues(q, D, E, target)
 *     immune = [0..q-1] \ fbd
 *     wheel := CRT-extend(wheel, m, q, immune); m *= q
 * ============================================================================ */
typedef struct {
    u64       primorial;
    u32       n_admissible;
    u64      *offsets;
    u64       fnv1a64;
} v16_wheel;

static int wheel_cmp_u64(const void *a, const void *b) {
    u64 x = *(const u64 *)a, y = *(const u64 *)b;
    if (x < y) return -1; if (x > y) return 1; return 0;
}

static u64 fnv1a64_u64s(const u64 *arr, size_t n) {
    u64 h = 0xcbf29ce484222325ULL;
    const u64 prime = 0x100000001b3ULL;
    for (size_t i = 0; i < n; i++) {
        u64 v = arr[i];
        for (int b = 0; b < 8; b++) { h ^= (u8)(v >> (8 * b)); h *= prime; }
    }
    return h;
}

static int build_wheel(const v16_cfg *c, v16_wheel *out) {
    out->primorial = 1;
    out->n_admissible = 0;
    out->offsets = NULL;
    out->fnv1a64 = 0;

    if (c->n_wheel_primes == 0) {
        out->offsets = (u64 *)malloc(sizeof(u64));
        if (!out->offsets) return -4;
        out->offsets[0] = 0;
        out->n_admissible = 1;
        return 0;
    }

    /* Compute final size up-front for an exact allocation. */
    size_t size_after[V16_MAX_WHEEL + 1];
    size_after[0] = 1;
    __uint128_t m128 = 1;
    for (int i = 0; i < c->n_wheel_primes; i++) {
        u32 q = c->wheel_primes[i];
        u32 fbd[V16_MAX_TARGET + 4];
        int nf = v16_forbidden_residues(q, c->D, c->exp_start, c->sieve_depth,
                                        fbd, (int)(sizeof(fbd) / sizeof(fbd[0])));
        if (nf < 0) {
            fprintf(stderr, "ERROR: wheel prime %u divides D — internal bug\n", q);
            return -1;
        }
        if (nf >= (int)q) {
            fprintf(stderr, "ERROR: family inadmissible at wheel prime %u (|F_q|=%d=q)\n",
                    q, nf);
            return -2;
        }
        size_after[i + 1] = size_after[i] * (size_t)(q - (u32)nf);
        __uint128_t mm = m128 * (__uint128_t)q;
        if (mm > (__uint128_t)UINT64_MAX) {
            fprintf(stderr, "ERROR: coefficient primorial overflows u64 (after q=%u)\n", q);
            return -3;
        }
        m128 = mm;
    }
    size_t final_size = size_after[c->n_wheel_primes];
    if (final_size > (size_t)V16_MAX_WHEEL_OFFSETS) {
        fprintf(stderr, "ERROR: wheel size %zu exceeds cap %u — reduce --wheel-max\n",
                final_size, V16_MAX_WHEEL_OFFSETS);
        return -1;
    }

    u64 *cur = (u64 *)malloc(sizeof(u64));
    if (!cur) return -4;
    cur[0] = 0;
    size_t cur_n = 1;
    u64 m = 1;

    for (int p = 0; p < c->n_wheel_primes; p++) {
        u32 q = c->wheel_primes[p];
        u32 fbd[V16_MAX_TARGET + 4];
        int nf = v16_forbidden_residues(q, c->D, c->exp_start, c->sieve_depth,
                                        fbd, (int)(sizeof(fbd) / sizeof(fbd[0])));
        u32 immune[64];
        int n_imm = 0;
        for (u32 r = 0; r < q; r++) {
            int isf = 0;
            for (int k = 0; k < nf; k++) if (fbd[k] == r) { isf = 1; break; }
            if (!isf) immune[n_imm++] = r;
        }

        u32 m_mod_q = (u32)(m % (u64)q);
        if (m_mod_q == 0) { free(cur); return -2; }
        u32 m_inv = mod_inv_u32(m_mod_q, q);

        size_t next_cap = size_after[p + 1];
        u64 *next = (u64 *)malloc(next_cap * sizeof(u64));
        if (!next) { free(cur); return -4; }

        size_t w = 0;
        for (size_t i = 0; i < cur_n; i++) {
            u64 a = cur[i];
            u32 a_mod_q = (u32)(a % (u64)q);
            for (int j = 0; j < n_imm; j++) {
                u32 b = immune[j];
                u32 bma = (b + q - a_mod_q) % q;
                u64 kk = ((u64)bma * (u64)m_inv) % (u64)q;
                next[w++] = a + m * kk;
            }
        }
        qsort(next, w, sizeof(u64), wheel_cmp_u64);
        size_t uw = 0;
        for (size_t i = 0; i < w; i++) {
            if (i == 0 || next[i] != next[i - 1]) next[uw++] = next[i];
        }
        free(cur);
        cur = next;
        cur_n = uw;
        m *= q;
    }

    out->primorial = m;
    out->n_admissible = (u32)cur_n;
    out->offsets = cur;
    out->fnv1a64 = fnv1a64_u64s(cur, cur_n);
    return 0;
}

static void wheel_free(v16_wheel *w) {
    if (!w) return;
    free(w->offsets);
    w->offsets = NULL;
    w->n_admissible = 0;
}

/* ============================================================================
 * Wheel/sieve prime selection: collect small primes coprime to D, partition
 * into wheel (<= wheel_prime_max, with primorial <= u64) and sieve (rest).
 * ============================================================================ */
static void select_wheel_sieve(v16_cfg *c) {
    /* Primes coprime to D, ascending. */
    u32 coprime[V16_MAX_SIEVE + V16_MAX_WHEEL + 16];
    int n_co = 0;
    for (size_t i = 0;
         i < V16_SMALL_PRIME_COUNT &&
         n_co < (int)(sizeof(coprime) / sizeof(coprime[0]));
         i++) {
        u32 q = V16_SMALL_PRIMES[i];
        if (q > c->sieve_prime_max && q > c->wheel_prime_max) break;
        int divides = 0;
        for (int j = 0; j < c->n_D_primes; j++) {
            if (c->D_primes[j] == q) { divides = 1; break; }
        }
        if (!divides) coprime[n_co++] = q;
    }

    /* Wheel = coprime primes <= wheel_prime_max, with primorial <= u64. */
    c->n_wheel_primes = 0;
    __uint128_t m = 1;
    for (int i = 0; i < n_co && c->n_wheel_primes < V16_MAX_WHEEL; i++) {
        u32 q = coprime[i];
        if (q > c->wheel_prime_max) break;
        __uint128_t mm = m * (__uint128_t)q;
        if (mm > (__uint128_t)UINT64_MAX) break;
        m = mm;
        c->wheel_primes[c->n_wheel_primes++] = q;
    }

    /* Sieve = coprime primes > wheel_prime_max, up to sieve_prime_max. */
    c->n_sieve_primes = 0;
    for (int i = 0; i < n_co && c->n_sieve_primes < V16_MAX_SIEVE; i++) {
        u32 q = coprime[i];
        if (q <= c->wheel_prime_max) continue;
        if (q > c->sieve_prime_max) break;
        if (q > V16_PACKED_WORDS_MAX * 64u) break;
        c->sieve_primes[c->n_sieve_primes++] = q;
    }
}

/* ============================================================================
 * Sieve packed forbidden masks (host).
 *   layout: mask[p_idx * V16_PACKED_WORDS_MAX + word_idx]
 *   bit r set iff residue r forbidden for prime p_idx.
 * ============================================================================ */
static int build_sieve_masks(const v16_cfg *c, u64 *masks /* sized n_sieve_primes * V16_PACKED_WORDS_MAX */) {
    memset(masks, 0,
           (size_t)c->n_sieve_primes * V16_PACKED_WORDS_MAX * sizeof(u64));
    for (int i = 0; i < c->n_sieve_primes; i++) {
        u32 q = c->sieve_primes[i];
        u32 fbd[V16_MAX_TARGET + 4];
        int nf = v16_forbidden_residues(q, c->D, c->exp_start, c->sieve_depth,
                                        fbd, (int)(sizeof(fbd) / sizeof(fbd[0])));
        if (nf < 0) continue;                   /* shouldn't happen — coprime list */
        if (nf >= (int)q) {
            fprintf(stderr, "ERROR: family inadmissible at sieve prime %u\n", q);
            return -1;
        }
        u64 *row = &masks[(size_t)i * V16_PACKED_WORDS_MAX];
        for (int j = 0; j < nf; j++) {
            u32 r = fbd[j];
            row[r >> 6] |= (1ULL << (r & 63u));
        }
    }
    return 0;
}

/* ============================================================================
 * GPU constants and state.
 * ============================================================================ */
__constant__ u64  d_D_lo;
__constant__ u64  d_D_hi;
__constant__ u64  d_coeff_primorial;
__constant__ int  d_exp_start;
__constant__ int  d_target_len;
__constant__ int  d_target_bits;
__constant__ int  d_prove_forward;
__constant__ int  d_min_report_len;
__constant__ int  d_n_sieve_primes;
__constant__ u32  d_n_admissible_const;        /* matches host wh.n_admissible */
__constant__ u32  d_sieve_primes[V16_MAX_SIEVE];
/* coeff_primorial mod q[i] — used by per-tile residue setup. q ≤ 211 < 256,
 * so u8 would fit, but u32 keeps the multiply path simple and aligned. */
__constant__ u32  d_cp_mod_q[V16_MAX_SIEVE];
/* Trial-division primes for the prove path. */
__constant__ u32  d_trial_primes[V16_TRIAL_PRIMES_MAX];
__constant__ int  d_n_trial_primes;

/* ============================================================================
 * GPU u128 (ported from v15, trimmed to what v16 needs).
 * ============================================================================ */
typedef struct { u64 lo, hi; } gu128;

__device__ __forceinline__ gu128 gu128_from_u64(u64 v)         { gu128 r; r.lo=v; r.hi=0; return r; }
__device__ __forceinline__ gu128 gu128_add(gu128 a, gu128 b)   { gu128 r; r.lo=a.lo+b.lo; r.hi=a.hi+b.hi+(r.lo<a.lo); return r; }
__device__ __forceinline__ gu128 gu128_sub(gu128 a, gu128 b)   { gu128 r; r.lo=a.lo-b.lo; r.hi=a.hi-b.hi-(a.lo<b.lo); return r; }
__device__ __forceinline__ int   gu128_lt (gu128 a, gu128 b)   { return (a.hi<b.hi)||(a.hi==b.hi&&a.lo<b.lo); }
__device__ __forceinline__ int   gu128_eq (gu128 a, gu128 b)   { return a.lo==b.lo && a.hi==b.hi; }
__device__ __forceinline__ int   gu128_gte(gu128 a, gu128 b)   { return !gu128_lt(a,b); }
__device__ __forceinline__ int   gu128_is0(gu128 a)            { return a.lo==0 && a.hi==0; }
__device__ __forceinline__ gu128 gu128_shl1(gu128 a)           { gu128 r; r.hi=(a.hi<<1)|(a.lo>>63); r.lo=a.lo<<1; return r; }
__device__ __forceinline__ gu128 gu128_shr1(gu128 a)           { gu128 r; r.lo=(a.lo>>1)|(a.hi<<63); r.hi=a.hi>>1; return r; }
__device__ __forceinline__ gu128 gu128_shl_n(gu128 a, int n) {
    gu128 r;
    if (n == 0) return a;
    if (n >= 64) { r.hi = a.lo << (n - 64); r.lo = 0; }
    else         { r.hi = (a.hi << n) | (a.lo >> (64 - n)); r.lo = a.lo << n; }
    return r;
}
__device__ __forceinline__ int gu128_clz(gu128 a) {
    if (a.hi) return __clzll(a.hi);
    if (a.lo) return 64 + __clzll(a.lo);
    return 128;
}
__device__ __forceinline__ gu128 gu128_mul64(u64 a, u64 b) {
    gu128 r; r.lo = a * b;
    asm("mul.hi.u64 %0, %1, %2;" : "=l"(r.hi) : "l"(a), "l"(b));
    return r;
}
__device__ __forceinline__ gu128 gu128_mul_u64(gu128 a, u64 b) {
    gu128 lo = gu128_mul64(a.lo, b);
    u64 hi_part = a.hi * b;
    lo.hi += hi_part;
    return lo;
}
/* Low-128-bit u128 x u128 multiply for K1.5. p_top fits in
 * u128 by construction, and this helper is only used by
 * K1.5 where caller has pre-checked target_len + bits_max < 128.
 *
 *   (a.hi*2^64 + a.lo) * (b.hi*2^64 + b.lo)
 *     = a.lo*b.lo                        (full 128 bits, kept)
 *     + (a.lo*b.hi + a.hi*b.lo) * 2^64   (low 64 bits of each cross-product
 *                                          shifted into r.hi)
 *     + a.hi*b.hi * 2^128                (overflows u128, discarded)
 *
 * Carry handling: r.lo from gu128_mul64(a.lo, b.lo) is already correct;
 * r.hi accumulates 3 contributions (the high 64 of a.lo*b.lo, plus the low
 * 64 of a.lo*b.hi and a.hi*b.lo). u64 wrap is the intended truncation —
 * anything that would carry into the (hypothetical) 129th bit is the
 * discarded high half. Self-tested at boot against a Python reference at
 * Q ≈ 2^64, 2^96, 2^120, and 2^127 boundaries. */
__device__ __forceinline__ gu128 gu128_mul_u128(gu128 a, gu128 b) {
    gu128 r = gu128_mul64(a.lo, b.lo);
    r.hi += a.lo * b.hi;  /* discard upper 64 of this 128-bit cross-product */
    r.hi += a.hi * b.lo;  /* same */
    return r;
}
__device__ gu128 gu128_mod(gu128 a, gu128 m) {
    if (gu128_lt(a, m)) return a;
    if (gu128_is0(m))   return a;
    int shift = gu128_clz(m) - gu128_clz(a);
    if (shift < 0) return a;
    gu128 div = m;
    #pragma unroll 1
    for (int i = 0; i < shift; i++) div = gu128_shl1(div);
    for (int i = shift; i >= 0; i--) {
        if (gu128_gte(a, div)) a = gu128_sub(a, div);
        div = gu128_shr1(div);
    }
    return a;
}
__device__ gu128 gu128_mulmod(gu128 a, gu128 b, gu128 m) {
    a = gu128_mod(a, m);
    gu128 r = gu128_from_u64(0);
    int top = 127 - gu128_clz(b);
    for (int bit = top; bit >= 0; bit--) {
        r = gu128_shl1(r);
        if (gu128_gte(r, m)) r = gu128_sub(r, m);
        u64 w = (bit >= 64) ? b.hi : b.lo;
        if ((w >> (bit & 63)) & 1) {
            r = gu128_add(r, a);
            if (gu128_gte(r, m)) r = gu128_sub(r, m);
        }
    }
    return r;
}
__device__ gu128 gu128_powmod(gu128 base, gu128 exp, gu128 mod) {
    gu128 r = gu128_from_u64(1);
    base = gu128_mod(base, mod);
    int top = 127 - gu128_clz(exp);
    for (int bit = top; bit >= 0; bit--) {
        r = gu128_mulmod(r, r, mod);
        u64 w = (bit >= 64) ? exp.hi : exp.lo;
        if ((w >> (bit & 63)) & 1) r = gu128_mulmod(r, base, mod);
    }
    return r;
}
__device__ int gu128_mr_witness(gu128 n, u64 w) {
    gu128 nm1 = gu128_sub(n, gu128_from_u64(1));
    gu128 d = nm1;
    int s = 0;
    while ((d.lo & 1) == 0) { d = gu128_shr1(d); s++; }
    gu128 x = gu128_powmod(gu128_from_u64(w), d, n);
    gu128 one = gu128_from_u64(1);
    if (gu128_eq(x, one) || gu128_eq(x, nm1)) return 1;
    for (int i = 0; i < s - 1; i++) {
        x = gu128_mulmod(x, x, n);
        if (gu128_eq(x, nm1)) return 1;
        if (gu128_eq(x, one)) return 0;
    }
    return 0;
}
__device__ int gu128_is_prime(gu128 n) {
    if (n.hi == 0) {
        if (n.lo < 2) return 0;
        if (n.lo == 2 || n.lo == 3) return 1;
        if ((n.lo & 1) == 0) return 0;
        if (n.lo < 9) return 1;
    } else if ((n.lo & 1) == 0) return 0;
    /* Small-prime quick reject. */
    const u64 sp[] = {3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71};
    #pragma unroll
    for (int i = 0; i < 19; i++) {
        gu128 r = gu128_mod(n, gu128_from_u64(sp[i]));
        if (gu128_is0(r))
            return (n.hi == 0 && n.lo == sp[i]) ? 1 : 0;
    }
    const u64 witnesses[] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71};
    #pragma unroll
    for (int i = 0; i < 20; i++) {
        if (!gu128_mr_witness(n, witnesses[i])) return 0;
    }
    return 1;
}

/* Chain position value: chain[pos] = n*2^pos + (2^pos - 1). */
__device__ __forceinline__ gu128 v16_chain_pos(gu128 n, int pos) {
    if (pos == 0) return n;
    gu128 shifted = gu128_shl_n(n, pos);
    gu128 mask = gu128_from_u64((1ULL << pos) - 1);
    return gu128_add(shifted, mask);
}

/* ============================================================================
 * GPU kernel: search the immune family.
 *
 *   Each thread handles one (A_tile, woff_idx) pair.
 *     A = (A_tile_base + tile_in_launch) * coeff_primorial + offsets[woff_idx]
 *     p = A * D * 2^E - 1     (computed lazily after sieve passes)
 *
 *   Grid layout (mirrors kt_filter_v8_swapped):
 *     gridDim.x = woff chunks = ceil(n_admissible / V16_THREADS)
 *     gridDim.y = tiles_per_launch
 *     blockDim.x = V16_THREADS
 *     tile_in_launch = blockIdx.y
 *     woff_idx = blockIdx.x * blockDim.x + threadIdx.x
 *
 *   For each thread:
 *     1. Compute A as u128.
 *     2. For each sieve prime q, compute A mod q (host-precomputed mods
 *        are unused — A varies per thread, so we just take a quick u128
 *        mod-by-u32 via Granlund-style or naive long-division).
 *        Check packed forbidden mask. If forbidden, return.
 *     3. Compute p = A * D * 2^E - 1.
 *     4. Top-pos MR (target-1). If composite, return.
 *     5. Forward chain-follow MR. Emit result if chain_len >= min_report_len.
 *
 * ============================================================================ */

/* ============================================================================
 * K1 — v16 sieve kernel (LUT version, v15-style hot path).
 *
 *   Each thread = one (A_tile, woff_idx) pair. Threads that survive every
 *   sieve prime emit a V16Survivor record to a global list with warp-
 *   aggregated atomics (one atomic per warp instead of one per surviving
 *   thread).
 *
 *   Grid:
 *     gridDim.x = ceil(n_admissible / V16_THREADS)
 *     gridDim.y = tiles_in_launch
 *     blockDim.x = V16_THREADS = 256
 *
 *   LUT layout: d_kmod_lut[i * d_n_admissible + woff_idx] (u8 per cell,
 *   consecutive threads -> consecutive bytes -> coalesced).
 *
 *   Per block: cooperative shared-mem residue setup
 *     s_tile[i] = (tile_base_r[i] + tile_in_launch * cp_mod_q[i]) mod q
 *   filled by threads 0..n_sieve_primes-1, then __syncthreads. Every thread
 *   then does
 *     r = s_tile[i] + kmod_lut[i][w]; if (r >= q) r -= q;
 *     forbidden? -> return.
 * ============================================================================ */
__global__ __launch_bounds__(V16_THREADS, 4)
void v16_sieve_kernel(
    u32                       d_n_admissible,
    const u32 * __restrict__ d_tile_base_r,
    u32                       tiles_in_launch,
    u32                       tile_stride,
    const u8  * __restrict__ d_kmod_lut,
    const u64 * __restrict__ d_sieve_masks,
    V16Survivor              *d_survivors,
    u32                      *d_survivor_count,
    u32                       d_max_survivors,
    unsigned long long       *d_cand_total,
    unsigned long long       *d_pass_sieve,
    /* New multi-D-in-kernel params. d_n_variants==1 keeps legacy fast/slow
     * paths byte-identical; d_n_variants>1 takes the variants path below.
     * d_variant_masks: per-(variant, sieve prime) packed forbidden mask,
     * layout [v * n_sieve_primes * V16_PACKED_WORDS_MAX +
     *         i * V16_PACKED_WORDS_MAX + word]. Unused when n_variants<=1.
     * d_variant_skip: per-(variant, sieve prime) byte, 1 iff q_i | D_V
     * (variant immune at q_i — mask row is all zeros, so skip the load).
     * Layout [v * n_sieve_primes + i]. Unused when n_variants<=1. */
    u32                       d_n_variants,
    const u64 * __restrict__ d_variant_masks,
    const u8  * __restrict__ d_variant_skip,
    unsigned long long       *d_variant_sieve_pass
) {
    u32 woff_idx = blockIdx.x * blockDim.x + threadIdx.x;

    /* Single-writer atomic for cand_total. */
    if (threadIdx.x == 0 && blockIdx.x == 0 && blockIdx.y == 0 && d_cand_total) {
        atomicAdd(d_cand_total,
                  (unsigned long long)tiles_in_launch *
                  (unsigned long long)d_n_admissible);
    }

    __shared__ u32 s_tile[V16_MAX_SIEVE];

    /* Fast path: one tile per block. Matches v15's K1 launch pattern, used
     * whenever n_admissible is large enough to fill the block (light-D).
     * Single-variant only (d_n_variants <= 1). */
    if (tile_stride == 1u && d_n_variants == 0u) {
        u32 tile_in_launch = blockIdx.y;
        int tx = (int)threadIdx.x;
        if (tx < d_n_sieve_primes) {
            u32 q   = d_sieve_primes[tx];
            u32 cpq = d_cp_mod_q[tx];
            u32 base = d_tile_base_r[tx];
            u32 r = base + (u32)((u64)tile_in_launch * (u64)cpq % (u64)q);
            if (r >= q) r -= q;
            s_tile[tx] = r;
        }
        __syncthreads();

        if (woff_idx >= d_n_admissible) return;

        int passed = 1;
        #pragma unroll 1
        for (int i = 0; i < d_n_sieve_primes; i++) {
            u32 q  = d_sieve_primes[i];
            u32 km = (u32)__ldg(&d_kmod_lut[(size_t)i * d_n_admissible + woff_idx]);
            u32 r  = s_tile[i] + km;
            if (r >= q) r -= q;
            const u64 *row = &d_sieve_masks[(size_t)i * V16_PACKED_WORDS_MAX];
            if ((row[r >> 6] >> (r & 63u)) & 1ULL) { passed = 0; break; }
        }
        if (!passed) return;

        unsigned active_mask = __activemask();
        unsigned pass_mask = __ballot_sync(active_mask, passed);
        if (pass_mask == 0) return;
        int lane = (int)threadIdx.x & 31;
        int leader = __ffs((int)pass_mask) - 1;
        u32 warp_base = 0;
        if (lane == leader) {
            warp_base = atomicAdd(d_survivor_count, (u32)__popc(pass_mask));
            if (d_pass_sieve)
                atomicAdd(d_pass_sieve, (unsigned long long)__popc(pass_mask));
        }
        warp_base = __shfl_sync(pass_mask, warp_base, leader);
        u32 lane_prefix = (u32)__popc(pass_mask & ((1u << lane) - 1u));
        u32 slot = warp_base + lane_prefix;
        if (slot < d_max_survivors) {
            d_survivors[slot].tile_in_launch = (u16)tile_in_launch;
            d_survivors[slot].variant_idx    = 0;
            d_survivors[slot].woff_idx       = woff_idx;
        }
        return;
    }

    /* Variants path: multi-D evaluation per (tile, woff). Each thread computes
     * the per-prime residue r_i once, then loops over variants v and walks the
     * per-variant packed forbidden mask. tile_stride==1 only here for clarity;
     * heavy-D + variants is uncommon (variants give the parallelism). */
    if (d_n_variants >= 1u && tile_stride == 1u) {
        u32 tile_in_launch = blockIdx.y;
        int tx = (int)threadIdx.x;
        if (tx < d_n_sieve_primes) {
            u32 q   = d_sieve_primes[tx];
            u32 cpq = d_cp_mod_q[tx];
            u32 base = d_tile_base_r[tx];
            u32 r = base + (u32)((u64)tile_in_launch * (u64)cpq % (u64)q);
            if (r >= q) r -= q;
            s_tile[tx] = r;
        }
        __syncthreads();

        if (woff_idx >= d_n_admissible) return;

        /* Compute r_i[] (one per sieve prime, kept in a thread-local array).
         * V16_MAX_SIEVE=256 is the upper bound; in practice <=60. The compiler
         * may spill to local memory if register pressure is too high. */
        u8 r_i[V16_MAX_SIEVE];
        #pragma unroll 1
        for (int i = 0; i < d_n_sieve_primes; i++) {
            u32 q  = d_sieve_primes[i];
            u32 km = (u32)__ldg(&d_kmod_lut[(size_t)i * d_n_admissible + woff_idx]);
            u32 r  = s_tile[i] + km;
            if (r >= q) r -= q;
            r_i[i] = (u8)r;
        }

        /* Loop over variants. Each surviving (tile, woff, v) is emitted as a
         * separate survivor record via warp-aggregated atomic (same pattern
         * as the fast path). */
        #pragma unroll 1
        for (u32 v = 0; v < d_n_variants; v++) {
            const u64 *vbase = &d_variant_masks[
                (size_t)v * (size_t)d_n_sieve_primes * V16_PACKED_WORDS_MAX];
            const u8 *vskip = &d_variant_skip[
                (size_t)v * (size_t)d_n_sieve_primes];
            int passed = 1;
            #pragma unroll 1
            for (int i = 0; i < d_n_sieve_primes; i++) {
                /* q_i | D_V -> V immune at q_i, mask row is all zeros.
                 * Skip the load entirely. */
                if (vskip[i]) continue;
                u32 r = (u32)r_i[i];
                const u64 *row = &vbase[(size_t)i * V16_PACKED_WORDS_MAX];
                if ((row[r >> 6] >> (r & 63u)) & 1ULL) { passed = 0; break; }
            }

            unsigned active_mask = __activemask();
            unsigned pass_mask = __ballot_sync(active_mask, passed);
            if (pass_mask == 0) continue;
            int lane = (int)threadIdx.x & 31;
            int leader = __ffs((int)pass_mask) - 1;
            u32 warp_base = 0;
            if (lane == leader) {
                warp_base = atomicAdd(d_survivor_count, (u32)__popc(pass_mask));
                if (d_pass_sieve)
                    atomicAdd(d_pass_sieve,
                              (unsigned long long)__popc(pass_mask));
                if (d_variant_sieve_pass)
                    atomicAdd(&d_variant_sieve_pass[v],
                              (unsigned long long)__popc(pass_mask));
            }
            warp_base = __shfl_sync(active_mask, warp_base, leader);
            if (passed) {
                u32 lane_prefix =
                    (u32)__popc(pass_mask & ((1u << lane) - 1u));
                u32 slot = warp_base + lane_prefix;
                if (slot < d_max_survivors) {
                    d_survivors[slot].tile_in_launch = (u16)tile_in_launch;
                    d_survivors[slot].variant_idx    = (u16)v;
                    d_survivors[slot].woff_idx       = woff_idx;
                }
            }
        }
        return;
    }

    /* Variants + tile_stride>1: walk tile_stride tiles per block, then per
     * tile compute r_i[] and loop variants. Used when block_size shrinks for
     * heavy-D (small n_admissible) while --sieve-variants is also ON. */
    if (d_n_variants >= 1u) {
        u32 tile_lo = blockIdx.y * tile_stride;
        u32 tile_hi = tile_lo + tile_stride;
        if (tile_hi > tiles_in_launch) tile_hi = tiles_in_launch;
        for (u32 tile_in_launch = tile_lo; tile_in_launch < tile_hi;
             tile_in_launch++) {
            for (u32 i = threadIdx.x; i < (u32)d_n_sieve_primes;
                 i += blockDim.x) {
                u32 q   = d_sieve_primes[i];
                u32 cpq = d_cp_mod_q[i];
                u32 base = d_tile_base_r[i];
                u32 r = base + (u32)((u64)tile_in_launch * (u64)cpq % (u64)q);
                if (r >= q) r -= q;
                s_tile[i] = r;
            }
            __syncthreads();

            int have_woff = (woff_idx < d_n_admissible);
            u8 r_i[V16_MAX_SIEVE];
            if (have_woff) {
                #pragma unroll 1
                for (int i = 0; i < d_n_sieve_primes; i++) {
                    u32 q  = d_sieve_primes[i];
                    u32 km = (u32)__ldg(
                        &d_kmod_lut[(size_t)i * d_n_admissible + woff_idx]);
                    u32 r  = s_tile[i] + km;
                    if (r >= q) r -= q;
                    r_i[i] = (u8)r;
                }
            }
            for (u32 v = 0; v < d_n_variants; v++) {
                const u64 *vbase = &d_variant_masks[
                    (size_t)v * (size_t)d_n_sieve_primes * V16_PACKED_WORDS_MAX];
                const u8 *vskip = &d_variant_skip[
                    (size_t)v * (size_t)d_n_sieve_primes];
                int passed = 0;
                if (have_woff) {
                    passed = 1;
                    #pragma unroll 1
                    for (int i = 0; i < d_n_sieve_primes; i++) {
                        /* q_i | D_V -> V immune at q_i, mask row all zeros. */
                        if (vskip[i]) continue;
                        u32 r = (u32)r_i[i];
                        const u64 *row = &vbase[(size_t)i * V16_PACKED_WORDS_MAX];
                        if ((row[r >> 6] >> (r & 63u)) & 1ULL) {
                            passed = 0; break;
                        }
                    }
                }
                unsigned pass_mask = __ballot_sync(0xFFFFFFFFu, passed);
                if (pass_mask) {
                    int lane = (int)threadIdx.x & 31;
                    int leader = __ffs((int)pass_mask) - 1;
                    u32 warp_base = 0;
                    if (lane == leader) {
                        warp_base = atomicAdd(d_survivor_count,
                                              (u32)__popc(pass_mask));
                        if (d_pass_sieve)
                            atomicAdd(d_pass_sieve,
                                      (unsigned long long)__popc(pass_mask));
                        if (d_variant_sieve_pass)
                            atomicAdd(&d_variant_sieve_pass[v],
                                      (unsigned long long)__popc(pass_mask));
                    }
                    warp_base = __shfl_sync(0xFFFFFFFFu, warp_base, leader);
                    if (passed) {
                        u32 lane_prefix =
                            (u32)__popc(pass_mask & ((1u << lane) - 1u));
                        u32 slot = warp_base + lane_prefix;
                        if (slot < d_max_survivors) {
                            d_survivors[slot].tile_in_launch =
                                (u16)tile_in_launch;
                            d_survivors[slot].variant_idx = (u16)v;
                            d_survivors[slot].woff_idx    = woff_idx;
                        }
                    }
                }
            }
            __syncthreads();
        }
        return;
    }

    /* Slow path: tile_stride > 1. Each block walks tile_stride consecutive
     * tiles. Used for wheel-collapse (heavy-D / sequential) mode where the
     * block is small enough that we want to keep it resident on the SM for
     * many tiles to amortize launch overhead. */
    u32 tile_lo = blockIdx.y * tile_stride;
    u32 tile_hi = tile_lo + tile_stride;
    if (tile_hi > tiles_in_launch) tile_hi = tiles_in_launch;

    for (u32 tile_in_launch = tile_lo; tile_in_launch < tile_hi; tile_in_launch++) {
        /* Strided s_tile init so block size is decoupled from n_sieve_primes. */
        for (u32 i = threadIdx.x; i < (u32)d_n_sieve_primes; i += blockDim.x) {
            u32 q   = d_sieve_primes[i];
            u32 cpq = d_cp_mod_q[i];
            u32 base = d_tile_base_r[i];
            u32 r = base + (u32)((u64)tile_in_launch * (u64)cpq % (u64)q);
            if (r >= q) r -= q;
            s_tile[i] = r;
        }
        __syncthreads();

        int passed = 0;
        if (woff_idx < d_n_admissible) {
            passed = 1;
            #pragma unroll 1
            for (int i = 0; i < d_n_sieve_primes; i++) {
                u32 q  = d_sieve_primes[i];
                u32 km = (u32)__ldg(&d_kmod_lut[(size_t)i * d_n_admissible + woff_idx]);
                u32 r  = s_tile[i] + km;
                if (r >= q) r -= q;
                const u64 *row = &d_sieve_masks[(size_t)i * V16_PACKED_WORDS_MAX];
                if ((row[r >> 6] >> (r & 63u)) & 1ULL) { passed = 0; break; }
            }
        }

        unsigned pass_mask = __ballot_sync(0xFFFFFFFFu, passed);
        if (pass_mask) {
            int lane = (int)threadIdx.x & 31;
            int leader = __ffs((int)pass_mask) - 1;
            u32 warp_base = 0;
            if (lane == leader) {
                warp_base = atomicAdd(d_survivor_count, (u32)__popc(pass_mask));
                if (d_pass_sieve)
                    atomicAdd(d_pass_sieve, (unsigned long long)__popc(pass_mask));
            }
            warp_base = __shfl_sync(0xFFFFFFFFu, warp_base, leader);
            if (passed) {
                u32 lane_prefix = (u32)__popc(pass_mask & ((1u << lane) - 1u));
                u32 slot = warp_base + lane_prefix;
                if (slot < d_max_survivors) {
                    d_survivors[slot].tile_in_launch = (u16)tile_in_launch;
                    d_survivors[slot].variant_idx    = 0;
                    d_survivors[slot].woff_idx       = woff_idx;
                }
            }
        }
        __syncthreads();
    }
}

/* ============================================================================
 * K1 (Q-iter): v16 Q-iteration kernel. Q is represented as a u128 host
 * base plus u32 thread_offset.
 *
 *   One launch = one variant V; threads sweep Q values for that V. The Q
 *   coordinate is split into a u128 host base (per launch, not seen by kernel)
 *   and a u32 thread_offset (kernel-side).
 *
 *   Bitmask forbid. For each prime q (<= 503), the SET of forbidden
 *     residues across all depth-j positions is a small subset of [0,q).
 *     The host precomputes a per-(V,q) bitmask `forbid_mask[V][q]`
 *     (V16_Q_MASK_WORDS = 8 u64 words = up to 512 bits). Per launch the
 *     host uploads only `Q_base_mod_q[q]` (n_sp u32 ~344 B) instead of the
 *     full depth×n_sp offset slab. The kernel does:
 *         idx = (om + Q_base_mod_q) mod q     (one cmp-sub, no divide)
 *         if (forbid_mask[V][q].words[idx>>6] & (1ULL<<(idx&63))) reject
 *     The depth-j inner loop is gone. Cost per (Q, prime) drops from O(depth)
 *     compares to O(1) shift+AND.
 *
 *   Barrett magic for `om = thread_offset mod q` (no HW divide).
 *
 *   Loop nesting: q outside j (was) — with the bitmask there is no j-loop, just one
 *   bitmask test per prime. Small q still reject early, warp-aligned.
 *
 *   Grid: gridDim.x = ceil(n_Q / (V16_THREADS · Q_per_thread)); blockDim.x=256.
 *   Each thread sweeps d_Q_per_thread consecutive thread_offset values from
 *   tid*Q_per_thread.
 *
 *   __launch_bounds__(256, 4) — Blackwell SM has ~100 KiB shared; bitmasking drops
 *   per-block shared from ~7.2 KiB to ~6.0 KiB, so 4 resident blocks fits
 *   easily.
 * ============================================================================ */
__global__ __launch_bounds__(V16_THREADS, 4)
void v16_q_iter_kernel(
    u32                       d_n_Q,             /* number of Q values in launch */
    u32                       d_Q_per_thread,
    u32                       d_variant_idx,
    u32                       d_q_n_sieve_primes,
    const u32 * __restrict__  d_q_sieve_primes,
    const BarrettQ * __restrict__ d_q_barrett,    /* per-prime (q, magic) */
    const u64 * __restrict__  d_q_forbid_mask,    /* V's [n_sp][WORDS] u64 */
    const u32 * __restrict__  d_q_base_mod,       /* per-launch [n_sp] u32 */
    u32                       d_q_depth,
    QSurvivor                *d_q_survivors,
    u32                      *d_q_survivor_count,
    u32                       d_q_max_survivors,
    unsigned long long       *d_q_total,
    unsigned long long       *d_q_pass,
    /* Per-V active-prime list. When d_active_primes
     * is non-NULL the inner loop iterates only n_active entries, each a u8
     * index into (s_qb, s_base_mod, s_forbid_mask). Skips primes where
     * q | D_V (mask is all-zero). NULL → legacy "all primes" path. */
    const u8 * __restrict__   d_active_primes,
    u32                       d_n_active,
    /* Overflow guard: set to 1 by K1 when warp_base + popc exceeds
     * d_q_max_survivors. K1.5 / host inspect this flag to know the next
     * count read may overshoot the buffer cap. Per-slot, reset to 0 each
     * launch by the memset at the top of the dispatch loop. */
    u32                      *d_q_overflow
) {
    /* Shared-mem copy: Barrett pairs, bitmask slab, and base_mod. */
    __shared__ BarrettQ s_qb[V16_MAX_SIEVE];                 /* (q, magic) */
    __shared__ u64 s_forbid_mask[V16_MAX_SIEVE * V16_Q_MASK_WORDS];
    __shared__ u32 s_base_mod[V16_MAX_SIEVE];                /* per-launch */
    /* Per-V active-prime list (u8 packed). */
    __shared__ u8  s_active_primes[V16_MAX_SIEVE];

    const u32 n_sp = d_q_n_sieve_primes;
    const u32 mask_count = n_sp * V16_Q_MASK_WORDS;
    const u32 use_active = (d_active_primes != NULL);
    const u32 n_iter = use_active ? d_n_active : n_sp;
    (void)d_q_sieve_primes;  /* kernel reads q from s_qb[i].q */
    (void)d_q_depth;         /* depth subsumed into bitmask */

    /* Cooperative copy. */
    for (u32 i = threadIdx.x; i < n_sp; i += blockDim.x) {
        s_qb[i]       = d_q_barrett[i];
        s_base_mod[i] = d_q_base_mod[i];
    }
    for (u32 i = threadIdx.x; i < mask_count; i += blockDim.x) {
        s_forbid_mask[i] = d_q_forbid_mask[i];
    }
    /* Cooperative copy of active-prime list (u8 per entry). */
    if (use_active) {
        for (u32 i = threadIdx.x; i < d_n_active; i += blockDim.x) {
            s_active_primes[i] = d_active_primes[i];
        }
    }
    __syncthreads();

    /* Single-writer Q_total atomic (one thread in the entire grid). */
    if (threadIdx.x == 0 && blockIdx.x == 0 && d_q_total) {
        atomicAdd(d_q_total, (unsigned long long)d_n_Q);
    }

    u32 tid = blockIdx.x * blockDim.x + threadIdx.x;
    u32 offset_base = tid * d_Q_per_thread;

    #pragma unroll 1
    for (u32 qi = 0; qi < d_Q_per_thread; qi++) {
        u32 thread_offset = offset_base + qi;
        if (thread_offset >= d_n_Q) break;

        int passed = 1;
        /* Outer: primes (ASC; early-reject at small q is warp-aligned).
         *
         * When use_active, iterate only the indices
         * in s_active_primes[] (skipping primes where q | D_V — mask all-
         * zero). Otherwise iterate all n_sp primes. The active list
         * preserves ascending-q order for warp-aligned early reject. */
        #pragma unroll 1
        for (u32 k = 0; k < n_iter; k++) {
            u32 i = use_active ? (u32)s_active_primes[k] : k;
            u32 q     = s_qb[i].q;
            u32 magic = s_qb[i].magic;
            /* Barrett: om = thread_offset mod q. */
            u32 hi = (u32)(((u64)thread_offset * (u64)magic) >> 32);
            u32 om = thread_offset - hi * q;
            if (om >= q) om -= q;
            /* idx = (om + Q_base_mod_q) mod q (one cmp-sub, no divide). */
            u32 idx = om + s_base_mod[i];
            if (idx >= q) idx -= q;
            /* Bitmask probe. base = i * WORDS. */
            const u64 *mask = &s_forbid_mask[(size_t)i * V16_Q_MASK_WORDS];
            if (mask[idx >> 6] & (1ULL << (idx & 63u))) {
                passed = 0;
                break;
            }
        }

        /* Warp-aggregated emit (same pattern as v16_sieve_kernel fast path). */
        unsigned active_mask = __activemask();
        unsigned pass_mask = __ballot_sync(active_mask, passed);
        if (pass_mask == 0) continue;
        int lane = (int)threadIdx.x & 31;
        int leader = __ffs((int)pass_mask) - 1;
        u32 warp_base = 0;
        if (lane == leader) {
            u32 popc = (u32)__popc(pass_mask);
            warp_base = atomicAdd(d_q_survivor_count, popc);
            if (d_q_pass)
                atomicAdd(d_q_pass, (unsigned long long)popc);
            /* Flag if this warp's emit range crosses
             * d_q_max_survivors. K1.5 / host will clamp/abort accordingly. */
            if (d_q_overflow && warp_base + popc > d_q_max_survivors)
                atomicExch(d_q_overflow, 1u);
        }
        warp_base = __shfl_sync(active_mask, warp_base, leader);
        if (passed) {
            u32 lane_prefix = (u32)__popc(pass_mask & ((1u << lane) - 1u));
            u32 slot = warp_base + lane_prefix;
            if (slot < d_q_max_survivors) {
                d_q_survivors[slot].thread_offset = thread_offset;
                d_q_survivors[slot].seed_idx      = (u16)d_variant_idx;
                d_q_survivors[slot]._pad          = 0;
            }
        }
    }
}

/* ============================================================================
 * K1.5: v16 top-position Miller-Rabin pre-filter.
 *
 *   Inserted between K1 (sieve) and the host D2H. Computes
 *     Q_full = Q_base + (u128)thread_offset
 *     p_top  = Q_full · D_V · 2^(E + top_pos) − 1
 *   and runs `mr_reps` Miller-Rabin witnesses on p_top. Survivors emit a
 *   K15Survivor (same wire layout as QSurvivor) so the CPU prove worker can
 *   skip its top-MR block.
 *
 *   Witness selection:
 *     reps=2 -> {2, 3}    (matches mpz_probab_prime_p(topv, 2))
 *     reps=3 -> {2, 3, 5} (tighter drop, marginally costlier kernel)
 *
 *   Grid: sized to the K1 survivor count from the previous launch on this
 *   slot. Kernel self-limits via *d_in_count,
 *   so the grid is just an upper bound on parallelism. Stream: launched on
 *   the same q_streams[cur_slot] as K1 so K1 → K1.5 → count-D2H is one
 *   serialized chain with no host sync. */
__global__ __launch_bounds__(V16_THREADS, 4)
void v16_k15_topmr_kernel(
    const QSurvivor * __restrict__ d_in_survivors,
    const u32       * __restrict__ d_in_count,
    const u64       * __restrict__ d_dv_lo,
    const u64       * __restrict__ d_dv_hi,
    u64                              q_base_lo,
    u64                              q_base_hi,
    int                              exp_start,
    int                              top_pos,
    int                              mr_reps,
    K15Survivor                     *d_out_survivors,
    u32                             *d_out_count,
    u32                              d_out_cap,
    unsigned long long              *d_k15_pass,
    /* K1's survivor buffer cap. K1's atomicAdd can
     * push *d_in_count past this when the prior K1 hit the in-buffer wall,
     * so K1.5 must clamp before indexing d_in_survivors. d_overflow lets
     * K1.5 also raise the flag (defense in depth) when it observes the
     * overshoot — host treats either source as a hard abort. */
    u32                              d_in_cap,
    u32                             *d_overflow)
{
    const u32 idx0   = blockIdx.x * blockDim.x + threadIdx.x;
    const u32 stride = gridDim.x  * blockDim.x;
    const u32 n_raw  = *d_in_count;
    const u32 n      = (n_raw > d_in_cap) ? d_in_cap : n_raw;
    if (n_raw > d_in_cap && idx0 == 0 && d_overflow)
        atomicExch(d_overflow, 1u);

    /* Bug 1a fix (fix/k15-safety): grid-stride loop. The launcher sizes the
     * grid from prev_k1_count[cur_slot] (the PREVIOUS launch's K1 survivor
     * count). The current launch's *d_in_count can be larger, in which case
     * threads beyond gridDim.x*blockDim.x would never execute and survivors
     * would be silently dropped. Iterate with a stride so every survivor is
     * processed. Whole-warp early-out when idx0 >= n keeps inactive warps
     * from entering the loop; surviving warps step uniformly (i += stride)
     * so the warp-collective emit below remains in lockstep. */
    if (idx0 >= n) return;
    for (u32 i = idx0; i < n; i += stride) {
        /* Reset per-iteration state. `passed` MUST be 0 by default so that
         * __ballot_sync excludes any lane whose witness chain rejected. */
        int passed = 0;
        QSurvivor qs;
        qs.thread_offset = 0;
        qs.seed_idx = 0;
        qs._pad = 0;

        qs = d_in_survivors[i];

        /* Q = Q_base + thread_offset (u128 add). */
        gu128 Q;
        Q.lo = q_base_lo + (u64)qs.thread_offset;
        Q.hi = q_base_hi + (Q.lo < q_base_lo ? 1ULL : 0ULL);

        /* D_V from per-seed table (global memory, read via __ldg for L2/RO
         * caching since the same V repeats across many threads in a
         * launch). */
        gu128 Dv;
        Dv.lo = __ldg(&d_dv_lo[qs.seed_idx]);
        Dv.hi = __ldg(&d_dv_hi[qs.seed_idx]);

        /* QD = Q * D_V (low 128 bits). */
        gu128 QD = gu128_mul_u128(Q, Dv);

        /* p_top = (QD << (exp_start + top_pos)) - 1. Host guards against
         * the >127-bit overflow case at startup (see config validation). */
        gu128 p_top = gu128_sub(gu128_shl_n(QD, exp_start + top_pos),
                                gu128_from_u64(1));

        /* Miller-Rabin witnesses. Strictly no-stronger-than CPU's reps=2. */
        passed = 1;
        const u64 witnesses[3] = {2ULL, 3ULL, 5ULL};
        #pragma unroll 1
        for (int w = 0; w < mr_reps; w++) {
            if (!gu128_mr_witness(p_top, witnesses[w])) {
                passed = 0;
                break;
            }
        }

        /* Warp-aggregated emit per stride iteration. The vote runs across
         * the warp's currently-active lanes; in a partial-final iteration
         * some lanes may already have exited via i<n, so __activemask()
         * shrinks correspondingly. The popc/prefix arithmetic still yields
         * a contiguous slot range that won't collide with other warps. */
        unsigned active_mask = __activemask();
        unsigned pass_mask = __ballot_sync(active_mask, passed);
        if (pass_mask != 0) {
            int lane = (int)threadIdx.x & 31;
            int leader = __ffs((int)pass_mask) - 1;
            u32 warp_base = 0;
            if (lane == leader) {
                warp_base = atomicAdd(d_out_count, (u32)__popc(pass_mask));
                if (d_k15_pass)
                    atomicAdd(d_k15_pass,
                              (unsigned long long)__popc(pass_mask));
            }
            warp_base = __shfl_sync(active_mask, warp_base, leader);
            if (passed) {
                u32 lane_prefix =
                    (u32)__popc(pass_mask & ((1u << lane) - 1u));
                u32 slot = warp_base + lane_prefix;
                if (slot < d_out_cap) {
                    d_out_survivors[slot].thread_offset = qs.thread_offset;
                    d_out_survivors[slot].seed_idx      = qs.seed_idx;
                    d_out_survivors[slot]._pad          = 0;
                }
            }
        }
    }
}

/* gu128_mul_u128 boot-time self-test kernel: boundary cases at Q ≈ 2^64,
 * 2^96, 2^120 + a few edges. Host calls this
 * with 7 (a, b) pairs and compares to the u128 reference computed on the
 * CPU. Hard-error on any mismatch. */
__global__ void v16_k15_mul_u128_selftest_kernel(
    const u64 *a_lo, const u64 *a_hi,
    const u64 *b_lo, const u64 *b_hi,
    u64 *out_lo, u64 *out_hi, int n)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    gu128 a, b;
    a.lo = a_lo[i]; a.hi = a_hi[i];
    b.lo = b_lo[i]; b.hi = b_hi[i];
    gu128 r = gu128_mul_u128(a, b);
    out_lo[i] = r.lo;
    out_hi[i] = r.hi;
}

/* ============================================================================
 * K2 — v16 prove kernel.
 *
 *   One thread per K1 survivor. Threads are densely packed (each warp has
 *   ~32 active candidates), so the expensive MR / chain-follow work no
 *   longer wastes 31/32 of the warp on tail-idle threads.
 *
 *   Per thread:
 *     1. Reconstruct A from (tile_in_launch, woff_idx).
 *     2. Build AD = A * D, p_top = (AD << (E+target-1)) - 1.
 *     3. Reverse-depth top-position MR.
 *     4. Trial divide / root MR / root check.
 *     5. Forward chain follow.
 *     6. Emit ChainResult if chain_len >= min_report_len.
 * ============================================================================ */
__global__ __launch_bounds__(V16_THREADS, 4)
void v16_prove_kernel(
    const u64 * __restrict__ d_offsets,
    const V16Survivor * __restrict__ d_survivors,
    const u32        * __restrict__ d_survivor_count_ptr,
    u64                       a_tile_base_lo,
    u32                      *d_chain_count,
    u64                      *d_chain_a_lo,
    u64                      *d_chain_a_hi,
    u32                      *d_chain_len,
    u32                      *d_chain_top_pass,
    u32                      *d_chain_root_prime,
    unsigned long long       *d_pass_top,
    unsigned long long       *d_chains_emitted
) {
    u32 idx = blockIdx.x * blockDim.x + threadIdx.x;
    u32 n = *d_survivor_count_ptr;
    if (idx >= n) return;

    V16Survivor s = d_survivors[idx];

    gu128 a_tile;
    a_tile.lo = a_tile_base_lo + s.tile_in_launch;
    a_tile.hi = (a_tile.lo < a_tile_base_lo) ? 1 : 0;
    gu128 A = gu128_mul_u64(a_tile, d_coeff_primorial);
    A = gu128_add(A, gu128_from_u64(__ldg(&d_offsets[s.woff_idx])));

    gu128 D; D.lo = d_D_lo; D.hi = d_D_hi;
    gu128 AD;
    if (D.hi == 0) {
        AD = gu128_mul_u64(A, D.lo);
    } else {
        AD = gu128_mul_u64(A, D.lo);
        gu128 hi_term = gu128_mul_u64(gu128_from_u64(A.hi), D.lo);
        gu128 sh; sh.hi = hi_term.lo; sh.lo = 0;
        AD = gu128_add(AD, sh);
        gu128 hi_term2 = gu128_mul_u64(A, D.hi);
        gu128 sh2; sh2.hi = hi_term2.lo; sh2.lo = 0;
        AD = gu128_add(AD, sh2);
    }

    int top_shift = d_exp_start + (d_target_len - 1);
    gu128 p_top = gu128_sub(gu128_shl_n(AD, top_shift), gu128_from_u64(1));

    if (!d_prove_forward) {
        if (!gu128_is_prime(p_top)) return;
        if (d_pass_top) {
            unsigned m = __activemask();
            if (((int)threadIdx.x & 31) == __ffs((int)m) - 1)
                atomicAdd(d_pass_top, (unsigned long long)__popc(m));
        }
    }

    gu128 p_lo = gu128_sub(gu128_shl_n(AD, d_exp_start), gu128_from_u64(1));

    if (d_target_bits > 0) {
        int p_bits = 128 - gu128_clz(p_lo);
        if (p_bits != d_target_bits) return;
    }

    #pragma unroll 1
    for (int i = 0; i < d_n_trial_primes; i++) {
        u32 q = d_trial_primes[i];
        if (gu128_is0(gu128_mod(p_lo, gu128_from_u64((u64)q)))) return;
    }

    if (!gu128_is_prime(p_lo)) return;

    {
        gu128 pred = gu128_shr1(gu128_sub(p_lo, gu128_from_u64(1)));
        if (pred.hi > 0 || pred.lo >= 2) {
            if (gu128_is_prime(pred)) return;
        }
    }
    if (d_chain_root_prime) atomicAdd(d_chain_root_prime, 1u);

    int chain_len = 1;
    gu128 current = p_lo;
    gu128 one = gu128_from_u64(1);
    while (chain_len < V16_MAX_TARGET * 2) {
        gu128 next = gu128_add(gu128_shl1(current), one);
        if (!gu128_is_prime(next)) break;
        chain_len++;
        current = next;
    }

    if (chain_len < d_min_report_len) return;

    u32 slot = atomicAdd(d_chain_count, 1u);
    if (slot < V16_MAX_RESULTS) {
        d_chain_a_lo[slot]    = A.lo;
        d_chain_a_hi[slot]    = A.hi;
        d_chain_len[slot]     = (u32)chain_len;
        d_chain_top_pass[slot] = (chain_len >= d_target_len) ? 1u : 0u;
    }
    if (d_chains_emitted) atomicAdd(d_chains_emitted, 1ULL);
}

/* ============================================================================
 * CPU prove pool — mirrors v15's production path. When --cpu-prove is on,
 * K1 emits to the survivor list (GPU 100% on sieve), the survivor list is
 * copied D->H, and a pthread pool runs GMP-based MR + chain follow.
 * GMP's MR for ~120-bit numbers is much faster than the branchy gu128 path,
 * which dominated K2 in early benches.
 * ============================================================================ */
typedef struct {
    int                       tid;
    int                       n_threads;
    const V16Survivor        *survivors;
    u32                       n_surv;
    u64                       a_tile_base;
    u64                       coeff_primorial;
    const u64                *wheel_offsets;
    const v16_cfg            *cfg;
    u64                      *pass_top;
    u64                      *root_prime;
    u64                      *chains_emitted;
    u64                      *chains_emit_dup;
    HostChainResult          *hits;
    u32                       hits_cap;
    u32                      *hits_count;
    pthread_mutex_t          *log_mutex;
    FILE                     *log_fp;
    const char               *seed_name;
    v16_pool_seed            *pool_seeds;    /* shared, mutate ->hits with atomic */
    int                       n_pool_seeds;
    u64                      *chain_hist;    /* shared histogram, indexed by chain_len */
    int                       chain_hist_max;
    /* Q-iter additions. When q_iter_mode != 0, q_survivors is non-NULL and
     * `survivors`/`wheel_offsets`/`a_tile_base`/`coeff_primorial` are
     * unused. Mutually exclusive paths in the worker body.
     *
     * Each survivor batch is tagged with the u128 Q_base
     * for the launch that produced it. Q_full = Q_base + (u128)thread_offset.
     * Carried as two u64 halves to keep the struct trivially copyable. */
    const QSurvivor          *q_survivors;
    int                       q_iter_mode;
    u64                       q_base_lo;
    u64                       q_base_hi;
    /* K1.5 already proved top-position primality on the GPU
     * with witnesses identical (and a superset) to what mpz_probab_prime_p
     * tests at reps=2. When k15_passed != 0, the CPU worker skips its top-MR
     * block — the work has already been done. */
    int                       k15_passed;
} ProveWorker;

/* ============================================================================
 * Emit-dedup on p_0.
 *
 * Same p_0 can be observed via:
 *   (a) random Q-order samples with replacement (xorshift64* per variant has
 *       no visited-Q memory; narrow Q-ranges → birthday collisions).
 *   (b) structural multiplicity: same p_0+1 has multiple V-factorizations in
 *       a primorial-quotient pool, so several (Q,V) framings discover the
 *       same chain.
 *
 * Both are intentional under the current design, but flood the hit.log with
 * redundant verbose lines. Dedup keeps the first emit verbose; suppresses
 * later ones to a compact HIT-DUP line.
 *
 * Storage is a single growable u64 array of fingerprints, guarded by its own
 * mutex (caller already holds log_mutex, so this never contends in practice).
 * A 64-bit fingerprint over hits in a single run keeps collision probability
 * <10^-9 at 10^4 hits — far above realistic counts.
 * ============================================================================ */
static pthread_mutex_t g_v16_dedup_mu = PTHREAD_MUTEX_INITIALIZER;
static u64           *g_v16_dedup_fps = NULL;
static size_t         g_v16_dedup_n   = 0;
static size_t         g_v16_dedup_cap = 0;

/* Fold mpz_t p to a 64-bit fingerprint: xor of limbs, then a Murmur-style
 * mixer to scatter bits so similar magnitudes don't cluster. */
static u64 v16_p0_fingerprint(const mpz_t p) {
    u64 h = 0;
    size_t n = mpz_size(p);
    for (size_t i = 0; i < n; i++) h ^= (u64)mpz_getlimbn(p, i);
    h ^= h >> 33; h *= 0xff51afd7ed558ccdULL;
    h ^= h >> 33; h *= 0xc4ceb9fe1a85ec53ULL;
    h ^= h >> 33;
    return h;
}

/* Returns 1 if fp was newly inserted (first sighting), 0 if duplicate.
 * On OOM the function fails-safe to "newly inserted" so the verbose emit
 * still happens — better to log redundantly than to silently drop a hit. */
static int v16_dedup_check_insert(u64 fp) {
    pthread_mutex_lock(&g_v16_dedup_mu);
    for (size_t i = 0; i < g_v16_dedup_n; i++) {
        if (g_v16_dedup_fps[i] == fp) {
            pthread_mutex_unlock(&g_v16_dedup_mu);
            return 0;
        }
    }
    if (g_v16_dedup_n >= g_v16_dedup_cap) {
        size_t nc = g_v16_dedup_cap ? g_v16_dedup_cap * 2 : 256;
        u64 *nb = (u64 *)realloc(g_v16_dedup_fps, nc * sizeof(u64));
        if (!nb) {
            pthread_mutex_unlock(&g_v16_dedup_mu);
            return 1;
        }
        g_v16_dedup_fps = nb;
        g_v16_dedup_cap = nc;
    }
    g_v16_dedup_fps[g_v16_dedup_n++] = fp;
    pthread_mutex_unlock(&g_v16_dedup_mu);
    return 1;
}

static void *cpu_prove_worker(void *arg) {
    ProveWorker *w = (ProveWorker *)arg;
    const v16_cfg *cfg = w->cfg;
    mpz_t A, D, p, topv, pred, cur;
    mpz_inits(A, D, p, topv, pred, cur, NULL);

    /* Set D from u128 (sieve D). When --sieve-variants is ON we override per
     * survivor below using the variant's D_V; otherwise D stays as cfg->D. */
    {
        u64 hi = (u64)(cfg->D >> 64);
        u64 lo = (u64)cfg->D;
        mpz_set_ui(D, (unsigned long)hi);
        mpz_mul_2exp(D, D, 64);
        mpz_add_ui(D, D, (unsigned long)lo);
    }
    const int target    = cfg->target_len;
    const int E         = cfg->exp_start;
    const int top_pos   = target - 1;
    const int prove_fwd = cfg->prove_forward;
    const int min_rep   = cfg->min_report_len;
    const u32 hits_cap  = w->hits_cap;
    const int use_variants = cfg->sieve_variants;
    const int q_iter_mode  = w->q_iter_mode;
    const int k15_passed   = w->k15_passed;

    u64 loc_top = 0, loc_root = 0, loc_emit = 0, loc_emit_dup = 0;

    for (u32 i = (u32)w->tid; i < w->n_surv; i += (u32)w->n_threads) {
        u128 Au = 0;       /* A for HIT-log fingerprint (q-iter: Q itself) */

        if (q_iter_mode) {
            /* Q-iter survivor: reconstruct
             *   Q_full = Q_base + (u128)thread_offset
             *   p_0   = Q_full · D_V · 2^E − 1
             * The target_bits == strict-equality block is skipped here.
             * Bit-size is informational only.
             * Q_base is carried per-batch in w->q_base_{lo,hi} (set by the
             * host launch loop alongside the survivor handoff). */
            QSurvivor qs = w->q_survivors[i];
            u16 v = qs.seed_idx;
            if ((int)v >= cfg->n_pool_seeds) continue;
            u128 Dv = cfg->pool_seeds[v].D;
            u64 dhi = (u64)(Dv >> 64);
            u64 dlo = (u64)Dv;
            mpz_set_ui(D, (unsigned long)dhi);
            mpz_mul_2exp(D, D, 64);
            mpz_add_ui(D, D, (unsigned long)dlo);

            /* Build mpz A = (q_base_hi << 64) | q_base_lo, then + thread_offset.
             * Same u128→mpz split pattern as the D upload above and as v15
             * (cc18_filter_cuda_CpC_v15.cu :3107-3117). */
            mpz_set_ui(A, (unsigned long)w->q_base_hi);
            mpz_mul_2exp(A, A, 64);
            mpz_add_ui(A, A, (unsigned long)w->q_base_lo);
            mpz_add_ui(A, A, (unsigned long)qs.thread_offset);

            /* Au tracks the low 128 bits of Q for the HIT log. */
            u128 Q_full = ((u128)w->q_base_hi << 64) | (u128)w->q_base_lo;
            Q_full += (u128)qs.thread_offset;
            Au = Q_full;

            /* Extended-sieve filter: check Q's residue mod each ext prime
             * against the precomputed forbid bitmask. Bails out before any
             * GMP work. Math is identical to the GPU forbid_mask, applied
             * to primes that exceed the 512-bit GPU mask cap. */
            if (g_ext_enabled) {
                int v_idx = (int)v;
                const u64 *base = &g_ext_forbid_table[
                    ((size_t)v_idx * V16_EXT_PRIME_COUNT) * V16_EXT_MASK_WORDS];
                int ext_rejected = 0;
                for (int ei = 0; ei < (int)V16_EXT_PRIME_COUNT; ei++) {
                    u32 q = V16_EXT_PRIMES[ei];
                    u32 r = (u32)(Q_full % (u128)q);
                    const u64 *mask = base + (size_t)ei * V16_EXT_MASK_WORDS;
                    if ((mask[r >> 6] >> (r & 63)) & 1ULL) {
                        ext_rejected = 1;
                        break;
                    }
                }
                if (ext_rejected) {
                    __atomic_fetch_add(&g_ext_hits, 1, __ATOMIC_RELAXED);
                    continue;
                }
            }

            mpz_mul(p, A, D);
            if (E > 0) mpz_mul_2exp(p, p, (unsigned long)E);
            mpz_sub_ui(p, p, 1);
            /* Skip target_bits filter intentionally (q-iter). */
        } else {
            V16Survivor s = w->survivors[i];
            u64 a_tile = w->a_tile_base + s.tile_in_launch;
            u64 woff   = w->wheel_offsets[s.woff_idx];

            /* --sieve-variants: each survivor names which variant V's sieve set
             * it passed. Compute p_0 with V's D_V (not cfg->D). When OFF,
             * variant_idx is always 0 and we use cfg->D. */
            if (use_variants) {
                u16 v = s.variant_idx;
                if ((int)v >= cfg->n_pool_seeds) continue;
                u128 Dv = cfg->pool_seeds[v].D;
                u64 hi = (u64)(Dv >> 64);
                u64 lo = (u64)Dv;
                mpz_set_ui(D, (unsigned long)hi);
                mpz_mul_2exp(D, D, 64);
                mpz_add_ui(D, D, (unsigned long)lo);
            }

            mpz_set_ui(A, a_tile);
            mpz_mul_ui(A, A, (unsigned long)w->coeff_primorial);
            mpz_add_ui(A, A, (unsigned long)woff);

            Au = (u128)a_tile * (u128)w->coeff_primorial + (u128)woff;

            mpz_mul(p, A, D);
            if (E > 0) mpz_mul_2exp(p, p, (unsigned long)E);
            mpz_sub_ui(p, p, 1);

            if (cfg->target_bits > 0) {
                size_t pb = mpz_sizeinbase(p, 2);
                if ((int)pb != cfg->target_bits) continue;
            }
        }

        /* Top-pos test (reverse-depth pre-filter). reps=2 is plenty: GMP's
         * probab_prime_p does BPSW + reps Miller-Rabin, so reps=2 already
         * gives a vanishingly small false-positive rate. False positives are
         * caught downstream at chain follow. Trims top-MR cost ~5-7x at
         * scale.
         *
         * When k15_passed != 0, the GPU's K1.5 stage already ran
         * Miller-Rabin with witnesses {2,3} or {2,3,5} on p_top — strictly
         * no-weaker than mpz_probab_prime_p(topv, 2) — so we skip the CPU
         * call entirely. Count K1.5 survivors into loc_top so downstream
         * accounting stays consistent. */
        if (!prove_fwd) {
            if (k15_passed) {
                loc_top++;
            } else {
                mpz_mul_2exp(topv, p, (unsigned long)top_pos);
                mpz_set_ui(pred, 1);
                mpz_mul_2exp(pred, pred, (unsigned long)top_pos);
                mpz_sub_ui(pred, pred, 1);
                mpz_add(topv, topv, pred);
                if (mpz_probab_prime_p(topv, 2) == 0) continue;
                loc_top++;
            }
        }

        /* Root MR (5 GMP reps; same false-pos tolerance as top — chain-follow
         * catches drift). */
        if (mpz_probab_prime_p(p, 5) == 0) continue;

        /* Root check: (p-1)/2 must NOT be prime. Birth-certificate fast path
         * (v34 OPT-I): if p ≡ 1 (mod q) for any small prime q > 2, then
         * q | (p-1)/2, so (p-1)/2 is composite → p IS the chain root and we
         * skip the expensive MR. Use sieve primes 67..503 (V16_SMALL_PRIMES
         * indices 18..95 — never in D_V since palette tops at 61). Covers
         * ~70% of roots; the remaining ~30% fall through to MR. */
        int birth_certified = 0;
        for (int bi = 18; bi < (int)V16_SMALL_PRIME_COUNT; bi++) {
            if (mpz_fdiv_ui(p, V16_SMALL_PRIMES[bi]) == 1) {
                birth_certified = 1;
                break;
            }
        }
        if (!birth_certified) {
            mpz_sub_ui(pred, p, 1);
            mpz_fdiv_q_2exp(pred, pred, 1);
            if (mpz_cmp_ui(pred, 2) >= 0 && mpz_probab_prime_p(pred, 2) != 0) continue;
        }
        loc_root++;

        int chain_len = 1;
        mpz_set(cur, p);
        while (chain_len < 2 * V16_MAX_TARGET) {
            mpz_mul_2exp(cur, cur, 1);
            mpz_add_ui(cur, cur, 1);
            if (mpz_probab_prime_p(cur, 2) == 0) break;
            chain_len++;
        }
        /* Reps=2 is fine for chain follow too: a false positive at step k
         * just inflates the recorded chain_len by 1; subsequent steps will
         * catch it. If a hit reaches min_report_len, the next block does a
         * high-reps reverification. */
        if (chain_len >= min_rep) {
            /* Reverify with high reps to weed out any reps=2 false positives. */
            mpz_set(cur, p);
            int verified = 1;
            for (int j = 1; j < chain_len; j++) {
                mpz_mul_2exp(cur, cur, 1);
                mpz_add_ui(cur, cur, 1);
                if (mpz_probab_prime_p(cur, 25) == 0) { verified = 0; chain_len = j; break; }
            }
            if (!verified && chain_len < min_rep) continue;
            /* Also reverify root with high reps. */
            if (mpz_probab_prime_p(p, 25) == 0) continue;
        }
        /* Always bin into the histogram so we see the full distribution. */
        if (chain_len >= 1) {
            int bucket = chain_len < w->chain_hist_max
                       ? chain_len : (w->chain_hist_max - 1);
            __atomic_add_fetch(&w->chain_hist[bucket], 1ULL, __ATOMIC_RELAXED);
        }
        if (chain_len >= min_rep) {
            u32 slot = __atomic_fetch_add(w->hits_count, 1u, __ATOMIC_RELAXED);
            if (slot < hits_cap) {
                w->hits[slot].A_lo      = (u64)Au;
                w->hits[slot].A_hi      = (u64)(Au >> 64);
                w->hits[slot].chain_len = (u32)chain_len;
                w->hits[slot].top_pass  = (chain_len >= target) ? 1u : 0u;
            }

            /* Bounded append helper: snprintf returns the would-be length, so
             * naive `wpos += snprintf(...)` lets wpos climb past sizeof(buf).
             * On the next call sizeof(buf) - wpos wraps to ~SIZE_MAX and the
             * subsequent snprintf writes off the stack frame → SIGSEGV. The
             * macro caps wpos at sizeof(buf)-1 and signals truncation. */
            #define V16_APPEND(buf, wpos, truncated, ...) do {              \
                if (!(truncated) && (size_t)(wpos) + 1u < sizeof(buf)) {    \
                    int _n = snprintf((buf) + (wpos),                       \
                                      sizeof(buf) - (size_t)(wpos),         \
                                      __VA_ARGS__);                         \
                    if (_n < 0 ||                                           \
                        (size_t)(wpos) + (size_t)_n >= sizeof(buf)) {       \
                        (truncated) = 1;                                    \
                        (wpos) = (int)sizeof(buf) - 1;                      \
                        if (sizeof(buf) >= 5) {                             \
                            memcpy((buf) + sizeof(buf) - 5, "...", 4);      \
                        }                                                   \
                    } else {                                                \
                        (wpos) += _n;                                       \
                    }                                                       \
                }                                                           \
            } while (0)

            /* Build fingerprint string: which pool primes divide (p+1)? */
            char fp_buf[1024] = "";
            if (cfg->n_pool_primes > 0) {
                mpz_t p1;
                mpz_init(p1);
                mpz_add_ui(p1, p, 1);
                int wpos = 0;
                int trunc = 0;
                V16_APPEND(fp_buf, wpos, trunc, "fp={");
                int first = 1;
                for (int k = 0; k < cfg->n_pool_primes; k++) {
                    u32 q = cfg->pool_primes[k];
                    if (mpz_divisible_ui_p(p1, q)) {
                        V16_APPEND(fp_buf, wpos, trunc, "%s%u",
                                   first ? "" : ",", q);
                        first = 0;
                    }
                }
                V16_APPEND(fp_buf, wpos, trunc, "}");
                mpz_clear(p1);
            }

            /* Attribute hit to any pool seed whose D divides (p+1).
             * Buffer is sized for typical long-run hits; if more pool
             * seeds match than fit, append is bounded and a "..." marker is
             * dropped at the tail (no segfault). Per-seed hit counters
             * (ps->hits) are still incremented for every match. */
            char pool_buf[8192] = "";
            if (w->n_pool_seeds > 0) {
                mpz_t p1, Dz;
                mpz_inits(p1, Dz, NULL);
                mpz_add_ui(p1, p, 1);
                int wpos = 0;
                int trunc = 0;
                V16_APPEND(pool_buf, wpos, trunc, " pool=[");
                int first = 1;
                for (int s = 0; s < w->n_pool_seeds; s++) {
                    v16_pool_seed *ps = &w->pool_seeds[s];
                    u64 hi = (u64)(ps->D >> 64);
                    u64 lo = (u64)ps->D;
                    mpz_set_ui(Dz, (unsigned long)hi);
                    mpz_mul_2exp(Dz, Dz, 64);
                    mpz_add_ui(Dz, Dz, (unsigned long)lo);
                    if (mpz_divisible_p(p1, Dz)) {
                        __atomic_add_fetch(&ps->hits, 1ULL, __ATOMIC_RELAXED);
                        V16_APPEND(pool_buf, wpos, trunc, "%s%s",
                                   first ? "" : ",", ps->name);
                        first = 0;
                    }
                }
                V16_APPEND(pool_buf, wpos, trunc, "]");
                mpz_clears(p1, Dz, NULL);
            }
            #undef V16_APPEND

            /* Live print every hit to stderr so it lands in the tailed
             * stdout.log (fluent-bit captures it). The log file is the
             * canonical record; this is for live visibility.
             *
             * When dedup_p0 is on, the first emit of a given p_0 stays
             * verbose (full fp/pool buffers); later emits collapse to
             * compact HIT-DUP lines so the log still shows the duplicate
             * arrived but doesn't repeat the multi-V framing. */
            int p0_is_new = 1;
            if (cfg->dedup_p0) {
                u64 p0_fp = v16_p0_fingerprint(p);
                p0_is_new = v16_dedup_check_insert(p0_fp);
            }
            pthread_mutex_lock(w->log_mutex);
            if (p0_is_new) {
                loc_emit++;
                gmp_fprintf(stderr,
                    "*** HIT CC%d seed=%s%s%s A=%Zd p_0=0x%Zx ***\n",
                    chain_len, w->seed_name,
                    fp_buf, pool_buf, A, p);
                fflush(stderr);
                if (w->log_fp) {
                    gmp_fprintf(w->log_fp,
                        "HIT seed=%s len=%d %s%s A=%Zd p_0=0x%Zx\n",
                        w->seed_name, chain_len,
                        fp_buf, pool_buf, A, p);
                    fflush(w->log_fp);
                }
            } else {
                loc_emit_dup++;
                gmp_fprintf(stderr,
                    "    HIT-DUP seed=%s len=%d A=%Zd p_0=0x%Zx\n",
                    w->seed_name, chain_len, A, p);
                fflush(stderr);
                if (w->log_fp) {
                    gmp_fprintf(w->log_fp,
                        "HIT-DUP seed=%s len=%d A=%Zd p_0=0x%Zx\n",
                        w->seed_name, chain_len, A, p);
                    fflush(w->log_fp);
                }
            }
            pthread_mutex_unlock(w->log_mutex);
        }
    }
    __atomic_add_fetch(w->pass_top,        loc_top,      __ATOMIC_RELAXED);
    __atomic_add_fetch(w->root_prime,      loc_root,     __ATOMIC_RELAXED);
    __atomic_add_fetch(w->chains_emitted,  loc_emit,     __ATOMIC_RELAXED);
    __atomic_add_fetch(w->chains_emit_dup, loc_emit_dup, __ATOMIC_RELAXED);
    mpz_clears(A, D, p, topv, pred, cur, NULL);
    return NULL;
}

static void run_cpu_prove(
    const V16Survivor *h_survivors, u32 n_surv,
    u64 a_tile_base, u64 coeff_primorial,
    const u64 *wheel_offsets,
    const v16_cfg *cfg,
    u64 *pass_top, u64 *root_prime, u64 *chains_emitted, u64 *chains_emit_dup,
    HostChainResult *hits, u32 hits_cap, u32 *hits_count,
    pthread_mutex_t *log_mutex, FILE *log_fp, const char *seed_name,
    u64 *chain_hist, int chain_hist_max,
    /* Q-iter args: when q_survivors != NULL, h_survivors must be NULL. The
     * worker dispatches on q_iter_mode internally. q_base_lo/hi carry the
     * per-batch u128 Q base for Q reconstruction.
     *
     * k15_passed != 0 tells the worker the top-pos MR has
     * already been done on the GPU (K1.5 emitted these survivors). The CPU
     * skips its mpz_probab_prime_p(topv, 2) call. */
    const QSurvivor *q_survivors, int q_iter_mode,
    u64 q_base_lo, u64 q_base_hi, int k15_passed) {
    int N = cfg->prove_threads > 0 ? cfg->prove_threads : 8;
    if (N > 64) N = 64;
    pthread_t th[64];
    ProveWorker ctx[64];
    for (int t = 0; t < N; t++) {
        ctx[t].tid             = t;
        ctx[t].n_threads       = N;
        ctx[t].survivors       = h_survivors;
        ctx[t].n_surv          = n_surv;
        ctx[t].a_tile_base     = a_tile_base;
        ctx[t].coeff_primorial = coeff_primorial;
        ctx[t].wheel_offsets   = wheel_offsets;
        ctx[t].cfg             = cfg;
        ctx[t].pass_top        = pass_top;
        ctx[t].root_prime      = root_prime;
        ctx[t].chains_emitted  = chains_emitted;
        ctx[t].chains_emit_dup = chains_emit_dup;
        ctx[t].hits            = hits;
        ctx[t].hits_cap        = hits_cap;
        ctx[t].hits_count      = hits_count;
        ctx[t].log_mutex       = log_mutex;
        ctx[t].log_fp          = log_fp;
        ctx[t].seed_name       = seed_name;
        ctx[t].pool_seeds      = (v16_pool_seed *)cfg->pool_seeds;
        ctx[t].n_pool_seeds    = cfg->n_pool_seeds;
        ctx[t].chain_hist      = chain_hist;
        ctx[t].chain_hist_max  = chain_hist_max;
        ctx[t].q_survivors     = q_survivors;
        ctx[t].q_iter_mode     = q_iter_mode;
        ctx[t].q_base_lo       = q_base_lo;
        ctx[t].q_base_hi       = q_base_hi;
        ctx[t].k15_passed      = k15_passed;
        pthread_create(&th[t], NULL, cpu_prove_worker, &ctx[t]);
    }
    for (int t = 0; t < N; t++) pthread_join(th[t], NULL);
}

/* ============================================================================
 * Async prove queue.
 *
 * Replaces the synchronous run_cpu_prove() call in the dispatch loop with
 * an independent-worker producer/consumer pool. The dispatch thread
 * submits a ProveJob and returns immediately so the next K1 can be queued
 * before the previous batch's CPU prove finishes.
 *
 * Architecture:
 *   - cfg->prove_threads pthread workers spawned once at run start.
 *   - Each worker independently fetches jobs from a bounded FIFO queue
 *     and runs cpu_prove_worker on the whole batch alone (n_threads=1).
 *     cpu_prove_worker accumulates output via __atomic_add_fetch, so
 *     multiple in-flight jobs are race-free.
 *   - Queue is bounded by prove_queue_depth. Each slot owns one
 *     pre-allocated pinned host survivor buffer (cudaMallocHost'd at
 *     init). Producer blocks on submit when full (backpressure).
 *   - Per-V drained-cursor advance happens inside the worker on job
 *     completion. Mutex-protected, monotonic-max.
 *   - prove_pool_drain_all() blocks until queue empty AND no job in
 *     flight. Called before every checkpoint and at run end.
 * ============================================================================ */

typedef struct ProveJob {
    int      buf_idx;
    u32      n_surv;
    int      q_iter_mode;
    int      k15_passed;
    u64      q_base_lo;
    u64      q_base_hi;
    int      v_idx;
    u64      drained_target_lo;
    u64      drained_target_hi;
    u64      a_tile_base;
    u64      coeff_primorial;
    const u64 *wheel_offsets;
} ProveJob;

#define PROVE_MAX_QUEUE 64

typedef struct ProvePool {
    const v16_cfg        *cfg;
    int                   n_workers;
    int                   queue_depth;
    u64                  *pass_top;
    u64                  *root_prime;
    u64                  *chains_emitted;
    u64                  *chains_emit_dup;
    HostChainResult      *hits;
    u32                   hits_cap;
    u32                  *hits_count;
    pthread_mutex_t      *log_mutex;
    FILE                 *log_fp;
    const char           *seed_name;
    u64                  *chain_hist;
    int                   chain_hist_max;
    u64                  *v_cursor_drained_lo;
    u64                  *v_cursor_drained_hi;
    pthread_mutex_t       drained_mutex;
    void                **bufs;
    size_t                buf_stride;
    u32                   buf_cap;
    int                  *free_stack;
    int                   free_top;
    pthread_mutex_t       free_mutex;
    pthread_cond_t        free_cond;
    ProveJob             *queue;
    int                   head, tail;
    int                   n_pending;
    int                   n_in_flight;
    pthread_mutex_t       queue_mutex;
    pthread_cond_t        job_avail;
    pthread_cond_t        queue_drained;
    pthread_cond_t        queue_space;
    pthread_t            *worker_ths;
    void                 *worker_args;
    int                   shutdown;
    double                total_prove_ms;
    pthread_mutex_t       stats_mutex;
} ProvePool;

struct PoolWorkerArg { ProvePool *pool; int tid; };

static int prove_pool_acquire_buf(ProvePool *pool) {
    pthread_mutex_lock(&pool->free_mutex);
    while (pool->free_top == 0) {
        pthread_cond_wait(&pool->free_cond, &pool->free_mutex);
    }
    int idx = pool->free_stack[--pool->free_top];
    pthread_mutex_unlock(&pool->free_mutex);
    return idx;
}

static void prove_pool_release_buf(ProvePool *pool, int idx) {
    pthread_mutex_lock(&pool->free_mutex);
    pool->free_stack[pool->free_top++] = idx;
    pthread_cond_signal(&pool->free_cond);
    pthread_mutex_unlock(&pool->free_mutex);
}

static void *prove_pool_buf_ptr(ProvePool *pool, int idx) {
    return pool->bufs[idx];
}

static void prove_pool_exec_job(ProvePool *pool, ProveJob *job) {
    ProveWorker w;
    memset(&w, 0, sizeof(w));
    w.tid             = 0;
    w.n_threads       = 1;
    w.n_surv          = job->n_surv;
    w.cfg             = pool->cfg;
    w.pass_top        = pool->pass_top;
    w.root_prime      = pool->root_prime;
    w.chains_emitted  = pool->chains_emitted;
    w.chains_emit_dup = pool->chains_emit_dup;
    w.hits            = pool->hits;
    w.hits_cap        = pool->hits_cap;
    w.hits_count      = pool->hits_count;
    w.log_mutex       = pool->log_mutex;
    w.log_fp          = pool->log_fp;
    w.seed_name       = pool->seed_name;
    w.pool_seeds      = (v16_pool_seed *)pool->cfg->pool_seeds;
    w.n_pool_seeds    = pool->cfg->n_pool_seeds;
    w.chain_hist      = pool->chain_hist;
    w.chain_hist_max  = pool->chain_hist_max;
    w.q_iter_mode     = job->q_iter_mode;
    w.k15_passed      = job->k15_passed;
    if (job->q_iter_mode) {
        w.q_survivors  = (const QSurvivor *)pool->bufs[job->buf_idx];
        w.q_base_lo    = job->q_base_lo;
        w.q_base_hi    = job->q_base_hi;
    } else {
        w.survivors       = (const V16Survivor *)pool->bufs[job->buf_idx];
        w.a_tile_base     = job->a_tile_base;
        w.coeff_primorial = job->coeff_primorial;
        w.wheel_offsets   = job->wheel_offsets;
    }
    cpu_prove_worker(&w);
}

static void *prove_pool_worker_main(void *arg) {
    struct PoolWorkerArg *a = (struct PoolWorkerArg *)arg;
    ProvePool *pool = a->pool;
    while (1) {
        pthread_mutex_lock(&pool->queue_mutex);
        while (!pool->shutdown && pool->n_pending == 0) {
            pthread_cond_wait(&pool->job_avail, &pool->queue_mutex);
        }
        if (pool->shutdown && pool->n_pending == 0) {
            pthread_mutex_unlock(&pool->queue_mutex);
            break;
        }
        ProveJob job = pool->queue[pool->head];
        pool->head = (pool->head + 1) % pool->queue_depth;
        pool->n_pending--;
        pool->n_in_flight++;
        pthread_cond_signal(&pool->queue_space);
        pthread_mutex_unlock(&pool->queue_mutex);

        struct timespec t_a, t_b;
        clock_gettime(CLOCK_MONOTONIC, &t_a);
        prove_pool_exec_job(pool, &job);
        clock_gettime(CLOCK_MONOTONIC, &t_b);
        double ms = (t_b.tv_sec - t_a.tv_sec) * 1000.0 +
                    (t_b.tv_nsec - t_a.tv_nsec) / 1e6;

        if (job.v_idx >= 0 && pool->v_cursor_drained_lo) {
            pthread_mutex_lock(&pool->drained_mutex);
            u128 cur = ((u128)pool->v_cursor_drained_hi[job.v_idx] << 64) |
                       (u128)pool->v_cursor_drained_lo[job.v_idx];
            u128 want = ((u128)job.drained_target_hi << 64) |
                        (u128)job.drained_target_lo;
            if (want > cur) {
                pool->v_cursor_drained_lo[job.v_idx] = job.drained_target_lo;
                pool->v_cursor_drained_hi[job.v_idx] = job.drained_target_hi;
            }
            pthread_mutex_unlock(&pool->drained_mutex);
        }

        pthread_mutex_lock(&pool->stats_mutex);
        pool->total_prove_ms += ms;
        pthread_mutex_unlock(&pool->stats_mutex);

        prove_pool_release_buf(pool, job.buf_idx);

        pthread_mutex_lock(&pool->queue_mutex);
        pool->n_in_flight--;
        if (pool->n_pending == 0 && pool->n_in_flight == 0) {
            pthread_cond_broadcast(&pool->queue_drained);
        }
        pthread_mutex_unlock(&pool->queue_mutex);
    }
    return NULL;
}

static int prove_pool_init(ProvePool *pool,
                           const v16_cfg *cfg,
                           u64 *pass_top, u64 *root_prime, u64 *chains_emitted,
                           u64 *chains_emit_dup,
                           HostChainResult *hits, u32 hits_cap, u32 *hits_count,
                           pthread_mutex_t *log_mutex, FILE *log_fp,
                           const char *seed_name,
                           u64 *chain_hist, int chain_hist_max,
                           u64 *v_cursor_drained_lo, u64 *v_cursor_drained_hi,
                           size_t buf_stride, u32 buf_cap) {
    memset(pool, 0, sizeof(*pool));
    pool->cfg          = cfg;
    pool->n_workers    = cfg->prove_threads > 0 ? cfg->prove_threads : 8;
    if (pool->n_workers > 64) pool->n_workers = 64;
    pool->queue_depth  = cfg->prove_queue_depth > 0 ? cfg->prove_queue_depth : 8;
    if (pool->queue_depth > PROVE_MAX_QUEUE) pool->queue_depth = PROVE_MAX_QUEUE;
    pool->pass_top       = pass_top;
    pool->root_prime     = root_prime;
    pool->chains_emitted  = chains_emitted;
    pool->chains_emit_dup = chains_emit_dup;
    pool->hits            = hits;
    pool->hits_cap       = hits_cap;
    pool->hits_count     = hits_count;
    pool->log_mutex      = log_mutex;
    pool->log_fp         = log_fp;
    pool->seed_name      = seed_name;
    pool->chain_hist     = chain_hist;
    pool->chain_hist_max = chain_hist_max;
    pool->v_cursor_drained_lo = v_cursor_drained_lo;
    pool->v_cursor_drained_hi = v_cursor_drained_hi;
    pool->buf_stride     = buf_stride;
    pool->buf_cap        = buf_cap;

    pthread_mutex_init(&pool->free_mutex, NULL);
    pthread_cond_init(&pool->free_cond, NULL);
    pthread_mutex_init(&pool->queue_mutex, NULL);
    pthread_cond_init(&pool->job_avail, NULL);
    pthread_cond_init(&pool->queue_drained, NULL);
    pthread_cond_init(&pool->queue_space, NULL);
    pthread_mutex_init(&pool->drained_mutex, NULL);
    pthread_mutex_init(&pool->stats_mutex, NULL);

    pool->bufs       = (void **)calloc((size_t)pool->queue_depth, sizeof(void *));
    pool->free_stack = (int *)calloc((size_t)pool->queue_depth, sizeof(int));
    pool->queue      = (ProveJob *)calloc((size_t)pool->queue_depth, sizeof(ProveJob));
    if (!pool->bufs || !pool->free_stack || !pool->queue) {
        fprintf(stderr, "ERROR: prove_pool_init: calloc OOM\n");
        return -1;
    }
    for (int i = 0; i < pool->queue_depth; i++) {
        cudaError_t e = cudaMallocHost(&pool->bufs[i],
                                       (size_t)buf_cap * buf_stride);
        if (e != cudaSuccess) {
            fprintf(stderr,
                "ERROR: prove_pool_init: cudaMallocHost buf %d (%zu bytes): %s\n",
                i, (size_t)buf_cap * buf_stride, cudaGetErrorString(e));
            return -1;
        }
        pool->free_stack[pool->free_top++] = i;
    }

    pool->worker_ths = (pthread_t *)calloc((size_t)pool->n_workers,
                                           sizeof(pthread_t));
    struct PoolWorkerArg *args =
        (struct PoolWorkerArg *)calloc((size_t)pool->n_workers, sizeof(*args));
    if (!pool->worker_ths || !args) {
        fprintf(stderr, "ERROR: prove_pool_init: worker calloc OOM\n");
        return -1;
    }
    pool->worker_args = args;
    for (int t = 0; t < pool->n_workers; t++) {
        args[t].pool = pool;
        args[t].tid  = t;
        if (pthread_create(&pool->worker_ths[t], NULL,
                           prove_pool_worker_main, &args[t]) != 0) {
            fprintf(stderr,
                "ERROR: prove_pool_init: pthread_create worker %d failed\n", t);
            return -1;
        }
    }
    printf("#   prove-pool : ON  workers=%d  queue_depth=%d  "
           "buf=%u survs × %zu B (%.2f MiB pinned total)\n",
           pool->n_workers, pool->queue_depth, buf_cap, buf_stride,
           (double)pool->queue_depth * (double)buf_cap * (double)buf_stride /
               (1024.0 * 1024.0));
    fflush(stdout);
    return 0;
}

static void prove_pool_submit(ProvePool *pool, ProveJob job) {
    pthread_mutex_lock(&pool->queue_mutex);
    while (pool->n_pending >= pool->queue_depth) {
        pthread_cond_wait(&pool->queue_space, &pool->queue_mutex);
    }
    pool->queue[pool->tail] = job;
    pool->tail = (pool->tail + 1) % pool->queue_depth;
    pool->n_pending++;
    pthread_cond_signal(&pool->job_avail);
    pthread_mutex_unlock(&pool->queue_mutex);
}

static void prove_pool_drain_all(ProvePool *pool) {
    pthread_mutex_lock(&pool->queue_mutex);
    while (pool->n_pending > 0 || pool->n_in_flight > 0) {
        pthread_cond_wait(&pool->queue_drained, &pool->queue_mutex);
    }
    pthread_mutex_unlock(&pool->queue_mutex);
}

static double prove_pool_total_prove_ms(ProvePool *pool) {
    pthread_mutex_lock(&pool->stats_mutex);
    double v = pool->total_prove_ms;
    pthread_mutex_unlock(&pool->stats_mutex);
    return v;
}

static void prove_pool_destroy(ProvePool *pool) {
    if (!pool->worker_ths) return;
    prove_pool_drain_all(pool);
    pthread_mutex_lock(&pool->queue_mutex);
    pool->shutdown = 1;
    pthread_cond_broadcast(&pool->job_avail);
    pthread_mutex_unlock(&pool->queue_mutex);
    for (int t = 0; t < pool->n_workers; t++) {
        pthread_join(pool->worker_ths[t], NULL);
    }
    free(pool->worker_ths);
    pool->worker_ths = NULL;
    free(pool->worker_args);
    pool->worker_args = NULL;
    if (pool->bufs) {
        for (int i = 0; i < pool->queue_depth; i++) {
            if (pool->bufs[i]) cudaFreeHost(pool->bufs[i]);
        }
        free(pool->bufs);
    }
    free(pool->free_stack);
    free(pool->queue);
    pthread_mutex_destroy(&pool->free_mutex);
    pthread_cond_destroy(&pool->free_cond);
    pthread_mutex_destroy(&pool->queue_mutex);
    pthread_cond_destroy(&pool->job_avail);
    pthread_cond_destroy(&pool->queue_drained);
    pthread_cond_destroy(&pool->queue_space);
    pthread_mutex_destroy(&pool->drained_mutex);
    pthread_mutex_destroy(&pool->stats_mutex);
}

/* ============================================================================
 * Host glue: CLI / validation / search driver.
 * ============================================================================ */
static void usage(const char *argv0) {
    fprintf(stderr,
        "Usage: %s [options]\n"
        "  --target N             chain length to report (default 18)\n"
        "  --depth N              sieve forbidden-residue positions per prime.\n"
        "                         Independent of --target (depth < target, ==,\n"
        "                         and > target are all legal, design §15).\n"
        "                         Default = 12 (silent default when --depth\n"
        "                         unset). depth < target logs more survivors\n"
        "                         (profiling / SW-verification use case).\n"
        "                         depth == target is the production sieve.\n"
        "                         depth > target prescreens deeper than the\n"
        "                         reported chain length; WARN: at depth>target\n"
        "                         the sieve will reject exact-length-target\n"
        "                         chains whose p_target is small-prime composite\n"
        "                         (the sieve cannot distinguish 'chain ends\n"
        "                         here' from 'chain continues with composite').\n"
        "  --bits N               root bit size to filter to (default 0 = no filter).\n"
        "                         For --q-iter: collapses the band to [N,N]\n"
        "                         (back-compat). Prefer --bits-min/--bits-max.\n"
        "  --bits-min N1          q-iter bit-band lower edge (default 90).\n"
        "  --bits-max N2          q-iter bit-band upper edge (default 127).\n"
        "  --bits-rotate-step D   rotate through sub-bands of width D bits\n"
        "                         across [bits-min, bits-max]. Default 0 (off).\n"
        "                         Counters the volume-bias toward bits-max — each\n"
        "                         sub-band gets equal GPU time instead of half\n"
        "                         the time burning bits-max-1 only. Try D=3 at\n"
        "                         band [90,107] (6 sub-bands of 3 bits).\n"
        "  --bits-rotate-period S seconds per sub-band before rotating to the\n"
        "                         next stripe (default 30). Active only with\n"
        "                         --bits-rotate-step > 0.\n"
        "  --immune-prime P       D = product of primes <= P\n"
        "  --omit-prime Q         omit Q from --immune-prime product (repeatable)\n"
        "  --immune-factor-set L  explicit comma list, e.g. 2,3,5,7,11,13\n"
        "  --extra-prime Q        add Q to D (repeatable)\n"
        "  --exp-start E          E in p = A*D*2^E - 1 (default 0)\n"
        "  --seed-name NAME       label written to checkpoint/log\n"
        "  --validate-seed        diagnostic mode (no GPU run)\n"
        "  --validate-A N         decimal A to validate (--validate-seed)\n"
        "  --wheel-max P          max prime in coefficient wheel (default 53)\n"
        "  --sieve-max P          max prime in line sieve (default 211)\n"
        "  --tiles N              tiles per launch (default 1024)\n"
        "  --max-tiles N          stop after N tiles (default = no limit)\n"
        "  --min-report-len N     emit chains of length >= N (default 8)\n"
        "  --prove-order ORDER    reverse (default) | forward\n"
        "  --cpu-prove            drain K1 survivors to a CPU GMP pool (skips K2)\n"
        "  --prove-threads N      CPU workers when --cpu-prove (default 8)\n"
        "  --prove-queue-depth N  bounded prove-job queue depth (default 8). [1,64]\n"
        "  --pool-primes P,Q,...  primes to record in each hit's fingerprint\n"
        "  --seed-pool L          comma list of pool seed specs. Each spec is either\n"
        "                         a primorial Z#/X/Y/... (drop primes from Z#),\n"
        "                         '*'-separated prime list P*Q*R (with optional P^k),\n"
        "                         or 'auto:Z#:Kmin-Kmax' (combinatorial expansion of\n"
        "                         primorial Z# with K primes dropped from the\n"
        "                         drop-eligible set = primes<=Z minus the sieve seed).\n"
        "                         e.g. '41#/31,41#,3*5*11*13*19,auto:53#:1-3'.\n"
        "                         A hit is attributed to every variant V where D_V|(p+1).\n"
        "  --seed-pool-file PATH  load pool specs from a file (one spec per line,\n"
        "                         same syntax as --seed-pool tokens; '#' starts a\n"
        "                         comment). Additive with --seed-pool. Use case:\n"
        "                         historical fingerprint pool from pools/.\n"
        "  --sieve-variants       run pool seeds as sieve variants (multi-D K1).\n"
        "                         Requires --cpu-prove and at least one --seed-pool\n"
        "                         entry. Default OFF (post-prove attribution only).\n"
        "  --q-iter               enable Q-iteration multi-seed kernel. Requires\n"
        "                         --cpu-prove and --seed-pool. Replaces the CRT\n"
        "                         wheel; --wheel-max ignored. Incompatible with\n"
        "                         --sieve-variants.\n"
        "  --q-order ORDER        sequential (default) | random. Random samples\n"
        "                         Q_base uniformly in [Q_min(V), Q_max(V)] per\n"
        "                         launch — covers the full [bits_min, bits_max]\n"
        "                         band instead of clustering at bits_min.\n"
        "                         Recommended for wide bit bands where Q_max(V)\n"
        "                         − Q_min(V) exceeds practical sequential\n"
        "                         coverage (v15-style sampling).\n"
        "  --q-seed N             PRNG seed for --q-order random. Default 0 =\n"
        "                         seed from /dev/urandom (machine-unique).\n"
        "                         NOTE: not reproducible — every launch xors\n"
        "                         CLOCK_MONOTONIC nsec into the rng state so\n"
        "                         independent processes diverge. Checkpoint\n"
        "                         resume does NOT restore rng state. Treat\n"
        "                         --q-seed as a banner tag, not a re-runnable\n"
        "                         identity.\n"
        "  --q-chunk N            Q values per kernel launch (default 1<<24).\n"
        "  --q-band-mode MODE     fixed (default) | exhaustive. In exhaustive\n"
        "                         mode the engine picks PER VARIANT: sequential\n"
        "                         when (Q_max - Q_min) bit-width <=\n"
        "                         --exhaustive-max-q-bits, random otherwise.\n"
        "                         Useful for high-V primorial pools whose A-range\n"
        "                         can be exhausted in seconds-to-hours; smaller-V\n"
        "                         variants in the same pool stay on random. The\n"
        "                         global --q-order is informational only in this\n"
        "                         mode.\n"
        "  --exhaustive-max-q-bits N\n"
        "                         Q-range bit threshold for exhaustive mode\n"
        "                         (default 40). Variants with bit-width(Q_max -\n"
        "                         Q_min + 1) <= N use sequential; rest use random.\n"
        "  --dedup-p0 | --no-dedup-p0\n"
        "                         Collapse duplicate HIT emits keyed on p_0.\n"
        "                         Default ON. Duplicates arise from random Q-order\n"
        "                         with-replacement sampling and from same-p_0\n"
        "                         framings via different (Q,V) factorizations.\n"
        "                         First emit is full; later emits become HIT-DUP.\n"
        "  --k15 | --no-k15       enable/disable K1.5 GPU top-MR pre-filter.\n"
        "                         Default OFF. ON inserts an MR test\n"
        "                         on p_top between K1 and the host D2H, dropping\n"
        "                         ~7-12x of K1 survivors before CPU prove.\n"
        "  --k15-mr-reps N        2 (witnesses {2,3}; matches CPU reps=2) or 3\n"
        "                         (witnesses {2,3,5}). Default 2 (matches CPU,\n"
        "                         7.82× drop ratio at smoke; reps=3 ~marginal\n"
        "                         additional drop at >1.5× kernel cost).\n"
        "  --active-primes        (default ON) skip K1 sieve primes whose\n"
        "                         forbid_mask is all-zero (q | D_V); ~10%% K1 throughput.\n"
        "  --no-active-primes     disable active-prime list (legacy path,\n"
        "                         iterates all n_sp primes per Q).\n"
        "  --streams N            N-deep GPU stream pool (default 4, range\n"
        "                         [2,%d]). Higher N hides per-launch host\n"
        "                         overhead at cost of N × 32 MiB pinned RAM\n"
        "                         per slot.\n"
        "  --test                 (default ON) startup engine-integrity self-test:\n"
        "                         11 classic CC roots (CC5..CC17) + 5 q-iter\n"
        "                         Q*D_V=p_0+1 round-trip vectors. ~30s. Aborts\n"
        "                         the run if any vector fails.\n"
        "  --no-test              skip the startup self-test.\n"
        "  --gpu N                CUDA device id\n"
        "  --time SEC             wallclock budget\n"
        "  --report SEC           periodic report interval\n"
        "  --checkpoint FILE      write resume state\n"
        "  --resume FILE          read resume state\n"
        "  --log FILE             append-only hit log\n",
        argv0, V16_MAX_STREAMS);
}

static int parse_args(int argc, char **argv, v16_cfg *c) {
    memset(c, 0, sizeof(*c));
    c->target_len = 18;
    c->sieve_depth = 0;          /* 0 sentinel → silent default applied after
                                  * parse (12, design §15). Explicit --depth
                                  * overrides; depth < target is legal. */
    c->target_bits = 0;
    /* Bit-band defaults. When --bits N is passed
     * we collapse to [N,N] for back-compat. --bits-min/--bits-max override. */
    c->bits_min = 90;
    c->bits_max = 127;
    c->bits_rotate_step = 0;          /* rotation off by default */
    c->bits_rotate_period_sec = 30;   /* 30s per sub-band when rotation on */
    c->exp_start = 0;
    c->wheel_prime_max = 53;
    c->sieve_prime_max = 211;
    c->tiles_per_launch = 1024;
    c->max_tiles = 0;
    c->min_report_len = 8;
    c->prove_forward = 0;
    c->cpu_prove = 0;
    c->prove_threads = 8;
    c->prove_queue_depth = 8;
    c->gpu_id = 0;
    c->time_limit_sec = 0;
    c->report_sec = 5;
    c->q_iter = 0;
    c->q_order_random = 0;
    c->q_seed = 0;
    c->q_chunk = (1u << 24);   /* 16M Q per launch */
    c->q_band_mode_exhaustive = 0;   /* default fixed (global q-order) */
    c->exhaustive_max_q_bits  = 40;  /* Q-range ≤ 2^40 → sequential per-V */
    c->dedup_p0               = 1;   /* dedup on by default */
    /* K1.5 GPU top-MR pre-filter. OFF by default; reps=2
     * (witnesses {2,3}) matches CPU mpz_probab_prime_p(topv, 2) exactly. */
    c->k15_enabled = 0;
    c->k15_mr_reps = 2;
    /* Active-prime list ON by default. */
    c->active_primes = 1;
    /* N-deep stream pool: default 4 slots in the ring. */
    c->n_streams = 4;
    /* Startup engine-integrity self-test: default ON. --no-test skips. */
    c->run_self_test = 1;
    /* Extended CPU sieve (primes 509..1009 depth-aware): default ON. */
    c->use_ext_sieve = 1;
    strncpy(c->seed_name, "unnamed", sizeof(c->seed_name) - 1);

    static const struct option long_opts[] = {
        {"target",             required_argument, 0, 't'},
        {"depth",              required_argument, 0, 'd'},
        {"bits",               required_argument, 0, 'b'},
        {"immune-prime",       required_argument, 0, 'p'},
        {"omit-prime",         required_argument, 0, 'o'},
        {"immune-factor-set",  required_argument, 0, 'I'},
        {"extra-prime",        required_argument, 0, 'x'},
        {"exp-start",          required_argument, 0, 'e'},
        {"seed-name",          required_argument, 0, 'n'},
        {"validate-seed",      no_argument,       0, 'V'},
        {"validate-A",         required_argument, 0, 'A'},
        {"wheel-max",          required_argument, 0, 'W'},
        {"sieve-max",          required_argument, 0, 'S'},
        {"tiles",              required_argument, 0, 'T'},
        {"max-tiles",          required_argument, 0, 'M'},
        {"start-tile",         required_argument, 0, 's'},
        {"end-tile",           required_argument, 0, 'E'},
        {"min-report-len",     required_argument, 0, 'm'},
        {"prove-order",        required_argument, 0, 'O'},
        {"cpu-prove",          no_argument,       0, 'P'},
        {"prove-threads",      required_argument, 0, 'q'},
        {"prove-queue-depth",  required_argument, 0, 1030},
        {"pool-primes",        required_argument, 0, 'u'},
        {"seed-pool",          required_argument, 0, 'Q'},
        {"seed-pool-file",     required_argument, 0, 1005},
        {"sieve-variants",     no_argument,       0, 'Z'},
        /* Q-iter (opt-in). Long-only options; getopt assigns synthetic chars. */
        {"q-iter",             no_argument,       0, 1001},
        {"q-order",            required_argument, 0, 1002},
        {"q-seed",             required_argument, 0, 1003},
        {"q-chunk",            required_argument, 0, 1004},
        /* Per-V auto-pick: sequential when Q-range bits ≤ threshold, else random. */
        {"q-band-mode",           required_argument, 0, 1060},
        {"exhaustive-max-q-bits", required_argument, 0, 1061},
        /* Emit-dedup on p_0: collapse repeat-p_0 emits to HIT-DUP. */
        {"dedup-p0",           no_argument,       0, 1062},
        {"no-dedup-p0",        no_argument,       0, 1063},
        {"bits-min",           required_argument, 0, 1006},
        {"bits-max",           required_argument, 0, 1007},
        {"bits-rotate-step",   required_argument, 0, 1008},
        {"bits-rotate-period", required_argument, 0, 1009},
        /* K1.5 GPU top-MR pre-filter. */
        {"k15",                no_argument,       0, 1010},
        {"no-k15",             no_argument,       0, 1011},
        {"k15-mr-reps",        required_argument, 0, 1012},
        /* Per-V active-prime list (default ON). */
        {"active-primes",      no_argument,       0, 1020},
        {"no-active-primes",   no_argument,       0, 1021},
        /* N-deep stream pool. Default 4. */
        {"streams",            required_argument, 0, 1022},
        /* Startup engine-integrity self-test. Default ON. */
        {"test",               no_argument,       0, 1040},
        {"no-test",            no_argument,       0, 1041},
        /* Extended CPU sieve (509..1009, depth-aware, pre-MR). */
        {"ext-sieve",          no_argument,       0, 1050},
        {"no-ext-sieve",       no_argument,       0, 1051},
        {"gpu",                required_argument, 0, 'g'},
        {"time",               required_argument, 0, 'L'},
        {"report",             required_argument, 0, 'r'},
        {"checkpoint",         required_argument, 0, 'C'},
        {"resume",             required_argument, 0, 'R'},
        {"log",                required_argument, 0, 'l'},
        {"help",               no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };

    int ch;
    while ((ch = getopt_long(argc, argv, "", long_opts, NULL)) != -1) {
        switch (ch) {
            case 't': c->target_len = atoi(optarg); break;
            case 'd': c->sieve_depth = atoi(optarg); break;
            case 'b': {
                /* --bits N collapses the band to [N,N] (back-compat). */
                int n = atoi(optarg);
                c->target_bits = n;
                c->bits_min = n;
                c->bits_max = n;
                break;
            }
            case 'p': c->immune_prime_cap = (u32)strtoul(optarg, NULL, 10); break;
            case 'o':
                if (c->n_omit < V16_MAX_IMMUNE)
                    c->omit_primes[c->n_omit++] = (u32)strtoul(optarg, NULL, 10);
                break;
            case 'I': {
                int n = parse_csv_primes(optarg, c->immune_set, V16_MAX_IMMUNE);
                if (n < 0) return -1;
                c->n_immune_set = n;
                break;
            }
            case 'x':
                if (c->n_extra < V16_MAX_IMMUNE)
                    c->extra_primes[c->n_extra++] = (u32)strtoul(optarg, NULL, 10);
                break;
            case 'e': c->exp_start = atoi(optarg); break;
            case 'n':
                strncpy(c->seed_name, optarg, sizeof(c->seed_name) - 1);
                c->seed_name[sizeof(c->seed_name) - 1] = '\0';
                break;
            case 'V': c->validate_seed = 1; break;
            case 'A':
                strncpy(c->validate_A_str, optarg, sizeof(c->validate_A_str) - 1);
                c->validate_A_str[sizeof(c->validate_A_str) - 1] = '\0';
                break;
            case 'W': c->wheel_prime_max = (u32)strtoul(optarg, NULL, 10); break;
            case 'S': c->sieve_prime_max = (u32)strtoul(optarg, NULL, 10); break;
            case 'T':
                c->tiles_per_launch = strtoull(optarg, NULL, 10);
                c->tiles_explicit = 1;
                break;
            case 'M': c->max_tiles = strtoull(optarg, NULL, 10); break;
            case 's': c->start_tile = strtoull(optarg, NULL, 10); break;
            case 'E': c->end_tile   = strtoull(optarg, NULL, 10); break;
            case 'm': c->min_report_len = atoi(optarg); break;
            case 'O':
                if (strcmp(optarg, "forward") == 0) c->prove_forward = 1;
                else if (strcmp(optarg, "reverse") == 0) c->prove_forward = 0;
                else { fprintf(stderr, "ERROR: --prove-order must be forward|reverse\n"); return -1; }
                break;
            case 'P': c->cpu_prove = 1; break;
            case 'q': c->prove_threads = atoi(optarg); break;
            case 1030: {
                int v = atoi(optarg);
                if (v < 1 || v > 64) {
                    fprintf(stderr,
                        "ERROR: --prove-queue-depth must be in [1,64] (got %d)\n", v);
                    return -1;
                }
                c->prove_queue_depth = v;
                break;
            }
            case 'u': {
                int n = parse_csv_primes(optarg, c->pool_primes, V16_MAX_IMMUNE);
                if (n < 0) return -1;
                c->n_pool_primes = n;
                break;
            }
            case 'Q':
                if (parse_seed_pool(optarg, c) != 0) return -1;
                break;
            case 1005:
                if (parse_seed_pool_file(optarg, c) != 0) return -1;
                break;
            case 'Z': c->sieve_variants = 1; break;
            case 1001: c->q_iter = 1; break;
            case 1002:
                if (strcmp(optarg, "sequential") == 0) c->q_order_random = 0;
                else if (strcmp(optarg, "random") == 0) c->q_order_random = 1;
                else {
                    fprintf(stderr,
                        "ERROR: --q-order must be sequential|random\n");
                    return -1;
                }
                break;
            case 1003: c->q_seed = strtoull(optarg, NULL, 10); break;
            case 1004: {
                unsigned long long v = strtoull(optarg, NULL, 10);
                if (v == 0 || v > (unsigned long long)UINT32_MAX) {
                    fprintf(stderr,
                        "ERROR: --q-chunk must be in (0, 2^32)\n");
                    return -1;
                }
                c->q_chunk = (u32)v;
                break;
            }
            case 1060:   /* --q-band-mode */
                if (strcmp(optarg, "fixed") == 0) c->q_band_mode_exhaustive = 0;
                else if (strcmp(optarg, "exhaustive") == 0) c->q_band_mode_exhaustive = 1;
                else {
                    fprintf(stderr,
                        "ERROR: --q-band-mode must be fixed|exhaustive\n");
                    return -1;
                }
                break;
            case 1061: {  /* --exhaustive-max-q-bits */
                int n = atoi(optarg);
                if (n < 1 || n > 100) {
                    fprintf(stderr,
                        "ERROR: --exhaustive-max-q-bits must be in [1,100] (got %d)\n", n);
                    return -1;
                }
                c->exhaustive_max_q_bits = n;
                break;
            }
            case 1062: c->dedup_p0 = 1; break;  /* --dedup-p0 (default) */
            case 1063: c->dedup_p0 = 0; break;  /* --no-dedup-p0 */
            case 1006: c->bits_min = atoi(optarg); break;
            case 1007: c->bits_max = atoi(optarg); break;
            case 1008: c->bits_rotate_step       = atoi(optarg); break;
            case 1009: c->bits_rotate_period_sec = atoi(optarg); break;
            case 1010: c->k15_enabled = 1; break;
            case 1011: c->k15_enabled = 0; break;
            case 1012: {
                int r = atoi(optarg);
                if (r != 2 && r != 3) {
                    fprintf(stderr,
                        "ERROR: --k15-mr-reps must be 2 or 3 (got %d)\n", r);
                    return -1;
                }
                c->k15_mr_reps = r;
                break;
            }
            case 1020: c->active_primes = 1; break;
            case 1021: c->active_primes = 0; break;
            case 1022: {
                int n = atoi(optarg);
                if (n < 2 || n > V16_MAX_STREAMS) {
                    fprintf(stderr,
                        "ERROR: --streams=%d out of range [2,%d]\n",
                        n, V16_MAX_STREAMS);
                    return -1;
                }
                c->n_streams = n;
                break;
            }
            case 1040: c->run_self_test = 1; break;  /* --test (explicit opt-IN; default) */
            case 1041: c->run_self_test = 0; break;  /* --no-test (opt-OUT) */
            case 1050: c->use_ext_sieve = 1; break;  /* --ext-sieve (explicit opt-IN; default) */
            case 1051: c->use_ext_sieve = 0; break;  /* --no-ext-sieve (opt-OUT) */
            case 'g': c->gpu_id = atoi(optarg); break;
            case 'L': c->time_limit_sec = atoi(optarg); break;
            case 'r': c->report_sec = atoi(optarg); break;
            case 'C': strncpy(c->checkpoint_file, optarg, sizeof(c->checkpoint_file) - 1); break;
            case 'R': strncpy(c->resume_file, optarg, sizeof(c->resume_file) - 1); break;
            case 'l': strncpy(c->log_file, optarg, sizeof(c->log_file) - 1); break;
            case 'h': usage(argv[0]); exit(0);
            default: usage(argv[0]); return -1;
        }
    }

    /* Seed-name aliases. */
    if (c->immune_prime_cap == 0 && c->n_immune_set == 0) {
        if (strcmp(c->seed_name, "41hash_over_31") == 0) {
            c->immune_prime_cap = 41;
            c->omit_primes[c->n_omit++] = 31;
        } else if (strcmp(c->seed_name, "custom_664870") == 0) {
            u32 set[] = {2,3,5,7,11,13,17,19,29,37,229};
            int n = (int)(sizeof(set) / sizeof(set[0]));
            for (int i = 0; i < n; i++) c->immune_set[c->n_immune_set++] = set[i];
        } else if (strcmp(c->seed_name, "nm_v15") == 0) {
            /* v15 lattice: M = 5*7*11*13*17*19 = 1,616,615.
             * 2 is always immune in first-kind; include explicitly. 3 was
             * a CRT base in v15 (not part of M) but is similarly always
             * immune in practice (chain density forces 3 | p+1 for any
             * chain length where ord(2,3)=2 ≤ target). Including {2,3} in
             * D matches v15's effective constraint and lets the coeff
             * wheel cover {23,29,31,...} like v15's WHEEL_PERIOD.
             * Yield baseline: matches v15's "CC10+ every 10-15s" rate. */
            u32 set[] = {2, 3, 5, 7, 11, 13, 17, 19};
            int n = (int)(sizeof(set) / sizeof(set[0]));
            for (int i = 0; i < n; i++) c->immune_set[c->n_immune_set++] = set[i];
        } else if (strcmp(c->seed_name, "nm_cc10_base") == 0 ||
                   strcmp(c->seed_name, "unnamed") == 0) {
            /* Data-driven default: 52% of CC10+ roots have fingerprint
             * {3,5,11,13,19}; 98.9% have it as a subset. (2 is always immune
             * in first-kind since p+1 is even — include it explicitly so the
             * candidate generator doesn't have to constrain parity.)
             * See docs/v16_design_note.md and scripts/verify_fingerprints.py.
             * Default when no --immune-* is given. */
            u32 set[] = {2, 3, 5, 11, 13, 19};
            int n = (int)(sizeof(set) / sizeof(set[0]));
            for (int i = 0; i < n; i++) c->immune_set[c->n_immune_set++] = set[i];
            if (c->seed_name[0] == 'u') {        /* "unnamed" → relabel */
                strncpy(c->seed_name, "nm_cc10_base", sizeof(c->seed_name) - 1);
                c->seed_name[sizeof(c->seed_name) - 1] = '\0';
            }
        }
    }

    if (c->target_len < 1 || c->target_len > V16_MAX_TARGET) {
        fprintf(stderr, "ERROR: --target must be in [1,%d]\n", V16_MAX_TARGET);
        return -1;
    }
    /* depth is independent of target (design §15). The 0 sentinel applies
     * the silent default = 12, chosen so the common case is not silently
     * over-filtering at high target (e.g. target=18). All three orderings
     * depth < target, ==, > target are legal use cases — explicitly do not
     * reject depth < target here. The depth>target startup WARN fires after
     * banner. Both q-iter and legacy --sieve-variants paths fall through
     * this single check (q-iter validation is below at line 2256+). */
    if (c->sieve_depth == 0) c->sieve_depth = 12;
    if (c->sieve_depth < 1) {
        fprintf(stderr, "ERROR: --depth must be >= 1 (got %d)\n",
                c->sieve_depth);
        return -1;
    }
    if (c->sieve_depth > V16_MAX_TARGET) {
        fprintf(stderr, "ERROR: --depth must be <= %d\n", V16_MAX_TARGET);
        return -1;
    }
    /* --sieve-variants requires --cpu-prove (GPU K2 stays single-D for now)
     * and at least one --seed-pool entry. The seed-pool count is rechecked
     * after expand_auto_pools(); here we only enforce the presence of at
     * least one spec on the CLI. */
    if (c->sieve_variants) {
        if (!c->cpu_prove) {
            fprintf(stderr,
                "ERROR: --sieve-variants requires --cpu-prove\n");
            return -1;
        }
        if (c->n_pool_seeds == 0 && c->n_pool_auto_specs == 0) {
            fprintf(stderr,
                "ERROR: --sieve-variants requires at least one --seed-pool entry\n");
            return -1;
        }
    }
    if (c->q_iter) {
        if (c->sieve_variants) {
            fprintf(stderr,
                "ERROR: --q-iter is incompatible with --sieve-variants\n");
            return -1;
        }
        if (!c->cpu_prove) {
            fprintf(stderr,
                "ERROR: --q-iter requires --cpu-prove\n");
            return -1;
        }
        if (c->n_pool_seeds == 0 && c->n_pool_auto_specs == 0) {
            fprintf(stderr,
                "ERROR: --q-iter requires at least one --seed-pool entry\n");
            return -1;
        }
        /* --bits N collapses to [N,N]; otherwise band defaults [90,127]
         * apply (design §17 record zone). target_bits stays for back-compat
         * checkpoint writes; treat it as bits_max when set. */
        if (c->bits_max > 127) {
            fprintf(stderr,
                "ERROR: --bits-max %d exceeds u128 ceiling 127 "
                "(design doc §17.2)\n", c->bits_max);
            return -1;
        }
        if (c->bits_min > c->bits_max) {
            fprintf(stderr,
                "ERROR: --bits-min %d must be <= --bits-max %d\n",
                c->bits_min, c->bits_max);
            return -1;
        }
        if (c->bits_min < 8) {
            fprintf(stderr,
                "WARN: --bits-min %d < 8 produces tiny candidates\n",
                c->bits_min);
        }
        /* --q-order random: PRNG seeded later (after gpu_id parsed). */
    }
    /* K1.5 computes
     *   p_top = Q · D_V · 2^(E + top_pos) − 1
     * in u128 (gu128_mul_u128 truncates above 128 bits). The precise
     * overflow guard runs after pool expansion in main() once max(D_V)
     * is known — see "K1.5 overflow guard (Bug 1b)" below. Only the
     * "k15 requires q-iter" sanity check stays here. */
    if (c->k15_enabled && !c->q_iter) {
        fprintf(stderr,
            "WARN: --k15 has no effect without --q-iter; ignoring.\n");
        c->k15_enabled = 0;
    }
    /* --sieve-max used to silently clamp to the largest
     * prime in V16_SMALL_PRIMES[]. The table is static and finite, so passing
     * a cap above the table's last entry just truncates at table end with no
     * indication. Make this an explicit, immediate failure. */
    {
        u32 table_max = V16_SMALL_PRIMES[V16_SMALL_PRIME_COUNT - 1];
        if (c->sieve_prime_max > table_max) {
            fprintf(stderr,
                "ERROR: --sieve-max=%u exceeds compiled small-prime table "
                "max (last prime = %u). Reduce --sieve-max or extend "
                "V16_SMALL_PRIMES.\n",
                c->sieve_prime_max, table_max);
            return -1;
        }
    }
    return 0;
}

/* ============================================================================
 * Host-side validation mode.
 *   --validate-seed [--seed-name NAME] [--target N] [--exp-start E] [--validate-A A]
 * ============================================================================ */
static int validate_seed(v16_cfg *c) {
    printf("=== v16 validate-seed ===\n");
    printf("  seed_name = %s\n", c->seed_name);
    printf("  target    = %d\n", c->target_len);
    printf("  exp_start = %d\n", c->exp_start);
    printf("  D primes  = [");
    for (int i = 0; i < c->n_D_primes; i++)
        printf("%s%u", i ? "," : "", c->D_primes[i]);
    printf("]\n");
    printf("  D         = %s\n", c->D_dec);

    /* Default A by seed name. */
    char a_buf[80] = "";
    if (c->validate_A_str[0]) {
        memcpy(a_buf, c->validate_A_str, sizeof(a_buf) - 1);
        a_buf[sizeof(a_buf) - 1] = '\0';
    } else if (strcmp(c->seed_name, "41hash_over_31") == 0) {
        strcpy(a_buf, "2600776728061305571");
    } else if (strcmp(c->seed_name, "custom_664870") == 0) {
        strcpy(a_buf, "139480563822237");        /* derived in design note */
    }
    if (!a_buf[0]) {
        printf("  (no default A for this seed; pass --validate-A NUM to test primality)\n");
        return 0;
    }
    printf("  A         = %s\n", a_buf);

    mpz_t A, D, p, two, top, predec, qi;
    mpz_inits(A, D, p, two, top, predec, qi, NULL);
    mpz_set_str(A, a_buf, 10);
    char D_dec[80];
    u128_to_dec(c->D, D_dec, sizeof(D_dec));
    mpz_set_str(D, D_dec, 10);
    mpz_set_ui(two, 2);

    /* Custom-seed extra check: verify divisibility of
     * 664870017491412573885064019 + 1 by the immune set. */
    if (strcmp(c->seed_name, "custom_664870") == 0) {
        mpz_t known, k_p1;
        mpz_inits(known, k_p1, NULL);
        mpz_set_str(known, "664870017491412573885064019", 10);
        mpz_add_ui(k_p1, known, 1);
        printf("  custom_664870 known root + 1 = ");
        mpz_out_str(stdout, 10, k_p1); printf("\n");
        int all_ok = 1;
        for (int i = 0; i < c->n_D_primes; i++) {
            int divides = mpz_divisible_ui_p(k_p1, c->D_primes[i]);
            printf("    p+1 %% %u = %s\n", c->D_primes[i], divides ? "0 [OK]" : "NONZERO [FAIL]");
            if (!divides) all_ok = 0;
        }
        mpz_clears(known, k_p1, NULL);
        if (!all_ok) { mpz_clears(A, D, p, two, top, predec, qi, NULL); return -1; }
    }

    /* Walk the chain. */
    int chain_actual = 0;
    for (int j = 0; j < c->target_len; j++) {
        /* p_j = A * D * 2^(E+j) - 1 */
        mpz_set(p, A);
        mpz_mul(p, p, D);
        mpz_mul_2exp(p, p, (unsigned long)(c->exp_start + j));
        mpz_sub_ui(p, p, 1);
        int is_prime = mpz_probab_prime_p(p, 25) > 0;
        size_t nbits = mpz_sizeinbase(p, 2);
        printf("    j=%2d  bits=%3zu  prime=%s  hex=0x", j, nbits, is_prime ? "Y" : "N");
        mpz_out_str(stdout, 16, p); printf("\n");
        if (is_prime) chain_actual = j + 1;
        if (!is_prime) break;
    }

    /* Root check: (A*D*2^E - 1 - 1) / 2 = A*D*2^(E-1) - 1 if E>=1 should be composite to be root. */
    if (chain_actual >= 1) {
        mpz_set(p, A);
        mpz_mul(p, p, D);
        if (c->exp_start >= 1) mpz_mul_2exp(p, p, (unsigned long)(c->exp_start - 1));
        mpz_sub_ui(p, p, 1);
        if (c->exp_start >= 1 && mpz_cmp_ui(p, 2) >= 0) {
            int pre = mpz_probab_prime_p(p, 25) > 0;
            printf("    predecessor (E-1=%d) is %s — root %s\n",
                   c->exp_start - 1,
                   pre ? "PRIME" : "composite",
                   pre ? "EXTENDS backward (not minimal)" : "OK (minimal start)");
        }
    }

    printf("  chain_actual = %d / target %d  %s\n",
           chain_actual, c->target_len,
           chain_actual >= c->target_len ? "[PASS]" : "[INCOMPLETE]");
    mpz_clears(A, D, p, two, top, predec, qi, NULL);
    return (chain_actual >= c->target_len) ? 0 : -1;
}

/* ============================================================================
 * Checkpoint format (text, key=value lines).
 * ============================================================================ */
static int save_checkpoint(const v16_cfg *c, u64 a_tile_cursor, u64 wheel_hash) {
    if (!c->checkpoint_file[0]) return 0;
    /* Atomic write: same temp+rename pattern as save_q_iter_checkpoint. */
    char tmp_path[sizeof(c->checkpoint_file) + 8];
    {
        size_t n = strlen(c->checkpoint_file);
        if (n + 5 >= sizeof(tmp_path)) {
            fprintf(stderr,
                "ERROR: checkpoint path too long for atomic-write tmp\n");
            return -1;
        }
        memcpy(tmp_path, c->checkpoint_file, n);
        memcpy(tmp_path + n, ".tmp", 5);
    }
    FILE *f = fopen(tmp_path, "w");
    if (!f) { perror("checkpoint write (tmp)"); return -1; }
    int werr = 0;
    #define V16_CKPT_LEGACY_FPRINTF(...) do { \
        if (!werr && fprintf(f, __VA_ARGS__) < 0) werr = 1; \
    } while (0)
    V16_CKPT_LEGACY_FPRINTF("version=%s\n", V16_VERSION);
    V16_CKPT_LEGACY_FPRINTF("target_len=%d\n", c->target_len);
    V16_CKPT_LEGACY_FPRINTF("sieve_depth=%d\n", c->sieve_depth);
    V16_CKPT_LEGACY_FPRINTF("target_bits=%d\n", c->target_bits);
    V16_CKPT_LEGACY_FPRINTF("seed_name=%s\n", c->seed_name);
    V16_CKPT_LEGACY_FPRINTF("D_dec=%s\n", c->D_dec);
    V16_CKPT_LEGACY_FPRINTF("D_primes=");
    for (int i = 0; i < c->n_D_primes; i++)
        V16_CKPT_LEGACY_FPRINTF("%s%u", i ? "," : "", c->D_primes[i]);
    V16_CKPT_LEGACY_FPRINTF("\n");
    V16_CKPT_LEGACY_FPRINTF("exp_start=%d\n", c->exp_start);
    V16_CKPT_LEGACY_FPRINTF("wheel_prime_max=%u\n", c->wheel_prime_max);
    V16_CKPT_LEGACY_FPRINTF("sieve_prime_max=%u\n", c->sieve_prime_max);
    V16_CKPT_LEGACY_FPRINTF("wheel_fnv1a64=0x%016llx\n",
                            (unsigned long long)wheel_hash);
    V16_CKPT_LEGACY_FPRINTF("a_tile_cursor=%llu\n",
                            (unsigned long long)a_tile_cursor);
    #undef V16_CKPT_LEGACY_FPRINTF
    if (werr || fflush(f) != 0) {
        fclose(f);
        unlink(tmp_path);
        fprintf(stderr,
                "ERROR: checkpoint write failed (tmp=%s) — original "
                "checkpoint untouched\n", tmp_path);
        return -1;
    }
    int fd = fileno(f);
    if (fd >= 0) (void)fsync(fd);
    if (fclose(f) != 0) {
        unlink(tmp_path);
        fprintf(stderr,
                "ERROR: checkpoint fclose failed — original untouched\n");
        return -1;
    }
    if (rename(tmp_path, c->checkpoint_file) != 0) {
        perror("checkpoint rename");
        unlink(tmp_path);
        return -1;
    }
    return 0;
}

static int load_checkpoint_cursor(const v16_cfg *c, u64 wheel_hash,
                                  u64 *cursor_out) {
    if (!c->resume_file[0]) { *cursor_out = 0; return 0; }
    FILE *f = fopen(c->resume_file, "r");
    if (!f) {
        fprintf(stderr, "WARN: resume file %s not found — starting at 0\n",
                c->resume_file);
        *cursor_out = 0;
        return 0;
    }
    char line[256];
    char rd_seed[128] = "";
    char rd_D[80] = "";
    char rd_mode[32] = "";
    int rd_target = -1, rd_depth = -1, rd_E = -1;
    u32 rd_wmax = 0, rd_smax = 0;
    u64 rd_hash = 0, rd_cursor = 0;
    while (fgets(line, sizeof(line), f)) {
        char *eq = strchr(line, '=');
        if (!eq) continue;
        *eq = '\0';
        char *val = eq + 1;
        size_t vl = strlen(val);
        if (vl && val[vl - 1] == '\n') val[vl - 1] = '\0';
        if (!strcmp(line, "seed_name")) strncpy(rd_seed, val, sizeof(rd_seed) - 1);
        else if (!strcmp(line, "D_dec")) strncpy(rd_D, val, sizeof(rd_D) - 1);
        else if (!strcmp(line, "mode")) strncpy(rd_mode, val, sizeof(rd_mode) - 1);
        else if (!strcmp(line, "target_len")) rd_target = atoi(val);
        else if (!strcmp(line, "sieve_depth")) rd_depth = atoi(val);
        else if (!strcmp(line, "exp_start"))  rd_E = atoi(val);
        else if (!strcmp(line, "wheel_prime_max")) rd_wmax = (u32)strtoul(val, NULL, 10);
        else if (!strcmp(line, "sieve_prime_max")) rd_smax = (u32)strtoul(val, NULL, 10);
        else if (!strcmp(line, "wheel_fnv1a64"))   rd_hash = strtoull(val, NULL, 0);
        else if (!strcmp(line, "a_tile_cursor"))   rd_cursor = strtoull(val, NULL, 10);
    }
    fclose(f);
    /* Mode discriminator: refuse to load q-iter checkpoint into the
     * non-q-iter path (and vice versa, handled in load_q_iter_checkpoint). */
    if (rd_mode[0] && strcmp(rd_mode, "q-iter") == 0) {
        fprintf(stderr,
            "ERROR: incompatible checkpoint mode (file=q-iter, requested=non-q-iter). "
            "Add --q-iter or use a different resume file.\n");
        return -1;
    }
    if (strcmp(rd_seed, c->seed_name) ||
        strcmp(rd_D, c->D_dec) ||
        rd_target != c->target_len ||
        (rd_depth > 0 && rd_depth != c->sieve_depth) ||
        rd_E != c->exp_start ||
        rd_wmax != c->wheel_prime_max ||
        rd_smax != c->sieve_prime_max ||
        rd_hash != wheel_hash) {
        fprintf(stderr, "ERROR: resume identity mismatch — refusing to resume\n");
        fprintf(stderr, "       seed=%s want=%s\n", rd_seed, c->seed_name);
        fprintf(stderr, "       D=%s want=%s\n", rd_D, c->D_dec);
        fprintf(stderr, "       target=%d want=%d  E=%d want=%d\n",
                rd_target, c->target_len, rd_E, c->exp_start);
        fprintf(stderr, "       wmax=%u want=%u  smax=%u want=%u\n",
                rd_wmax, c->wheel_prime_max, rd_smax, c->sieve_prime_max);
        fprintf(stderr, "       wheel_hash=0x%016llx want=0x%016llx\n",
                (unsigned long long)rd_hash, (unsigned long long)wheel_hash);
        return -1;
    }
    *cursor_out = rd_cursor;
    return 0;
}

/* FNV-1a 64-bit over an arbitrary byte buffer (used by Q-iter checkpoint
 * identity fields: sieve_primes hash and pool_seeds hash). */
static u64 fnv1a64_bytes(const void *buf, size_t n) {
    u64 h = 0xcbf29ce484222325ULL;
    const u64 prime = 0x100000001b3ULL;
    const u8 *p = (const u8 *)buf;
    for (size_t i = 0; i < n; i++) {
        h ^= p[i];
        h *= prime;
    }
    return h;
}

/* FNV-1a 64-bit hash of the sieve-prime list (u32 array). */
static u64 q_iter_sieve_fnv1a64(const u32 *primes, int n) {
    return fnv1a64_bytes(primes, (size_t)n * sizeof(u32));
}

/* FNV-1a 64-bit hash of the pool-seed list. Hashes the concatenation of
 * each seed's D_dec string + a comma separator. Deterministic across runs
 * with the same pool spec. */
static u64 q_iter_pool_fnv1a64(const v16_pool_seed *seeds, int n) {
    u64 h = 0xcbf29ce484222325ULL;
    const u64 prime = 0x100000001b3ULL;
    for (int i = 0; i < n; i++) {
        const char *s = seeds[i].D_dec;
        for (size_t k = 0; s[k]; k++) { h ^= (u8)s[k]; h *= prime; }
        h ^= (u8)','; h *= prime;
    }
    return h;
}

/* Q-iter checkpoint save. Sparse per-V cursor emit: only V's with cursor
 * advanced past Q_min are written.
 *
 * NOTE: target_bits IS written (informational), but NOT identity-checked
 * on load (user may resume with a different --bits; pool fingerprint
 * is the coverage discriminator).
 *
 * v_cursor is u128 (carried as paired lo/hi u64). Serialized as a
 * decimal string per V so the format handles Q > 2^63 (full record zone).
 */
static int save_q_iter_checkpoint(const v16_cfg *c,
                                  const u64 *v_cursor_lo,
                                  const u64 *v_cursor_hi,
                                  const u64 *v_qmin_lo,
                                  const u64 *v_qmin_hi,
                                  int n_v,
                                  u64 sieve_hash,
                                  u64 pool_hash) {
    if (!c->checkpoint_file[0]) return 0;
    /* Atomic write. Old code fopen'd the live file
     * with "w" — a SIGKILL / OOM / GPU exception mid-write left the file
     * truncated. Strategy: write to "<file>.tmp", fflush+fsync+fclose, then
     * rename(2). rename(2) is atomic on POSIX filesystems for paths on the
     * same mount, so either the old file or the fully-written new file is
     * visible — never a half-written prefix. On any write failure we unlink
     * the tmp and leave the original untouched. */
    char tmp_path[sizeof(c->checkpoint_file) + 8];
    {
        size_t n = strlen(c->checkpoint_file);
        if (n + 5 >= sizeof(tmp_path)) {
            fprintf(stderr,
                "ERROR: checkpoint path too long for atomic-write tmp\n");
            return -1;
        }
        memcpy(tmp_path, c->checkpoint_file, n);
        memcpy(tmp_path + n, ".tmp", 5);  /* incl. trailing \0 */
    }
    FILE *f = fopen(tmp_path, "w");
    if (!f) { perror("checkpoint write (tmp)"); return -1; }
    /* Write-and-track: every fprintf return is checked; on any error we
     * unlink the tmp and bail without touching the live file. */
    int werr = 0;
    #define V16_CKPT_FPRINTF(...) do { \
        if (!werr && fprintf(f, __VA_ARGS__) < 0) werr = 1; \
    } while (0)
    V16_CKPT_FPRINTF("version=%s\n", V16_VERSION_Q_ITER);
    V16_CKPT_FPRINTF("mode=q-iter\n");
    V16_CKPT_FPRINTF("target_len=%d\n", c->target_len);
    V16_CKPT_FPRINTF("sieve_depth=%d\n", c->sieve_depth);
    V16_CKPT_FPRINTF("target_bits=%d\n", c->target_bits);
    V16_CKPT_FPRINTF("bits_min=%d\n", c->bits_min);
    V16_CKPT_FPRINTF("bits_max=%d\n", c->bits_max);
    V16_CKPT_FPRINTF("seed_name=%s\n", c->seed_name);
    V16_CKPT_FPRINTF("D_dec=%s\n", c->D_dec);
    V16_CKPT_FPRINTF("D_primes=");
    for (int i = 0; i < c->n_D_primes; i++)
        V16_CKPT_FPRINTF("%s%u", i ? "," : "", c->D_primes[i]);
    V16_CKPT_FPRINTF("\n");
    V16_CKPT_FPRINTF("exp_start=%d\n", c->exp_start);
    V16_CKPT_FPRINTF("sieve_prime_max=%u\n", c->sieve_prime_max);
    V16_CKPT_FPRINTF("sieve_primes=");
    for (int i = 0; i < c->n_sieve_primes; i++)
        V16_CKPT_FPRINTF("%s%u", i ? "," : "", c->sieve_primes[i]);
    V16_CKPT_FPRINTF("\n");
    V16_CKPT_FPRINTF("sieve_fnv1a64=0x%016llx\n",
                     (unsigned long long)sieve_hash);
    V16_CKPT_FPRINTF("pool_seeds_count=%d\n", n_v);
    V16_CKPT_FPRINTF("pool_seeds_fnv1a64=0x%016llx\n",
                     (unsigned long long)pool_hash);
    V16_CKPT_FPRINTF("q_order=%s\n",
                     c->q_order_random ? "random" : "sequential");
    if (c->q_order_random)
        V16_CKPT_FPRINTF("q_seed=%llu\n", (unsigned long long)c->q_seed);
    V16_CKPT_FPRINTF("q_chunk=%u\n", c->q_chunk);
    /* Persist q-band-mode so a resume can't silently switch modes
     * (fixed↔exhaustive flips the per-V coverage semantics). */
    V16_CKPT_FPRINTF("q_band_mode=%s\n",
                     c->q_band_mode_exhaustive ? "exhaustive" : "fixed");
    if (c->q_band_mode_exhaustive)
        V16_CKPT_FPRINTF("exhaustive_max_q_bits=%d\n", c->exhaustive_max_q_bits);
    /* Sparse drained-cursor emit: only V's that have advanced past Q_min(V).
     * Decimal u128.
     *
     * Field is now v_cursor_drained[v], the *rear*
     * edge — the largest Q for which all survivors have been CPU-proved.
     * The loader accepts the new field name; old-style "v_cursor[v]" is
     * rejected via the version bump (V16_VERSION_Q_ITER). */
    for (int v = 0; v < n_v; v++) {
        u128 cur = ((u128)v_cursor_hi[v] << 64) | (u128)v_cursor_lo[v];
        u128 lo  = ((u128)v_qmin_hi[v]   << 64) | (u128)v_qmin_lo[v];
        if (cur > lo) {
            char dec[64];
            u128_to_dec(cur, dec, sizeof(dec));
            V16_CKPT_FPRINTF("v_cursor_drained[%d]=%s\n", v, dec);
        }
    }
    #undef V16_CKPT_FPRINTF
    if (werr || fflush(f) != 0) {
        fclose(f);
        unlink(tmp_path);
        fprintf(stderr,
                "ERROR: checkpoint write failed (tmp=%s) — original "
                "checkpoint untouched\n", tmp_path);
        return -1;
    }
    /* fsync best-effort: durability of the tmp content before rename. Errors
     * here are non-fatal (filesystem may not support it, e.g. tmpfs). */
    int fd = fileno(f);
    if (fd >= 0) (void)fsync(fd);
    if (fclose(f) != 0) {
        unlink(tmp_path);
        fprintf(stderr,
                "ERROR: checkpoint fclose failed — original untouched\n");
        return -1;
    }
    if (rename(tmp_path, c->checkpoint_file) != 0) {
        perror("checkpoint rename");
        unlink(tmp_path);
        return -1;
    }
    return 0;
}

/* Parse a decimal u128 string into (lo, hi). Returns 0 on success. */
static int parse_u128_dec(const char *s, u64 *lo, u64 *hi) {
    u128 v = 0;
    for (; *s; s++) {
        if (*s < '0' || *s > '9') {
            if (*s == '\n' || *s == '\r' || *s == ' ' || *s == '\t') break;
            return -1;
        }
        v = v * 10 + (u128)(*s - '0');
    }
    *lo = (u64)v;
    *hi = (u64)(v >> 64);
    return 0;
}

/* Q-iter checkpoint load. Populates v_cursor[]; missing V's stay at
 * Q_min(V) (caller pre-fills with q_min before calling). Identity-checks
 * version/mode/target_len/sieve_depth/seed_name/D_dec/exp_start/
 * sieve_prime_max/sieve_fnv1a64/pool_seeds_count/pool_seeds_fnv1a64/q_order.
 * Returns 0 on success/no-resume-file, -1 on identity mismatch.
 * Cursors are u128 (paired lo/hi). */
static int load_q_iter_checkpoint(const v16_cfg *c,
                                  u64 *v_cursor_lo,
                                  u64 *v_cursor_hi,
                                  int n_v,
                                  u64 sieve_hash,
                                  u64 pool_hash) {
    if (!c->resume_file[0]) return 0;
    FILE *f = fopen(c->resume_file, "r");
    if (!f) {
        fprintf(stderr, "WARN: resume file %s not found — starting fresh\n",
                c->resume_file);
        return 0;
    }
    char line[1024];
    char rd_version[64] = "";
    char rd_mode[32] = "";
    char rd_seed[128] = "";
    char rd_D[80] = "";
    char rd_qorder[16] = "";
    char rd_qbandmode[16] = "";
    int rd_target = -1, rd_depth = -1, rd_E = -1;
    int rd_bits_min = -1, rd_bits_max = -1;
    int rd_pool_count = -1;
    int rd_exhaustive_max_q_bits = -1;
    u32 rd_smax = 0;
    u64 rd_sieve_hash = 0, rd_pool_hash = 0;
    while (fgets(line, sizeof(line), f)) {
        char *eq = strchr(line, '=');
        if (!eq) continue;
        *eq = '\0';
        char *val = eq + 1;
        size_t vl = strlen(val);
        if (vl && val[vl - 1] == '\n') val[vl - 1] = '\0';
        if (!strcmp(line, "version")) strncpy(rd_version, val,
                                              sizeof(rd_version) - 1);
        else if (!strcmp(line, "mode")) strncpy(rd_mode, val,
                                                sizeof(rd_mode) - 1);
        else if (!strcmp(line, "seed_name")) strncpy(rd_seed, val,
                                                     sizeof(rd_seed) - 1);
        else if (!strcmp(line, "D_dec")) strncpy(rd_D, val,
                                                 sizeof(rd_D) - 1);
        else if (!strcmp(line, "target_len")) rd_target = atoi(val);
        else if (!strcmp(line, "sieve_depth")) rd_depth = atoi(val);
        else if (!strcmp(line, "exp_start")) rd_E = atoi(val);
        else if (!strcmp(line, "bits_min"))  rd_bits_min = atoi(val);
        else if (!strcmp(line, "bits_max"))  rd_bits_max = atoi(val);
        else if (!strcmp(line, "sieve_prime_max"))
            rd_smax = (u32)strtoul(val, NULL, 10);
        else if (!strcmp(line, "sieve_fnv1a64"))
            rd_sieve_hash = strtoull(val, NULL, 0);
        else if (!strcmp(line, "pool_seeds_count")) rd_pool_count = atoi(val);
        else if (!strcmp(line, "pool_seeds_fnv1a64"))
            rd_pool_hash = strtoull(val, NULL, 0);
        else if (!strcmp(line, "q_order")) strncpy(rd_qorder, val,
                                                   sizeof(rd_qorder) - 1);
        else if (!strcmp(line, "q_band_mode"))
            strncpy(rd_qbandmode, val, sizeof(rd_qbandmode) - 1);
        else if (!strcmp(line, "exhaustive_max_q_bits"))
            rd_exhaustive_max_q_bits = atoi(val);
        else if (!strncmp(line, "v_cursor_drained[", 17)) {
            /* Field renamed from v_cursor[] to mark
             * that this is the *drained* (CPU-proved) rear edge, not the
             * launched front edge. Old format is rejected by version
             * check above. */
            int idx = atoi(line + 17);
            u64 cur_lo = 0, cur_hi = 0;
            if (parse_u128_dec(val, &cur_lo, &cur_hi) != 0) {
                fprintf(stderr,
                    "WARN: q-iter checkpoint: bad v_cursor_drained[%d]=%s\n",
                    idx, val);
                continue;
            }
            if (idx >= 0 && idx < n_v) {
                v_cursor_lo[idx] = cur_lo;
                v_cursor_hi[idx] = cur_hi;
            }
        }
        /* q_chunk, q_seed (when sequential), target_bits: NOT
         * identity-checked (per design §8). */
    }
    fclose(f);
    /* Mode discriminator. */
    if (strcmp(rd_mode, "q-iter") != 0) {
        fprintf(stderr,
            "ERROR: incompatible checkpoint mode (file=%s, requested=q-iter). "
            "Drop --q-iter or use a different resume file.\n",
            rd_mode[0] ? rd_mode : "non-q-iter");
        return -1;
    }
    if (strcmp(rd_version, V16_VERSION_Q_ITER) != 0) {
        fprintf(stderr,
            "ERROR: Q-iter checkpoint version mismatch (file=%s want=%s)\n",
            rd_version, V16_VERSION_Q_ITER);
        return -1;
    }
    /* Bit-band identity check. Mismatch on either edge rejects resume.
     * Allow resume from older checkpoints (rd_bits_* == -1). */
    int bits_band_mismatch =
        (rd_bits_min >= 0 && rd_bits_min != c->bits_min) ||
        (rd_bits_max >= 0 && rd_bits_max != c->bits_max);
    if (strcmp(rd_seed, c->seed_name) ||
        strcmp(rd_D, c->D_dec) ||
        rd_target != c->target_len ||
        (rd_depth > 0 && rd_depth != c->sieve_depth) ||
        rd_E != c->exp_start ||
        bits_band_mismatch ||
        rd_smax != c->sieve_prime_max ||
        rd_sieve_hash != sieve_hash ||
        rd_pool_count != n_v ||
        rd_pool_hash != pool_hash ||
        strcmp(rd_qorder,
               c->q_order_random ? "random" : "sequential") ||
        /* q-band-mode mismatch rejects resume. Empty rd_qbandmode
         * means an older checkpoint pre-dating the field — accept that. */
        (rd_qbandmode[0] && strcmp(rd_qbandmode,
               c->q_band_mode_exhaustive ? "exhaustive" : "fixed")) ||
        (c->q_band_mode_exhaustive &&
         rd_exhaustive_max_q_bits > 0 &&
         rd_exhaustive_max_q_bits != c->exhaustive_max_q_bits)) {
        fprintf(stderr,
            "ERROR: Q-iter resume identity mismatch — refusing to resume\n");
        fprintf(stderr, "       seed=%s want=%s\n", rd_seed, c->seed_name);
        fprintf(stderr, "       D=%s want=%s\n", rd_D, c->D_dec);
        fprintf(stderr, "       target=%d want=%d  depth=%d want=%d  E=%d want=%d\n",
                rd_target, c->target_len, rd_depth, c->sieve_depth,
                rd_E, c->exp_start);
        fprintf(stderr, "       bits=[%d,%d] want=[%d,%d]\n",
                rd_bits_min, rd_bits_max, c->bits_min, c->bits_max);
        fprintf(stderr, "       smax=%u want=%u\n",
                rd_smax, c->sieve_prime_max);
        fprintf(stderr, "       sieve_hash=0x%016llx want=0x%016llx\n",
                (unsigned long long)rd_sieve_hash,
                (unsigned long long)sieve_hash);
        fprintf(stderr, "       pool_count=%d want=%d  pool_hash=0x%016llx want=0x%016llx\n",
                rd_pool_count, n_v,
                (unsigned long long)rd_pool_hash,
                (unsigned long long)pool_hash);
        fprintf(stderr, "       q_order=%s want=%s\n",
                rd_qorder, c->q_order_random ? "random" : "sequential");
        fprintf(stderr, "       q_band_mode=%s want=%s\n",
                rd_qbandmode[0] ? rd_qbandmode : "(absent)",
                c->q_band_mode_exhaustive ? "exhaustive" : "fixed");
        if (c->q_band_mode_exhaustive)
            fprintf(stderr, "       exhaustive_max_q_bits=%d want=%d\n",
                    rd_exhaustive_max_q_bits, c->exhaustive_max_q_bits);
        return -1;
    }
    return 0;
}

/* ============================================================================
 * Per-V Q-window recomputation for an arbitrary [N1, N2] sub-band.
 *
 * Used by:
 *   - The startup q-iter init pass (mark_exhausted=1, runs once at full band).
 *   - The sub-band rotation tick in the main launch loop (mark_exhausted=0;
 *     transient empty windows must not become permanent).
 *
 * Updates q_min_lo/hi[v], q_max_lo/hi[v] in place. Variants with empty Q
 * windows get sentinel values (q_min=~0, q_max=0) so the dispatch-loop check
 * `qmx_local < qmn` skips them naturally. Returns the count of variants with
 * non-empty windows at this band, or -1 on Q_max > 2^128 (caller should bail).
 * ============================================================================ */
static int compute_q_windows_subband(int n_v_total, int N1, int N2, int E,
                                     v16_pool_seed *pool_seeds,
                                     u64 *q_min_lo, u64 *q_min_hi,
                                     u64 *q_max_lo, u64 *q_max_hi,
                                     u8  *v_exhausted,
                                     int  mark_exhausted)
{
    mpz_t z_lower, z_upper, z_DvE, z_qmin, z_qmax, z_Dv, z_tmp;
    mpz_inits(z_lower, z_upper, z_DvE, z_qmin, z_qmax, z_Dv, z_tmp, NULL);

    mpz_set_ui(z_lower, 1);
    mpz_mul_2exp(z_lower, z_lower, (unsigned long)(N1 - 1));
    mpz_set_ui(z_upper, 1);
    mpz_mul_2exp(z_upper, z_upper, (unsigned long)N2);

    int n_active = 0;
    for (int v = 0; v < n_v_total; v++) {
        if (v_exhausted && v_exhausted[v]) {
            /* Already permanently dead — leave sentinel. */
            q_min_lo[v] = ~(u64)0; q_min_hi[v] = ~(u64)0;
            q_max_lo[v] = 0;       q_max_hi[v] = 0;
            continue;
        }
        v16_pool_seed *ps = &pool_seeds[v];
        mpz_set_str(z_Dv, ps->D_dec, 10);
        mpz_set(z_DvE, z_Dv);
        if (E > 0) mpz_mul_2exp(z_DvE, z_DvE, (unsigned long)E);

        mpz_cdiv_q(z_qmin, z_lower, z_DvE);
        {
            mpz_t z_upper_m1; mpz_init(z_upper_m1);
            mpz_sub_ui(z_upper_m1, z_upper, 1);
            mpz_fdiv_q(z_qmax, z_upper_m1, z_DvE);
            mpz_clear(z_upper_m1);
        }
        if (mpz_sizeinbase(z_qmax, 2) > 128) {
            fprintf(stderr,
                "ERROR: variant %s Q_max bits=%zu exceeds u128. "
                "Raise D_V or lower --bits-max.\n",
                ps->name, mpz_sizeinbase(z_qmax, 2));
            mpz_clears(z_lower, z_upper, z_DvE, z_qmin, z_qmax, z_Dv, z_tmp, NULL);
            return -1;
        }
        if (mpz_cmp(z_qmin, z_qmax) > 0) {
            /* Empty window — sentinel. */
            q_min_lo[v] = ~(u64)0; q_min_hi[v] = ~(u64)0;
            q_max_lo[v] = 0;       q_max_hi[v] = 0;
            if (mark_exhausted && v_exhausted) v_exhausted[v] = 1;
            continue;
        }
        /* Pack u128 lo/hi. */
        mpz_set(z_tmp, z_qmin);
        q_min_lo[v] = (u64)mpz_get_ui(z_tmp);
        mpz_fdiv_q_2exp(z_tmp, z_tmp, 64);
        q_min_hi[v] = (u64)mpz_get_ui(z_tmp);
        mpz_set(z_tmp, z_qmax);
        q_max_lo[v] = (u64)mpz_get_ui(z_tmp);
        mpz_fdiv_q_2exp(z_tmp, z_tmp, 64);
        q_max_hi[v] = (u64)mpz_get_ui(z_tmp);
        n_active++;
    }
    mpz_clears(z_lower, z_upper, z_DvE, z_qmin, z_qmax, z_Dv, z_tmp, NULL);
    return n_active;
}

/* ============================================================================
 * Main.
 * ============================================================================ */
int main(int argc, char **argv) {
    signal(SIGINT, on_sigint);
    if (parse_args(argc, argv, &g_cfg) != 0) return 1;
    if (build_D(&g_cfg) != 0) return 1;
    /* Expand any 'auto:Z#:K' specs now that c->D_primes is set. */
    if (expand_auto_pools(&g_cfg) != 0) return 1;

    /* Bug 1b fix (fix/k15-safety): K1.5 overflow guard.
     *
     * p_top = Q · D_V · 2^(exp_start + top_pos) − 1 is computed in
     * gu128 (128-bit) inside v16_k15_topmr_kernel. If the bit budget
     * exceeds 127, gu128_shl_n / gu128_mul_u128 silently truncate and
     * the Miller-Rabin witnesses run on garbage — survivors get wrongly
     * accepted/rejected.
     *
     * Bit budget: bits(p_top) ≈ bits(Q) + bits(D_V) + exp_start + top_pos,
     * where top_pos = target_len − 1 (K1.5 probes the top-of-chain
     * position) and bits(Q) ≤ bits_max in the q-iter band. Require the
     * sum ≤ 127 (one bit of headroom for the "− 1" at the end).
     *
     * The check is post-pool because max(D_V) is only known after
     * expand_auto_pools(). Triggered ⇒ force --no-k15 with a clear WARN. */
    if (g_cfg.k15_enabled && g_cfg.q_iter && g_cfg.n_pool_seeds > 0) {
        int max_dv_bits = 0;
        for (int v = 0; v < g_cfg.n_pool_seeds; v++) {
            int b = u128_bits(g_cfg.pool_seeds[v].D);
            if (b > max_dv_bits) max_dv_bits = b;
        }
        int top_pos = g_cfg.target_len - 1;
        int budget  = g_cfg.bits_max + max_dv_bits + g_cfg.exp_start
                                                   + top_pos;
        if (budget > 127) {
            fprintf(stderr,
                "WARN: K1.5 disabled — p_top requires %d bits at this "
                "config (bits_max=%d top_pos=%d exp_start=%d "
                "max_Dv≈2^%d, gu128 caps at 127). Re-enable when "
                "bits_max + (target-1) + exp_start + log2(max_Dv) <= 127.\n",
                budget,
                g_cfg.bits_max, top_pos, g_cfg.exp_start,
                max_dv_bits);
            g_cfg.k15_enabled = 0;
        }
        /* K1.5 guard #2 (Test A.3/A.4 finding): K1.5 emits only survivors
         * where p_top = p_{target-1} is MR-prime. If --min-report-len < target,
         * the engine accepts shorter chains (e.g. CC8 when hunting CC14),
         * which can have a composite p_top. K1.5 would silently drop those.
         * Force --no-k15 in this regime; re-enable only when
         * min_report_len == target (the engine cares only about length-target
         * chains). */
        if (g_cfg.k15_enabled && g_cfg.min_report_len > 0 &&
            g_cfg.min_report_len < g_cfg.target_len) {
            fprintf(stderr,
                "WARN: K1.5 disabled — min_report_len=%d < target=%d would "
                "drop valid CC%d..CC%d chains whose p_top is composite. "
                "Re-enable only when min_report_len == target.\n",
                g_cfg.min_report_len, g_cfg.target_len,
                g_cfg.min_report_len, g_cfg.target_len - 1);
            g_cfg.k15_enabled = 0;
        }
    }

    printf("# cc20_first_kind_immune_v16  %s\n", V16_VERSION);
    printf("#   seed       : %s\n", g_cfg.seed_name);
    printf("#   D primes   : ");
    for (int i = 0; i < g_cfg.n_D_primes; i++)
        printf("%s%u", i ? "," : "", g_cfg.D_primes[i]);
    printf("\n");
    printf("#   D          : %s (%d bits)\n",
           g_cfg.D_dec, u128_bits(g_cfg.D));
    printf("#   target     : %d\n", g_cfg.target_len);
    printf("#   sieve_depth: %d\n", g_cfg.sieve_depth);
    printf("#   exp_start  : %d\n", g_cfg.exp_start);
    /* Depth/target relation warnings (design §15). Covers both q-iter and
     * legacy --sieve-variants paths — banner runs before either dispatch. */
    if (g_cfg.sieve_depth > g_cfg.target_len) {
        fprintf(stderr,
            "WARN: depth=%d > target=%d -- chains of exact length %d will be\n"
            "      filtered out if p_%d is small-prime-composite (sieve cannot\n"
            "      distinguish 'chain ends here' from 'chain continues with\n"
            "      composite p_%d'). See the depth/target discussion in DESIGN.md.\n",
            g_cfg.sieve_depth, g_cfg.target_len, g_cfg.target_len,
            g_cfg.target_len, g_cfg.target_len);
    } else if (g_cfg.sieve_depth < g_cfg.target_len) {
        fprintf(stderr,
            "INFO: depth=%d < target=%d -- weaker sieve, higher survivor rate.\n"
            "      Intended use: profiling / SW-verification (design §15).\n"
            "      Production sieve sets depth == target.\n",
            g_cfg.sieve_depth, g_cfg.target_len);
    }

    if (g_cfg.validate_seed) {
        return validate_seed(&g_cfg) == 0 ? 0 : 2;
    }

    select_wheel_sieve(&g_cfg);
    printf("#   wheel primes (%d): ", g_cfg.n_wheel_primes);
    for (int i = 0; i < g_cfg.n_wheel_primes; i++)
        printf("%s%u", i ? "," : "", g_cfg.wheel_primes[i]);
    printf("\n");
    printf("#   sieve primes (%d): ", g_cfg.n_sieve_primes);
    for (int i = 0; i < g_cfg.n_sieve_primes; i++)
        printf("%s%u", i ? "," : "", g_cfg.sieve_primes[i]);
    printf("\n");
    /* Explicit "active sieve" diagnostic. Surfaces the
     * largest prime actually in use so a silent clamp from a too-large
     * --sieve-max (now caught at parse time) or from D-coprime filtering is
     * immediately visible in the banner. */
    {
        u32 active_max =
            (g_cfg.n_sieve_primes > 0)
                ? g_cfg.sieve_primes[g_cfg.n_sieve_primes - 1]
                : ((g_cfg.n_wheel_primes > 0)
                       ? g_cfg.wheel_primes[g_cfg.n_wheel_primes - 1]
                       : 0u);
        u32 active_min = (g_cfg.n_wheel_primes > 0)
                             ? g_cfg.wheel_primes[0]
                             : ((g_cfg.n_sieve_primes > 0)
                                    ? g_cfg.sieve_primes[0]
                                    : 0u);
        int total = g_cfg.n_wheel_primes + g_cfg.n_sieve_primes;
        fprintf(stderr,
            "INFO: sieve uses %d primes [%u..%u] (cap --sieve-max=%u)\n",
            total, active_min, active_max, g_cfg.sieve_prime_max);
    }

    /* --sieve-variants: validate every pool seed obeys the kernel's hard
     * precondition (V's D is a superset of sieve D, and V adds no primes
     * that collide with active wheel/sieve primes). */
    if (g_cfg.sieve_variants && validate_sieve_variants(&g_cfg) != 0) {
        return 1;
    }
    if (g_cfg.sieve_variants && g_cfg.n_pool_seeds > V16_MAX_POOL_SEEDS) {
        fprintf(stderr,
            "ERROR: --sieve-variants: pool seed count %d exceeds u16 cap\n",
            g_cfg.n_pool_seeds);
        return 1;
    }

    /* Extended-sieve forbid table: depth-aware CPU pre-MR filter for primes
     * 509..1009 (above the V16_Q_MASK_WORDS=8 cap of the GPU sieve). Built
     * once here, drained on demand inside cpu_prove_worker. ~43 MiB at the
     * 4248-variant pool. --no-ext-sieve to disable. */
    if (build_ext_forbid_table(&g_cfg) != 0) return 2;

    /* Startup engine-integrity self-test (default ON; --no-test to skip).
     * Catches the class of bug discovered 2026-05-14 where a structural
     * pool error silently blinded the engine to 76% of CC10+ chain shapes.
     * CPU/GMP only; ~30s. Failure exits before any GPU init. */
    if (g_cfg.run_self_test) {
        int rc = run_v16_self_test(&g_cfg);
        if (rc != 0) return rc;
    } else {
        fprintf(stderr, "INFO: self-test skipped (--no-test)\n");
    }

    /* ============================================================================
     * Q-iter mode entrypoint. When --q-iter is on we take a wholly separate
     * code path (no wheel, no kmod LUT, no variant masks); on completion we
     * return directly from main(). The non-q-iter path below is byte-
     * identical to the pre-q-iter baseline.
     * ============================================================================ */
    if (g_cfg.q_iter) {
        /* ---- 0) Echo Q-iter config + pool. ---- */
        printf("#   q-iter     : ON  (target ≥ 10 G (Q,V)/sec; design §11)\n");
        /* q-order banner deferred to the seeded-init site below (shows seed). */
        printf("#   q-chunk    : %u Q/launch\n", g_cfg.q_chunk);
        printf("#   bits range : [%d, %d]\n",
               g_cfg.bits_min, g_cfg.bits_max);
        if (g_cfg.k15_enabled) {
            printf("#   k15        : ON  (reps=%d)\n", g_cfg.k15_mr_reps);
        } else {
            printf("#   k15        : off (legacy K1 → CPU prove path)\n");
        }
        if (g_cfg.n_pool_seeds == 0) {
            fprintf(stderr,
                "ERROR: --q-iter requires --seed-pool (expanded to >=1 V)\n");
            return 1;
        }
        if (g_cfg.n_pool_seeds > V16_MAX_POOL_SEEDS) {
            fprintf(stderr,
                "ERROR: pool seed count %d exceeds V16_MAX_POOL_SEEDS=%d\n",
                g_cfg.n_pool_seeds, V16_MAX_POOL_SEEDS);
            return 1;
        }
        if (g_cfg.n_pool_seeds > (1 << 16)) {
            fprintf(stderr,
                "ERROR: pool seed count %d exceeds u16 (QSurvivor.seed_idx)\n",
                g_cfg.n_pool_seeds);
            return 1;
        }
        if (g_cfg.n_pool_primes > 0) {
            printf("#   pool-primes: ");
            for (int i = 0; i < g_cfg.n_pool_primes; i++)
                printf("%s%u", i ? "," : "", g_cfg.pool_primes[i]);
            printf("\n");
        }
        printf("#   pool-seeds: %d\n", g_cfg.n_pool_seeds);
        for (int i = 0; i < g_cfg.n_pool_seeds && i < 8; i++) {
            v16_pool_seed *ps = &g_cfg.pool_seeds[i];
            printf("#     %s  D=%s  primes=", ps->name, ps->D_dec);
            for (int j = 0; j < ps->n_primes; j++)
                printf("%s%u", j ? "*" : "", ps->primes[j]);
            printf("\n");
        }
        if (g_cfg.n_pool_seeds > 8)
            printf("#     ... (%d more)\n", g_cfg.n_pool_seeds - 8);

        /* ---- 1) Build the Q-iter sieve prime list. Reuse select_wheel_sieve
         * output: under Q-iter, wheel+sieve are merged into one ascending
         * list (no CRT wheel exists). */
        u32 q_primes[V16_MAX_SIEVE + V16_MAX_WHEEL];
        int n_qp = 0;
        for (int i = 0; i < g_cfg.n_wheel_primes &&
                       n_qp < (int)(sizeof(q_primes)/sizeof(q_primes[0]));
             i++) {
            q_primes[n_qp++] = g_cfg.wheel_primes[i];
        }
        for (int i = 0; i < g_cfg.n_sieve_primes &&
                       n_qp < (int)(sizeof(q_primes)/sizeof(q_primes[0]));
             i++) {
            q_primes[n_qp++] = g_cfg.sieve_primes[i];
        }
        /* Sort ascending (both lists are already ascending and disjoint by
         * cfg.wheel_prime_max). Just merge then assert ascending. */
        for (int i = 0; i < n_qp; i++) {
            for (int j = i + 1; j < n_qp; j++) {
                if (q_primes[j] < q_primes[i]) {
                    u32 t = q_primes[i]; q_primes[i] = q_primes[j];
                    q_primes[j] = t;
                }
            }
        }
        for (int i = 1; i < n_qp; i++) {
            if (q_primes[i] <= q_primes[i - 1]) {
                fprintf(stderr,
                    "ERROR: Q-iter sieve primes not strictly ascending "
                    "at i=%d (%u <= %u)\n",
                    i, q_primes[i], q_primes[i - 1]);
                return 1;
            }
        }
        if (n_qp > V16_MAX_SIEVE) {
            fprintf(stderr,
                "ERROR: Q-iter sieve prime count %d exceeds V16_MAX_SIEVE=%d\n",
                n_qp, V16_MAX_SIEVE);
            return 1;
        }
        if (n_qp == 0) {
            fprintf(stderr,
                "ERROR: Q-iter has 0 sieve primes — every Q would pass. "
                "Raise --sieve-max or check --immune-prime.\n");
            return 1;
        }
        printf("#   q-iter sieve primes (%d): ", n_qp);
        for (int i = 0; i < n_qp; i++)
            printf("%s%u", i ? "," : "", q_primes[i]);
        printf("\n");

        /* Barrett magic numbers for every q-iter prime.
         * Self-test against HW mod over a thorough sweep at boot. */
        BarrettQ *q_barrett = (BarrettQ *)malloc((size_t)n_qp * sizeof(BarrettQ));
        if (!q_barrett) {
            fprintf(stderr, "OOM BarrettQ table\n");
            return 1;
        }
        for (int i = 0; i < n_qp; i++) {
            q_barrett[i].q     = q_primes[i];
            q_barrett[i].magic = barrett_magic_u32(q_primes[i]);
        }
        {
            size_t barrett_ok = 0, barrett_fail = 0;
            static const u32 probes[] = {
                0u, 1u, 2u, 3u, 4u, 7u, 8u, 15u, 16u, 31u, 32u, 63u, 64u,
                127u, 128u, 255u, 256u, 463u, 503u, 1023u, 1024u, 65535u,
                65536u, 1048575u, 1048576u, 16777215u, 16777216u,
                0x7fffffffu, 0x80000000u, 0x80000001u, 0xc0000000u,
                0xffffff00u, 0xfffffffeu, 0xffffffffu
            };
            for (int i = 0; i < n_qp; i++) {
                u32 q = q_barrett[i].q, m = q_barrett[i].magic;
                for (size_t k = 0; k < sizeof(probes)/sizeof(probes[0]); k++) {
                    u32 x = probes[k];
                    u32 ref = x % q;
                    u32 got = barrett_mod_u32(x, q, m);
                    if (ref == got) barrett_ok++;
                    else            barrett_fail++;
                }
                /* Dense sweep across [0, 4q) and around 2^k boundaries
                 * — covers worst-case under-estimate cases. */
                for (u32 x = 0; x < 4u * q; x++) {
                    u32 ref = x % q;
                    u32 got = barrett_mod_u32(x, q, m);
                    if (ref == got) barrett_ok++;
                    else            barrett_fail++;
                }
                static const u32 corners[] = {
                    0x7fffff00u, 0x80000000u, 0x80000100u,
                    0xc0000000u, 0xfffff000u, 0xffffff00u, 0xfffffff0u,
                    0xfffffffcu, 0xfffffffdu, 0xfffffffeu, 0xffffffffu
                };
                for (size_t k = 0; k < sizeof(corners)/sizeof(corners[0]); k++) {
                    for (int d = -3; d <= 3; d++) {
                        u32 x = corners[k] + (u32)d;
                        u32 ref = x % q;
                        u32 got = barrett_mod_u32(x, q, m);
                        if (ref == got) barrett_ok++;
                        else            barrett_fail++;
                    }
                }
            }
            if (barrett_fail) {
                /* Dump first ~10 mismatches for triage. */
                int dumped = 0;
                for (int i = 0; i < n_qp && dumped < 10; i++) {
                    u32 q = q_barrett[i].q, m = q_barrett[i].magic;
                    for (size_t k = 0; k < sizeof(probes)/sizeof(probes[0]) && dumped < 10; k++) {
                        u32 x = probes[k];
                        u32 ref = x % q, got = barrett_mod_u32(x, q, m);
                        if (ref != got) {
                            fprintf(stderr,
                                "  barrett mismatch q=%u magic=0x%08x x=0x%08x ref=%u got=%u\n",
                                q, m, x, ref, got);
                            dumped++;
                        }
                    }
                }
                fprintf(stderr,
                    "ERROR: Barrett self-test FAILED (%zu/%zu probes) -- abort.\n",
                    barrett_fail, barrett_ok + barrett_fail);
                return 1;
            }
            printf("#   q-iter Barrett self-test: %zu/%zu probes OK\n",
                   barrett_ok, barrett_ok + barrett_fail);
        }

        /* Boot-time self-test for gu128_mul_u128. Carry handling is
         * the bug-prone part; we test
         * against a u128 reference at boundary cases Q ≈ 2^64, 2^96, 2^120,
         * 2^127, plus a full-width and two small cases. Hard-error on
         * mismatch. The kernel + helper are stateless, so this fires once
         * per binary launch regardless of --k15. */
        if (g_cfg.k15_enabled) {
            const int N_ST = 7;
            const u64 sa_lo[7] = {
                0x8000000000000000ULL,
                0ULL,
                0xFFFFFFFFFFFFFFFFULL,
                0xFFFFFFFFFFFFFFFFULL,
                0xDEADBEEFCAFEBABEULL,
                0x0000000000000123ULL,
                0x0000000100000000ULL,
            };
            const u64 sa_hi[7] = {
                0ULL,
                0x0000000100000000ULL,
                0x00FFFFFFFFFFFFFFULL,
                0x7FFFFFFFFFFFFFFFULL,
                0ULL,
                0ULL,
                0ULL,
            };
            const u64 sb_lo[7] = {
                2ULL,
                0xFFFFFFFFULL,
                17ULL,
                3ULL,
                0x123456789ABCDEF0ULL,
                0x0000000000000456ULL,
                0x0000000100000000ULL,
            };
            const u64 sb_hi[7] = {0,0,0,0,0,0,0};
            u64 *d_a_lo, *d_a_hi, *d_b_lo, *d_b_hi, *d_o_lo, *d_o_hi;
            CUDA_CHECK(cudaMalloc(&d_a_lo, N_ST * sizeof(u64)));
            CUDA_CHECK(cudaMalloc(&d_a_hi, N_ST * sizeof(u64)));
            CUDA_CHECK(cudaMalloc(&d_b_lo, N_ST * sizeof(u64)));
            CUDA_CHECK(cudaMalloc(&d_b_hi, N_ST * sizeof(u64)));
            CUDA_CHECK(cudaMalloc(&d_o_lo, N_ST * sizeof(u64)));
            CUDA_CHECK(cudaMalloc(&d_o_hi, N_ST * sizeof(u64)));
            CUDA_CHECK(cudaMemcpy(d_a_lo, sa_lo, N_ST * sizeof(u64),
                                  cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy(d_a_hi, sa_hi, N_ST * sizeof(u64),
                                  cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy(d_b_lo, sb_lo, N_ST * sizeof(u64),
                                  cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy(d_b_hi, sb_hi, N_ST * sizeof(u64),
                                  cudaMemcpyHostToDevice));
            v16_k15_mul_u128_selftest_kernel<<<1, 32>>>(
                d_a_lo, d_a_hi, d_b_lo, d_b_hi, d_o_lo, d_o_hi, N_ST);
            CUDA_CHECK(cudaGetLastError());
            CUDA_CHECK(cudaDeviceSynchronize());
            u64 got_lo[7], got_hi[7];
            CUDA_CHECK(cudaMemcpy(got_lo, d_o_lo, N_ST * sizeof(u64),
                                  cudaMemcpyDeviceToHost));
            CUDA_CHECK(cudaMemcpy(got_hi, d_o_hi, N_ST * sizeof(u64),
                                  cudaMemcpyDeviceToHost));
            int n_pass = 0;
            for (int i = 0; i < N_ST; i++) {
                u128 a = ((u128)sa_hi[i] << 64) | (u128)sa_lo[i];
                u128 b = ((u128)sb_hi[i] << 64) | (u128)sb_lo[i];
                u128 ref = a * b;
                u64 ref_lo = (u64)ref;
                u64 ref_hi = (u64)(ref >> 64);
                if (got_lo[i] == ref_lo && got_hi[i] == ref_hi) {
                    n_pass++;
                } else {
                    fprintf(stderr,
                        "K1.5 gu128_mul_u128 self-test FAIL case %d:\n"
                        "  a    = 0x%016llx%016llx\n"
                        "  b    = 0x%016llx%016llx\n"
                        "  got  = 0x%016llx%016llx\n"
                        "  want = 0x%016llx%016llx\n",
                        i,
                        (unsigned long long)sa_hi[i],
                        (unsigned long long)sa_lo[i],
                        (unsigned long long)sb_hi[i],
                        (unsigned long long)sb_lo[i],
                        (unsigned long long)got_hi[i],
                        (unsigned long long)got_lo[i],
                        (unsigned long long)ref_hi,
                        (unsigned long long)ref_lo);
                }
            }
            cudaFree(d_a_lo); cudaFree(d_a_hi);
            cudaFree(d_b_lo); cudaFree(d_b_hi);
            cudaFree(d_o_lo); cudaFree(d_o_hi);
            if (n_pass != N_ST) {
                fprintf(stderr,
                    "K1.5 gu128_mul_u128 self-test: FAIL (%d/%d passed)\n",
                    n_pass, N_ST);
                return 1;
            }
            printf("# K1.5 gu128_mul_u128 self-test: PASS (%d cases)\n",
                   n_pass);
        }

        const int depth = g_cfg.sieve_depth;
        const int E = g_cfg.exp_start;
        /* Bit-band [N1, N2]. */
        const int N1 = g_cfg.bits_min;
        const int N2 = g_cfg.bits_max;

        /* ---- 2) Shared-mem budget check:
         *   - s_qb (BarrettQ pairs) : n_sp * 8 B
         *   - s_forbid_mask         : n_sp * V16_Q_MASK_WORDS * 8 B
         *   - s_base_mod            : n_sp * 4 B
         *   Independent of --depth because the bitmask collapses the j-loop. */
        size_t needed_shared =
            (size_t)n_qp * sizeof(BarrettQ) +
            (size_t)n_qp * (size_t)V16_Q_MASK_WORDS * sizeof(u64) +
            (size_t)n_qp * sizeof(u32);
        if (needed_shared > 48u * 1024u) {
            fprintf(stderr,
                "ERROR: Q-iter shared mem need %zu B exceeds 48 KiB "
                "(n_sp=%d, mask_words=%d). Lower --sieve-max or shrink mask.\n",
                needed_shared, n_qp, V16_Q_MASK_WORDS);
            return 1;
        }
        printf("#   q-iter shared mem: %zu B (n_sp=%d × (8B Barrett + %dB mask + 4B base))\n",
               needed_shared, n_qp, V16_Q_MASK_WORDS * 8);

        /* ---- 3) Per-V Q windows (Q_min, Q_max) via GMP. Mark V's with
         * empty windows as exhausted.
         *
         * Q_min/Q_max/Q_cursor are u128 (carried as paired u64 lo/hi).
         * Keep u128 on the host and pass only per-launch residues to the
         * kernel. Hard cap: full u128, i.e. < 2^128. */
        int n_v_total = g_cfg.n_pool_seeds;
        u64 *q_min_lo = (u64 *)calloc((size_t)n_v_total, sizeof(u64));
        u64 *q_min_hi = (u64 *)calloc((size_t)n_v_total, sizeof(u64));
        u64 *q_max_lo = (u64 *)calloc((size_t)n_v_total, sizeof(u64));
        u64 *q_max_hi = (u64 *)calloc((size_t)n_v_total, sizeof(u64));
        u8  *v_exhausted = (u8 *)calloc((size_t)n_v_total, 1);
        u64 *v_cursor_lo = (u64 *)calloc((size_t)n_v_total, sizeof(u64));
        u64 *v_cursor_hi = (u64 *)calloc((size_t)n_v_total, sizeof(u64));
        /* Split cursor into launched / drained.
         *   launched (v_cursor_lo/hi)  = front edge — what to launch next.
         *   drained  (v_cursor_drained) = rear edge  — last Q whose survivors
         *                                  are CPU-proved. Only the drained
         *                                  cursor is checkpointed, so a
         *                                  kill between launch and prove
         *                                  doesn't lose chains in that range. */
        u64 *v_cursor_drained_lo = (u64 *)calloc((size_t)n_v_total, sizeof(u64));
        u64 *v_cursor_drained_hi = (u64 *)calloc((size_t)n_v_total, sizeof(u64));
        if (!q_min_lo || !q_min_hi || !q_max_lo || !q_max_hi ||
            !v_exhausted || !v_cursor_lo || !v_cursor_hi ||
            !v_cursor_drained_lo || !v_cursor_drained_hi) {
            fprintf(stderr, "OOM Q-iter per-V arrays\n");
            return 1;
        }
        /* Pool D_V bit-width range (across the whole loaded pool, before any
         * band filtering). Cheap scan; used by the startup-band warning. */
        size_t pool_dv_bits_min = (size_t)-1;
        size_t pool_dv_bits_max = 0;
        {
            mpz_t z_Dv; mpz_init(z_Dv);
            for (int v = 0; v < n_v_total; v++) {
                mpz_set_str(z_Dv, g_cfg.pool_seeds[v].D_dec, 10);
                size_t b = mpz_sizeinbase(z_Dv, 2);
                if (b < pool_dv_bits_min) pool_dv_bits_min = b;
                if (b > pool_dv_bits_max) pool_dv_bits_max = b;
            }
            mpz_clear(z_Dv);
        }

        /* Band [N1, N2] = "p_0 has between N1 and N2 bits inclusive",
         * i.e. p_0 in [2^(N1-1), 2^N2). N1==N2 reproduces exact-bit
         * filter. Helper marks empty-window variants permanently exhausted. */
        int n_active = compute_q_windows_subband(
            n_v_total, N1, N2, E, g_cfg.pool_seeds,
            q_min_lo, q_min_hi, q_max_lo, q_max_hi,
            v_exhausted, /*mark_exhausted=*/1);
        if (n_active < 0) return 1;   /* Q_max >= 2^128 — message already printed. */

        /* Cursor init for sequential mode: front edge = drained edge = q_min.
         * (Random mode ignores cursors.) */
        for (int v = 0; v < n_v_total; v++) {
            if (v_exhausted[v]) continue;
            v_cursor_lo[v]         = q_min_lo[v];
            v_cursor_hi[v]         = q_min_hi[v];
            v_cursor_drained_lo[v] = q_min_lo[v];
            v_cursor_drained_hi[v] = q_min_hi[v];
        }

        /* Per-V random/sequential flag.
         *   fixed mode      → all V's follow g_cfg.q_order_random (legacy).
         *   exhaustive mode → per-V: sequential when bit-width(Q_max-Q_min+1)
         *                     <= exhaustive_max_q_bits, random otherwise. */
        u8 *v_use_random = (u8 *)calloc((size_t)n_v_total, 1);
        if (!v_use_random) {
            fprintf(stderr, "OOM v_use_random[]\n");
            return 1;
        }
        int n_v_sequential = 0, n_v_random = 0;
        for (int v = 0; v < n_v_total; v++) {
            if (v_exhausted[v]) continue;
            if (g_cfg.q_band_mode_exhaustive) {
                u128 qmn = ((u128)q_min_hi[v] << 64) | (u128)q_min_lo[v];
                u128 qmx = ((u128)q_max_hi[v] << 64) | (u128)q_max_lo[v];
                u128 span = (qmx >= qmn) ? (qmx - qmn) + (u128)1 : (u128)0;
                /* u128 bit-length. */
                int span_bits;
                u64 span_hi = (u64)(span >> 64);
                u64 span_lo = (u64)span;
                if (span_hi)      span_bits = 128 - __builtin_clzll(span_hi);
                else if (span_lo) span_bits = 64  - __builtin_clzll(span_lo);
                else              span_bits = 0;
                v_use_random[v] = (span_bits > g_cfg.exhaustive_max_q_bits) ? 1 : 0;
            } else {
                v_use_random[v] = (u8)(g_cfg.q_order_random ? 1 : 0);
            }
            if (v_use_random[v]) n_v_random++; else n_v_sequential++;
        }

        /* Startup banner. Whether or not anything got dropped, print
         * the configured-vs-actual coverage so users can see at a glance how
         * much of their loaded pool the engine will actually scan. If variants
         * got dropped, scream about it and point at --bits-rotate-step. */
        int n_dropped = n_v_total - n_active;
        printf("#   q-iter band      : bits ∈ [%d, %d]  (delta = %d bits)\n",
               N1, N2, N2 - N1 + 1);
        printf("#   pool D_V bits    : [%zu, %zu]  across %d loaded variants\n",
               pool_dv_bits_min, pool_dv_bits_max, n_v_total);
        printf("#   pool active here : %d / %d  (%d dropped: D_V too wide to fit band)\n",
               n_active, n_v_total, n_dropped);
        if (n_dropped > 0) {
            fprintf(stderr,
                "\n"
                "WARN: %d / %d pool variants have empty Q-windows at band [%d, %d]\n"
                "      and will NOT be scanned.\n"
                "      Cause: their D_V values (or D_V·2^E with --exp-start) leave\n"
                "      no room for Q in the configured bit-band.\n"
                "      To use all variants, EITHER:\n"
                "        (a) widen the band: lower --bits-min / raise --bits-max\n"
                "        (b) rotate sub-bands: --bits-rotate-step D (e.g. 3) — cycles\n"
                "            through [N1, N1+D], [N1+D, N1+2D], ... so each variant\n"
                "            sees a band it fits, instead of being permanently dropped\n"
                "        (c) use a lighter pool with smaller D_V values\n"
                "      Pool D_V bits=[%zu, %zu], band delta=%d bits.\n\n",
                n_dropped, n_v_total, N1, N2,
                pool_dv_bits_min, pool_dv_bits_max, N2 - N1 + 1);
        }
        if (n_active == 0) {
            fprintf(stderr,
                "ERROR: every pool variant has an empty Q-window at "
                "bits=[%d,%d]. Widen --bits-min/--bits-max or use a "
                "lighter pool spec.\n", N1, N2);
            return 1;
        }

        /* Validate rotation knobs and emit summary if enabled. */
        int rotate_step   = g_cfg.bits_rotate_step;
        int rotate_period = g_cfg.bits_rotate_period_sec;
        if (rotate_step < 0) {
            fprintf(stderr, "ERROR: --bits-rotate-step must be >= 0 (got %d)\n", rotate_step);
            return 1;
        }
        if (rotate_step > 0 && rotate_step > (N2 - N1 + 1)) {
            fprintf(stderr,
                "WARN: --bits-rotate-step %d exceeds band width %d; disabling rotation.\n",
                rotate_step, N2 - N1 + 1);
            rotate_step = 0;
        }
        if (rotate_step > 0 && rotate_period < 1) {
            fprintf(stderr, "ERROR: --bits-rotate-period must be >= 1 (got %d)\n", rotate_period);
            return 1;
        }
        if (rotate_step > 0) {
            int n_subbands = (N2 - N1 + 1 + rotate_step - 1) / rotate_step;
            printf("#   bit-band rotate  : step=%d period=%ds → %d sub-bands cycling [%d, %d]\n",
                   rotate_step, rotate_period, n_subbands, N1, N2);
        }

        /* --q-band-mode exhaustive classifies every V's sequential/random
         * choice from the *initial* full-band Q-range. Rotation later carves
         * sub-bands and recomputes Q_min/Q_max per sub-band, but the
         * sequential/random flag is not refreshed — a V that's "too wide"
         * over the full band would stay random even if a rotated sub-band
         * is small enough for deterministic coverage. Refuse rather than
         * lie about coverage. */
        if (g_cfg.q_band_mode_exhaustive && rotate_step > 0) {
            fprintf(stderr,
                "ERROR: --q-band-mode exhaustive is incompatible with "
                "--bits-rotate-step > 0.\n"
                "       The per-V sequential/random classification is "
                "computed once from the\n"
                "       initial Q-range and would not adapt to rotated "
                "sub-bands. Drop one\n"
                "       of the two flags.\n");
            return 1;
        }

        /* ---- 4) Build forbid_Q table (transient) and per-(V,q) bitmask.
         *
         *   forbid_Q[V][j][i] is the single forbidden residue at depth j for
         *   prime q[i] in variant V. We OR all (j) bits into a per-(V,q)
         *   bitmask `forbid_mask[V][q]` so the kernel can probe with one
         *   shift+AND instead of a depth-long compare loop.
         *
         *   Sentinel (q | D_V) -> no forbidden residues at q, mask all zero.
         *
         *   forbid_Q itself is freed after construction + round-trip self-test.
         */
        if (n_qp > V16_MAX_SIEVE) {
            fprintf(stderr, "ERROR: n_qp %d > V16_MAX_SIEVE %d\n",
                    n_qp, V16_MAX_SIEVE); return 1;
        }
        /* Sanity: all primes fit V16_Q_MASK_WORDS u64s. */
        for (int i = 0; i < n_qp; i++) {
            if (q_primes[i] > (u32)(V16_Q_MASK_WORDS * 64)) {
                fprintf(stderr,
                    "ERROR: sieve prime %u exceeds V16_Q_MASK_WORDS*64=%d "
                    "-- bump V16_Q_MASK_WORDS.\n",
                    q_primes[i], V16_Q_MASK_WORDS * 64);
                return 1;
            }
        }
        size_t fb_slab_u32 = (size_t)depth * (size_t)n_qp;
        size_t fb_total_u32 = (size_t)n_v_total * fb_slab_u32;
        u32 *h_forbid_all = (u32 *)malloc(fb_total_u32 * sizeof(u32));
        if (!h_forbid_all) {
            fprintf(stderr, "OOM forbid_Q (%.2f MiB)\n",
                    fb_total_u32 * sizeof(u32) / (1024.0 * 1024.0));
            return 1;
        }
        for (int v = 0; v < n_v_total; v++) {
            v16_pool_seed *ps = &g_cfg.pool_seeds[v];
            u32 *slab = h_forbid_all + (size_t)v * fb_slab_u32;
            for (int j = 0; j < depth; j++) {
                for (int i = 0; i < n_qp; i++) {
                    u32 q = q_primes[i];
                    u32 D_mod_q = (u32)(ps->D % (u128)q);
                    u32 *cell = &slab[(size_t)j * n_qp + i];
                    if (D_mod_q == 0) {
                        /* q | D_V → V immune at q, no forbidden Q-residue. */
                        *cell = V16_Q_FORBID_SENTINEL;
                        continue;
                    }
                    u32 p2 = pow2_mod_u32(E + j, q);
                    /* inv(D_V · 2^(E+j) mod q) */
                    u32 m = (u32)(((u64)D_mod_q * (u64)p2) % (u64)q);
                    u32 inv = mod_inv_u32(m, q);
                    *cell = inv;
                }
            }
        }
        /* Build bitmask slab: one bit per residue, OR'd across depths. */
        size_t mask_slab_u64  = (size_t)n_qp * (size_t)V16_Q_MASK_WORDS;
        size_t mask_total_u64 = (size_t)n_v_total * mask_slab_u64;
        u64 *h_forbid_mask_all = (u64 *)calloc(mask_total_u64, sizeof(u64));
        if (!h_forbid_mask_all) {
            fprintf(stderr, "OOM forbid_mask (%.2f MiB)\n",
                    mask_total_u64 * sizeof(u64) / (1024.0 * 1024.0));
            return 1;
        }
        for (int v = 0; v < n_v_total; v++) {
            u32 *src   = h_forbid_all       + (size_t)v * fb_slab_u32;
            u64 *mslab = h_forbid_mask_all  + (size_t)v * mask_slab_u64;
            for (int i = 0; i < n_qp; i++) {
                u64 *mask = &mslab[(size_t)i * V16_Q_MASK_WORDS];
                for (int j = 0; j < depth; j++) {
                    u32 fb = src[(size_t)j * n_qp + i];
                    if (fb == V16_Q_FORBID_SENTINEL) continue;
                    /* fb < q < V16_Q_MASK_WORDS*64; safe. */
                    mask[fb >> 6] |= (1ULL << (fb & 63u));
                }
            }
        }
        /* Round-trip self-test: enumerate set bits of mask, confirm each
         * appears as a forbid_Q entry for that (V,q); confirm every
         * non-sentinel forbid_Q entry has its bit set in mask. */
        {
            size_t rt_ok = 0, rt_fail = 0;
            for (int v = 0; v < n_v_total; v++) {
                u32 *src   = h_forbid_all      + (size_t)v * fb_slab_u32;
                u64 *mslab = h_forbid_mask_all + (size_t)v * mask_slab_u64;
                for (int i = 0; i < n_qp; i++) {
                    u32 q = q_primes[i];
                    u64 *mask = &mslab[(size_t)i * V16_Q_MASK_WORDS];
                    /* (a) Every non-sentinel forbid_Q value -> bit set. */
                    for (int j = 0; j < depth; j++) {
                        u32 fb = src[(size_t)j * n_qp + i];
                        if (fb == V16_Q_FORBID_SENTINEL) continue;
                        int got = (mask[fb >> 6] >> (fb & 63u)) & 1ULL;
                        if (got) rt_ok++; else rt_fail++;
                    }
                    /* (b) Every set bit b < q -> b appears in src. */
                    for (u32 b = 0; b < q; b++) {
                        int set = (mask[b >> 6] >> (b & 63u)) & 1ULL;
                        if (!set) continue;
                        int found = 0;
                        for (int j = 0; j < depth; j++) {
                            if (src[(size_t)j * n_qp + i] == b) { found = 1; break; }
                        }
                        if (found) rt_ok++; else rt_fail++;
                    }
                    /* (c) No bit set at position >= q. */
                    for (u32 b = q; b < (u32)V16_Q_MASK_WORDS * 64; b++) {
                        int set = (mask[b >> 6] >> (b & 63u)) & 1ULL;
                        if (!set) rt_ok++; else rt_fail++;
                    }
                }
            }
            if (rt_fail) {
                fprintf(stderr,
                    "ERROR: forbid_mask round-trip FAILED (%zu/%zu) -- abort.\n",
                    rt_fail, rt_ok + rt_fail);
                return 1;
            }
            printf("#   q-iter forbid_mask round-trip: %zu/%zu OK\n",
                   rt_ok, rt_ok + rt_fail);
        }
        printf("#   q-iter forbid_Q (transient): %.2f MiB (built then freed)\n",
               fb_total_u32 * sizeof(u32) / (1024.0 * 1024.0));
        printf("#   q-iter forbid_mask: %.2f MiB (%d V × %d primes × %d words × 8 B)\n",
               mask_total_u64 * sizeof(u64) / (1024.0 * 1024.0),
               n_v_total, n_qp, V16_Q_MASK_WORDS);
        free(h_forbid_all);
        h_forbid_all = NULL;
        (void)fb_total_u32;

        /* ---- 4b) Per-V active-prime list.
         *
         * For each V, walk q_primes[0..n_qp) and emit indices i where the
         * forbid_mask for (V, i) is NOT all-zero. Primes where q | D_V have
         * an all-zero mask (no Q-residue can make p ≡ 0 mod q, since q
         * already divides D), so the K1 inner loop wastes ~8 IMAD per such
         * prime per Q. Skipping them buys ~10% K1 throughput on the
         * 3433V exp_d pool.
         *
         * Sanity: when q | D_V the mask MUST be all-zero (forbid_Q hit the
         * sentinel branch above). When q does not divide D_V the mask
         * SHOULD have at least one bit set; all-zero is a mask-builder bug
         * — warn and bail. */
        u8 *h_active_primes_all = NULL;
        u8 *h_n_active          = NULL;
        if (g_cfg.active_primes) {
            if (n_qp > 255) {
                fprintf(stderr,
                    "ERROR: --active-primes requires n_sp <= 255 (got %d)\n",
                    n_qp);
                return 1;
            }
            size_t ap_total = (size_t)n_v_total * (size_t)n_qp;
            h_active_primes_all = (u8 *)calloc(ap_total, 1);
            h_n_active          = (u8 *)calloc((size_t)n_v_total, 1);
            if (!h_active_primes_all || !h_n_active) {
                fprintf(stderr, "OOM active_primes (%zu B)\n",
                        ap_total + (size_t)n_v_total);
                return 1;
            }
            u64 sum_active = 0, sum_total = 0;
            u32 min_active = (u32)n_qp, max_active = 0;
            u32 max_pruned = 0;
            for (int v = 0; v < n_v_total; v++) {
                v16_pool_seed *ps = &g_cfg.pool_seeds[v];
                u64 *mslab = h_forbid_mask_all + (size_t)v * mask_slab_u64;
                u8  *aslab = h_active_primes_all + (size_t)v * (size_t)n_qp;
                u32 n_act = 0;
                for (int i = 0; i < n_qp; i++) {
                    u32 q = q_primes[i];
                    u64 *mask = &mslab[(size_t)i * V16_Q_MASK_WORDS];
                    int any = 0;
                    for (int w = 0; w < V16_Q_MASK_WORDS; w++) {
                        if (mask[w]) { any = 1; break; }
                    }
                    u32 D_mod_q = (u32)(ps->D % (u128)q);
                    if (!any) {
                        if (D_mod_q != 0) {
                            fprintf(stderr,
                                "ERROR: forbid_mask all-zero for V=%d "
                                "q=%u but q does not divide D_V "
                                "-- mask builder bug\n", v, q);
                            return 1;
                        }
                        /* q | D_V → prune. */
                        continue;
                    }
                    aslab[n_act++] = (u8)i;
                }
                h_n_active[v] = (u8)n_act;
                sum_active += n_act;
                sum_total  += (u64)n_qp;
                if (n_act < min_active) min_active = n_act;
                if (n_act > max_active) max_active = n_act;
                u32 pruned = (u32)n_qp - n_act;
                if (pruned > max_pruned) max_pruned = pruned;
            }
            double mean_pruned =
                (double)(sum_total - sum_active) / (double)n_v_total;
            printf("#   q-iter active_primes: %.2f KiB "
                   "(%d V x %d u8) -- pruned mean=%.2f max=%u "
                   "(min_active=%u max_active=%u of %d)\n",
                   (ap_total + (size_t)n_v_total) / 1024.0,
                   n_v_total, n_qp,
                   mean_pruned, max_pruned,
                   min_active, max_active, n_qp);
        } else {
            printf("#   q-iter active_primes: DISABLED "
                   "(--no-active-primes) -- kernel iterates all %d primes\n",
                   n_qp);
        }

        /* ---- 5) Identity hashes for checkpoint. */
        u64 sieve_hash = q_iter_sieve_fnv1a64(q_primes, n_qp);
        u64 pool_hash  = q_iter_pool_fnv1a64(g_cfg.pool_seeds, n_v_total);
        printf("#   q-iter sieve fnv1a64 = 0x%016llx\n",
               (unsigned long long)sieve_hash);
        printf("#   q-iter pool  fnv1a64 = 0x%016llx\n",
               (unsigned long long)pool_hash);

        /* ---- 6) Resume cursors (sparse v_cursor). v_cursor[v] already
         * == q_min[v] from window-computation above; load_q_iter_checkpoint
         * overwrites only the V's present in the file.
         *
         * Checkpoint stores the *drained* cursor (rear edge of proved
         * coverage). At resume time launched and drained are reset to that
         * common value — we re-launch from the rear edge, so anything that
         * was launched-but-not-drained on the previous run is recomputed
         * (no chain loss). */
        if (g_cfg.resume_file[0]) {
            if (load_q_iter_checkpoint(&g_cfg, v_cursor_lo, v_cursor_hi,
                                       n_v_total,
                                       sieve_hash, pool_hash) != 0) {
                return 1;
            }
            /* Mirror loaded launched-cursor into drained-cursor. */
            for (int v = 0; v < n_v_total; v++) {
                v_cursor_drained_lo[v] = v_cursor_lo[v];
                v_cursor_drained_hi[v] = v_cursor_hi[v];
            }
            /* Any V advanced past Q_max during resume -> mark exhausted. */
            for (int v = 0; v < n_v_total; v++) {
                if (v_exhausted[v]) continue;
                u128 cur = ((u128)v_cursor_hi[v] << 64) | (u128)v_cursor_lo[v];
                u128 qmx = ((u128)q_max_hi[v]    << 64) | (u128)q_max_lo[v];
                if (cur > qmx) v_exhausted[v] = 1;
            }
            printf("#   q-iter resume = OK (per-V cursors loaded)\n");
        }

        /* ---- 7) GPU setup. */
        CUDA_CHECK(cudaSetDevice(g_cfg.gpu_id));

        /* Device-side sieve primes (one upload). */
        u32 *d_q_sieve_primes;
        CUDA_CHECK(cudaMalloc(&d_q_sieve_primes,
                              (size_t)n_qp * sizeof(u32)));
        CUDA_CHECK(cudaMemcpy(d_q_sieve_primes, q_primes,
                              (size_t)n_qp * sizeof(u32),
                              cudaMemcpyHostToDevice));

        /* Device-side Barrett table (one upload). */
        BarrettQ *d_q_barrett;
        CUDA_CHECK(cudaMalloc(&d_q_barrett,
                              (size_t)n_qp * sizeof(BarrettQ)));
        CUDA_CHECK(cudaMemcpy(d_q_barrett, q_barrett,
                              (size_t)n_qp * sizeof(BarrettQ),
                              cudaMemcpyHostToDevice));

        /* All-V forbid_mask resident on device (one upload at boot). */
        u64 *d_q_forbid_mask_all;
        CUDA_CHECK(cudaMalloc(&d_q_forbid_mask_all,
                              mask_total_u64 * sizeof(u64)));
        CUDA_CHECK(cudaMemcpy(d_q_forbid_mask_all, h_forbid_mask_all,
                              mask_total_u64 * sizeof(u64),
                              cudaMemcpyHostToDevice));

        /* All-V active-prime list resident on
         * device (one upload at boot). Per-launch we just re-aim the
         * kernel at the V's slab. NULL when disabled — kernel falls back
         * to the all-primes path. */
        u8 *d_active_primes_all = NULL;
        if (g_cfg.active_primes) {
            size_t ap_total = (size_t)n_v_total * (size_t)n_qp;
            CUDA_CHECK(cudaMalloc(&d_active_primes_all, ap_total));
            CUDA_CHECK(cudaMemcpy(d_active_primes_all, h_active_primes_all,
                                  ap_total, cudaMemcpyHostToDevice));
        }

        /* N-deep stream pool. N slots in a ring; we drain a slot
         * lazily — only when about to re-launch on that slot (i.e., the ring
         * wraps around). With N=4 (default), up to 4 K1 launches are queued
         * to the GPU before the oldest needs to drain, which hides per-launch
         * host overhead that capped gpu= at ~85% on the previous 2-slot path.
         * Arrays are statically sized to V16_MAX_STREAMS; runtime loops use
         * g_cfg.n_streams. */
        const int N_STREAMS = g_cfg.n_streams;
        cudaStream_t q_streams[V16_MAX_STREAMS];
        /* Per-slot cudaEvent pair times the
         * full GPU work on that slot (K1, plus K1.5 when --k15). Recorded
         * pre-launch and post-last-GPU-op (count D2H). cudaEventElapsedTime
         * is called only AFTER cudaStreamSynchronize on the slot, so it
         * never host-blocks the dispatch loop. Mirrors v15
         * (cc18_filter_cuda_CpC_v15.cu:386, 1337-1338, 4067, 4080, 4124). */
        cudaEvent_t k1_start_evt[V16_MAX_STREAMS], k1_end_evt[V16_MAX_STREAMS];
        u32 *d_q_base_mod[V16_MAX_STREAMS];
        QSurvivor *d_q_survivors[V16_MAX_STREAMS];
        u32 *d_q_survivor_count[V16_MAX_STREAMS];
        u32 *d_q_overflow[V16_MAX_STREAMS] = {NULL};
        u32  h_q_overflow[V16_MAX_STREAMS] = {0};
        QSurvivor *h_q_survivors_pinned[V16_MAX_STREAMS] = {NULL};
        u32 *h_q_base_mod_pinned[V16_MAX_STREAMS] = {NULL};
        u32 h_q_survivor_count[V16_MAX_STREAMS] = {0};
        int q_slot_active[V16_MAX_STREAMS] = {0};
        u16 q_slot_v[V16_MAX_STREAMS] = {0};
        /* Per-slot Q_base (u128 paired); passed to prove worker so it can
         * reconstruct Q_full = Q_base + thread_offset. */
        u64 q_slot_base_lo[V16_MAX_STREAMS] = {0};
        u64 q_slot_base_hi[V16_MAX_STREAMS] = {0};
        /* Per-slot Q chunk size. Recorded at launch so
         * the drain side can advance the drained-cursor for q_slot_v[slot]
         * to (q_slot_base + q_slot_chunk) without recomputing. */
        u32 q_slot_chunk[V16_MAX_STREAMS] = {0};
        /* K1.5 GPU top-MR pre-filter slots. Sized identically
         * to the K1 cap (V16_Q_MAX_SURVIVORS) since worst case is "K1.5 passes
         * everything K1 emitted". */
        K15Survivor *d_k15_survivors[V16_MAX_STREAMS] = {NULL};
        u32 *d_k15_count[V16_MAX_STREAMS] = {NULL};
        K15Survivor *h_k15_survivors_pinned[V16_MAX_STREAMS] = {NULL};
        u32 h_k15_count[V16_MAX_STREAMS] = {0};
        /* Track per-slot K1 survivor count
         * from the previous launch on this slot. The kernel self-limits via
         * *d_in_count, so over-/under-shooting is safe — just wasted launch
         * tail. Bootstrap = 256 blocks for the very first launch on a slot. */
        u32 prev_k1_count[V16_MAX_STREAMS] = {0};
        const u32 K15_GRID_BOOTSTRAP_BLOCKS = 256u;
        const u32 K15_GRID_CAP_BLOCKS =
            (u32)(((u64)V16_Q_MAX_SURVIVORS + V16_THREADS - 1) / V16_THREADS);
        printf("#   stream pool: N=%d × (32 MiB device survivor + 32 MiB pinned host + small)\n",
               N_STREAMS);
        /* --q-order random: seed PRNG. /dev/urandom by default; --q-seed N
         * forces a reproducible seed (debug). In exhaustive mode the seed is
         * still needed for the random-half of variants; the global
         * q_order_random is ignored for mode selection but still controls
         * --q-seed semantics. */
        int prng_needed = g_cfg.q_order_random ||
                          (g_cfg.q_band_mode_exhaustive && n_v_random > 0);
        if (prng_needed) {
            u64 used_seed;
            if (g_cfg.q_seed) {
                rng_seed(g_cfg.q_seed);
                used_seed = g_cfg.q_seed;
                printf("#   q-order   : random  seed=0x%016llx (source=--q-seed)\n",
                       (unsigned long long)used_seed);
            } else {
                used_seed = rng_seed_entropy(g_cfg.gpu_id);
                printf("#   q-order   : random  seed=0x%016llx (source=urandom)\n",
                       (unsigned long long)used_seed);
            }
        } else {
            printf("#   q-order   : sequential (every variant starts at Q_min, walks up)\n");
        }
        if (g_cfg.q_band_mode_exhaustive) {
            printf("#   q-band-mode: exhaustive  threshold=2^%d Q-range  "
                   "(sequential V's=%d  random V's=%d)\n",
                   g_cfg.exhaustive_max_q_bits, n_v_sequential, n_v_random);
        }
        printf("#   dedup-p0  : %s  (collapse repeat-p_0 emits to HIT-DUP)\n",
               g_cfg.dedup_p0 ? "on" : "off");
        for (int s = 0; s < N_STREAMS; s++) {
            CUDA_CHECK(cudaStreamCreate(&q_streams[s]));
            /* Per-slot K1 timing events. */
            CUDA_CHECK(cudaEventCreate(&k1_start_evt[s]));
            CUDA_CHECK(cudaEventCreate(&k1_end_evt[s]));
            CUDA_CHECK(cudaMalloc(&d_q_base_mod[s],
                                  (size_t)n_qp * sizeof(u32)));
            CUDA_CHECK(cudaMalloc(&d_q_survivors[s],
                                  (size_t)V16_Q_MAX_SURVIVORS * sizeof(QSurvivor)));
            CUDA_CHECK(cudaMalloc(&d_q_survivor_count[s], sizeof(u32)));
            /* Per-slot K1/K1.5 overflow flag. Memset
             * to 0 each launch alongside the survivor count; read after the
             * stream sync; abort on set. */
            CUDA_CHECK(cudaMalloc(&d_q_overflow[s], sizeof(u32)));
            CUDA_CHECK(cudaMallocHost((void **)&h_q_survivors_pinned[s],
                                      (size_t)V16_Q_MAX_SURVIVORS *
                                      sizeof(QSurvivor)));
            CUDA_CHECK(cudaMallocHost((void **)&h_q_base_mod_pinned[s],
                                      (size_t)n_qp * sizeof(u32)));
            if (g_cfg.k15_enabled) {
                CUDA_CHECK(cudaMalloc(&d_k15_survivors[s],
                                      (size_t)V16_Q_MAX_SURVIVORS *
                                      sizeof(K15Survivor)));
                CUDA_CHECK(cudaMalloc(&d_k15_count[s], sizeof(u32)));
                CUDA_CHECK(cudaMallocHost(
                    (void **)&h_k15_survivors_pinned[s],
                    (size_t)V16_Q_MAX_SURVIVORS * sizeof(K15Survivor)));
            }
        }
        unsigned long long *d_q_total, *d_q_pass;
        CUDA_CHECK(cudaMalloc(&d_q_total, sizeof(unsigned long long)));
        CUDA_CHECK(cudaMalloc(&d_q_pass,  sizeof(unsigned long long)));
        CUDA_CHECK(cudaMemset(d_q_total, 0, sizeof(unsigned long long)));
        CUDA_CHECK(cudaMemset(d_q_pass,  0, sizeof(unsigned long long)));
        /* K1.5 device-side stat counter + per-seed D_V table. */
        unsigned long long *d_k15_pass = NULL;
        u64 *d_k15_dv_lo = NULL;
        u64 *d_k15_dv_hi = NULL;
        if (g_cfg.k15_enabled) {
            CUDA_CHECK(cudaMalloc(&d_k15_pass, sizeof(unsigned long long)));
            CUDA_CHECK(cudaMemset(d_k15_pass, 0, sizeof(unsigned long long)));
            int npool = g_cfg.n_pool_seeds;
            u64 *h_dv_lo = (u64 *)calloc((size_t)npool, sizeof(u64));
            u64 *h_dv_hi = (u64 *)calloc((size_t)npool, sizeof(u64));
            if (!h_dv_lo || !h_dv_hi) {
                fprintf(stderr, "OOM K1.5 D_V table\n");
                return 1;
            }
            for (int v = 0; v < npool; v++) {
                u128 Dv = g_cfg.pool_seeds[v].D;
                h_dv_lo[v] = (u64)Dv;
                h_dv_hi[v] = (u64)(Dv >> 64);
            }
            CUDA_CHECK(cudaMalloc(&d_k15_dv_lo, (size_t)npool * sizeof(u64)));
            CUDA_CHECK(cudaMalloc(&d_k15_dv_hi, (size_t)npool * sizeof(u64)));
            CUDA_CHECK(cudaMemcpy(d_k15_dv_lo, h_dv_lo,
                                  (size_t)npool * sizeof(u64),
                                  cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy(d_k15_dv_hi, h_dv_hi,
                                  (size_t)npool * sizeof(u64),
                                  cudaMemcpyHostToDevice));
            free(h_dv_lo); free(h_dv_hi);
        }

        /* CPU prove output. */
        HostChainResult *h_chain_hits = (HostChainResult *)calloc(
            V16_MAX_RESULTS, sizeof(HostChainResult));
        u32 h_chain_hits_count = 0;
        u64 h_pass_top_cpu = 0, h_root_prime_cpu = 0, h_chains_emitted_cpu = 0;
        u64 h_chains_emit_dup_cpu = 0;
        pthread_mutex_t log_mutex = PTHREAD_MUTEX_INITIALIZER;
        u64 chain_hist[V16_HIST_MAX] = {0};

        FILE *log_fp = NULL;
        if (g_cfg.log_file[0]) {
            log_fp = fopen(g_cfg.log_file, "a");
            if (log_fp) fprintf(log_fp,
                "# v16 q-iter run start seed=%s D=%s target=%d E=%d bits=[%d,%d] pool=%d\n",
                g_cfg.seed_name, g_cfg.D_dec, g_cfg.target_len,
                g_cfg.exp_start, g_cfg.bits_min, g_cfg.bits_max, n_v_total);
        }

        /* Persistent prove worker pool. */
        ProvePool prove_pool;
        size_t pool_buf_stride =
            g_cfg.k15_enabled ? sizeof(K15Survivor) : sizeof(QSurvivor);
        if (prove_pool_init(&prove_pool, &g_cfg,
                            &h_pass_top_cpu, &h_root_prime_cpu,
                            &h_chains_emitted_cpu,
                            &h_chains_emit_dup_cpu,
                            h_chain_hits, V16_MAX_RESULTS, &h_chain_hits_count,
                            &log_mutex, log_fp, g_cfg.seed_name,
                            chain_hist, V16_HIST_MAX,
                            v_cursor_drained_lo, v_cursor_drained_hi,
                            pool_buf_stride, V16_Q_MAX_SURVIVORS) != 0) {
            fprintf(stderr, "FATAL: prove_pool_init failed\n");
            return 1;
        }

        struct timeval t_start, t_now, t_last_report;
        gettimeofday(&t_start, NULL);
        t_last_report = t_start;

        /* Cumulative timers + per-interval
         * baselines for the gpu= and prove/k1 fields in the periodic line.
         *   total_k1_ms     -- sum of cudaEventElapsedTime over drained
         *                      slots (full GPU work: K1 + K1.5 when on).
         *   total_prove_ms  -- sum of run_cpu_prove wall-clock; since
         *                      run_cpu_prove joins all worker threads,
         *                      its wall time == the parallel CPU bottleneck.
         *                      No per-worker mutex needed here.
         *   last_report_*_ms -- snapshot at the last [t=...] tick; the
         *                      reporter uses (current - last) to compute
         *                      *interval* gpu_util, and tracks the lowest
         *                      observed interval as min_gpu_util_pct.
         * Mirrors v15 lines 3997-3998, 4302-4309 and kt v8 §237-238. */
        double total_k1_ms = 0.0;
        double total_prove_ms = 0.0;
        double last_report_k1_ms = 0.0;
        double min_gpu_util_pct = 100.0;
        int    min_gpu_util_seen = 0;  /* set true once we've measured 1 interval */

        /* ---- 8) Main launch loop. Round-robin over active V's; per launch
         * emits one (V, Q_chunk). N-deep stream pool: cur_slot rotates through
         * N streams. A slot is drained lazily — only when the ring wraps
         * around and that slot already has work in flight. */
        const u32 Q_per_thread = 16;
        int cur_v = 0;
        u32 q_launches = 0;

        /* Rotation state. n_subbands>0 enables rotation; cur_subband_idx=-1
         * forces immediate rotation to sub-band 0 on first iteration. */
        struct timeval t_last_rotate;
        gettimeofday(&t_last_rotate, NULL);
        int n_subbands = (rotate_step > 0)
            ? (N2 - N1 + 1 + rotate_step - 1) / rotate_step
            : 0;
        int cur_subband_idx = -1;

        while (!g_stop) {
            /* Rotation tick. Cheap (one gettimeofday + compare) when off
             * or between rotations; ~ms GMP cost on transitions (every 30s
             * by default — negligible vs ~5h hit cadence). */
            if (rotate_step > 0) {
                struct timeval t_now_rot;
                gettimeofday(&t_now_rot, NULL);
                double dt_rot = (t_now_rot.tv_sec - t_last_rotate.tv_sec)
                              + (t_now_rot.tv_usec - t_last_rotate.tv_usec) * 1e-6;
                int do_rotate = (cur_subband_idx < 0) || (dt_rot >= (double)rotate_period);
                if (do_rotate) {
                    cur_subband_idx = (cur_subband_idx + 1) % n_subbands;
                    t_last_rotate = t_now_rot;
                    int sub_N1 = N1 + cur_subband_idx * rotate_step;
                    int sub_N2 = sub_N1 + rotate_step - 1;
                    if (sub_N2 > N2) sub_N2 = N2;
                    int sub_active = compute_q_windows_subband(
                        n_v_total, sub_N1, sub_N2, E, g_cfg.pool_seeds,
                        q_min_lo, q_min_hi, q_max_lo, q_max_hi,
                        v_exhausted, /*mark_exhausted=*/0);
                    if (sub_active < 0) { g_stop = 1; break; }
                    fprintf(stderr,
                        "rotate: sub-band [%d, %d]  "
                        "(active V: %d / %d, sub-band %d/%d)\n",
                        sub_N1, sub_N2, sub_active, n_v_total,
                        cur_subband_idx + 1, n_subbands);
                    /* Reset sequential cursors to fresh q_min. (Random ignores.) */
                    for (int v = 0; v < n_v_total; v++) {
                        if (v_exhausted[v]) continue;
                        v_cursor_lo[v]         = q_min_lo[v];
                        v_cursor_hi[v]         = q_min_hi[v];
                        v_cursor_drained_lo[v] = q_min_lo[v];
                        v_cursor_drained_hi[v] = q_min_hi[v];
                    }
                }
            }

            int v_found = -1;
            u128 cur = 0;        /* sequential cursor (unused for random V's) */
            u128 qmx = 0;
            u32  this_chunk = 0;
            u128 Q_base_pick = 0;
            /* Per-V mode branch. Pick next non-exhausted V common to both
             * modes, then dispatch on v_use_random[v] for the Q-walker.
             * In fixed mode v_use_random[] is uniform on g_cfg.q_order_random,
             * so this loop reproduces the legacy global-mode behaviour. */
            {
                int found_any = 0;
                for (int step = 0; step < n_v_total; step++) {
                    int v = (cur_v + step) % n_v_total;
                    if (v_exhausted[v]) continue;
                    u128 qmn       = ((u128)q_min_hi[v] << 64) | (u128)q_min_lo[v];
                    u128 qmx_local = ((u128)q_max_hi[v] << 64) | (u128)q_max_lo[v];
                    if (qmx_local < qmn) {
                        /* Empty window. Skip transient in rotation mode. */
                        if (rotate_step == 0) v_exhausted[v] = 1;
                        continue;
                    }
                    if (v_use_random[v]) {
                        qmx = qmx_local;
                        u128 span = (qmx_local - qmn) + (u128)1;
                        this_chunk = g_cfg.q_chunk;
                        if ((u128)this_chunk > span) this_chunk = (u32)span;
                        u128 range = (span - (u128)this_chunk) + (u128)1;
                        struct timespec _ts;
                        clock_gettime(CLOCK_MONOTONIC, &_ts);
                        rng_state ^= (u64)_ts.tv_nsec;
                        u128 rand128 = ((u128)rng_next() << 64) | (u128)rng_next();
                        Q_base_pick = qmn + (rand128 % range);
                    } else {
                        cur = ((u128)v_cursor_hi[v] << 64) | (u128)v_cursor_lo[v];
                        qmx = qmx_local;
                        if (cur > qmx) {
                            /* Sequential cursor past qmax: exhausted
                             * (or sub-band transient under rotation). */
                            if (rotate_step == 0) v_exhausted[v] = 1;
                            continue;
                        }
                        u128 remaining = (qmx - cur) + (u128)1;
                        this_chunk = g_cfg.q_chunk;
                        if (remaining < (u128)this_chunk) this_chunk = (u32)remaining;
                        if (this_chunk == 0) {
                            if (rotate_step == 0) v_exhausted[v] = 1;
                            continue;
                        }
                        Q_base_pick = cur;
                    }
                    v_found = v;
                    found_any = 1;
                    break;
                }
                if (!found_any) break;   /* all variants exhausted */
                cur_v = (v_found + 1) % n_v_total;
            }
            int v = v_found;

            int cur_slot = (int)(q_launches % (u32)N_STREAMS);
            u128 Q_base = Q_base_pick;

            /* N-deep drain-on-reuse: if this slot already has work pending
             * (the ring wrapped), drain it before re-launching. With N=4,
             * the slot's K1 has been running for ~N-1 iterations, so the
             * cudaStreamSynchronize below is typically near-immediate. */
            if (q_slot_active[cur_slot]) {
                int prev_slot = cur_slot;  /* local alias for the existing drain block */
                CUDA_CHECK(cudaStreamSynchronize(q_streams[prev_slot]));
                {
                    float slot_k1_ms = 0.0f;
                    cudaEventElapsedTime(&slot_k1_ms,
                                         k1_start_evt[prev_slot],
                                         k1_end_evt[prev_slot]);
                    total_k1_ms += (double)slot_k1_ms;
                }
                u32 n_surv_prev = h_q_survivor_count[prev_slot];
                /* K1.5 OOB guard. Either K1 or K1.5
                 * raised the flag → the survivor buffer is at capacity and
                 * K1.5 may already have skipped tail records. Abort before
                 * submitting a partial prove job. */
                if (h_q_overflow[prev_slot]) {
                    fprintf(stderr,
                        "ERROR: Q-iter survivor cap reached (slot=%d, "
                        "n_surv=%u, cap=%u). Aborting — lower --q-chunk or "
                        "raise V16_Q_MAX_SURVIVORS.\n",
                        prev_slot, n_surv_prev, V16_Q_MAX_SURVIVORS);
                    return 1;
                }
                if (n_surv_prev > V16_Q_MAX_SURVIVORS) {
                    fprintf(stderr,
                        "ERROR: Q-iter survivor overflow (%u > %u) — bug "
                        "(design §5 caps at 4M; expected ~2.2)\n",
                        n_surv_prev, V16_Q_MAX_SURVIVORS);
                    return 1;
                }
                prev_k1_count[prev_slot] = n_surv_prev;

                u128 _q_base = ((u128)q_slot_base_hi[prev_slot] << 64) |
                               (u128)q_slot_base_lo[prev_slot];
                u128 _drained_next = _q_base + (u128)q_slot_chunk[prev_slot];

                if (g_cfg.k15_enabled) {
                    u32 n_k15_prev = h_k15_count[prev_slot];
                    if (n_k15_prev > V16_Q_MAX_SURVIVORS) {
                        fprintf(stderr,
                            "ERROR: K1.5 survivor overflow (%u > %u)\n",
                            n_k15_prev, V16_Q_MAX_SURVIVORS);
                        return 1;
                    }
                    int buf_idx = prove_pool_acquire_buf(&prove_pool);
                    K15Survivor *buf = (K15Survivor *)
                        prove_pool_buf_ptr(&prove_pool, buf_idx);
                    if (n_k15_prev > 0) {
                        CUDA_CHECK(cudaMemcpyAsync(
                            buf, d_k15_survivors[prev_slot],
                            (size_t)n_k15_prev * sizeof(K15Survivor),
                            cudaMemcpyDeviceToHost, q_streams[prev_slot]));
                        CUDA_CHECK(cudaStreamSynchronize(q_streams[prev_slot]));
                    }
                    ProveJob job;
                    memset(&job, 0, sizeof(job));
                    job.buf_idx           = buf_idx;
                    job.n_surv            = n_k15_prev;
                    job.q_iter_mode       = 1;
                    job.k15_passed        = 1;
                    job.q_base_lo         = q_slot_base_lo[prev_slot];
                    job.q_base_hi         = q_slot_base_hi[prev_slot];
                    job.v_idx             = (int)q_slot_v[prev_slot];
                    job.drained_target_lo = (u64)_drained_next;
                    job.drained_target_hi = (u64)(_drained_next >> 64);
                    prove_pool_submit(&prove_pool, job);
                } else {
                    int buf_idx = prove_pool_acquire_buf(&prove_pool);
                    QSurvivor *buf = (QSurvivor *)
                        prove_pool_buf_ptr(&prove_pool, buf_idx);
                    if (n_surv_prev > 0) {
                        CUDA_CHECK(cudaMemcpyAsync(
                            buf, d_q_survivors[prev_slot],
                            (size_t)n_surv_prev * sizeof(QSurvivor),
                            cudaMemcpyDeviceToHost, q_streams[prev_slot]));
                        CUDA_CHECK(cudaStreamSynchronize(q_streams[prev_slot]));
                    }
                    ProveJob job;
                    memset(&job, 0, sizeof(job));
                    job.buf_idx           = buf_idx;
                    job.n_surv            = n_surv_prev;
                    job.q_iter_mode       = 1;
                    job.k15_passed        = 0;
                    job.q_base_lo         = q_slot_base_lo[prev_slot];
                    job.q_base_hi         = q_slot_base_hi[prev_slot];
                    job.v_idx             = (int)q_slot_v[prev_slot];
                    job.drained_target_lo = (u64)_drained_next;
                    job.drained_target_hi = (u64)(_drained_next >> 64);
                    prove_pool_submit(&prove_pool, job);
                }
                q_slot_active[prev_slot] = 0;
            }

            /* Per-launch precompute: Q_base_mod_q[i] = Q_base mod q[i].
             * Only n_qp u128%u32 ops (~86 total) -- the offset-form slab and
             * its depth*n_qp upload are gone. */
            u32 *base_mod_dst = h_q_base_mod_pinned[cur_slot];
            for (int i = 0; i < n_qp; i++) {
                base_mod_dst[i] = (u32)(Q_base % (u128)q_primes[i]);
            }

            /* Upload the small base_mod vector to cur_slot (async). */
            CUDA_CHECK(cudaMemcpyAsync(
                d_q_base_mod[cur_slot], base_mod_dst,
                (size_t)n_qp * sizeof(u32),
                cudaMemcpyHostToDevice, q_streams[cur_slot]));
            CUDA_CHECK(cudaMemsetAsync(d_q_survivor_count[cur_slot], 0,
                                       sizeof(u32),
                                       q_streams[cur_slot]));
            /* Reset the K1/K1.5 overflow flag for
             * this slot on the same stream so the memset precedes K1's
             * atomicExch. */
            CUDA_CHECK(cudaMemsetAsync(d_q_overflow[cur_slot], 0,
                                       sizeof(u32),
                                       q_streams[cur_slot]));
            if (g_cfg.k15_enabled) {
                /* Reset K1.5 output count. Same stream so the memset
                 * deterministically precedes K1.5's atomicAdd's. */
                CUDA_CHECK(cudaMemsetAsync(d_k15_count[cur_slot], 0,
                                           sizeof(u32),
                                           q_streams[cur_slot]));
            }

            /* Launch grid: one thread per Q_per_thread Q values. */
            u64 threads_needed =
                ((u64)this_chunk + Q_per_thread - 1) / Q_per_thread;
            u64 blocks_needed =
                (threads_needed + V16_THREADS - 1) / V16_THREADS;
            if (blocks_needed > (u64)0x7fffffffULL) {
                fprintf(stderr,
                    "ERROR: grid_x %llu exceeds 2^31 -- lower --q-chunk\n",
                    (unsigned long long)blocks_needed);
                return 1;
            }
            dim3 grid((u32)blocks_needed, 1, 1);
            dim3 block(V16_THREADS, 1, 1);
            /* Aim kernel at this V's mask slab inside the resident all-V mask. */
            const u64 *d_v_mask =
                d_q_forbid_mask_all + (size_t)v * mask_slab_u64;
            /* Aim at V's active-prime slab. NULL
             * pointer + n_active=0 selects the legacy all-primes path. */
            const u8 *d_v_active = NULL;
            u32       v_n_active = 0;
            if (g_cfg.active_primes) {
                d_v_active = d_active_primes_all + (size_t)v * (size_t)n_qp;
                v_n_active = (u32)h_n_active[v];
            }
            /* Mark start of GPU work for this
             * slot. Recorded on the slot's stream so it is ordered between
             * the small memsets above and the K1 kernel that follows -- but
             * since events are async, the dispatch loop is not blocked. */
            CUDA_CHECK(cudaEventRecord(k1_start_evt[cur_slot],
                                       q_streams[cur_slot]));
            v16_q_iter_kernel<<<grid, block, 0, q_streams[cur_slot]>>>(
                this_chunk, Q_per_thread, (u32)v,
                (u32)n_qp, d_q_sieve_primes, d_q_barrett,
                d_v_mask, d_q_base_mod[cur_slot],
                (u32)depth,
                d_q_survivors[cur_slot], d_q_survivor_count[cur_slot],
                V16_Q_MAX_SURVIVORS,
                d_q_total, d_q_pass,
                d_v_active, v_n_active,
                d_q_overflow[cur_slot]);
            CUDA_CHECK(cudaGetLastError());

            if (g_cfg.k15_enabled) {
                /* K1.5 launch: grid sized to
                 * the K1 survivor count observed on this slot's previous
                 * launch (or a small bootstrap if we haven't seen one yet).
                 * The kernel reads *d_in_count and self-limits, so the grid
                 * is just an upper bound on parallelism. Same stream as K1,
                 * so the count load is deterministic without a host sync
                 * the same stream order. */
                u32 k15_blocks_u32;
                if (prev_k1_count[cur_slot] == 0) {
                    k15_blocks_u32 = K15_GRID_BOOTSTRAP_BLOCKS;
                } else {
                    u64 want = ((u64)prev_k1_count[cur_slot] +
                                V16_THREADS - 1) / V16_THREADS;
                    if (want > (u64)K15_GRID_CAP_BLOCKS)
                        want = (u64)K15_GRID_CAP_BLOCKS;
                    if (want < 1) want = 1;
                    k15_blocks_u32 = (u32)want;
                }
                dim3 k15_grid(k15_blocks_u32, 1, 1);
                v16_k15_topmr_kernel<<<k15_grid, block, 0,
                                       q_streams[cur_slot]>>>(
                    d_q_survivors[cur_slot], d_q_survivor_count[cur_slot],
                    d_k15_dv_lo, d_k15_dv_hi,
                    (u64)Q_base, (u64)(Q_base >> 64),
                    g_cfg.exp_start, g_cfg.target_len - 1,
                    g_cfg.k15_mr_reps,
                    d_k15_survivors[cur_slot], d_k15_count[cur_slot],
                    V16_Q_MAX_SURVIVORS,
                    d_k15_pass,
                    V16_Q_MAX_SURVIVORS,
                    d_q_overflow[cur_slot]);
                CUDA_CHECK(cudaGetLastError());

                /* Copy K1.5 count back (K1 count is still useful as a
                 * stat; we copy it too — both are cheap async). */
                CUDA_CHECK(cudaMemcpyAsync(&h_k15_count[cur_slot],
                                           d_k15_count[cur_slot],
                                           sizeof(u32),
                                           cudaMemcpyDeviceToHost,
                                           q_streams[cur_slot]));
            }

            /* Async copy survivor count back. */
            CUDA_CHECK(cudaMemcpyAsync(&h_q_survivor_count[cur_slot],
                                       d_q_survivor_count[cur_slot],
                                       sizeof(u32),
                                       cudaMemcpyDeviceToHost,
                                       q_streams[cur_slot]));
            /* Async copy the overflow flag back on
             * the same stream so the host sees it after the K1+K1.5 chain
             * has fully serialized. Read at slot-drain time below. */
            CUDA_CHECK(cudaMemcpyAsync(&h_q_overflow[cur_slot],
                                       d_q_overflow[cur_slot],
                                       sizeof(u32),
                                       cudaMemcpyDeviceToHost,
                                       q_streams[cur_slot]));
            /* Mark end of GPU work for this
             * slot, AFTER the survivor-count D2H so the elapsed delta covers
             * the full K1 (+ K1.5 + count copy) latency. Async record on
             * the slot's stream. cudaEventElapsedTime is only invoked once
             * we've synchronized this slot (in the drain branch below), so
             * it never blocks the launch pipeline. */
            CUDA_CHECK(cudaEventRecord(k1_end_evt[cur_slot],
                                       q_streams[cur_slot]));

            q_slot_active[cur_slot] = 1;
            q_slot_v[cur_slot] = (u16)v;
            q_slot_base_lo[cur_slot] = (u64)Q_base;
            q_slot_base_hi[cur_slot] = (u64)(Q_base >> 64);
            /* Fix 3: record chunk so drained-cursor advance is exact. */
            q_slot_chunk[cur_slot] = this_chunk;

            /* Advance v_cursor by this_chunk (sequential V's only). Random V's
             * don't track per-V coverage. Per-V flag, not the global q_order. */
            if (!v_use_random[v]) {
                u128 new_cur = cur + (u128)this_chunk;
                v_cursor_lo[v] = (u64)new_cur;
                v_cursor_hi[v] = (u64)(new_cur >> 64);
                if (new_cur > qmx) v_exhausted[v] = 1;
            }

            /* N-deep: no per-iteration drain after launch. The drain happens
             * at the TOP of the next iteration that maps cur_slot back to
             * this slot (after N-1 other launches have queued ahead). */
            q_launches++;

            /* Periodic reporting + checkpoint. */
            gettimeofday(&t_now, NULL);
            double now_dt = (t_now.tv_sec - t_last_report.tv_sec) +
                            (t_now.tv_usec - t_last_report.tv_usec) / 1e6;
            if (now_dt >= g_cfg.report_sec) {
                unsigned long long qt = 0, qp = 0;
                CUDA_CHECK(cudaMemcpy(&qt, d_q_total,
                                      sizeof(unsigned long long),
                                      cudaMemcpyDeviceToHost));
                CUDA_CHECK(cudaMemcpy(&qp, d_q_pass,
                                      sizeof(unsigned long long),
                                      cudaMemcpyDeviceToHost));
                /* Read K1.5 device-side pass counter (if enabled). */
                unsigned long long k15p = 0;
                if (g_cfg.k15_enabled && d_k15_pass) {
                    CUDA_CHECK(cudaMemcpy(&k15p, d_k15_pass,
                                          sizeof(unsigned long long),
                                          cudaMemcpyDeviceToHost));
                }
                double elapsed = (t_now.tv_sec - t_start.tv_sec) +
                                 (t_now.tv_usec - t_start.tv_usec) / 1e6;
                /* Active V count for live report. */
                int n_active_now = 0;
                for (int vi = 0; vi < n_v_total; vi++)
                    if (!v_exhausted[vi]) n_active_now++;
                /* Compute per-interval GPU
                 * utilization from the slot-event deltas accumulated since
                 * the last report. Track running min across intervals.
                 * utilization from slot-event timing. */
                double interval_k1_ms    = total_k1_ms    - last_report_k1_ms;
                double interval_wall_ms  = now_dt * 1000.0;
                int    gpu_known = (interval_k1_ms > 0.0 && interval_wall_ms > 0.0);
                double interval_gpu_pct = 0.0;
                if (gpu_known) {
                    interval_gpu_pct = interval_k1_ms / interval_wall_ms * 100.0;
                    if (interval_gpu_pct > 100.0) interval_gpu_pct = 100.0;
                    if (!min_gpu_util_seen || interval_gpu_pct < min_gpu_util_pct) {
                        min_gpu_util_pct  = interval_gpu_pct;
                        min_gpu_util_seen = 1;
                    }
                }
                /* gpu=NN%(min=MM%) string. */
                char gpu_str[64];
                if (gpu_known) {
                    snprintf(gpu_str, sizeof(gpu_str),
                             "gpu=%.0f%%(min=%.0f%%)",
                             interval_gpu_pct, min_gpu_util_pct);
                } else {
                    snprintf(gpu_str, sizeof(gpu_str),
                             "gpu=N/A(min=N/A)");
                }
                /* prove/k1=X.Xx string (only when --cpu-prove).
                 * Pool accumulates prove_ms under
                 * its own mutex; sample it for the live report. */
                total_prove_ms = prove_pool_total_prove_ms(&prove_pool);
                char prove_str[64];
                prove_str[0] = '\0';
                if (g_cfg.cpu_prove) {
                    if (total_k1_ms > 0.0) {
                        snprintf(prove_str, sizeof(prove_str),
                                 "  prove/k1=%.1fx",
                                 total_prove_ms / total_k1_ms);
                    } else {
                        snprintf(prove_str, sizeof(prove_str),
                                 "  prove/k1=N/A");
                    }
                }
                /* s/c=X.XXe-N string. */
                char sc_str[48];
                if (qt > 0) {
                    snprintf(sc_str, sizeof(sc_str),
                             "s/c=%.2e", (double)qp / (double)qt);
                } else {
                    snprintf(sc_str, sizeof(sc_str), "s/c=N/A");
                }
                if (g_cfg.k15_enabled) {
                    printf("[t=%.2fs] launches=%u  Q_eval=%llu (%.2fG (Q,V)/s)  "
                           "K1=%llu K15=%llu (drop=%.2fx)  %s%s  %s  "
                           "hits>=%d=%u  emitted=%llu  dups=%llu  active_V=%d/%d\n",
                           elapsed, q_launches, qt,
                           qt / 1e9 / (elapsed > 0 ? elapsed : 1),
                           qp, k15p,
                           k15p > 0 ? (double)qp / (double)k15p : 0.0,
                           gpu_str, prove_str, sc_str,
                           g_cfg.min_report_len, h_chain_hits_count,
                           h_chains_emitted_cpu, h_chains_emit_dup_cpu,
                           n_active_now, n_v_total);
                } else {
                    printf("[t=%.2fs] launches=%u  Q_eval=%llu (%.2fG (Q,V)/s)  "
                           "survivors=%llu  %s%s  %s  hits>=%d=%u  "
                           "emitted=%llu  dups=%llu  active_V=%d/%d\n",
                           elapsed, q_launches, qt,
                           qt / 1e9 / (elapsed > 0 ? elapsed : 1),
                           qp,
                           gpu_str, prove_str, sc_str,
                           g_cfg.min_report_len, h_chain_hits_count,
                           h_chains_emitted_cpu, h_chains_emit_dup_cpu,
                           n_active_now, n_v_total);
                }
                last_report_k1_ms    = total_k1_ms;
                t_last_report = t_now;
                if (g_cfg.checkpoint_file[0]) {
                    /* Checkpoint the *drained* cursor.
                     * Drain pool first so every
                     * in-flight job's drained-cursor advance has landed. */
                    prove_pool_drain_all(&prove_pool);
                    save_q_iter_checkpoint(&g_cfg,
                                           v_cursor_drained_lo,
                                           v_cursor_drained_hi,
                                           q_min_lo, q_min_hi,
                                           n_v_total,
                                           sieve_hash, pool_hash);
                }
            }
            if (g_cfg.time_limit_sec) {
                double el = (t_now.tv_sec - t_start.tv_sec) +
                            (t_now.tv_usec - t_start.tv_usec) / 1e6;
                if (el >= g_cfg.time_limit_sec) break;
            }
        }

        /* Drain ALL still-active slots after the main loop exits. With N-deep
         * streams there can be up to N-1 in-flight slots that never got
         * re-launched (the loop ended before wrapping). Iterate the entire
         * ring; skip the inactive ones. */
        for (int prev_slot = 0; prev_slot < N_STREAMS; prev_slot++) {
            if (!q_slot_active[prev_slot]) continue;
            CUDA_CHECK(cudaStreamSynchronize(q_streams[prev_slot]));
            /* Collect GPU time for this slot. */
            {
                float slot_k1_ms = 0.0f;
                cudaEventElapsedTime(&slot_k1_ms,
                                     k1_start_evt[prev_slot],
                                     k1_end_evt[prev_slot]);
                total_k1_ms += (double)slot_k1_ms;
            }
            u32 n_surv_prev = h_q_survivor_count[prev_slot];
            /* K1.5 OOB guard (final-drain mirror). */
            if (h_q_overflow[prev_slot]) {
                fprintf(stderr,
                    "ERROR: Q-iter survivor cap reached at final drain "
                    "(slot=%d, n_surv=%u, cap=%u)\n",
                    prev_slot, n_surv_prev, V16_Q_MAX_SURVIVORS);
                return 1;
            }
            if (n_surv_prev > V16_Q_MAX_SURVIVORS) {
                fprintf(stderr,
                    "ERROR: Q-iter survivor overflow (%u > %u) — bug "
                    "(design §5 caps at 4M; expected ~2.2)\n",
                    n_surv_prev, V16_Q_MAX_SURVIVORS);
                return 1;
            }
            /* Final-drain mirrors in-loop submit. */
            u128 _q_base = ((u128)q_slot_base_hi[prev_slot] << 64) |
                           (u128)q_slot_base_lo[prev_slot];
            u128 _drained_next = _q_base + (u128)q_slot_chunk[prev_slot];

            if (g_cfg.k15_enabled) {
                u32 n_k15_prev = h_k15_count[prev_slot];
                if (n_k15_prev > V16_Q_MAX_SURVIVORS) {
                    fprintf(stderr,
                        "ERROR: K1.5 survivor overflow (%u > %u)\n",
                        n_k15_prev, V16_Q_MAX_SURVIVORS);
                    return 1;
                }
                int buf_idx = prove_pool_acquire_buf(&prove_pool);
                K15Survivor *buf = (K15Survivor *)
                    prove_pool_buf_ptr(&prove_pool, buf_idx);
                if (n_k15_prev > 0) {
                    CUDA_CHECK(cudaMemcpyAsync(
                        buf, d_k15_survivors[prev_slot],
                        (size_t)n_k15_prev * sizeof(K15Survivor),
                        cudaMemcpyDeviceToHost, q_streams[prev_slot]));
                    CUDA_CHECK(cudaStreamSynchronize(q_streams[prev_slot]));
                }
                ProveJob job;
                memset(&job, 0, sizeof(job));
                job.buf_idx           = buf_idx;
                job.n_surv            = n_k15_prev;
                job.q_iter_mode       = 1;
                job.k15_passed        = 1;
                job.q_base_lo         = q_slot_base_lo[prev_slot];
                job.q_base_hi         = q_slot_base_hi[prev_slot];
                job.v_idx             = (int)q_slot_v[prev_slot];
                job.drained_target_lo = (u64)_drained_next;
                job.drained_target_hi = (u64)(_drained_next >> 64);
                prove_pool_submit(&prove_pool, job);
            } else {
                int buf_idx = prove_pool_acquire_buf(&prove_pool);
                QSurvivor *buf = (QSurvivor *)
                    prove_pool_buf_ptr(&prove_pool, buf_idx);
                if (n_surv_prev > 0) {
                    CUDA_CHECK(cudaMemcpyAsync(
                        buf, d_q_survivors[prev_slot],
                        (size_t)n_surv_prev * sizeof(QSurvivor),
                        cudaMemcpyDeviceToHost, q_streams[prev_slot]));
                    CUDA_CHECK(cudaStreamSynchronize(q_streams[prev_slot]));
                }
                ProveJob job;
                memset(&job, 0, sizeof(job));
                job.buf_idx           = buf_idx;
                job.n_surv            = n_surv_prev;
                job.q_iter_mode       = 1;
                job.k15_passed        = 0;
                job.q_base_lo         = q_slot_base_lo[prev_slot];
                job.q_base_hi         = q_slot_base_hi[prev_slot];
                job.v_idx             = (int)q_slot_v[prev_slot];
                job.drained_target_lo = (u64)_drained_next;
                job.drained_target_hi = (u64)(_drained_next >> 64);
                prove_pool_submit(&prove_pool, job);
            }
            q_slot_active[prev_slot] = 0;
        }

        /* Drain everything before final summary. */
        prove_pool_drain_all(&prove_pool);
        total_prove_ms = prove_pool_total_prove_ms(&prove_pool);

        /* Final summary. */
        {
            unsigned long long qt = 0, qp = 0;
            CUDA_CHECK(cudaMemcpy(&qt, d_q_total,
                                  sizeof(unsigned long long),
                                  cudaMemcpyDeviceToHost));
            CUDA_CHECK(cudaMemcpy(&qp, d_q_pass,
                                  sizeof(unsigned long long),
                                  cudaMemcpyDeviceToHost));
            gettimeofday(&t_now, NULL);
            double elapsed = (t_now.tv_sec - t_start.tv_sec) +
                             (t_now.tv_usec - t_start.tv_usec) / 1e6;
            printf("\n=== SUMMARY (q-iter) ===\n");
            printf("  elapsed       : %.2fs\n", elapsed);
            printf("  launches      : %u\n", q_launches);
            printf("  Q evaluated   : %llu\n", qt);
            printf("  (Q,V) rate    : %.2f G/sec\n",
                   qt / 1e9 / (elapsed > 0 ? elapsed : 1));
            printf("  sieve_pass    : %llu\n", qp);
            if (g_cfg.k15_enabled && d_k15_pass) {
                unsigned long long k15p = 0;
                CUDA_CHECK(cudaMemcpy(&k15p, d_k15_pass,
                                      sizeof(unsigned long long),
                                      cudaMemcpyDeviceToHost));
                printf("  k15_pass      : %llu (drop=%.2fx, reps=%d)\n",
                       k15p,
                       k15p > 0 ? (double)qp / (double)k15p : 0.0,
                       g_cfg.k15_mr_reps);
            }
            printf("  top_pass      : %llu\n", h_pass_top_cpu);
            printf("  root_prime    : %llu\n", h_root_prime_cpu);
            printf("  chains>=%d    : %u\n",
                   g_cfg.min_report_len, h_chain_hits_count);
            printf("  chains_emitted: %llu (unique p_0)\n", h_chains_emitted_cpu);
            printf("  chains_dup    : %llu (HIT-DUP collapses)\n", h_chains_emit_dup_cpu);
            if (g_cfg.n_pool_seeds > 0) {
                printf("  pool-seed hits (D | (p+1)):\n");
                for (int i = 0; i < g_cfg.n_pool_seeds; i++) {
                    v16_pool_seed *ps = &g_cfg.pool_seeds[i];
                    if (ps->hits == 0) continue;
                    printf("    %-32s  %10llu\n", ps->name,
                           (unsigned long long)ps->hits);
                }
            }
            printf("  chain-len histogram (CC1+):\n");
            int max_seen = 0;
            for (int i = 1; i < V16_HIST_MAX; i++)
                if (chain_hist[i]) max_seen = i;
            for (int i = 1; i <= max_seen; i++) {
                if (chain_hist[i] == 0) continue;
                u64 cum_high = 0;
                for (int j = i; j < V16_HIST_MAX; j++) cum_high += chain_hist[j];
                printf("    CC%-3d  count=%-10llu  >=CC%d=%llu\n",
                       i, (unsigned long long)chain_hist[i],
                       i, (unsigned long long)cum_high);
            }
            /* Final summary timing block. */
            double wall_s   = elapsed;
            double k1_s     = total_k1_ms / 1000.0;
            double prove_s  = total_prove_ms / 1000.0;
            double mean_gpu = (wall_s > 0.0)
                                ? (total_k1_ms / (wall_s * 1000.0) * 100.0)
                                : 0.0;
            if (mean_gpu > 100.0) mean_gpu = 100.0;
            double sc_ovr = (qt > 0) ? (double)qp / (double)qt : 0.0;
            printf("\n=== run summary ===\n");
            printf("  total wall:        %.1f s\n", wall_s);
            printf("  total K1 kernel:   %.1f s\n", k1_s);
            if (g_cfg.cpu_prove)
                printf("  total CPU prove:   %.1f s\n", prove_s);
            printf("  gpu_util_pct:      %.2f (mean over run)\n", mean_gpu);
            if (min_gpu_util_seen)
                printf("  gpu_util_pct_min:  %.2f (worst interval)\n",
                       min_gpu_util_pct);
            else
                printf("  gpu_util_pct_min:  N/A (no full interval)\n");
            if (g_cfg.cpu_prove && total_k1_ms > 0.0)
                printf("  prove_to_k1:       %.2fx\n",
                       total_prove_ms / total_k1_ms);
            if (qt > 0)
                printf("  s/c overall:       %.2e\n", sc_ovr);
            else
                printf("  s/c overall:       N/A\n");
            printf("===================\n");
        }

        /* Final checkpoint. */
        if (g_cfg.checkpoint_file[0]) {
            /* Final checkpoint also writes the drained cursor.
             * After the final-drain block above, drained == launched for
             * every slot's V, so this checkpoint records the true tail
             * of proved coverage. */
            save_q_iter_checkpoint(&g_cfg,
                                   v_cursor_drained_lo,
                                   v_cursor_drained_hi,
                                   q_min_lo, q_min_hi,
                                   n_v_total,
                                   sieve_hash, pool_hash);
        }
        if (log_fp) fclose(log_fp);

        /* Cleanup — N-deep: free all N slots. */
        for (int s = 0; s < N_STREAMS; s++) {
            cudaFree(d_q_base_mod[s]);
            cudaFree(d_q_survivors[s]);
            cudaFree(d_q_survivor_count[s]);
            if (d_q_overflow[s]) cudaFree(d_q_overflow[s]);
            if (h_q_survivors_pinned[s]) cudaFreeHost(h_q_survivors_pinned[s]);
            if (h_q_base_mod_pinned[s])
                cudaFreeHost(h_q_base_mod_pinned[s]);
            /* K1.5 buffers (only allocated when k15_enabled). */
            if (d_k15_survivors[s]) cudaFree(d_k15_survivors[s]);
            if (d_k15_count[s])     cudaFree(d_k15_count[s]);
            if (h_k15_survivors_pinned[s])
                cudaFreeHost(h_k15_survivors_pinned[s]);
            /* Tear down timing events. */
            cudaEventDestroy(k1_start_evt[s]);
            cudaEventDestroy(k1_end_evt[s]);
            cudaStreamDestroy(q_streams[s]);
        }
        cudaFree(d_q_sieve_primes);
        cudaFree(d_q_barrett);
        cudaFree(d_q_forbid_mask_all);
        cudaFree(d_q_total);
        cudaFree(d_q_pass);
        if (d_k15_pass) cudaFree(d_k15_pass);
        if (d_k15_dv_lo) cudaFree(d_k15_dv_lo);
        if (d_k15_dv_hi) cudaFree(d_k15_dv_hi);
        /* Free per-V active-prime list. */
        if (d_active_primes_all) cudaFree(d_active_primes_all);
        free(h_forbid_mask_all);
        if (h_active_primes_all) free(h_active_primes_all);
        if (h_n_active)          free(h_n_active);
        free(q_barrett);
        free(q_min_lo); free(q_min_hi);
        free(q_max_lo); free(q_max_hi);
        free(v_exhausted);
        free(v_cursor_lo); free(v_cursor_hi);
        /* Destroy the pool before freeing
         * v_cursor_drained_* (workers reference these). */
        prove_pool_destroy(&prove_pool);
        free(v_cursor_drained_lo); free(v_cursor_drained_hi);
        free(h_chain_hits);
        pthread_mutex_destroy(&log_mutex);
        return 0;
    }

    v16_wheel wh;
    int rc = build_wheel(&g_cfg, &wh);
    if (rc != 0) {
        fprintf(stderr, "ERROR: build_wheel failed rc=%d\n", rc);
        return 1;
    }
    printf("#   coeff primorial = %llu\n", (unsigned long long)wh.primorial);
    printf("#   wheel n_admissible = %u\n", wh.n_admissible);
    if (g_cfg.n_pool_primes > 0) {
        printf("#   pool-primes: ");
        for (int i = 0; i < g_cfg.n_pool_primes; i++)
            printf("%s%u", i ? "," : "", g_cfg.pool_primes[i]);
        printf("\n");
    }
    if (g_cfg.n_pool_seeds > 0) {
        printf("#   pool-seeds: %d\n", g_cfg.n_pool_seeds);
        for (int i = 0; i < g_cfg.n_pool_seeds; i++) {
            v16_pool_seed *ps = &g_cfg.pool_seeds[i];
            printf("#     %s  D=%s  primes=", ps->name, ps->D_dec);
            for (int j = 0; j < ps->n_primes; j++)
                printf("%s%u", j ? "*" : "", ps->primes[j]);
            printf("\n");
        }
    }
    printf("#   wheel fnv1a64 = 0x%016llx\n", (unsigned long long)wh.fnv1a64);

    if (g_cfg.target_bits == 0 && g_cfg.max_tiles == 0) {
        fprintf(stderr,
            "WARN: neither --bits nor --max-tiles given; running until SIGINT.\n");
    }

    /* Build sieve masks (host). */
    size_t mask_words = (size_t)g_cfg.n_sieve_primes * V16_PACKED_WORDS_MAX;
    u64 *h_masks = (u64 *)calloc(mask_words, sizeof(u64));
    if (!h_masks) { fprintf(stderr, "OOM masks\n"); return 1; }
    if (build_sieve_masks(&g_cfg, h_masks) != 0) {
        free(h_masks); wheel_free(&wh); return 1;
    }

    /* Upload to GPU. */
    CUDA_CHECK(cudaSetDevice(g_cfg.gpu_id));
    u64 D_lo = (u64)g_cfg.D;
    u64 D_hi = (u64)(g_cfg.D >> 64);
    CUDA_CHECK(cudaMemcpyToSymbol(d_D_lo,            &D_lo,                 sizeof(u64)));
    CUDA_CHECK(cudaMemcpyToSymbol(d_D_hi,            &D_hi,                 sizeof(u64)));
    CUDA_CHECK(cudaMemcpyToSymbol(d_coeff_primorial, &wh.primorial,         sizeof(u64)));
    CUDA_CHECK(cudaMemcpyToSymbol(d_exp_start,       &g_cfg.exp_start,      sizeof(int)));
    CUDA_CHECK(cudaMemcpyToSymbol(d_target_len,      &g_cfg.target_len,     sizeof(int)));
    CUDA_CHECK(cudaMemcpyToSymbol(d_target_bits,     &g_cfg.target_bits,    sizeof(int)));
    CUDA_CHECK(cudaMemcpyToSymbol(d_prove_forward,   &g_cfg.prove_forward,  sizeof(int)));
    CUDA_CHECK(cudaMemcpyToSymbol(d_min_report_len,  &g_cfg.min_report_len, sizeof(int)));
    CUDA_CHECK(cudaMemcpyToSymbol(d_n_sieve_primes,  &g_cfg.n_sieve_primes, sizeof(int)));
    CUDA_CHECK(cudaMemcpyToSymbol(d_n_admissible_const, &wh.n_admissible,   sizeof(u32)));
    CUDA_CHECK(cudaMemcpyToSymbol(d_sieve_primes,    g_cfg.sieve_primes,
                                  (size_t)g_cfg.n_sieve_primes * sizeof(u32)));

    /* coeff_primorial mod q for each sieve prime — used for per-tile residue
     * setup and per-launch tile_base_r computation. */
    u32 h_cp_mod_q[V16_MAX_SIEVE] = {0};
    for (int i = 0; i < g_cfg.n_sieve_primes; i++) {
        h_cp_mod_q[i] = (u32)(wh.primorial % g_cfg.sieve_primes[i]);
    }
    CUDA_CHECK(cudaMemcpyToSymbol(d_cp_mod_q, h_cp_mod_q,
                                  (size_t)g_cfg.n_sieve_primes * sizeof(u32)));

    /* Trial primes for prove kernel = sieve primes union {3..71} (small primes
     * already filtered by the wheel for D≠q, but we keep an explicit list to
     * harden the prove path against numerical drift). */
    u32 trial[V16_TRIAL_PRIMES_MAX];
    int n_trial = 0;
    for (size_t i = 0; i < V16_SMALL_PRIME_COUNT && n_trial < V16_TRIAL_PRIMES_MAX; i++) {
        u32 q = V16_SMALL_PRIMES[i];
        if (q > 71) break;
        int dividesD = 0;
        for (int j = 0; j < g_cfg.n_D_primes; j++)
            if (g_cfg.D_primes[j] == q) { dividesD = 1; break; }
        if (!dividesD) trial[n_trial++] = q;
    }
    CUDA_CHECK(cudaMemcpyToSymbol(d_trial_primes, trial, (size_t)n_trial * sizeof(u32)));
    CUDA_CHECK(cudaMemcpyToSymbol(d_n_trial_primes, &n_trial, sizeof(int)));

    u64 *d_offsets;
    u64 *d_masks;
    CUDA_CHECK(cudaMalloc(&d_offsets, (size_t)wh.n_admissible * sizeof(u64)));
    CUDA_CHECK(cudaMemcpy(d_offsets, wh.offsets,
                          (size_t)wh.n_admissible * sizeof(u64),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMalloc(&d_masks, mask_words * sizeof(u64)));
    CUDA_CHECK(cudaMemcpy(d_masks, h_masks,
                          mask_words * sizeof(u64), cudaMemcpyHostToDevice));

    /* --sieve-variants: per-variant packed forbidden masks. For each pool
     * seed V and each sieve prime q_i, recompute F_q_i(D_V) directly using
     * the variant's full D_V (instead of bit-rotating sieve_mask). Same
     * code path as build_sieve_masks; just substitutes D_V for cfg->D. */
    u32   n_variants    = g_cfg.sieve_variants ? (u32)g_cfg.n_pool_seeds : 0u;
    u64  *d_variant_masks = NULL;
    u8   *d_variant_skip  = NULL;
    unsigned long long *d_variant_sieve_pass = NULL;
    if (g_cfg.sieve_variants) {
        size_t vrow = (size_t)g_cfg.n_sieve_primes * V16_PACKED_WORDS_MAX;
        size_t vwords = (size_t)n_variants * vrow;
        u64 *h_vmasks = (u64 *)calloc(vwords, sizeof(u64));
        if (!h_vmasks) {
            fprintf(stderr, "OOM variant masks (%.2f MiB)\n",
                    vwords * sizeof(u64) / (1024.0 * 1024.0));
            return 1;
        }
        /* Per-(variant, sieve prime) skip flag: 1 iff q_i | D_V (variant
         * immune at q_i). Used by the kernel to short-circuit the per-prime
         * mask load when the mask row is all zeros. <=50 KB at 1024x50. */
        size_t vskip_bytes = (size_t)n_variants * (size_t)g_cfg.n_sieve_primes;
        u8 *h_vskip = (u8 *)calloc(vskip_bytes, 1);
        if (!h_vskip) {
            fprintf(stderr, "OOM variant skip flags (%.2f KiB)\n",
                    vskip_bytes / 1024.0);
            free(h_vmasks); return 1;
        }
        for (u32 v = 0; v < n_variants; v++) {
            v16_pool_seed *ps = &g_cfg.pool_seeds[v];
            for (int i = 0; i < g_cfg.n_sieve_primes; i++) {
                u32 q = g_cfg.sieve_primes[i];
                if ((ps->D % (u128)q) == 0) {
                    h_vskip[(size_t)v * (size_t)g_cfg.n_sieve_primes
                            + (size_t)i] = 1;
                }
                u32 fbd[V16_MAX_TARGET + 4];
                int nf = v16_forbidden_residues(
                    q, ps->D, g_cfg.exp_start, g_cfg.sieve_depth,
                    fbd, (int)(sizeof(fbd) / sizeof(fbd[0])));
                if (nf < 0) {
                    /* q | D_V — variant is fully immune at q. Mask row
                     * stays all zeros; skip flag set above lets the kernel
                     * short-circuit the load. */
                    continue;
                }
                if (nf >= (int)q) {
                    fprintf(stderr,
                        "ERROR: variant %s inadmissible at sieve prime %u "
                        "(|F_q|=%d=q)\n", ps->name, q, nf);
                    free(h_vmasks); free(h_vskip); return 1;
                }
                u64 *row = &h_vmasks[(size_t)v * vrow
                                   + (size_t)i * V16_PACKED_WORDS_MAX];
                for (int j = 0; j < nf; j++) {
                    u32 r = fbd[j];
                    row[r >> 6] |= (1ULL << (r & 63u));
                }
            }
        }
        CUDA_CHECK(cudaMalloc(&d_variant_masks, vwords * sizeof(u64)));
        CUDA_CHECK(cudaMemcpy(d_variant_masks, h_vmasks,
                              vwords * sizeof(u64), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMalloc(&d_variant_skip, vskip_bytes));
        CUDA_CHECK(cudaMemcpy(d_variant_skip, h_vskip,
                              vskip_bytes, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMalloc(&d_variant_sieve_pass,
                              (size_t)n_variants * sizeof(unsigned long long)));
        CUDA_CHECK(cudaMemset(d_variant_sieve_pass, 0,
                              (size_t)n_variants * sizeof(unsigned long long)));
        printf("#   variant masks: %.2f MiB (%u variants x %d primes x %d words)\n",
               vwords * sizeof(u64) / (1024.0 * 1024.0),
               n_variants, g_cfg.n_sieve_primes, V16_PACKED_WORDS_MAX);
        printf("#   variant skip : %.2f KiB (%u variants x %d primes)\n",
               vskip_bytes / 1024.0,
               n_variants, g_cfg.n_sieve_primes);
        free(h_vmasks);
        free(h_vskip);
    }

    /* kmod LUT: d_kmod_lut[i * n_admissible + w] = wh.offsets[w] mod sieve_primes[i].
     * Layout chosen so consecutive threads (consecutive w) read consecutive bytes
     * for the same prime row — coalesced. u8 fits since every sieve prime <= 211. */
    size_t kmod_bytes = (size_t)g_cfg.n_sieve_primes * (size_t)wh.n_admissible;
    u8 *h_kmod = (u8 *)malloc(kmod_bytes);
    if (!h_kmod) { fprintf(stderr, "OOM kmod_lut (%.2f MiB)\n",
                           kmod_bytes / (1024.0 * 1024.0)); return 1; }
    printf("#   kmod LUT   : %.2f MiB (%d primes x %u offsets)\n",
           kmod_bytes / (1024.0 * 1024.0),
           g_cfg.n_sieve_primes, wh.n_admissible);
    for (int i = 0; i < g_cfg.n_sieve_primes; i++) {
        u32 q = g_cfg.sieve_primes[i];
        u8 *row = h_kmod + (size_t)i * (size_t)wh.n_admissible;
        for (u32 w = 0; w < wh.n_admissible; w++) {
            row[w] = (u8)(wh.offsets[w] % (u64)q);
        }
    }
    u8 *d_kmod_lut;
    CUDA_CHECK(cudaMalloc(&d_kmod_lut, kmod_bytes));
    CUDA_CHECK(cudaMemcpy(d_kmod_lut, h_kmod, kmod_bytes, cudaMemcpyHostToDevice));
    free(h_kmod);

    /* Double-buffered slots: two CUDA streams, two survivor buffers, two
     * pinned host buffers. Used to overlap K1 on the GPU with CPU GMP prove
     * (mirrors v15's dual-stream pipeline). The synchronous path uses
     * slot[0] only. */
    cudaStream_t streams[2];
    u32 *d_tile_base_r[2];
    V16Survivor *d_survivors[2];
    u32 *d_survivor_count[2];
    V16Survivor *h_survivors_pinned[2] = {NULL, NULL};
    u32 h_survivor_count[2] = {0, 0};
    u64 a_tile_for_slot[2] = {0, 0};
    int slot_active[2] = {0, 0};            /* 1 if slot has a pending K1 batch */
    for (int s = 0; s < 2; s++) {
        CUDA_CHECK(cudaStreamCreate(&streams[s]));
        CUDA_CHECK(cudaMalloc(&d_tile_base_r[s], V16_MAX_SIEVE * sizeof(u32)));
        CUDA_CHECK(cudaMalloc(&d_survivors[s],
                              V16_MAX_SURVIVORS * sizeof(V16Survivor)));
        CUDA_CHECK(cudaMalloc(&d_survivor_count[s], sizeof(u32)));
    }

    /* L2 cache persistence for the kmod LUT. RTX 5090 has 96 MiB L2 and
     * our LUT is typically 24-48 MiB; pinning keeps every warp's kmod
     * reads out of HBM (saves dozens of GB/s of L2 fill traffic). Apply
     * per stream so both K1 launches benefit. */
    {
        cudaDeviceProp prop;
        CUDA_CHECK(cudaGetDeviceProperties(&prop, g_cfg.gpu_id));
        size_t persistL2 = prop.persistingL2CacheMaxSize;
        if (persistL2 > 0 && kmod_bytes <= persistL2) {
            CUDA_CHECK(cudaDeviceSetLimit(cudaLimitPersistingL2CacheSize,
                                          kmod_bytes));
            for (int s = 0; s < 2; s++) {
                cudaLaunchAttributeValue attr = {};
                attr.accessPolicyWindow.base_ptr  = (void *)d_kmod_lut;
                attr.accessPolicyWindow.num_bytes = kmod_bytes;
                attr.accessPolicyWindow.hitRatio  = 1.0f;
                attr.accessPolicyWindow.hitProp   = cudaAccessPropertyPersisting;
                attr.accessPolicyWindow.missProp  = cudaAccessPropertyStreaming;
                CUDA_CHECK(cudaStreamSetAttribute(
                    streams[s], cudaLaunchAttributeAccessPolicyWindow, &attr));
            }
            printf("#   L2 pin     : %.2f MiB pinned in L2 (cap %.0f MiB)\n",
                   kmod_bytes / (1024.0 * 1024.0),
                   persistL2 / (1024.0 * 1024.0));
        }
    }
    printf("#   survivors  : 2 x %u slots (2 x %.0f MiB device, dual-stream)\n",
           V16_MAX_SURVIVORS,
           (double)V16_MAX_SURVIVORS * sizeof(V16Survivor) / (1024.0 * 1024.0));

    /* CPU prove state (only used when --cpu-prove). */
    HostChainResult *h_chain_hits = NULL;
    u32 h_chain_hits_count = 0;
    u64 h_pass_top_cpu = 0, h_root_prime_cpu = 0, h_chains_emitted_cpu = 0;
    u64 h_chains_emit_dup_cpu = 0;
    pthread_mutex_t log_mutex = PTHREAD_MUTEX_INITIALIZER;
    /* Chain-length histogram (CPU prove). Indexed by chain_len. */
    u64 chain_hist[V16_HIST_MAX] = {0};
    if (g_cfg.cpu_prove) {
        for (int s = 0; s < 2; s++) {
            CUDA_CHECK(cudaMallocHost((void **)&h_survivors_pinned[s],
                                      (size_t)V16_MAX_SURVIVORS * sizeof(V16Survivor)));
        }
        h_chain_hits = (HostChainResult *)calloc(V16_MAX_RESULTS, sizeof(HostChainResult));
        if (!h_chain_hits) { fprintf(stderr, "OOM chain hits\n"); return 1; }
        printf("#   cpu-prove  : ON  threads=%d (dual-stream pipeline)\n",
               g_cfg.prove_threads);
    }
    /* Pool struct lives across the legacy main
     * loop. Initialized only when --cpu-prove; the GPU-prove path keeps
     * its v16_prove_kernel chain unchanged. */
    ProvePool prove_pool;
    memset(&prove_pool, 0, sizeof(prove_pool));

    /* Output buffers. */
    u32 *d_chain_count;
    u64 *d_chain_a_lo, *d_chain_a_hi;
    u32 *d_chain_len_out, *d_chain_top_pass, *d_chain_root_prime;
    unsigned long long *d_cand_total, *d_pass_sieve, *d_pass_top, *d_chains_emitted;
    CUDA_CHECK(cudaMalloc(&d_chain_count, sizeof(u32)));
    CUDA_CHECK(cudaMalloc(&d_chain_a_lo,  V16_MAX_RESULTS * sizeof(u64)));
    CUDA_CHECK(cudaMalloc(&d_chain_a_hi,  V16_MAX_RESULTS * sizeof(u64)));
    CUDA_CHECK(cudaMalloc(&d_chain_len_out, V16_MAX_RESULTS * sizeof(u32)));
    CUDA_CHECK(cudaMalloc(&d_chain_top_pass, V16_MAX_RESULTS * sizeof(u32)));
    CUDA_CHECK(cudaMalloc(&d_chain_root_prime, sizeof(u32)));
    CUDA_CHECK(cudaMalloc(&d_cand_total, sizeof(unsigned long long)));
    CUDA_CHECK(cudaMalloc(&d_pass_sieve, sizeof(unsigned long long)));
    CUDA_CHECK(cudaMalloc(&d_pass_top,   sizeof(unsigned long long)));
    CUDA_CHECK(cudaMalloc(&d_chains_emitted, sizeof(unsigned long long)));

    /* Determine A range from --bits (if given).
     * Bit-count of a positive integer p is floor(log2 p) + 1, so a "bits=N"
     * search means p ∈ [2^(N-1), 2^N). For target_bits ≥ 128 we use GMP
     * since 2^(N-1) overflows u128; --cpu-prove is required because the
     * GPU prove kernel uses u128 internally. */
    u128 A_min = 1, A_max = 0;
    if (g_cfg.target_bits > 0) {
        if (g_cfg.target_bits >= 128 && !g_cfg.cpu_prove) {
            fprintf(stderr,
                "ERROR: --bits >= 128 requires --cpu-prove (GPU MR uses u128)\n");
            return 1;
        }
        if (g_cfg.target_bits >= 200) {
            fprintf(stderr, "ERROR: --bits must be < 200 (got %d)\n",
                    g_cfg.target_bits);
            return 1;
        }
        mpz_t z_lower, z_upper, z_D2E, z_amin, z_amax;
        mpz_inits(z_lower, z_upper, z_D2E, z_amin, z_amax, NULL);
        mpz_set_ui(z_lower, 1);
        mpz_mul_2exp(z_lower, z_lower, (unsigned long)(g_cfg.target_bits - 1));
        mpz_set_ui(z_upper, 1);
        mpz_mul_2exp(z_upper, z_upper, (unsigned long)g_cfg.target_bits);
        mpz_set_str(z_D2E, g_cfg.D_dec, 10);
        if (g_cfg.exp_start > 0)
            mpz_mul_2exp(z_D2E, z_D2E, (unsigned long)g_cfg.exp_start);
        mpz_cdiv_q(z_amin, z_lower, z_D2E);
        mpz_fdiv_q(z_amax, z_upper, z_D2E);
        if (mpz_sizeinbase(z_amax, 2) > 128) {
            fprintf(stderr, "ERROR: A_max exceeds u128 — D is too small or"
                            " --bits too large\n");
            return 1;
        }
        /* Convert via hex roundtrip. */
        char *amin_hex = mpz_get_str(NULL, 16, z_amin);
        char *amax_hex = mpz_get_str(NULL, 16, z_amax);
        A_min = 0;
        for (char *p = amin_hex; *p; p++) {
            int d = (*p >= '0' && *p <= '9') ? *p - '0'
                  : (*p >= 'a' && *p <= 'f') ? *p - 'a' + 10
                  : (*p >= 'A' && *p <= 'F') ? *p - 'A' + 10 : -1;
            if (d < 0) break;
            A_min = (A_min << 4) | (u128)d;
        }
        A_max = 0;
        for (char *p = amax_hex; *p; p++) {
            int d = (*p >= '0' && *p <= '9') ? *p - '0'
                  : (*p >= 'a' && *p <= 'f') ? *p - 'a' + 10
                  : (*p >= 'A' && *p <= 'F') ? *p - 'A' + 10 : -1;
            if (d < 0) break;
            A_max = (A_max << 4) | (u128)d;
        }
        free(amin_hex);
        free(amax_hex);
        mpz_clears(z_lower, z_upper, z_D2E, z_amin, z_amax, NULL);
        if (A_max <= A_min) {
            fprintf(stderr, "ERROR: --bits %d leaves empty A range for this seed\n",
                    g_cfg.target_bits);
            return 1;
        }
        char buf1[80], buf2[80];
        u128_to_dec(A_min, buf1, sizeof(buf1));
        u128_to_dec(A_max, buf2, sizeof(buf2));
        printf("#   A range    = [%s, %s)\n", buf1, buf2);
    }

    /* Tile cursor: A_tile = A / coeff_primorial. */
    u64 a_tile_start = 0, a_tile_end = 0;
    if (g_cfg.target_bits > 0) {
        a_tile_start = (u64)(A_min / wh.primorial);
        u128 a_tile_end_128 = (A_max + wh.primorial - 1) / wh.primorial;
        if (a_tile_end_128 >> 64) {
            fprintf(stderr, "ERROR: A_tile range exceeds u64 — narrow --bits\n");
            return 1;
        }
        a_tile_end = (u64)a_tile_end_128;
    } else {
        a_tile_start = 0;
        a_tile_end = ~(u64)0;
    }
    /* --start-tile / --end-tile manual overrides (useful for known-A tests
     * and for resume-style sharding). Both clamp into the --bits range when
     * --bits was given. */
    if (g_cfg.start_tile > a_tile_start) a_tile_start = g_cfg.start_tile;
    if (g_cfg.end_tile && g_cfg.end_tile < a_tile_end) a_tile_end = g_cfg.end_tile;

    /* Resume. */
    u64 a_tile_cursor = a_tile_start;
    if (g_cfg.resume_file[0]) {
        u64 rcur = 0;
        if (load_checkpoint_cursor(&g_cfg, wh.fnv1a64, &rcur) != 0) return 1;
        if (rcur > a_tile_cursor) a_tile_cursor = rcur;
        printf("#   resume     = %llu\n", (unsigned long long)a_tile_cursor);
    }

    /* Run. */
    struct timeval t_start, t_now, t_last_report;
    gettimeofday(&t_start, NULL);
    t_last_report = t_start;

    u32 zero32 = 0; unsigned long long zero64 = 0;
    CUDA_CHECK(cudaMemcpy(d_chain_count, &zero32, sizeof(u32), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_chain_root_prime, &zero32, sizeof(u32), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_cand_total, &zero64, sizeof(unsigned long long), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_pass_sieve, &zero64, sizeof(unsigned long long), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_pass_top,   &zero64, sizeof(unsigned long long), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_chains_emitted, &zero64, sizeof(unsigned long long), cudaMemcpyHostToDevice));

    FILE *log_fp = NULL;
    if (g_cfg.log_file[0]) {
        log_fp = fopen(g_cfg.log_file, "a");
        if (log_fp) fprintf(log_fp, "# v16 run start seed=%s D=%s target=%d E=%d\n",
                            g_cfg.seed_name, g_cfg.D_dec, g_cfg.target_len, g_cfg.exp_start);
    }

    /* Init pool for legacy --cpu-prove path.
     * V16_MAX_SURVIVORS = 32M × 8 B = 256 MiB / buffer — clamp the
     * legacy queue depth to max 4 to avoid pinning >1 GiB. */
    if (g_cfg.cpu_prove) {
        int saved_depth = g_cfg.prove_queue_depth;
        v16_cfg *cfg_mut = (v16_cfg *)&g_cfg;
        if (saved_depth > 4) cfg_mut->prove_queue_depth = 4;
        if (prove_pool_init(&prove_pool, &g_cfg,
                            &h_pass_top_cpu, &h_root_prime_cpu,
                            &h_chains_emitted_cpu,
                            &h_chains_emit_dup_cpu,
                            h_chain_hits, V16_MAX_RESULTS, &h_chain_hits_count,
                            &log_mutex, log_fp, g_cfg.seed_name,
                            chain_hist, V16_HIST_MAX,
                            NULL, NULL,
                            sizeof(V16Survivor), V16_MAX_SURVIVORS) != 0) {
            fprintf(stderr, "FATAL: prove_pool_init (legacy) failed\n");
            return 1;
        }
        cfg_mut->prove_queue_depth = saved_depth;
    }

    u64 tiles_done = 0;
    u32 launch_tiles = (u32)g_cfg.tiles_per_launch;
    if (launch_tiles == 0) launch_tiles = 1024;
    /* --sieve-variants emits ~n_variants survivors per (tile, woff) baseline,
     * so we shrink the default tiles-per-launch to keep total survivors per
     * batch under V16_MAX_SURVIVORS. The user's explicit --tiles always wins. */
    if (g_cfg.sieve_variants && !g_cfg.tiles_explicit && n_variants > 0) {
        u32 divider = n_variants / 16u;
        if (divider < 1u) divider = 1u;
        u32 adj = launch_tiles / divider;
        if (adj < 64u) adj = 64u;
        if (adj < launch_tiles) {
            launch_tiles = adj;
            printf("#   tiles auto : --sieve-variants n=%u -> tiles/launch=%u\n",
                   n_variants, launch_tiles);
        }
    }
    if (launch_tiles > 65535) {
        /* CUDA gridDim.y caps at 65535. Larger requested batches get
         * clamped here; the main loop just runs more launches. */
        launch_tiles = 65535;
    }
    /* Adaptive block/tile-stride: when n_admissible is small (heavy-D / wheel
     * collapse), shrinking the block to a single warp keeps thread-occupancy
     * high, while a tile-stride loop inside the kernel pumps more work per
     * block so the GPU stays saturated. For large n_admissible the dispatch
     * degenerates to the original 256-thread / 1-tile-per-block layout. */
    u32 block_size = V16_THREADS;
    {
        u32 padded = (wh.n_admissible + 31u) & ~31u;   /* round up to warp */
        if (padded == 0) padded = 32u;
        if (padded < block_size) block_size = padded;
    }
    u32 tile_stride = V16_THREADS / block_size;
    if (tile_stride == 0) tile_stride = 1;
    u32 grid_x = (wh.n_admissible + block_size - 1) / block_size;
    u32 grid_y_cap = (launch_tiles + tile_stride - 1) / tile_stride;

    printf("#   launch     : grid=(%u, <=%u) block=%u  tile_stride=%u  tiles/launch=%u\n",
           grid_x, grid_y_cap, block_size, tile_stride, launch_tiles);

    while (!g_stop) {
        if (g_cfg.max_tiles && tiles_done >= g_cfg.max_tiles) break;
        if (a_tile_cursor >= a_tile_end) break;
        u64 remaining = a_tile_end - a_tile_cursor;
        u32 this_launch = (u32)((remaining < launch_tiles) ? remaining : launch_tiles);
        if (g_cfg.max_tiles) {
            u64 rem_total = g_cfg.max_tiles - tiles_done;
            if (rem_total < this_launch) this_launch = (u32)rem_total;
        }

        /* Pipelined launch on slot[cur]. Per-iter:
         *   1. Launch K1 on slot[cur] (async on its stream).
         *   2. If slot[prev] has unconsumed survivors, run CPU prove on them
         *      in parallel with the GPU K1 (CPU and GPU different hardware).
         *   3. After prove, sync slot[cur]'s stream and copy its survivors
         *      to slot[cur]'s pinned host buffer.
         *   4. Swap cur ↔ prev.
         * For the GPU-prove path we keep the original single-stream
         * behavior by reusing slot[0] only. */
        u32 h_tile_base_r[V16_MAX_SIEVE] = {0};
        for (int i = 0; i < g_cfg.n_sieve_primes; i++) {
            u32 q = g_cfg.sieve_primes[i];
            u32 cur_mod_q = (u32)(a_tile_cursor % (u64)q);
            h_tile_base_r[i] = (u32)(((u64)cur_mod_q * (u64)h_cp_mod_q[i]) % (u64)q);
        }

        /* Slot selection: alternate for CPU-prove, fixed 0 for GPU-prove. */
        int cur_slot  = g_cfg.cpu_prove ? (int)((tiles_done / launch_tiles) & 1) : 0;
        int prev_slot = 1 - cur_slot;

        CUDA_CHECK(cudaMemcpyAsync(d_tile_base_r[cur_slot], h_tile_base_r,
                                   (size_t)g_cfg.n_sieve_primes * sizeof(u32),
                                   cudaMemcpyHostToDevice, streams[cur_slot]));
        CUDA_CHECK(cudaMemsetAsync(d_survivor_count[cur_slot], 0, sizeof(u32),
                                   streams[cur_slot]));

        u32 this_grid_y = (this_launch + tile_stride - 1) / tile_stride;
        dim3 grid(grid_x, this_grid_y, 1);
        dim3 block(block_size, 1, 1);

        v16_sieve_kernel<<<grid, block, 0, streams[cur_slot]>>>(
            wh.n_admissible,
            d_tile_base_r[cur_slot],
            this_launch,
            tile_stride,
            d_kmod_lut,
            d_masks,
            d_survivors[cur_slot], d_survivor_count[cur_slot], V16_MAX_SURVIVORS,
            d_cand_total, d_pass_sieve,
            n_variants, d_variant_masks, d_variant_skip,
            d_variant_sieve_pass);
        CUDA_CHECK(cudaGetLastError());
        a_tile_for_slot[cur_slot] = a_tile_cursor;
        slot_active[cur_slot] = 1;

        /* Async copy survivor count back. We need it before we can size the
         * survivor-list copy; stream-sync below blocks on this implicitly. */
        CUDA_CHECK(cudaMemcpyAsync(&h_survivor_count[cur_slot],
                                   d_survivor_count[cur_slot], sizeof(u32),
                                   cudaMemcpyDeviceToHost, streams[cur_slot]));

        if (!g_cfg.cpu_prove) {
            /* Synchronous GPU prove path. */
            CUDA_CHECK(cudaStreamSynchronize(streams[cur_slot]));
            u32 n_surv = h_survivor_count[cur_slot];
            if (n_surv > V16_MAX_SURVIVORS) {
                fprintf(stderr, "WARN: survivor overflow (%u > %u)\n",
                        n_surv, V16_MAX_SURVIVORS);
                n_surv = V16_MAX_SURVIVORS;
            }
            if (n_surv > 0) {
                u32 prove_grid = (n_surv + V16_THREADS - 1) / V16_THREADS;
                v16_prove_kernel<<<prove_grid, V16_THREADS, 0, streams[cur_slot]>>>(
                    d_offsets,
                    d_survivors[cur_slot], d_survivor_count[cur_slot],
                    a_tile_for_slot[cur_slot],
                    d_chain_count,
                    d_chain_a_lo, d_chain_a_hi, d_chain_len_out, d_chain_top_pass,
                    d_chain_root_prime,
                    d_pass_top, d_chains_emitted);
                CUDA_CHECK(cudaGetLastError());
            }
            CUDA_CHECK(cudaStreamSynchronize(streams[cur_slot]));
            slot_active[cur_slot] = 0;
        } else {
            /* Sync stream to get the survivor
             * count, D2H into a pool-acquired buffer, then submit a
             * non-blocking ProveJob. Pool workers run concurrent with
             * the next K1 launch. */
            CUDA_CHECK(cudaStreamSynchronize(streams[cur_slot]));
            u32 n_surv = h_survivor_count[cur_slot];
            if (n_surv > V16_MAX_SURVIVORS) {
                fprintf(stderr, "WARN: survivor overflow (%u > %u)\n",
                        n_surv, V16_MAX_SURVIVORS);
                n_surv = V16_MAX_SURVIVORS;
                h_survivor_count[cur_slot] = n_surv;
            }
            int buf_idx = prove_pool_acquire_buf(&prove_pool);
            V16Survivor *buf = (V16Survivor *)
                prove_pool_buf_ptr(&prove_pool, buf_idx);
            if (n_surv > 0) {
                CUDA_CHECK(cudaMemcpyAsync(
                    buf, d_survivors[cur_slot],
                    (size_t)n_surv * sizeof(V16Survivor),
                    cudaMemcpyDeviceToHost, streams[cur_slot]));
                CUDA_CHECK(cudaStreamSynchronize(streams[cur_slot]));
            }
            ProveJob job;
            memset(&job, 0, sizeof(job));
            job.buf_idx         = buf_idx;
            job.n_surv          = n_surv;
            job.q_iter_mode     = 0;
            job.k15_passed      = 0;
            job.v_idx           = -1;
            job.a_tile_base     = a_tile_for_slot[cur_slot];
            job.coeff_primorial = wh.primorial;
            job.wheel_offsets   = wh.offsets;
            prove_pool_submit(&prove_pool, job);
            slot_active[cur_slot] = 0;
            (void)prev_slot;
        }

        a_tile_cursor += this_launch;
        tiles_done    += this_launch;

        gettimeofday(&t_now, NULL);
        double now_dt = (t_now.tv_sec - t_last_report.tv_sec) +
                        (t_now.tv_usec - t_last_report.tv_usec) / 1e6;
        if (now_dt >= g_cfg.report_sec) {
            unsigned long long cand=0, ps=0, pt=0, ce=0; u32 ccount=0, croot=0;
            CUDA_CHECK(cudaMemcpy(&cand, d_cand_total, sizeof(unsigned long long), cudaMemcpyDeviceToHost));
            CUDA_CHECK(cudaMemcpy(&ps,   d_pass_sieve, sizeof(unsigned long long), cudaMemcpyDeviceToHost));
            if (g_cfg.cpu_prove) {
                pt = h_pass_top_cpu;
                ce = h_chains_emitted_cpu;
                croot = (u32)h_root_prime_cpu;
                ccount = h_chain_hits_count;
            } else {
                CUDA_CHECK(cudaMemcpy(&pt,   d_pass_top,   sizeof(unsigned long long), cudaMemcpyDeviceToHost));
                CUDA_CHECK(cudaMemcpy(&ce,   d_chains_emitted, sizeof(unsigned long long), cudaMemcpyDeviceToHost));
                CUDA_CHECK(cudaMemcpy(&ccount, d_chain_count, sizeof(u32), cudaMemcpyDeviceToHost));
                CUDA_CHECK(cudaMemcpy(&croot, d_chain_root_prime, sizeof(u32), cudaMemcpyDeviceToHost));
            }
            double elapsed = (t_now.tv_sec - t_start.tv_sec) +
                             (t_now.tv_usec - t_start.tv_usec) / 1e6;
            printf("[t=%.2fs] tiles=%llu cand=%llu (%.2fM/s)  "
                   "sieve-pass=%llu  top-pass=%llu  root-prime=%u  chains>=%d=%u  emitted=%llu\n",
                   elapsed, (unsigned long long)tiles_done, cand,
                   cand / 1e6 / (elapsed > 0 ? elapsed : 1),
                   ps, pt, croot, g_cfg.min_report_len, ccount, ce);
            t_last_report = t_now;
            if (g_cfg.checkpoint_file[0]) {
                /* Drain so a_tile_cursor matches
                 * proved coverage (dispatch has already advanced it past
                 * in-flight jobs). */
                if (g_cfg.cpu_prove) prove_pool_drain_all(&prove_pool);
                save_checkpoint(&g_cfg, a_tile_cursor, wh.fnv1a64);
            }
        }
        if (g_cfg.time_limit_sec) {
            double el = (t_now.tv_sec - t_start.tv_sec) +
                        (t_now.tv_usec - t_start.tv_usec) / 1e6;
            if (el >= g_cfg.time_limit_sec) break;
        }
    }

    /* Drain pending + in-flight jobs at run end. */
    if (g_cfg.cpu_prove) {
        prove_pool_drain_all(&prove_pool);
    }

    /* Final summary. */
    {
        unsigned long long cand=0, ps=0, pt=0, ce=0;
        u32 ccount=0, croot=0;
        CUDA_CHECK(cudaMemcpy(&cand, d_cand_total, sizeof(unsigned long long), cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(&ps,   d_pass_sieve, sizeof(unsigned long long), cudaMemcpyDeviceToHost));
        if (g_cfg.cpu_prove) {
            pt    = h_pass_top_cpu;
            ce    = h_chains_emitted_cpu;
            croot = (u32)h_root_prime_cpu;
            ccount = h_chain_hits_count;
        } else {
            CUDA_CHECK(cudaMemcpy(&pt,   d_pass_top,   sizeof(unsigned long long), cudaMemcpyDeviceToHost));
            CUDA_CHECK(cudaMemcpy(&ce,   d_chains_emitted, sizeof(unsigned long long), cudaMemcpyDeviceToHost));
            CUDA_CHECK(cudaMemcpy(&ccount, d_chain_count, sizeof(u32), cudaMemcpyDeviceToHost));
            CUDA_CHECK(cudaMemcpy(&croot, d_chain_root_prime, sizeof(u32), cudaMemcpyDeviceToHost));
        }
        gettimeofday(&t_now, NULL);
        double elapsed = (t_now.tv_sec - t_start.tv_sec) +
                         (t_now.tv_usec - t_start.tv_usec) / 1e6;
        printf("\n=== SUMMARY ===\n");
        printf("  prove path    : %s\n", g_cfg.cpu_prove ? "CPU/GMP" : "GPU/u128");
        printf("  elapsed       : %.2fs\n", elapsed);
        printf("  tiles_done    : %llu\n", (unsigned long long)tiles_done);
        printf("  candidates    : %llu\n", cand);
        printf("  coeff/sec     : %.2fM\n", cand / 1e6 / (elapsed > 0 ? elapsed : 1));
        printf("  sieve_pass    : %llu\n", ps);
        printf("  top_pass      : %llu\n", pt);
        printf("  root_prime    : %u\n", croot);
        printf("  chains>=%d    : %u\n", g_cfg.min_report_len, ccount);
        printf("  chains_emitted: %llu\n", ce);
        if (g_cfg.sieve_variants && d_variant_sieve_pass && n_variants > 0) {
            unsigned long long *h_vps = (unsigned long long *)calloc(
                n_variants, sizeof(unsigned long long));
            if (h_vps) {
                CUDA_CHECK(cudaMemcpy(h_vps, d_variant_sieve_pass,
                    (size_t)n_variants * sizeof(unsigned long long),
                    cudaMemcpyDeviceToHost));
                printf("  pool-seed sieve survivors:\n");
                for (u32 v = 0; v < n_variants; v++) {
                    if (h_vps[v] == 0) continue;
                    printf("    %-32s  %10llu\n",
                           g_cfg.pool_seeds[v].name, h_vps[v]);
                }
                free(h_vps);
            }
        }
        if (g_cfg.n_pool_seeds > 0) {
            printf("  pool-seed hits (D | (p+1)):\n");
            for (int i = 0; i < g_cfg.n_pool_seeds; i++) {
                v16_pool_seed *ps = &g_cfg.pool_seeds[i];
                if (ps->hits == 0) continue;
                printf("    %-32s  %10llu\n", ps->name,
                       (unsigned long long)ps->hits);
            }
        }
        if (g_cfg.cpu_prove) {
            printf("  chain-len histogram (CC1+):\n");
            u64 cum_high = 0;
            int max_seen = 0;
            for (int i = 1; i < V16_HIST_MAX; i++)
                if (chain_hist[i]) max_seen = i;
            for (int i = 1; i <= max_seen; i++) {
                if (chain_hist[i] == 0) continue;
                cum_high = 0;
                for (int j = i; j < V16_HIST_MAX; j++) cum_high += chain_hist[j];
                printf("    CC%-3d  count=%-10llu  >=CC%d=%llu\n",
                       i, (unsigned long long)chain_hist[i],
                       i, (unsigned long long)cum_high);
            }
        }
    }

    /* Drain chain results. */
    if (g_cfg.cpu_prove) {
        u32 drain = h_chain_hits_count > V16_MAX_RESULTS
                  ? V16_MAX_RESULTS : h_chain_hits_count;
        if (drain > 0) {
            printf("\n=== CHAIN HITS (drain %u of %u, CPU prove) ===\n",
                   drain, h_chain_hits_count);
            for (u32 i = 0; i < drain; i++) {
                u128 A = ((u128)h_chain_hits[i].A_hi << 64) | h_chain_hits[i].A_lo;
                char buf[80]; u128_to_dec(A, buf, sizeof(buf));
                printf("  chain_len=%u  top_pass=%u  A=%s\n",
                       h_chain_hits[i].chain_len, h_chain_hits[i].top_pass, buf);
            }
        }
    } else {
        u32 ccount = 0;
        CUDA_CHECK(cudaMemcpy(&ccount, d_chain_count, sizeof(u32), cudaMemcpyDeviceToHost));
        u32 drain = ccount > V16_MAX_RESULTS ? V16_MAX_RESULTS : ccount;
        if (drain > 0) {
            u64 *h_a_lo = (u64 *)malloc(drain * sizeof(u64));
            u64 *h_a_hi = (u64 *)malloc(drain * sizeof(u64));
            u32 *h_len  = (u32 *)malloc(drain * sizeof(u32));
            u32 *h_top  = (u32 *)malloc(drain * sizeof(u32));
            CUDA_CHECK(cudaMemcpy(h_a_lo, d_chain_a_lo, drain * sizeof(u64), cudaMemcpyDeviceToHost));
            CUDA_CHECK(cudaMemcpy(h_a_hi, d_chain_a_hi, drain * sizeof(u64), cudaMemcpyDeviceToHost));
            CUDA_CHECK(cudaMemcpy(h_len,  d_chain_len_out, drain * sizeof(u32), cudaMemcpyDeviceToHost));
            CUDA_CHECK(cudaMemcpy(h_top,  d_chain_top_pass, drain * sizeof(u32), cudaMemcpyDeviceToHost));
            printf("\n=== CHAIN HITS (drain %u of %u, GPU prove) ===\n", drain, ccount);
            for (u32 i = 0; i < drain; i++) {
                u128 A = ((u128)h_a_hi[i] << 64) | h_a_lo[i];
                char buf[80]; u128_to_dec(A, buf, sizeof(buf));
                printf("  chain_len=%u  top_pass=%u  A=%s\n",
                       h_len[i], h_top[i], buf);
                if (log_fp) fprintf(log_fp, "HIT len=%u top_pass=%u A=%s seed=%s\n",
                                    h_len[i], h_top[i], buf, g_cfg.seed_name);
            }
            free(h_a_lo); free(h_a_hi); free(h_len); free(h_top);
        }
    }
    if (log_fp) fclose(log_fp);

    if (g_cfg.checkpoint_file[0])
        save_checkpoint(&g_cfg, a_tile_cursor, wh.fnv1a64);

    cudaFree(d_offsets); cudaFree(d_masks); cudaFree(d_kmod_lut);
    for (int s = 0; s < 2; s++) {
        cudaFree(d_tile_base_r[s]);
        cudaFree(d_survivors[s]);
        cudaFree(d_survivor_count[s]);
        if (h_survivors_pinned[s]) cudaFreeHost(h_survivors_pinned[s]);
        cudaStreamDestroy(streams[s]);
    }
    /* Shut down pool (drains + joins workers). */
    if (g_cfg.cpu_prove) prove_pool_destroy(&prove_pool);
    free(h_chain_hits);
    pthread_mutex_destroy(&log_mutex);
    cudaFree(d_chain_count); cudaFree(d_chain_a_lo); cudaFree(d_chain_a_hi);
    cudaFree(d_chain_len_out); cudaFree(d_chain_top_pass); cudaFree(d_chain_root_prime);
    cudaFree(d_cand_total); cudaFree(d_pass_sieve); cudaFree(d_pass_top); cudaFree(d_chains_emitted);
    if (d_variant_masks) cudaFree(d_variant_masks);
    if (d_variant_skip) cudaFree(d_variant_skip);
    if (d_variant_sieve_pass) cudaFree(d_variant_sieve_pass);
    free(h_masks); wheel_free(&wh);
    return 0;
}
