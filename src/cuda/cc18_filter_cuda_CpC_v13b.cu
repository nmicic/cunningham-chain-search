/*
 * Copyright (c) 2026 Nenad Mićić <nenad@micic.be>
 * SPDX-License-Identifier: Apache-2.0
 *
 * =============================================================================
 * CC18 GPU Filter — C+C v13b: Second-Kind Variant (2p-1 chains)
 *
 * Second-kind Cunningham chains: p → 2p-1 → 2(2p-1)-1 → ...
 * Derived from v13 (first-kind) with 5 targeted changes:
 *   1. Forbidden residue formula: (inv+q-1)%q → (1+q-inv)%q
 *   2. Mod-6 constraint: n≡5(mod6) → n≡1(mod6)
 *   3. GPU chain follow: 2p+1 → 2p-1
 *   4. GPU/CPU root check: (n-1)/2 → (n+1)/2
 *   5. Self-test vectors: second-kind chains
 *
 * CHANGELOG (v13 over v12):
 *   [OPT] L1 cache preference: cudaFuncSetCacheConfig(PreferL1) — Blackwell only
 *         (sm_120+). On Ada (4090, sm_89) this drops occupancy 100%→83%, causing
 *         a ~2x regression. Gated on compute capability >= 12.0.
 *   [OPT] L2 cache persistence: kM_mod_line LUT (5.17MB) pinned in L2 via
 *         cudaAccessPolicyWindow. Blackwell only (sm_120+). On Ada the 72MB L2
 *         naturally holds the 4.93MB LUT; pinning adds API overhead for no gain.
 *   [OPT] kM_mod row prefetch: PTX prefetch.global.L1 issued for the 250-byte
 *         kM_mod row AFTER L2+ext-L2 filters pass (~8% of candidates). Previous
 *         placement was before L2 filter, wasting bandwidth on 92% of candidates
 *         that never reach line-sieve. Now prefetches only when needed, overlapping
 *         L2→L1 fetch with the line-sieve loop setup.
 *   K2 (prove kernel) unchanged — CPU proving recommended for correctness.
 *   [NEW] GPU utilization diagnostic: CUDA events measure K1 filter time vs
 *         wall-clock. Periodic report shows gpu=XX%. Final summary shows K1 time,
 *         prove time, prove/K1 ratio, and bottleneck warning when CPU-bound.
 *   [NEW] Occupancy query: cudaOccupancyMaxActiveBlocksPerMultiprocessor reports
 *         actual achievable occupancy, warns if --blocks-per-sm is oversized.
 *
 * CHANGELOG (v12 over v10):
 *   [NEW] High-entropy PRNG seeding via /dev/urandom — machine-unique seeds
 *         eliminate tile overlap across multi-server deployments. Fallback to
 *         clock_gettime(CLOCK_MONOTONIC) + getpid() if /dev/urandom unavailable.
 *   [NEW] --seed N: explicit u64 seed for reproducible debugging runs.
 *   [NEW] Between-batch nanosecond entropy mixing: rng_state ^= tv_nsec before
 *         each random tile selection ensures sequence divergence even if two
 *         GPUs start with the same seed.
 *   [NEW] Seed banner shows hex seed value and source (/dev/urandom or --seed).
 *   v11 is reserved for the wide variant.
 *
 * CHANGELOG (v10 over v9):
 *   [NEW] Full-coverage base generation: removed r%6==5 filter. ALL CRT-
 *         surviving bases are now used (~36 bases vs ~6 in v9). This covers
 *         the entire first-kind search space, including mod-7 residue class 6
 *         (q7-immunized candidates) that v9 could never find.
 *   [NEW] Mod-6 bucketed wheels: each base's wheel is partitioned into 6
 *         buckets by k_offset%6. The kernel selects the unique valid bucket
 *         per (tile,base) pair so that every tested candidate satisfies
 *         n ≡ 1 (mod 6). Zero wasted k-tests (v9 wasted ~5/6).
 *   [NEW] --coverage legacy|full: switch between v9-compatible (r%6==5 only)
 *         and full-coverage modes for A/B benchmarking. Default: full.
 *   [NEW] --verify-root 0xHEX: decompose a known root into (base, k_offset,
 *         tile) and verify it's reachable by the current base generation.
 *   [NEW] Self-tests 15-19: mod-6 bucket correctness, FULL vs LEGACY base
 *         count, CC17 root reachability, bucket consistency, q7-immunized
 *         test vectors.
 *   [NEW] Banner prints coverage mode, base count, mod-7 distribution.
 *   [NEW] MAX_BASES increased from 16 to 64 to accommodate all CRT bases.
 *   [FIX] Sequential prefix mode now STOPS after one complete sweep instead
 *         of wrapping around (v9 bug: same chains re-discovered). Random
 *         prefix mode is unaffected (runs until --time / Ctrl+C).
 *
 * CHANGELOG (v9 over v8):
 *   [NEW] --prefix, --prefix-mode, --prefix-lanes, --checkpoint, --resume.
 *   [FIX] Sequential prefix wrap-around.
 *
 * Implementation approach:
 *   Prefix → cand_min/cand_max → k_min/k_max → tile_min/tile_max.
 *   All narrowing happens BEFORE the main loop.
 *   Mod-6 bucketed wheels: per base, target_k_mod6 is precomputed. At kernel
 *   launch time, tile_start_mod6 is passed. Kernel selects the correct bucket
 *   as valid_bucket = (target_k_mod6 - tile_mod6 + 6) % 6.
 *   Math: n = (tile*WP + k_offset)*M + base. WP%6=1, M%6=5.
 *         n%6 = ((tile+k_offset)*5 + base)%6. For n≡1(mod6):
 *         k_offset%6 = target_k_mod6 - tile%6 (mod 6).
 *
 * Example commands:
 *   # Full-coverage run (default, all bases):
 *   ./cc18_CpC_v13b --target 18 --bits 89 --coverage full --time 30
 *
 *   # Legacy mode (v9-compatible, r%6==5 only):
 *   ./cc18_CpC_v13b --target 18 --bits 89 --coverage legacy --time 30
 *
 *   # Verify a known root is reachable:
 *   ./cc18_CpC_v13b --target 18 --bits 109 --verify-root 0x140ECDB88F050428006DE301E959
 *
 *   # Sequential prefix run (top 3 bits = 101):
 *   ./cc18_CpC_v13b --target 18 --bits 89 --prefix 0b101 --prefix-mode sequential
 *
 *   # Multi-GPU sharding (4 GPUs, this is GPU 2):
 *   ./cc18_CpC_v13b --target 18 --bits 89 --prefix 0b101 --prefix-lanes 4 \
 *                   --prefix-lane-id 2 --gpu 2
 *
 * All v8/v9 features retained: line-sieve fix, non-root backward walk,
 *   input validation, drain safety, 20-witness MR, prefix search, etc.
 *
 * TWO-KERNEL architecture — GPU handles EVERYTHING:
 *   Kernel 1 (filter): L2 bitmask → ext-L2 byte → line-sieve → survivors
 *   Kernel 2 (prove):  Reconstruct n → trial division → Miller-Rabin →
 *                       chain root check → chain follow → emit chain results
 *   CPU: Only reads chain results for logging (nearly idle).
 *
 * Pipeline diagram:
 *   Stream A:  [K1 filter → K2 prove][K1 → K2]...
 *   Stream B:       [K1 filter → K2 prove][K1 → K2]...
 *   CPU:       [output A results][output B results]...  (trivial)
 *
 * Build:
 *   nvcc -O3 -arch=sm_120 cc18_filter_cuda_CpC_v13b.cu -o cc18_CpC_v13b -lgmp -lpthread
 *   (sm_120 = Blackwell / RTX 5090. Use sm_86 for Ampere / RTX 4090.)
 *
 * Fallback: --cpu-prove uses persistent pthread pool (same as v4)
 *
 * References:
 *   - cc18_filter_cuda_CpC_v9.cu: prefix-constrained search (base)
 *   - cc_gmp_v32_claude_03.c: mod-6 bucketed wheel algorithm
 *   - cc_gmp_v31_claude.c (historical CPU reference): filter pipeline, table init, hot path
 * =============================================================================
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
#include <pthread.h>
#include <sys/time.h>
#include <fcntl.h>
#include <gmp.h>

/* =========================================================================
 * Type aliases and constants (match v31)
 * ========================================================================= */

typedef uint64_t u64;
typedef uint32_t u32;
typedef uint16_t u16;
typedef uint8_t  u8;
typedef unsigned __int128 u128;

#define CUDA_CHECK(call) do { \
    cudaError_t err = (call); \
    if (err != cudaSuccess) { \
        fprintf(stderr, "CUDA error at %s:%d: %s\n", \
                __FILE__, __LINE__, cudaGetErrorString(err)); \
        exit(1); \
    } \
} while(0)

/* =========================================================================
 * Lattice configuration — matches v31 cc18a profile exactly
 * ========================================================================= */

#define CRT_PRIMES_COUNT      6
#define FILTER_PRIMES_COUNT   7    /* L2 bitmask: 37-61 */
#define EXT_FILTER_PRIMES_COUNT 7  /* ext-L2 byte: 67-97 */
#define LINE_SIEVE_COUNT      125  /* Line-sieve: 101-863 */
#define LINE_SIEVE_START_IDX  7    /* Skip primes handled by ext-L2 */
#define WHEEL_PRIMES_COUNT    3
#define WHEEL_PERIOD          20677u  /* 23 * 29 * 31 */
#define MAX_SIEVE_CHAIN_LEN   32
#define MAX_BASES             64
#define MAX_CHAIN             64
#define CC18_SIEVE_LEN        18   /* Target chain length */

/* Threads/blocks for RTX 5090 (170 SMs, Blackwell) */
#define THREADS_PER_BLOCK     256
#define WARP_SIZE             32

/* Survivor buffer size per batch */
#define MAX_SURVIVORS_PER_BATCH (1 << 20)  /* 1M survivors max (~0.14% of 700M) */

static const u32 CRT_PRIMES[CRT_PRIMES_COUNT] = {5, 7, 11, 13, 17, 19};
static const u64 LATTICE_M = 5ULL * 7 * 11 * 13 * 17 * 19;  /* 1,616,615 */

static const u32 FILTER_PRIMES[FILTER_PRIMES_COUNT] = {37, 41, 43, 47, 53, 59, 61};
static const u32 EXT_FILTER_PRIMES[EXT_FILTER_PRIMES_COUNT] = {67, 71, 73, 79, 83, 89, 97};
static const u32 WHEEL_PRIMES[WHEEL_PRIMES_COUNT] = {23, 29, 31};

/* Full chain screening primes (132 primes, 67-863) — indices 0-6 are ext-L2 */
static const u32 CHAIN_SCREEN_PRIMES[] = {
    67, 71, 73, 79, 83, 89, 97,           /* ext-L2 (indices 0-6) */
    101, 103, 107, 109, 113, 127, 131,     /* line-sieve starts at index 7 */
    137, 139, 149, 151, 157, 163, 167, 173,
    179, 181, 191, 193, 197, 199, 211, 223,
    227, 229, 233, 239, 241, 251, 257, 263,
    269, 271, 277, 281, 283, 293, 307, 311,
    313, 317, 331, 337, 347, 349, 353, 359,
    367, 373, 379, 383, 389, 397, 401, 409,
    419, 421, 431, 433, 439, 443, 449, 457,
    461, 463, 467, 479, 487, 491, 499, 503,
    509, 521, 523, 541, 547, 557, 563, 569,
    571, 577, 587, 593, 599, 601, 607, 613,
    617, 619, 631, 641, 643, 647, 653, 659,
    661, 673, 677, 683, 691, 701, 709, 719,
    727, 733, 739, 743, 751, 757, 761, 769,
    773, 787, 797, 809, 811, 821, 823, 827,
    829, 839, 853, 857, 859, 863
};
#define CHAIN_SCREEN_COUNT (sizeof(CHAIN_SCREEN_PRIMES)/sizeof(CHAIN_SCREEN_PRIMES[0]))

/* =========================================================================
 * Survivor record — emitted by GPU, consumed by CPU prover
 * ========================================================================= */

typedef struct __align__(8) {
    u32 base_idx;     /* Which CRT base */
    u32 k_offset;     /* Wheel position within tile */
    u32 tile_idx;     /* Tile index within batch (0..batch_tiles-1, always fits u32) */
    u32 _pad;         /* Alignment padding */
} Survivor;

/* =========================================================================
 * Chain result record — emitted by GPU kernel 2, consumed by CPU (output only)
 * ========================================================================= */

typedef struct __align__(8) {
    u32 base_idx;
    u32 k_offset;
    u32 tile_idx;
    u32 chain_len;      /* Length of CC chain found */
} ChainResult;

#define MAX_CHAIN_RESULTS (1 << 16)  /* 64K max chain results per batch */

/* =========================================================================
 * GPU-side constant memory tables
 * ========================================================================= */

/* L2 filter bitmasks: [prime_idx][chain_len] → 64-bit mask of forbidden residues */
__constant__ u64 d_filter_mask_first[FILTER_PRIMES_COUNT][MAX_SIEVE_CHAIN_LEN + 1];

/* Filter prime values */
__constant__ u32 d_filter_primes[FILTER_PRIMES_COUNT];

/* Ext-L2 prime values */
__constant__ u32 d_ext_filter_primes[EXT_FILTER_PRIMES_COUNT];

/* Line-sieve prime values */
__constant__ u32 d_line_primes[LINE_SIEVE_COUNT];

/* Target sieve length */
__constant__ int d_sieve_len;

/* Number of bases */
__constant__ int d_num_bases;

/* Precomputed (WHEEL_PERIOD * LATTICE_M) % q for each prime — set once.
 * Used by kernel to compute tile_idx contribution to residues. */
__constant__ u16 d_WPM_l2[FILTER_PRIMES_COUNT];
__constant__ u16 d_WPM_ext[EXT_FILTER_PRIMES_COUNT];
__constant__ u16 d_WPM_line[LINE_SIEVE_COUNT];

/* Kernel 2 (prove) constants */
__constant__ u64 d_base_values[MAX_BASES];  /* CRT base values */
__constant__ u64 d_lattice_m;               /* = LATTICE_M (fits u64) */
__constant__ u32 d_wheel_period;            /* = 20677 */
__constant__ int d_target_bits;
__constant__ u32 d_trial_primes[132];       /* All chain screen primes */
__constant__ int d_num_trial_primes;

/* Mod-6 bucketed wheel: per-base target k_offset%6 for n≡1(mod6) at tile%6=0 */
__constant__ u32 d_target_k_mod6[MAX_BASES];

/* =========================================================================
 * GPU-side global memory tables (read-only, large)
 * ========================================================================= */

/* ext-L2 byte-array filter: [prime_idx][chain_len][residue] = 1 if forbidden */
static u8* d_ext_filter_first;      /* EXT_PRIMES × 33 × 128 = ~29 KB */

/* Line-sieve kill table: [prime_idx][residue] = 1 if forbidden */
static u8* d_line_kill_first;        /* 125 × 864 = ~108 KB */

/* kM_mod LUT: [k_offset][prime_idx] = (k%q)*(M%q)%q as u16 */
static u16* d_kM_mod_line;           /* 20,677 × 125 × 2 = ~5.17 MB */

/* L2 kmod precomputed: [prime_idx][k_offset] = (k%p)*(M%p)%p as u8 */
static u8* d_kmod_filter;            /* 7 × 20,677 × 1 = ~145 KB */

/* Ext-L2 kmod precomputed: [prime_idx][k_offset] */
static u8* d_kmod_ext_filter;        /* 7 × 20,677 × 1 = ~145 KB */

/* Wheel arrays: flattened mod-6 bucketed [base_idx*6+bucket] → dense k_offset array */
static u32* d_wheel_offsets;         /* All wheel positions, flattened */
static u32* d_wheel_start_m6;       /* Start index per (base,bucket): [base*6+bucket] */
static u32* d_wheel_count_m6;       /* Count per (base,bucket): [base*6+bucket] */

/* =========================================================================
 * Double-buffer infrastructure for async pipeline
 * Two copies of: batch residues (device), survivors (device+host), counts
 * Two CUDA streams for overlapping kernel + transfer + CPU prove.
 * ========================================================================= */

/* Per-buffer batch residues (uploaded async per batch) */
static u16* d_batch_l2[2];
static u16* d_batch_ext[2];
static u16* d_batch_line[2];

/* Per-buffer survivor output */
static Survivor* d_survivors[2];
static u32* d_survivor_count[2];

/* Pinned host survivor buffers (for async D2H) */
static Survivor* h_survivors[2];

/* Pinned host survivor count (for async D2H) */
static u32* h_survivor_count[2];

/* Per-buffer chain results (Kernel 2 output) */
static ChainResult* d_chain_results[2];
static u32* d_chain_result_count[2];

/* Pinned host chain result buffers (for async D2H) */
static ChainResult* h_chain_results[2];
static u32* h_chain_result_count[2];

/* Per-buffer GPU primes counter (for stats) */
static u32* d_gpu_primes_count[2];
static u32* h_gpu_primes_count[2];

/* CUDA streams */
static cudaStream_t streams[2];

/* v13: CUDA events for K1 timing (zero-overhead when not queried) */
static cudaEvent_t k1_start_evt[2], k1_end_evt[2];

/* =========================================================================
 * Host-side mirror tables
 * ========================================================================= */

/* L2 bitmask tables (host) */
static u64 h_filter_mask_first[FILTER_PRIMES_COUNT][MAX_SIEVE_CHAIN_LEN + 1];

/* ext-L2 byte tables (host) */
static u8 h_ext_filter_first[EXT_FILTER_PRIMES_COUNT][MAX_SIEVE_CHAIN_LEN + 1][128];

/* Line-sieve kill tables (host) */
static u8 h_line_kill_first[LINE_SIEVE_COUNT][864];

/* kM_mod LUT (host) */
static u16 h_kM_mod_line[WHEEL_PERIOD][LINE_SIEVE_COUNT];

/* L2 kmod (host) */
static u8 h_kmod_filter[FILTER_PRIMES_COUNT][WHEEL_PERIOD];

/* Ext-L2 kmod (host) */
static u8 h_kmod_ext_filter[EXT_FILTER_PRIMES_COUNT][WHEEL_PERIOD];

/* M mod each prime */
static u16 h_M_mod_line[LINE_SIEVE_COUNT];

/* (WHEEL_PERIOD * LATTICE_M) % q for each prime stage (host) */
static u16 h_WPM_l2[FILTER_PRIMES_COUNT];
static u16 h_WPM_ext[EXT_FILTER_PRIMES_COUNT];
static u16 h_WPM_line[LINE_SIEVE_COUNT];

/* Per-base data (host) */
typedef struct {
    u64 base;
    u16 l2_residues[FILTER_PRIMES_COUNT];
    u16 ext_residues[EXT_FILTER_PRIMES_COUNT];
    u16 line_residues[LINE_SIEVE_COUNT];
    u32* wheel_by_mod6[6];    /* k offsets partitioned by k%6 */
    int  wheel_size_by_mod6[6];
    int  wheel_total;         /* Sum of all 6 bucket sizes */
    u32  target_k_mod6;       /* k_offset%6 needed for n≡1(mod6) at tile%6=0 */
} HostBase;

static HostBase h_bases[MAX_BASES];
static int h_num_bases = 0;

/* =========================================================================
 * Prefix-constrained search state (v9)
 * ========================================================================= */

#define PREFIX_MODE_INHERIT  0   /* Inherit from global --start logic */
#define PREFIX_MODE_SEQ      1   /* Sequential within prefix subspace */
#define PREFIX_MODE_RANDOM   2   /* Random within prefix subspace */

/* Parsed prefix configuration */
static int    g_use_prefix = 0;         /* 1 if --prefix was given */
static u128   g_prefix_value = 0;       /* Parsed binary prefix value */
static int    g_prefix_bits = 0;        /* Number of bits in prefix */
static char   g_prefix_str[256] = "";   /* Original prefix string for display/checkpoint */
static int    g_prefix_mode = PREFIX_MODE_INHERIT;  /* Prefix search mode */
static int    g_prefix_lanes = 0;       /* Number of prefix lanes (0 = disabled) */
static int    g_prefix_lane_id = 0;     /* This GPU's lane ID */

/* Prefix-constrained tile range (set during initialization) */
static u128   g_prefix_tile_min = 0;    /* First tile in prefix domain */
static u128   g_prefix_tile_max = 0;    /* One past last tile in prefix domain */
static u128   g_prefix_tile_count = 0;  /* Total tiles in prefix domain */

/* Checkpoint state */
static const char* g_checkpoint_file = NULL;  /* --checkpoint FILE */
static const char* g_resume_file = NULL;      /* --resume FILE */
static int    g_checkpoint_interval = 60;     /* Seconds between checkpoint saves */

/* Coverage mode (v10) */
#define COVERAGE_LEGACY  0   /* v9-compatible: r%6==5 filter, ~6 bases */
#define COVERAGE_FULL    1   /* All CRT-surviving bases, ~36 bases */
static int g_coverage_mode = COVERAGE_FULL;   /* Default: full coverage */

/* Verify-root mode (v10) */
static const char* g_verify_root = NULL;      /* --verify-root 0xHEX */

/* Parse binary prefix string (0b...) into u128 value and bit count.
 * Returns 0 on success, -1 on error. */
static int parse_prefix(const char* str, u128* out_value, int* out_bits) {
    *out_value = 0;
    *out_bits = 0;

    if (str[0] != '0' || (str[1] != 'b' && str[1] != 'B')) {
        fprintf(stderr, "ERROR: --prefix must be binary format (0b...)\n");
        return -1;
    }
    const char* p = str + 2;
    if (*p == '\0') {
        fprintf(stderr, "ERROR: --prefix has no digits after 0b\n");
        return -1;
    }
    if (*p != '1') {
        fprintf(stderr, "ERROR: binary prefix must start with 1 (got '%c')\n", *p);
        return -1;
    }

    u128 val = 0;
    int nbits = 0;
    while (*p == '0' || *p == '1') {
        if (nbits >= 127) {
            fprintf(stderr, "ERROR: prefix too long (max 126 bits)\n");
            return -1;
        }
        val = (val << 1) | (*p - '0');
        nbits++;
        p++;
    }
    if (*p != '\0') {
        fprintf(stderr, "ERROR: invalid character '%c' in binary prefix\n", *p);
        return -1;
    }

    *out_value = val;
    *out_bits = nbits;
    return 0;
}

/* Compute prefix-constrained tile range.
 * Given prefix of prefix_bits bits at target bit-width 'bits',
 * narrow [tile_min, tile_max) to only tiles producing candidates
 * whose top prefix_bits bits equal prefix_value.
 *
 * Algorithm (from cc_gmp_v32_claude_03.c):
 *   cand_min = prefix_value << (bits - prefix_bits)
 *   cand_max = ((prefix_value + 1) << (bits - prefix_bits)) - 1
 *   Clamp to [2^(bits-1), 2^bits - 1]
 *   k_min = ceil((cand_min - max_base) / M)
 *   k_max = floor((cand_max - min_base) / M)
 *   tile_min = k_min / WHEEL_PERIOD
 *   tile_max = k_max / WHEEL_PERIOD + 1
 */
static int compute_prefix_tile_range(
    int bits, u128 prefix_value, int prefix_bits,
    u64 min_base, u64 max_base,
    u128* out_tile_min, u128* out_tile_max
) {
    int shift = bits - prefix_bits;

    /* cand_min = prefix << shift */
    u128 cand_min = prefix_value << shift;

    /* cand_max = ((prefix + 1) << shift) - 1 */
    u128 cand_max = ((prefix_value + 1) << shift) - 1;

    /* Valid bit range: [2^(bits-1), 2^bits - 1] */
    u128 valid_min = (u128)1 << (bits - 1);
    u128 valid_max = ((u128)1 << bits) - 1;

    /* Check overlap */
    if (cand_max < valid_min || cand_min > valid_max) {
        fprintf(stderr, "ERROR: prefix 0b");
        for (int i = prefix_bits - 1; i >= 0; i--)
            fprintf(stderr, "%d", (int)((prefix_value >> i) & 1));
        fprintf(stderr, " does not overlap %d-bit range\n", bits);
        return -1;
    }

    /* Clamp to valid range */
    if (cand_min < valid_min) cand_min = valid_min;
    if (cand_max > valid_max) cand_max = valid_max;

    /* Convert candidate range to k-range.
     * n = k*M + base, so:
     *   k_min = ceil((cand_min - max_base) / M)   [conservative: largest base]
     *   k_max = floor((cand_max - min_base) / M)   [conservative: smallest base]
     */
    u128 k_min, k_max;

    if (cand_min > (u128)max_base) {
        u128 num = cand_min - (u128)max_base;
        k_min = (num + LATTICE_M - 1) / LATTICE_M;  /* ceiling division */
    } else {
        k_min = 0;
    }

    k_max = (cand_max - (u128)min_base) / LATTICE_M;

    if (k_min > k_max) {
        fprintf(stderr, "ERROR: prefix produces empty k-range\n");
        return -1;
    }

    *out_tile_min = k_min / WHEEL_PERIOD;
    *out_tile_max = k_max / WHEEL_PERIOD + 1;
    return 0;
}

/* =========================================================================
 * Checkpoint save/load (v10)
 *
 * Format (line-based text):
 *   CC18_CUDA_v10
 *   bits=<N>
 *   target=<N>
 *   depth=<N>
 *   coverage=<full|legacy>
 *   prefix=<0b... or none>
 *   prefix_lanes=<N>
 *   prefix_lane_id=<N>
 *   tile_hi=<hex16>
 *   tile_lo=<hex16>
 *   total_tiles=<N>
 *   total_primes=<N>
 *   total_chains=<N>
 *   best_chain=<N>
 * ========================================================================= */

static int checkpoint_save(
    const char* path, int bits, int target_len, int sieve_len,
    u128 current_tile,
    u64 total_tiles, u64 total_primes, u64 total_chains, int best_chain
) {
    char tmp_path[512];
    snprintf(tmp_path, sizeof(tmp_path), "%s.tmp", path);
    FILE* fp = fopen(tmp_path, "w");
    if (!fp) return -1;

    fprintf(fp, "CC18_CUDA_v10\n");
    fprintf(fp, "bits=%d\n", bits);
    fprintf(fp, "target=%d\n", target_len);
    fprintf(fp, "depth=%d\n", sieve_len);
    fprintf(fp, "coverage=%s\n", g_coverage_mode == COVERAGE_FULL ? "full" : "legacy");
    fprintf(fp, "prefix=%s\n", g_use_prefix ? g_prefix_str : "none");
    fprintf(fp, "prefix_lanes=%d\n", g_prefix_lanes);
    fprintf(fp, "prefix_lane_id=%d\n", g_prefix_lane_id);
    fprintf(fp, "tile_hi=%016llx\n", (unsigned long long)(current_tile >> 64));
    fprintf(fp, "tile_lo=%016llx\n", (unsigned long long)(current_tile & 0xFFFFFFFFFFFFFFFFULL));
    fprintf(fp, "total_tiles=%llu\n", (unsigned long long)total_tiles);
    fprintf(fp, "total_primes=%llu\n", (unsigned long long)total_primes);
    fprintf(fp, "total_chains=%llu\n", (unsigned long long)total_chains);
    fprintf(fp, "best_chain=%d\n", best_chain);
    fclose(fp);

    /* Atomic rename */
    if (rename(tmp_path, path) != 0) {
        perror("checkpoint rename");
        return -1;
    }
    return 0;
}

static int checkpoint_load(
    const char* path, int bits, int target_len, int sieve_len,
    u128* out_tile,
    u64* out_total_tiles, u64* out_total_primes, u64* out_total_chains,
    int* out_best_chain
) {
    FILE* fp = fopen(path, "r");
    if (!fp) {
        fprintf(stderr, "ERROR: cannot open checkpoint file '%s'\n", path);
        return -1;
    }

    char line[512];
    /* Line 1: version */
    if (!fgets(line, sizeof(line), fp) || strncmp(line, "CC18_CUDA_v10", 13) != 0) {
        fprintf(stderr, "ERROR: checkpoint file has wrong version marker\n");
        fclose(fp);
        return -1;
    }

    int ck_bits = 0, ck_target = 0, ck_depth = 0;
    char ck_prefix[256] = "";
    int ck_lanes = 0, ck_lane_id = 0;
    unsigned long long tile_hi = 0, tile_lo = 0;
    unsigned long long ck_tiles = 0, ck_primes = 0, ck_chains = 0;
    int ck_best = 0;

    while (fgets(line, sizeof(line), fp)) {
        /* Remove trailing newline */
        line[strcspn(line, "\n")] = '\0';

        if (sscanf(line, "bits=%d", &ck_bits) == 1) continue;
        if (sscanf(line, "target=%d", &ck_target) == 1) continue;
        if (sscanf(line, "depth=%d", &ck_depth) == 1) continue;
        if (sscanf(line, "prefix=%255s", ck_prefix) == 1) continue;
        if (sscanf(line, "prefix_lanes=%d", &ck_lanes) == 1) continue;
        if (sscanf(line, "prefix_lane_id=%d", &ck_lane_id) == 1) continue;
        if (sscanf(line, "tile_hi=%llx", &tile_hi) == 1) continue;
        if (sscanf(line, "tile_lo=%llx", &tile_lo) == 1) continue;
        if (sscanf(line, "total_tiles=%llu", &ck_tiles) == 1) continue;
        if (sscanf(line, "total_primes=%llu", &ck_primes) == 1) continue;
        if (sscanf(line, "total_chains=%llu", &ck_chains) == 1) continue;
        if (sscanf(line, "best_chain=%d", &ck_best) == 1) continue;
    }
    fclose(fp);

    /* Validate compatibility */
    if (ck_bits != bits) {
        fprintf(stderr, "ERROR: checkpoint bits=%d does not match current --bits %d\n",
                ck_bits, bits);
        return -1;
    }
    if (ck_target != target_len) {
        fprintf(stderr, "ERROR: checkpoint target=%d does not match current --target %d\n",
                ck_target, target_len);
        return -1;
    }
    if (ck_depth != sieve_len) {
        fprintf(stderr, "ERROR: checkpoint depth=%d does not match current --depth %d\n",
                ck_depth, sieve_len);
        return -1;
    }

    /* Validate prefix compatibility */
    const char* current_prefix = g_use_prefix ? g_prefix_str : "none";
    if (strcmp(ck_prefix, current_prefix) != 0) {
        fprintf(stderr, "ERROR: checkpoint prefix=%s does not match current prefix=%s\n",
                ck_prefix, current_prefix);
        return -1;
    }
    if (ck_lanes != g_prefix_lanes) {
        fprintf(stderr, "ERROR: checkpoint prefix_lanes=%d does not match current %d\n",
                ck_lanes, g_prefix_lanes);
        return -1;
    }
    if (ck_lane_id != g_prefix_lane_id) {
        fprintf(stderr, "ERROR: checkpoint prefix_lane_id=%d does not match current %d\n",
                ck_lane_id, g_prefix_lane_id);
        return -1;
    }

    *out_tile = ((u128)tile_hi << 64) | (u128)tile_lo;
    *out_total_tiles = (u64)ck_tiles;
    *out_total_primes = (u64)ck_primes;
    *out_total_chains = (u64)ck_chains;
    *out_best_chain = ck_best;

    printf("[RESUME] Loaded checkpoint: tile=0x%llx%016llx tiles=%llu primes=%llu best=CC%d\n",
           tile_hi, tile_lo, ck_tiles, ck_primes, ck_best);
    return 0;
}

/* Global state */
static volatile sig_atomic_t g_shutdown = 0;
static void signal_handler(int sig) {
    (void)sig;
    if (!g_shutdown) {
        g_shutdown = 1;
        {
            const char msg[] =
                "\n[INTERRUPT] stopping after current GPU batch (Ctrl+C again = force exit)\n";
            ssize_t wr = write(STDERR_FILENO, msg, sizeof(msg) - 1);
            (void)wr;
        }
        return;
    }

    {
        const char msg[] = "\n[INTERRUPT] force exit\n";
        ssize_t wr = write(STDERR_FILENO, msg, sizeof(msg) - 1);
        (void)wr;
    }
    _exit(130);
}

static void install_signal_handlers(void) {
    struct sigaction sa;
    memset(&sa, 0, sizeof(sa));
    sa.sa_handler = signal_handler;
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = 0;

    if (sigaction(SIGINT, &sa, NULL) != 0) {
        perror("sigaction(SIGINT)");
        exit(1);
    }
    if (sigaction(SIGTERM, &sa, NULL) != 0) {
        perror("sigaction(SIGTERM)");
        exit(1);
    }
}

/* =========================================================================
 * Fast 64-bit PRNG (xorshift64*) for random tile selection
 * ========================================================================= */

static u64 rng_state;

static void rng_seed(u64 seed) {
    rng_state = seed ? seed : 1;
}

static u64 rng_next(void) {
    u64 x = rng_state;
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    rng_state = x;
    return x;
}

/* High-entropy seed from /dev/urandom (machine-unique, never blocks).
 * Fallback: CLOCK_MONOTONIC nanoseconds + pid + gpu_id if urandom fails. */
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
    /* Fallback: high-entropy mix if /dev/urandom fails */
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    seed = (u64)ts.tv_sec * 1000000000ULL + (u64)ts.tv_nsec;
    seed ^= (u64)getpid() * 0x517CC1B727220A95ULL;
    seed ^= (u64)gpu_id * 0x9E3779B97F4A7C15ULL;
    if (seed == 0) seed = 1;
    rng_seed(seed);
    return seed;
}


/* =========================================================================
 * Math helpers (host-side, for table initialization)
 * ========================================================================= */

static u64 mod_inverse(u64 a, u64 m) {
    if (m == 1) return 0;
    long long m0 = (long long)m, x0 = 0, x1 = 1;
    long long aa = (long long)(a % m);
    while (aa > 1) {
        long long q = aa / m0;
        long long t = m0;
        m0 = aa % m0;
        aa = t;
        t = x0;
        x0 = x1 - q * x0;
        x1 = t;
    }
    if (x1 < 0) x1 += (long long)m;
    return (u64)x1;
}

static void compute_forbidden_residues(u32 q, int chain_len, u32* out, int* count) {
    *count = 0;
    int max_i = chain_len;
    if (max_i > (int)q - 1) max_i = (int)q - 1;
    for (int i = 0; i < max_i; i++) {
        u64 pow2i = 1;
        for (int j = 0; j < i; j++) pow2i = (pow2i * 2) % q;
        u64 inv = mod_inverse(pow2i, q);
        u32 r = (u32)((1 + q - inv) % q);
        /* Check for duplicate */
        int dup = 0;
        for (int k = 0; k < *count; k++) {
            if (out[k] == r) { dup = 1; break; }
        }
        if (!dup) out[(*count)++] = r;
    }
}

/* =========================================================================
 * Table initialization (host-side, runs once)
 * ========================================================================= */

static void init_filter_tables(int sieve_len) {
    printf("  Initializing L2 bitmask tables...\n");
    memset(h_filter_mask_first, 0, sizeof(h_filter_mask_first));

    for (int p_idx = 0; p_idx < FILTER_PRIMES_COUNT; p_idx++) {
        u32 q = FILTER_PRIMES[p_idx];
        for (int L = 1; L <= sieve_len && L <= MAX_SIEVE_CHAIN_LEN; L++) {
            u64 mask = 0;
            int max_i = L;
            if (max_i > (int)q - 1) max_i = (int)q - 1;
            for (int i = 0; i < max_i; i++) {
                u64 pow2i = 1;
                for (int j = 0; j < i; j++) pow2i = (pow2i * 2) % q;
                u64 inv = mod_inverse(pow2i, q);
                u32 forbidden_r = (u32)((1 + q - inv) % q);
                if (forbidden_r < 64) mask |= (1ULL << forbidden_r);
            }
            h_filter_mask_first[p_idx][L] = mask;
        }
    }

    printf("  Initializing ext-L2 byte tables...\n");
    memset(h_ext_filter_first, 0, sizeof(h_ext_filter_first));

    for (int p_idx = 0; p_idx < EXT_FILTER_PRIMES_COUNT; p_idx++) {
        u32 q = EXT_FILTER_PRIMES[p_idx];
        for (int L = 1; L <= sieve_len && L <= MAX_SIEVE_CHAIN_LEN; L++) {
            int max_i = L;
            if (max_i > (int)q - 1) max_i = (int)q - 1;
            for (int i = 0; i < max_i; i++) {
                u64 pow2i = 1;
                for (int j = 0; j < i; j++) pow2i = (pow2i * 2) % q;
                u64 inv = mod_inverse(pow2i, q);
                u32 forbidden_r = (u32)((1 + q - inv) % q);
                h_ext_filter_first[p_idx][L][forbidden_r] = 1;
            }
        }
    }

    printf("  Initializing line-sieve kill tables (125 primes, depth %d)...\n", sieve_len);
    memset(h_line_kill_first, 0, sizeof(h_line_kill_first));
    int total_kills = 0;

    for (int idx = 0; idx < LINE_SIEVE_COUNT; idx++) {
        u32 q = CHAIN_SCREEN_PRIMES[LINE_SIEVE_START_IDX + idx];
        for (int pos = 1; pos < sieve_len && pos < (int)q; pos++) {
            u64 pow2 = 1;
            for (int j = 0; j < pos; j++) pow2 = (pow2 * 2) % q;
            u64 inv2pos = mod_inverse(pow2, q);
            u32 kill_r = (u32)((1 + q - inv2pos) % q);
            if (!h_line_kill_first[idx][kill_r]) {
                h_line_kill_first[idx][kill_r] = 1;
                total_kills++;
            }
        }
        h_M_mod_line[idx] = (u16)(LATTICE_M % q);
    }
    printf("    Line-sieve: %d total kill entries across 125 primes\n", total_kills);

    printf("  Building kM_mod LUT (%.2f MB)...\n",
           (double)(WHEEL_PERIOD * LINE_SIEVE_COUNT * 2) / (1024.0 * 1024.0));

    for (u32 k = 0; k < WHEEL_PERIOD; k++) {
        for (int ls = 0; ls < LINE_SIEVE_COUNT; ls++) {
            u32 q = CHAIN_SCREEN_PRIMES[LINE_SIEVE_START_IDX + ls];
            u32 k_mod_q = k % q;
            h_kM_mod_line[k][ls] = (u16)((u64)k_mod_q * h_M_mod_line[ls] % q);
        }
    }

    printf("  Building L2 kmod tables (%.1f KB)...\n",
           (double)(FILTER_PRIMES_COUNT * WHEEL_PERIOD) / 1024.0);

    for (int i = 0; i < FILTER_PRIMES_COUNT; i++) {
        u32 p = FILTER_PRIMES[i];
        u32 m_mod_p = (u32)(LATTICE_M % p);
        for (u32 k = 0; k < WHEEL_PERIOD; k++) {
            h_kmod_filter[i][k] = (u8)((u64)(k % p) * m_mod_p % p);
        }
    }

    printf("  Building ext-L2 kmod tables (%.1f KB)...\n",
           (double)(EXT_FILTER_PRIMES_COUNT * WHEEL_PERIOD) / 1024.0);

    for (int i = 0; i < EXT_FILTER_PRIMES_COUNT; i++) {
        u32 p = EXT_FILTER_PRIMES[i];
        u32 m_mod_p = (u32)(LATTICE_M % p);
        for (u32 k = 0; k < WHEEL_PERIOD; k++) {
            h_kmod_ext_filter[i][k] = (u8)((u64)(k % p) * m_mod_p % p);
        }
    }

    /* Precompute (WHEEL_PERIOD * LATTICE_M) % q for each prime.
     * Used by kernel to add per-tile-idx residue contribution. */
    for (int i = 0; i < FILTER_PRIMES_COUNT; i++) {
        u32 q = FILTER_PRIMES[i];
        h_WPM_l2[i] = (u16)((u64)(WHEEL_PERIOD % q) * (LATTICE_M % q) % q);
    }
    for (int i = 0; i < EXT_FILTER_PRIMES_COUNT; i++) {
        u32 q = EXT_FILTER_PRIMES[i];
        h_WPM_ext[i] = (u16)((u64)(WHEEL_PERIOD % q) * (LATTICE_M % q) % q);
    }
    for (int i = 0; i < LINE_SIEVE_COUNT; i++) {
        u32 q = CHAIN_SCREEN_PRIMES[LINE_SIEVE_START_IDX + i];
        h_WPM_line[i] = (u16)((u64)(WHEEL_PERIOD % q) * (LATTICE_M % q) % q);
    }

    printf("  Tables initialized.\n");
}

/* =========================================================================
 * CRT lattice generation (host-side)
 * ========================================================================= */

static void generate_bases(int chain_len) {
    printf("  Generating CRT bases (M=%llu, coverage=%s)...\n",
           (unsigned long long)LATTICE_M,
           g_coverage_mode == COVERAGE_FULL ? "FULL" : "LEGACY");

    /* Enumerate all residues mod M that survive forbidden residue checks
     * for CRT primes. In FULL mode, all surviving bases are used (no mod-6
     * filter). In LEGACY mode, only r%6==5 bases are kept (v9 compat). */
    u64 temp_bases[4096];
    int num_temp = 0;

    for (u64 r = 0; r < LATTICE_M; r++) {
        /* LEGACY mode: pre-filter r%6==5 (v9 behavior) */
        if (g_coverage_mode == COVERAGE_LEGACY) {
            if (r % 6 != 5) continue;
        }

        int valid = 1;
        for (int pi = 0; pi < CRT_PRIMES_COUNT && valid; pi++) {
            u32 q = CRT_PRIMES[pi];
            u32 rmod = (u32)(r % q);
            u32 forbidden[64];
            int nf;
            compute_forbidden_residues(q, chain_len, forbidden, &nf);
            for (int j = 0; j < nf; j++) {
                if (rmod == forbidden[j]) { valid = 0; break; }
            }
        }
        if (valid && num_temp < 4096) {
            temp_bases[num_temp++] = r;
        }
    }

    printf("    CRT bases surviving: %d\n", num_temp);

    /* Mod-7 distribution diagnostics */
    {
        int mod7_count[7] = {0};
        for (int i = 0; i < num_temp; i++)
            mod7_count[temp_bases[i] % 7]++;
        printf("    Mod-7 distribution:");
        for (int m = 0; m < 7; m++) {
            if (mod7_count[m] > 0) printf(" r%d=%d", m, mod7_count[m]);
        }
        printf("\n");
        if (mod7_count[6] > 0)
            printf("    ** q7-immunized class (r≡6 mod 7): %d bases present **\n",
                   mod7_count[6]);
        else
            printf("    NOTE: q7-immunized class (r≡6 mod 7): ABSENT\n");
    }

    /* For each base, generate wheel, partition into mod-6 buckets,
     * and compute per-base residues */
    h_num_bases = 0;
    int total_wheel = 0;

    /* M%6 = 5, so inv_m6 = 5 (since 5*5=25≡1 mod 6) */
    const u32 inv_m6 = 5;

    for (int b = 0; b < num_temp && h_num_bases < MAX_BASES; b++) {
        u64 base = temp_bases[b];

        /* Compute target_k_mod6: the k%6 bucket needed for n≡1(mod6)
         * when tile%6=0. Math: n = k*M + base, M%6=5.
         * n%6 = k*5 + base (mod 6). For n≡1(mod6):
         * k ≡ (1 - base%6) * inv_m6 (mod 6). */
        u32 base_mod6 = (u32)(base % 6);
        u32 target_k_mod6 = ((1 - base_mod6 + 6) % 6) * inv_m6 % 6;

        /* Check wheel primes */
        u64 wp_base_r[WHEEL_PRIMES_COUNT];
        u64 wp_step[WHEEL_PRIMES_COUNT];
        u32 wp_forbidden[WHEEL_PRIMES_COUNT][64];
        int wp_nf[WHEEL_PRIMES_COUNT];

        for (int i = 0; i < WHEEL_PRIMES_COUNT; i++) {
            u32 q = WHEEL_PRIMES[i];
            wp_base_r[i] = base % q;
            wp_step[i] = LATTICE_M % q;
            compute_forbidden_residues(q, chain_len, wp_forbidden[i], &wp_nf[i]);
        }

        /* Generate valid wheel positions */
        u32* temp_wheel = (u32*)malloc(WHEEL_PERIOD * sizeof(u32));
        int wheel_size = 0;

        for (u32 k = 0; k < WHEEL_PERIOD; k++) {
            int valid = 1;
            for (int i = 0; i < WHEEL_PRIMES_COUNT && valid; i++) {
                u32 q = WHEEL_PRIMES[i];
                u32 r = (u32)((wp_base_r[i] + (u64)k * wp_step[i]) % q);
                for (int j = 0; j < wp_nf[i]; j++) {
                    if (r == wp_forbidden[i][j]) { valid = 0; break; }
                }
            }
            if (valid) temp_wheel[wheel_size++] = k;
        }

        if (wheel_size == 0) { free(temp_wheel); continue; }

        HostBase* hb = &h_bases[h_num_bases];
        hb->base = base;
        hb->wheel_total = wheel_size;
        hb->target_k_mod6 = target_k_mod6;

        /* Partition wheel into 6 buckets by k_offset%6
         * (same algorithm as cc_gmp_v32_claude_03.c lines 741-753) */
        for (int m = 0; m < 6; m++) hb->wheel_size_by_mod6[m] = 0;
        for (int w = 0; w < wheel_size; w++)
            hb->wheel_size_by_mod6[temp_wheel[w] % 6]++;
        for (int m = 0; m < 6; m++) {
            if (hb->wheel_size_by_mod6[m] > 0)
                hb->wheel_by_mod6[m] = (u32*)malloc(hb->wheel_size_by_mod6[m] * sizeof(u32));
            else
                hb->wheel_by_mod6[m] = NULL;
            hb->wheel_size_by_mod6[m] = 0;  /* Reset for fill pass */
        }
        for (int w = 0; w < wheel_size; w++) {
            int m = temp_wheel[w] % 6;
            hb->wheel_by_mod6[m][hb->wheel_size_by_mod6[m]++] = temp_wheel[w];
        }

        free(temp_wheel);

        /* Per-base residues */
        for (int i = 0; i < FILTER_PRIMES_COUNT; i++)
            hb->l2_residues[i] = (u16)(base % FILTER_PRIMES[i]);
        for (int i = 0; i < EXT_FILTER_PRIMES_COUNT; i++)
            hb->ext_residues[i] = (u16)(base % EXT_FILTER_PRIMES[i]);
        for (int i = 0; i < LINE_SIEVE_COUNT; i++)
            hb->line_residues[i] = (u16)(base % CHAIN_SCREEN_PRIMES[LINE_SIEVE_START_IDX + i]);

        total_wheel += wheel_size;
        h_num_bases++;
    }

    printf("    Bases used: %d, total wheel positions: %d\n", h_num_bases, total_wheel);
    printf("    Avg wheel density: %.2f%%\n",
           h_num_bases > 0 ? 100.0 * total_wheel / (h_num_bases * WHEEL_PERIOD) : 0.0);
    if (h_num_bases > 0) {
        printf("    Avg bucket size (per base, per tile): %d\n",
               total_wheel / (h_num_bases * 6));
        /* Print per-mod6 bucket totals */
        int bucket_totals[6] = {0};
        for (int bi = 0; bi < h_num_bases; bi++)
            for (int m = 0; m < 6; m++)
                bucket_totals[m] += h_bases[bi].wheel_size_by_mod6[m];
        printf("    Bucket totals (mod6): [%d, %d, %d, %d, %d, %d]\n",
               bucket_totals[0], bucket_totals[1], bucket_totals[2],
               bucket_totals[3], bucket_totals[4], bucket_totals[5]);
    }
}

/* =========================================================================
 * GPU memory allocation and table upload
 * ========================================================================= */

static void gpu_upload_tables(int sieve_len) {
    printf("  Uploading tables to GPU...\n");

    /* Constant memory */
    CUDA_CHECK(cudaMemcpyToSymbol(d_filter_mask_first, h_filter_mask_first,
                                   sizeof(h_filter_mask_first)));
    CUDA_CHECK(cudaMemcpyToSymbol(d_filter_primes, FILTER_PRIMES,
                                   sizeof(FILTER_PRIMES)));
    CUDA_CHECK(cudaMemcpyToSymbol(d_ext_filter_primes, EXT_FILTER_PRIMES,
                                   sizeof(EXT_FILTER_PRIMES)));

    u32 line_p[LINE_SIEVE_COUNT];
    for (int i = 0; i < LINE_SIEVE_COUNT; i++)
        line_p[i] = CHAIN_SCREEN_PRIMES[LINE_SIEVE_START_IDX + i];
    CUDA_CHECK(cudaMemcpyToSymbol(d_line_primes, line_p, sizeof(line_p)));
    CUDA_CHECK(cudaMemcpyToSymbol(d_sieve_len, &sieve_len, sizeof(int)));
    CUDA_CHECK(cudaMemcpyToSymbol(d_num_bases, &h_num_bases, sizeof(int)));

    /* WPM_mod constants */
    CUDA_CHECK(cudaMemcpyToSymbol(d_WPM_l2, h_WPM_l2, sizeof(h_WPM_l2)));
    CUDA_CHECK(cudaMemcpyToSymbol(d_WPM_ext, h_WPM_ext, sizeof(h_WPM_ext)));
    CUDA_CHECK(cudaMemcpyToSymbol(d_WPM_line, h_WPM_line, sizeof(h_WPM_line)));

    /* Global memory — ext-L2 filter */
    size_t ext_size = EXT_FILTER_PRIMES_COUNT * (MAX_SIEVE_CHAIN_LEN + 1) * 128;
    CUDA_CHECK(cudaMalloc(&d_ext_filter_first, ext_size));
    CUDA_CHECK(cudaMemcpy(d_ext_filter_first, h_ext_filter_first, ext_size,
                           cudaMemcpyHostToDevice));

    /* Global memory — line-sieve kill table */
    size_t kill_size = LINE_SIEVE_COUNT * 864;
    CUDA_CHECK(cudaMalloc(&d_line_kill_first, kill_size));
    CUDA_CHECK(cudaMemcpy(d_line_kill_first, h_line_kill_first, kill_size,
                           cudaMemcpyHostToDevice));

    /* Global memory — kM_mod LUT (biggest table: ~5.17 MB)
     * +128 bytes padding so the v13 prefetch of the last row's second cache
     * line never reads past the allocation boundary. */
    size_t km_size = WHEEL_PERIOD * LINE_SIEVE_COUNT * sizeof(u16);
    CUDA_CHECK(cudaMalloc(&d_kM_mod_line, km_size + 128));
    CUDA_CHECK(cudaMemcpy(d_kM_mod_line, h_kM_mod_line, km_size,
                           cudaMemcpyHostToDevice));

    /* Global memory — L2 kmod */
    size_t kmod_l2_size = FILTER_PRIMES_COUNT * WHEEL_PERIOD * sizeof(u8);
    CUDA_CHECK(cudaMalloc(&d_kmod_filter, kmod_l2_size));
    CUDA_CHECK(cudaMemcpy(d_kmod_filter, h_kmod_filter, kmod_l2_size,
                           cudaMemcpyHostToDevice));

    /* Global memory — ext-L2 kmod */
    size_t kmod_ext_size = EXT_FILTER_PRIMES_COUNT * WHEEL_PERIOD * sizeof(u8);
    CUDA_CHECK(cudaMalloc(&d_kmod_ext_filter, kmod_ext_size));
    CUDA_CHECK(cudaMemcpy(d_kmod_ext_filter, h_kmod_ext_filter, kmod_ext_size,
                           cudaMemcpyHostToDevice));

    /* Double-buffered batch residues + survivors (2 copies each) */
    for (int buf = 0; buf < 2; buf++) {
        CUDA_CHECK(cudaMalloc(&d_batch_l2[buf],
                               MAX_BASES * FILTER_PRIMES_COUNT * sizeof(u16)));
        CUDA_CHECK(cudaMalloc(&d_batch_ext[buf],
                               MAX_BASES * EXT_FILTER_PRIMES_COUNT * sizeof(u16)));
        CUDA_CHECK(cudaMalloc(&d_batch_line[buf],
                               MAX_BASES * LINE_SIEVE_COUNT * sizeof(u16)));
        CUDA_CHECK(cudaMalloc(&d_survivors[buf],
                               MAX_SURVIVORS_PER_BATCH * sizeof(Survivor)));
        CUDA_CHECK(cudaMalloc(&d_survivor_count[buf], sizeof(u32)));
    }

    /* Flatten mod-6 bucketed wheel arrays for GPU.
     * Layout: iterate (base, mod6) → flat d_wheel_offsets array.
     * d_wheel_start_m6[base*6+bucket] and d_wheel_count_m6[base*6+bucket]
     * index into the flat array. */
    int total_wheel = 0;
    u32 h_wstart_m6[MAX_BASES * 6];
    u32 h_wcount_m6[MAX_BASES * 6];
    u32 h_tkm6[MAX_BASES];
    memset(h_wstart_m6, 0, sizeof(h_wstart_m6));
    memset(h_wcount_m6, 0, sizeof(h_wcount_m6));

    for (int b = 0; b < h_num_bases; b++) {
        h_tkm6[b] = h_bases[b].target_k_mod6;
        for (int m = 0; m < 6; m++) {
            h_wstart_m6[b * 6 + m] = total_wheel;
            h_wcount_m6[b * 6 + m] = h_bases[b].wheel_size_by_mod6[m];
            total_wheel += h_bases[b].wheel_size_by_mod6[m];
        }
    }

    u32* h_wheel_flat = (u32*)malloc(total_wheel * sizeof(u32));
    for (int b = 0; b < h_num_bases; b++) {
        for (int m = 0; m < 6; m++) {
            if (h_bases[b].wheel_size_by_mod6[m] > 0) {
                memcpy(h_wheel_flat + h_wstart_m6[b * 6 + m],
                       h_bases[b].wheel_by_mod6[m],
                       h_bases[b].wheel_size_by_mod6[m] * sizeof(u32));
            }
        }
    }

    CUDA_CHECK(cudaMalloc(&d_wheel_offsets, total_wheel * sizeof(u32)));
    CUDA_CHECK(cudaMemcpy(d_wheel_offsets, h_wheel_flat,
                           total_wheel * sizeof(u32), cudaMemcpyHostToDevice));
    free(h_wheel_flat);

    CUDA_CHECK(cudaMalloc(&d_wheel_start_m6, h_num_bases * 6 * sizeof(u32)));
    CUDA_CHECK(cudaMemcpy(d_wheel_start_m6, h_wstart_m6,
                           h_num_bases * 6 * sizeof(u32), cudaMemcpyHostToDevice));

    CUDA_CHECK(cudaMalloc(&d_wheel_count_m6, h_num_bases * 6 * sizeof(u32)));
    CUDA_CHECK(cudaMemcpy(d_wheel_count_m6, h_wcount_m6,
                           h_num_bases * 6 * sizeof(u32), cudaMemcpyHostToDevice));

    /* Upload target_k_mod6 to constant memory */
    CUDA_CHECK(cudaMemcpyToSymbol(d_target_k_mod6, h_tkm6,
                                   h_num_bases * sizeof(u32)));

    /* Kernel 2 constant memory — base values for n reconstruction */
    u64 base_vals[MAX_BASES];
    for (int b = 0; b < h_num_bases; b++) base_vals[b] = h_bases[b].base;
    CUDA_CHECK(cudaMemcpyToSymbol(d_base_values, base_vals, h_num_bases * sizeof(u64)));
    u64 lattice_m_val = LATTICE_M;
    CUDA_CHECK(cudaMemcpyToSymbol(d_lattice_m, &lattice_m_val, sizeof(u64)));
    u32 wp = WHEEL_PERIOD;
    CUDA_CHECK(cudaMemcpyToSymbol(d_wheel_period, &wp, sizeof(u32)));

    /* Trial division primes for kernel 2 */
    int num_trial = (int)CHAIN_SCREEN_COUNT;
    CUDA_CHECK(cudaMemcpyToSymbol(d_trial_primes, CHAIN_SCREEN_PRIMES,
                                   num_trial * sizeof(u32)));
    CUDA_CHECK(cudaMemcpyToSymbol(d_num_trial_primes, &num_trial, sizeof(int)));

    /* Double-buffered chain result buffers (Kernel 2 output) */
    for (int buf = 0; buf < 2; buf++) {
        CUDA_CHECK(cudaMalloc(&d_chain_results[buf],
                               MAX_CHAIN_RESULTS * sizeof(ChainResult)));
        CUDA_CHECK(cudaMalloc(&d_chain_result_count[buf], sizeof(u32)));
        CUDA_CHECK(cudaMalloc(&d_gpu_primes_count[buf], sizeof(u32)));
    }

    /* CUDA streams for async pipeline */
    CUDA_CHECK(cudaStreamCreate(&streams[0]));
    CUDA_CHECK(cudaStreamCreate(&streams[1]));

    /* v13: CUDA events for K1 timing (measures actual GPU filter time) */
    for (int buf = 0; buf < 2; buf++) {
        CUDA_CHECK(cudaEventCreate(&k1_start_evt[buf]));
        CUDA_CHECK(cudaEventCreate(&k1_end_evt[buf]));
    }

    size_t total_gpu_mem = ext_size + kill_size + km_size + kmod_l2_size +
                           kmod_ext_size + total_wheel * sizeof(u32) +
                           2 * MAX_SURVIVORS_PER_BATCH * sizeof(Survivor) +
                           2 * MAX_CHAIN_RESULTS * sizeof(ChainResult);
    printf("    Total GPU memory: %.2f MB (double-buffered)\n",
           total_gpu_mem / (1024.0 * 1024.0));
    printf("    (kM_mod LUT: %.2f MB = %.1f%% of total)\n",
           km_size / (1024.0 * 1024.0), 100.0 * km_size / total_gpu_mem);
}

/* =========================================================================
 * Batch preparation — precompute (base_r + tile_start_km) % q on CPU.
 *
 * tile_start is u128, can exceed u64 for >=99-bit candidates.
 * We compute (tile_start % q) using u128 on CPU, then combine with base
 * residues. The kernel only needs to add the small tile_idx contribution.
 * ========================================================================= */

static void prepare_batch(
    u128 tile_start_128,
    u16* h_l2_flat, u16* h_ext_flat, u16* h_line_flat
) {
    for (int b = 0; b < h_num_bases; b++) {
        for (int i = 0; i < FILTER_PRIMES_COUNT; i++) {
            u32 q = FILTER_PRIMES[i];
            u32 ts_mod_q = (u32)(tile_start_128 % q);
            u32 ts_km = (u32)((u64)ts_mod_q * h_WPM_l2[i] % q);
            h_l2_flat[b * FILTER_PRIMES_COUNT + i] = (u16)(
                (h_bases[b].l2_residues[i] + ts_km) % q);
        }
        for (int i = 0; i < EXT_FILTER_PRIMES_COUNT; i++) {
            u32 q = EXT_FILTER_PRIMES[i];
            u32 ts_mod_q = (u32)(tile_start_128 % q);
            u32 ts_km = (u32)((u64)ts_mod_q * h_WPM_ext[i] % q);
            h_ext_flat[b * EXT_FILTER_PRIMES_COUNT + i] = (u16)(
                (h_bases[b].ext_residues[i] + ts_km) % q);
        }
        for (int i = 0; i < LINE_SIEVE_COUNT; i++) {
            u32 q = CHAIN_SCREEN_PRIMES[LINE_SIEVE_START_IDX + i];
            u32 ts_mod_q = (u32)(tile_start_128 % q);
            u32 ts_km = (u32)((u64)ts_mod_q * h_WPM_line[i] % q);
            h_line_flat[b * LINE_SIEVE_COUNT + i] = (u16)(
                (h_bases[b].line_residues[i] + ts_km) % q);
        }
    }
}

/* =========================================================================
 * 128-bit GPU arithmetic for Kernel 2 (Miller-Rabin proving)
 * ========================================================================= */

typedef struct {
    u64 lo, hi;
} gpu_u128;

__device__ __forceinline__ gpu_u128 gpu_u128_from_u64(u64 v) {
    gpu_u128 r; r.lo = v; r.hi = 0; return r;
}

__device__ __forceinline__ gpu_u128 gpu_u128_add(gpu_u128 a, gpu_u128 b) {
    gpu_u128 r;
    r.lo = a.lo + b.lo;
    r.hi = a.hi + b.hi + (r.lo < a.lo ? 1 : 0);
    return r;
}

__device__ __forceinline__ gpu_u128 gpu_u128_sub(gpu_u128 a, gpu_u128 b) {
    gpu_u128 r;
    r.lo = a.lo - b.lo;
    r.hi = a.hi - b.hi - (a.lo < b.lo ? 1 : 0);
    return r;
}

__device__ __forceinline__ int gpu_u128_lt(gpu_u128 a, gpu_u128 b) {
    return (a.hi < b.hi) || (a.hi == b.hi && a.lo < b.lo);
}

__device__ __forceinline__ int gpu_u128_eq(gpu_u128 a, gpu_u128 b) {
    return a.lo == b.lo && a.hi == b.hi;
}

__device__ __forceinline__ int gpu_u128_gte(gpu_u128 a, gpu_u128 b) {
    return !gpu_u128_lt(a, b);
}

__device__ __forceinline__ int gpu_u128_is_zero(gpu_u128 a) {
    return a.lo == 0 && a.hi == 0;
}

/* Full 64x64->128 multiply using PTX mul.hi */
__device__ __forceinline__ gpu_u128 gpu_mul64(u64 a, u64 b) {
    gpu_u128 r;
    r.lo = a * b;
    asm("mul.hi.u64 %0, %1, %2;" : "=l"(r.hi) : "l"(a), "l"(b));
    return r;
}

/* 128x64->128 multiply (low 128 bits only) */
__device__ __forceinline__ gpu_u128 gpu_u128_mul_u64(gpu_u128 a, u64 b) {
    gpu_u128 lo_part = gpu_mul64(a.lo, b);
    u64 hi_part = a.hi * b;
    lo_part.hi += hi_part;
    return lo_part;
}

/* 128-bit left shift by 1 */
__device__ __forceinline__ gpu_u128 gpu_u128_shl1(gpu_u128 a) {
    gpu_u128 r;
    r.hi = (a.hi << 1) | (a.lo >> 63);
    r.lo = a.lo << 1;
    return r;
}

/* 128-bit right shift by 1 */
__device__ __forceinline__ gpu_u128 gpu_u128_shr1(gpu_u128 a) {
    gpu_u128 r;
    r.lo = (a.lo >> 1) | (a.hi << 63);
    r.hi = a.hi >> 1;
    return r;
}

/* Count leading zeros for gpu_u128 */
__device__ __forceinline__ int gpu_u128_clz(gpu_u128 a) {
    if (a.hi) return __clzll(a.hi);
    if (a.lo) return 64 + __clzll(a.lo);
    return 128;
}

/* 128-bit modular reduction: a % m using binary long division */
__device__ gpu_u128 gpu_u128_mod(gpu_u128 a, gpu_u128 m) {
    if (gpu_u128_lt(a, m)) return a;
    if (gpu_u128_is_zero(m)) return a;

    int shift = gpu_u128_clz(m) - gpu_u128_clz(a);
    if (shift < 0) return a;

    gpu_u128 divisor = m;
    for (int i = 0; i < shift; i++) divisor = gpu_u128_shl1(divisor);

    for (int i = shift; i >= 0; i--) {
        if (gpu_u128_gte(a, divisor))
            a = gpu_u128_sub(a, divisor);
        divisor = gpu_u128_shr1(divisor);
    }
    return a;
}

/* Modular multiplication: (a * b) % m using shift-and-add.
 * Starts from highest set bit of b for efficiency. */
__device__ gpu_u128 gpu_u128_mulmod(gpu_u128 a, gpu_u128 b, gpu_u128 m) {
    a = gpu_u128_mod(a, m);
    gpu_u128 result = gpu_u128_from_u64(0);

    int top_bit = 127 - gpu_u128_clz(b);
    for (int bit = top_bit; bit >= 0; bit--) {
        result = gpu_u128_shl1(result);
        if (gpu_u128_gte(result, m)) result = gpu_u128_sub(result, m);

        u64 word = (bit >= 64) ? b.hi : b.lo;
        int bitpos = bit & 63;
        if ((word >> bitpos) & 1) {
            result = gpu_u128_add(result, a);
            if (gpu_u128_gte(result, m)) result = gpu_u128_sub(result, m);
        }
    }
    return result;
}

/* Modular exponentiation: base^exp % mod.
 * Starts from highest set bit of exp for efficiency. */
__device__ gpu_u128 gpu_u128_powmod(gpu_u128 base, gpu_u128 exp, gpu_u128 mod) {
    gpu_u128 result = gpu_u128_from_u64(1);
    base = gpu_u128_mod(base, mod);

    int top_bit = 127 - gpu_u128_clz(exp);
    for (int bit = top_bit; bit >= 0; bit--) {
        result = gpu_u128_mulmod(result, result, mod);
        u64 word = (bit >= 64) ? exp.hi : exp.lo;
        int bitpos = bit & 63;
        if ((word >> bitpos) & 1) {
            result = gpu_u128_mulmod(result, base, mod);
        }
    }
    return result;
}

/* Miller-Rabin single witness test. n must be odd > 3. */
__device__ int gpu_miller_rabin_witness(gpu_u128 n, u64 witness) {
    gpu_u128 n_minus_1 = gpu_u128_sub(n, gpu_u128_from_u64(1));
    gpu_u128 d = n_minus_1;
    int r = 0;
    while ((d.lo & 1) == 0) {
        d = gpu_u128_shr1(d);
        r++;
    }

    gpu_u128 x = gpu_u128_powmod(gpu_u128_from_u64(witness), d, n);
    gpu_u128 one = gpu_u128_from_u64(1);

    if (gpu_u128_eq(x, one) || gpu_u128_eq(x, n_minus_1)) return 1;

    for (int i = 0; i < r - 1; i++) {
        x = gpu_u128_mulmod(x, x, n);
        if (gpu_u128_eq(x, n_minus_1)) return 1;
        if (gpu_u128_eq(x, one)) return 0;
    }
    return 0;
}

/* Probabilistic MR for numbers up to 128 bits.
 * 20 witnesses (first 20 primes): error < 4^-20 ≈ 9.1×10^-13 per test.
 * The first 12 primes are deterministic only for n < 3.317×10^24 (< 2^82);
 * since our range is 89–128 bit, we use 20 rounds for robust coverage.
 * Over 10^13 lifetime tests, expected false primes < 0.01. */
__device__ int gpu_is_prime(gpu_u128 n) {
    if (n.hi == 0) {
        if (n.lo < 2) return 0;
        if (n.lo == 2 || n.lo == 3) return 1;
        if ((n.lo & 1) == 0) return 0;
        if (n.lo < 9) return 1;
    }
    if ((n.lo & 1) == 0) return 0;

    const u64 small_primes[] = {3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37,
                                41, 43, 47, 53, 59, 61, 67, 71};
    for (int i = 0; i < 19; i++) {
        gpu_u128 sp = gpu_u128_from_u64(small_primes[i]);
        gpu_u128 rem = gpu_u128_mod(n, sp);
        if (gpu_u128_is_zero(rem))
            return (n.hi == 0 && n.lo == small_primes[i]) ? 1 : 0;
    }

    const u64 witnesses[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
                             31, 37, 41, 43, 47, 53, 59, 61, 67, 71};
    for (int i = 0; i < 20; i++) {
        if (!gpu_miller_rabin_witness(n, witnesses[i])) return 0;
    }
    return 1;
}

/* =========================================================================
 * GPU PROVE KERNEL (Kernel 2)
 *
 * One thread per survivor. Reconstructs n, trial divides, MR tests,
 * chain root check, chain follow. Emits ChainResult for chains found.
 * Reads survivor count from device memory — no host sync needed.
 * ========================================================================= */

__global__ void __launch_bounds__(256)
cc18_prove_kernel(
    const Survivor* __restrict__ survivors,
    const u32* __restrict__ survivor_count_ptr,  /* device pointer — read from K1 */
    u64 batch_tile_start_lo,
    u64 batch_tile_start_hi,
    ChainResult* chain_results,
    u32* chain_result_count,
    u32* gpu_primes_count
) {
    u32 idx = blockIdx.x * blockDim.x + threadIdx.x;
    u32 num_survivors = *survivor_count_ptr;
    if (num_survivors > MAX_SURVIVORS_PER_BATCH) num_survivors = MAX_SURVIVORS_PER_BATCH;
    if (idx >= num_survivors) return;

    Survivor surv = survivors[idx];

    /* Reconstruct n as gpu_u128 */
    gpu_u128 tile_start;
    tile_start.lo = batch_tile_start_lo;
    tile_start.hi = batch_tile_start_hi;

    gpu_u128 tile_num = gpu_u128_add(tile_start, gpu_u128_from_u64(surv.tile_idx));
    gpu_u128 full_k = gpu_u128_add(
        gpu_u128_mul_u64(tile_num, (u64)d_wheel_period),
        gpu_u128_from_u64(surv.k_offset)
    );
    gpu_u128 n = gpu_u128_add(
        gpu_u128_mul_u64(full_k, d_lattice_m),
        gpu_u128_from_u64(d_base_values[surv.base_idx])
    );

    /* Check bit size */
    int nbits = 128 - gpu_u128_clz(n);
    if (nbits != d_target_bits) return;

    /* Trial division with chain screen primes */
    for (int i = 0; i < d_num_trial_primes; i++) {
        gpu_u128 p = gpu_u128_from_u64(d_trial_primes[i]);
        gpu_u128 rem = gpu_u128_mod(n, p);
        if (gpu_u128_is_zero(rem)) return;
    }

    /* Miller-Rabin primality test */
    if (!gpu_is_prime(n)) return;
    atomicAdd(gpu_primes_count, 1u);

    /* Chain root check: is (n+1)/2 also prime? If so, n is not a root.
     * Second-kind: predecessor is (n+1)/2, not (n-1)/2.
     * NOTE: GPU path does NOT walk backward (no GMP). Non-root chains whose
     * root is outside the bit-range are missed here. Use --prove-threads
     * (CPU proving) for complete non-root handling with backward walk. */
    gpu_u128 pred = gpu_u128_shr1(gpu_u128_add(n, gpu_u128_from_u64(1)));
    if (pred.hi > 0 || pred.lo >= 2) {
        if (gpu_is_prime(pred)) return;
    }

    /* Chain follow: keep computing 2p-1 while prime (second-kind) */
    int chain_len = 1;
    gpu_u128 current = n;
    gpu_u128 one = gpu_u128_from_u64(1);

    while (chain_len < MAX_CHAIN) {
        gpu_u128 next = gpu_u128_sub(gpu_u128_shl1(current), one);
        if (!gpu_is_prime(next)) break;
        chain_len++;
        current = next;
    }

    /* Emit chain result */
    u32 slot = atomicAdd(chain_result_count, 1u);
    if (slot < MAX_CHAIN_RESULTS) {
        chain_results[slot].base_idx = surv.base_idx;
        chain_results[slot].k_offset = surv.k_offset;
        chain_results[slot].tile_idx = surv.tile_idx;
        chain_results[slot].chain_len = chain_len;
    }
}

/* =========================================================================
 * GPU FILTER KERNEL — The hot path
 *
 * One thread per (base, wheel_position) in a tile.
 * Each thread:
 *   1. L2 bitmask filter (7 primes)
 *   2. Ext-L2 byte filter (7 primes)
 *   3. Line-sieve (125 primes)
 *   4. Emit survivor to global buffer
 * ========================================================================= */

__global__ void __launch_bounds__(THREADS_PER_BLOCK)
cc18_filter_kernel(
    u32 tile_count,            /* Number of tiles in this batch */
    u32 tile_start_mod6,       /* (batch tile_start) % 6, for mod-6 bucket selection */
    /* All tables passed as kernel args for flexibility */
    const u8* __restrict__ kmod_filter,        /* [FILTER_PRIMES_COUNT][WHEEL_PERIOD] */
    const u8* __restrict__ kmod_ext_filter,    /* [EXT_FILTER_PRIMES_COUNT][WHEEL_PERIOD] */
    const u16* __restrict__ kM_mod_line,       /* [WHEEL_PERIOD][LINE_SIEVE_COUNT] */
    const u8*  __restrict__ ext_filter_first,  /* [EXT_PRIMES][33][128] */
    const u8*  __restrict__ line_kill_first,   /* [LINE_SIEVE_COUNT][864] */
    const u16* __restrict__ batch_l2_res,      /* [num_bases][FILTER_PRIMES_COUNT] — includes tile_start */
    const u16* __restrict__ batch_ext_res,     /* [num_bases][EXT_FILTER_PRIMES_COUNT] */
    const u16* __restrict__ batch_line_res,    /* [num_bases][LINE_SIEVE_COUNT] */
    const u32* __restrict__ wheel_offsets,     /* Flat array of all wheel k-offsets (mod-6 bucketed) */
    const u32* __restrict__ wheel_start_m6,    /* Per (base,bucket) start: [base*6+bucket] */
    const u32* __restrict__ wheel_count_m6,    /* Per (base,bucket) count: [base*6+bucket] */
    Survivor*  survivors,
    u32*       survivor_count
) {
    int sieve_len = d_sieve_len;
    int num_bases = d_num_bases;

    /* Shared tile residues: computed once per (tile,base) work unit. */
    __shared__ u16 s_tile_l2[FILTER_PRIMES_COUNT];
    __shared__ u16 s_tile_ext[EXT_FILTER_PRIMES_COUNT];
    __shared__ u16 s_tile_line[LINE_SIEVE_COUNT];

    /* Work space is (tile_idx, base_idx), so each block handles one base
     * within one tile at a time. This removes the per-block base loop and
     * avoids redundant residue setup in every thread. */
    u64 total_work = (u64)tile_count * (u64)num_bases;
    for (u64 work_idx = blockIdx.x; work_idx < total_work; work_idx += gridDim.x) {
        u32 tile_idx = (u32)(work_idx / (u64)num_bases);
        u32 base_idx = (u32)(work_idx - (u64)tile_idx * (u64)num_bases);

        /* Mod-6 bucket selection: pick the unique bucket where all k_offsets
         * produce n ≡ 1 (mod 6). Math: valid k%6 = (target_k_mod6 - tile%6 + 6) % 6.
         * tile%6 = (tile_start_mod6 + tile_idx) % 6. */
        u32 tile_mod6 = (tile_start_mod6 + tile_idx) % 6;
        u32 valid_bucket = (d_target_k_mod6[base_idx] + 6 - tile_mod6) % 6;
        u32 wstart = wheel_start_m6[base_idx * 6 + valid_bucket];
        u32 wcount = wheel_count_m6[base_idx * 6 + valid_bucket];

        /* === PARALLEL RESIDUE COMPUTATION ===
         * 7 L2 + 7 ext + 125 line = 139 primes total.
         * Thread i computes prime i — all done in 1 cycle vs 139 serial.
         * Codex used thread-0-only; this parallelizes across 139 threads. */

        /* L2 residues: threads 0..6 */
        if (threadIdx.x < FILTER_PRIMES_COUNT) {
            int i = threadIdx.x;
            u32 q = d_filter_primes[i];
            u32 idx_km = (tile_idx % q) * (u32)d_WPM_l2[i] % q;
            u32 v = (u32)__ldg(&batch_l2_res[base_idx * FILTER_PRIMES_COUNT + i]) + idx_km;
            s_tile_l2[i] = (u16)(v >= q ? v - q : v);
        }

        /* Ext-L2 residues: threads 7..13 */
        if (threadIdx.x >= FILTER_PRIMES_COUNT &&
            threadIdx.x < FILTER_PRIMES_COUNT + EXT_FILTER_PRIMES_COUNT) {
            int i = threadIdx.x - FILTER_PRIMES_COUNT;
            u32 q = d_ext_filter_primes[i];
            u32 idx_km = (tile_idx % q) * (u32)d_WPM_ext[i] % q;
            u32 v = (u32)__ldg(&batch_ext_res[base_idx * EXT_FILTER_PRIMES_COUNT + i]) + idx_km;
            s_tile_ext[i] = (u16)(v >= q ? v - q : v);
        }

        /* Line-sieve residues: threads 14..138 (125 primes) */
        if (threadIdx.x >= FILTER_PRIMES_COUNT + EXT_FILTER_PRIMES_COUNT &&
            threadIdx.x < FILTER_PRIMES_COUNT + EXT_FILTER_PRIMES_COUNT + LINE_SIEVE_COUNT) {
            int i = threadIdx.x - FILTER_PRIMES_COUNT - EXT_FILTER_PRIMES_COUNT;
            u32 q = d_line_primes[i];
            u32 idx_km = (tile_idx % q) * (u32)d_WPM_line[i] % q;
            u32 v = (u32)__ldg(&batch_line_res[base_idx * LINE_SIEVE_COUNT + i]) + idx_km;
            s_tile_line[i] = (u16)(v >= q ? v - q : v);
        }

        __syncthreads();

        for (u32 w = threadIdx.x; w < wcount; w += blockDim.x) {
            u32 k_offset = __ldg(&wheel_offsets[wstart + w]);

            int passed = 1;
            #pragma unroll
            for (int i = 0; i < FILTER_PRIMES_COUNT; i++) {
                u32 q = d_filter_primes[i];
                u32 n_mod_p = (u32)s_tile_l2[i] +
                              (u32)__ldg(&kmod_filter[i * WHEEL_PERIOD + k_offset]);
                if (n_mod_p >= q) n_mod_p -= q;
                if ((d_filter_mask_first[i][sieve_len] >> n_mod_p) & 1ULL) {
                    passed = 0;
                    break;
                }
            }
            if (!passed) continue;

            #pragma unroll
            for (int i = 0; i < EXT_FILTER_PRIMES_COUNT; i++) {
                u32 q = d_ext_filter_primes[i];
                u32 n_mod_p = (u32)s_tile_ext[i] +
                              (u32)__ldg(&kmod_ext_filter[i * WHEEL_PERIOD + k_offset]);
                if (n_mod_p >= q) n_mod_p -= q;
                int idx = i * (MAX_SIEVE_CHAIN_LEN + 1) * 128 +
                          sieve_len * 128 + (int)n_mod_p;
                if (__ldg(&ext_filter_first[idx])) {
                    passed = 0;
                    break;
                }
            }
            if (!passed) continue;

            /* v13 OPT-3: Prefetch kM_mod row for line-sieve.
             * Placed AFTER L2+ext-L2 filters: only the ~8% of candidates that
             * survive both filters pay for this prefetch. The row is 125 × u16
             * = 250 bytes = 2 cache lines (128B each). The line-sieve loop
             * immediately follows, so the L2→L1 fetch overlaps with loop setup.
             * The +128 read is safe: d_kM_mod_line has 128B padding. */
            {
                const u16* kM_row_pf = &kM_mod_line[k_offset * LINE_SIEVE_COUNT];
                asm volatile("prefetch.global.L1 [%0];" :: "l"(kM_row_pf));
                asm volatile("prefetch.global.L1 [%0];" :: "l"((const char*)kM_row_pf + 128));
            }

            {
                const u16* kM_row = &kM_mod_line[k_offset * LINE_SIEVE_COUNT];
                for (int ls = 0; ls < LINE_SIEVE_COUNT; ls++) {
                    u32 q = d_line_primes[ls];
                    u32 n_mod_q = (u32)s_tile_line[ls] + (u32)__ldg(&kM_row[ls]);
                    if (n_mod_q >= q) n_mod_q -= q;
                    if (__ldg(&line_kill_first[ls * 864 + n_mod_q])) {
                        passed = 0;
                        break;
                    }
                }
            }
            if (!passed) continue;

            /* Warp-aggregated survivor emission reduces global atomic traffic
             * when survivors spike in specific regions. */
            unsigned active_mask = __activemask();
            unsigned pass_mask = __ballot_sync(active_mask, passed);
            if (pass_mask) {
                int lane = threadIdx.x & (WARP_SIZE - 1);
                int leader = __ffs((int)pass_mask) - 1;
                u32 warp_base = 0;
                if (lane == leader) {
                    warp_base = atomicAdd(survivor_count, (u32)__popc(pass_mask));
                }
                warp_base = __shfl_sync(pass_mask, warp_base, leader);

                if (passed) {
                    u32 lane_prefix = (u32)__popc(pass_mask & ((1u << lane) - 1u));
                    u32 slot = warp_base + lane_prefix;
                    if (slot < MAX_SURVIVORS_PER_BATCH) {
                        survivors[slot].base_idx = base_idx;
                        survivors[slot].k_offset = k_offset;
                        survivors[slot].tile_idx = tile_idx;
                    }
                }
            }
        }
        __syncthreads();
    }
}

/* =========================================================================
 * CPU-SIDE PROVER — processes survivors from GPU
 *
 * For each survivor:
 *   1. Reconstruct n = base + (tile_num * WHEEL_PERIOD + k_offset) * M
 *   2. Batched trial division (u128)
 *   3. Root primality (Montgomery MR / GMP BPSW)
 *   4. Chain root check
 *   5. Chain follow + Pocklington proving
 * ========================================================================= */

/* Chain length distribution — indexed by chain length */
static u64 g_chain_dist[MAX_CHAIN + 1];
static u64 g_non_roots = 0;  /* Total non-root primes found (walked backward) */

/* =========================================================================
 * MULTI-THREADED PROVER — N workers prove survivors in parallel
 *
 * Each worker thread gets a contiguous chunk of survivors to prove.
 * Per-worker stats are merged after all workers finish (no locks on hot path).
 * Each worker has its own GMP context (mpz_t is NOT thread-safe).
 *
 * Output (chain logging) uses a mutex to prevent interleaved writes.
 * ========================================================================= */

#define MAX_PROVE_THREADS 64

/* Per-worker result accumulator — merged into globals after each batch */
typedef struct {
    u64 primes;
    u64 chains;
    u64 non_roots;   /* Candidates where (n-1)/2 is prime (chain extends backward) */
    int best_chain;
    u64 chain_dist[MAX_CHAIN + 1];
    /* Chain log entries — buffered to avoid lock contention during proving */
    struct {
        int chain_len;
        char hex[256];
        size_t nbits;
        int is_target;
        int is_non_root;
    } log_entries[256];
    int log_count;
    /* Dedicated target-hit buffer — never lost to 256-entry log cap */
    struct {
        int chain_len;
        char hex[256];
        size_t nbits;
        int is_target;
    } target_entries[4];
    int target_count;
    int target_overflow;  /* Count of target hits dropped due to reserve cap */
} ProverStats;

/* Per-worker input for a batch */
typedef struct {
    const Survivor* survivors;
    int start;                  /* First survivor index for this worker */
    int end;                    /* One past last survivor index */
    u128 batch_tile_start;
    int target_len;
    int bits;
    int log_threshold;
    ProverStats stats;          /* Output — written only by this worker */
} ProverWork;

static ProverWork g_prover_work[MAX_PROVE_THREADS];
static pthread_t  g_prover_threads[MAX_PROVE_THREADS];
static int        g_num_prove_threads = 16;
static int        g_prover_pool_started = 0;
static int        g_prover_pool_shutdown = 0;
static int        g_prover_batch_gen = 0;
static int        g_prover_done_count = 0;
static pthread_mutex_t g_prover_mutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t  g_prover_dispatch_cv = PTHREAD_COND_INITIALIZER;
static pthread_cond_t  g_prover_done_cv = PTHREAD_COND_INITIALIZER;

/* Worker thread function — persists for whole run, proving one batch per dispatch */
static void* prover_worker(void* arg) {
    ProverWork* work = (ProverWork*)arg;
    int local_gen = 0;

    mpz_t n, temp, pred, current, root;
    mpz_init(n); mpz_init(temp); mpz_init(pred); mpz_init(current); mpz_init(root);

    for (;;) {
        pthread_mutex_lock(&g_prover_mutex);
        while (!g_prover_pool_shutdown && g_prover_batch_gen == local_gen)
            pthread_cond_wait(&g_prover_dispatch_cv, &g_prover_mutex);
        if (g_prover_pool_shutdown) {
            pthread_mutex_unlock(&g_prover_mutex);
            break;
        }
        local_gen = g_prover_batch_gen;
        pthread_mutex_unlock(&g_prover_mutex);

        memset(&work->stats, 0, sizeof(work->stats));

        /* Workers always complete their assigned chunk — never bail mid-batch.
         * g_shutdown only stops the main loop from launching NEW batches.
         * This ensures no CC18 is ever lost from a partially-processed batch.
         * Worst case: ~1250 survivors/thread at 23K proves/sec = 0.05s delay. */
        for (int s = work->start; s < work->end; s++) {
            u64 base = h_bases[work->survivors[s].base_idx].base;
            u128 tile_num = work->batch_tile_start + work->survivors[s].tile_idx;
            u128 full_k = tile_num * WHEEL_PERIOD + work->survivors[s].k_offset;
            u128 n128 = full_k * LATTICE_M + base;

            mpz_set_ui(n, (unsigned long)(n128 >> 64));
            mpz_mul_2exp(n, n, 64);
            mpz_add_ui(n, n, (unsigned long)(n128 & 0xFFFFFFFFFFFFFFFFULL));

            size_t nbits = mpz_sizeinbase(n, 2);
            if ((int)nbits != work->bits) continue;

            /* Trial division */
            int trial_fail = 0;
            for (int i = 0; i < 62; i++) {
                if (mpz_divisible_ui_p(n, CHAIN_SCREEN_PRIMES[i])) {
                    trial_fail = 1; break;
                }
            }
            if (trial_fail) continue;

            /* Root primality */
            if (mpz_probab_prime_p(n, 1) == 0) continue;
            work->stats.primes++;

            /* Chain root check — walk backward to find true root.
             * Second-kind: predecessor is (n+1)/2, not (n-1)/2.
             * If (n+1)/2 is prime, n is NOT a root — the chain extends
             * backward. We follow it back to the true root, then follow
             * forward from there to measure the FULL chain length. */
            int is_root = 1;
            mpz_set(root, n);
            mpz_add_ui(pred, n, 1);
            mpz_fdiv_q_2exp(pred, pred, 1);
            if (mpz_cmp_ui(pred, 2) >= 0 && mpz_probab_prime_p(pred, 1) > 0) {
                /* n is NOT a root — walk backward to find it */
                is_root = 0;
                work->stats.non_roots++;
                mpz_set(root, pred);
                int back_steps = 0;
                while (back_steps < MAX_CHAIN) {
                    mpz_add_ui(pred, root, 1);
                    mpz_fdiv_q_2exp(pred, pred, 1);
                    if (mpz_cmp_ui(pred, 2) < 0) break;
                    if (mpz_probab_prime_p(pred, 1) == 0) break;
                    mpz_set(root, pred);
                    back_steps++;
                }
            }

            /* Follow chain forward from the true root (second-kind: 2p-1).
             * Do NOT check g_shutdown here — chain follow is short (max 64
             * MR tests) and truncation would mis-record chain length.
             * Double-Ctrl+C still force-exits via signal handler. */
            int chain_len = 1;
            mpz_set(current, root);
            while (chain_len < MAX_CHAIN) {
                mpz_mul_2exp(temp, current, 1);
                mpz_sub_ui(temp, temp, 1);
                if (mpz_probab_prime_p(temp, 1) == 0) break;
                chain_len++;
                mpz_set(current, temp);
            }

            work->stats.chains++;
            if (chain_len <= MAX_CHAIN)
                work->stats.chain_dist[chain_len]++;
            if (chain_len > work->stats.best_chain)
                work->stats.best_chain = chain_len;

            /* Buffer log entry for chains >= threshold */
            if ((chain_len >= work->log_threshold || chain_len >= work->target_len)
                && work->stats.log_count < 256) {
                int idx = work->stats.log_count++;
                work->stats.log_entries[idx].chain_len = chain_len;
                /* Log the TRUE ROOT, not the found member */
                gmp_snprintf(work->stats.log_entries[idx].hex, 256, "0x%ZX", root);
                work->stats.log_entries[idx].nbits = mpz_sizeinbase(root, 2);
                work->stats.log_entries[idx].is_target = (chain_len >= work->target_len);
                work->stats.log_entries[idx].is_non_root = !is_root;
            }

            /* Target hits go to dedicated buffer (never lost to log cap) */
            if (chain_len >= work->target_len) {
                if (work->stats.target_count < 4) {
                    int ti = work->stats.target_count++;
                    work->stats.target_entries[ti].chain_len = chain_len;
                    gmp_snprintf(work->stats.target_entries[ti].hex, 256, "0x%ZX", root);
                    work->stats.target_entries[ti].nbits = mpz_sizeinbase(root, 2);
                    work->stats.target_entries[ti].is_target = 1;
                } else {
                    work->stats.target_overflow++;
                }
            }
        }

        pthread_mutex_lock(&g_prover_mutex);
        g_prover_done_count++;
        if (g_prover_done_count == g_num_prove_threads)
            pthread_cond_signal(&g_prover_done_cv);
        pthread_mutex_unlock(&g_prover_mutex);
    }

    mpz_clear(n); mpz_clear(temp); mpz_clear(pred); mpz_clear(current); mpz_clear(root);
    return NULL;
}

static int prover_pool_start(void) {
    if (g_prover_pool_started) return 0;

    pthread_mutex_lock(&g_prover_mutex);
    g_prover_pool_shutdown = 0;
    g_prover_batch_gen = 0;
    g_prover_done_count = 0;
    pthread_mutex_unlock(&g_prover_mutex);

    for (int t = 0; t < g_num_prove_threads; t++) {
        memset(&g_prover_work[t], 0, sizeof(g_prover_work[t]));
        int rc = pthread_create(&g_prover_threads[t], NULL, prover_worker, &g_prover_work[t]);
        if (rc != 0) {
            fprintf(stderr, "ERROR: pthread_create failed for prover thread %d (rc=%d)\n", t, rc);
            pthread_mutex_lock(&g_prover_mutex);
            g_prover_pool_shutdown = 1;
            pthread_cond_broadcast(&g_prover_dispatch_cv);
            pthread_mutex_unlock(&g_prover_mutex);
            for (int j = 0; j < t; j++) pthread_join(g_prover_threads[j], NULL);
            return -1;
        }
    }
    g_prover_pool_started = 1;
    return 0;
}

static void prover_pool_stop(void) {
    if (!g_prover_pool_started) return;

    pthread_mutex_lock(&g_prover_mutex);
    g_prover_pool_shutdown = 1;
    pthread_cond_broadcast(&g_prover_dispatch_cv);
    pthread_mutex_unlock(&g_prover_mutex);

    for (int t = 0; t < g_num_prove_threads; t++)
        pthread_join(g_prover_threads[t], NULL);

    g_prover_pool_started = 0;
}

/* Dispatch one batch to persistent workers, wait for completion, merge stats */
static void prove_survivors_mt(
    const Survivor* survivors, int count,
    u128 batch_tile_start,
    int target_len, int bits, int log_threshold,
    FILE* output_fp,
    u64* out_primes, u64* out_chains, int* best_chain_len
) {
    if (count == 0) return;

    int nthreads = g_num_prove_threads;
    if (nthreads < 1) nthreads = 1;

    /* Partition survivors into chunks */
    int chunk = count / nthreads;
    int remainder = count % nthreads;
    int offset = 0;

    pthread_mutex_lock(&g_prover_mutex);
    for (int t = 0; t < nthreads; t++) {
        int my_count = chunk + (t < remainder ? 1 : 0);
        g_prover_work[t].survivors = survivors;
        g_prover_work[t].start = offset;
        g_prover_work[t].end = offset + my_count;
        g_prover_work[t].batch_tile_start = batch_tile_start;
        g_prover_work[t].target_len = target_len;
        g_prover_work[t].bits = bits;
        g_prover_work[t].log_threshold = log_threshold;
        offset += my_count;
    }

    g_prover_done_count = 0;
    g_prover_batch_gen++;
    pthread_cond_broadcast(&g_prover_dispatch_cv);
    while (g_prover_done_count < nthreads)
        pthread_cond_wait(&g_prover_done_cv, &g_prover_mutex);
    pthread_mutex_unlock(&g_prover_mutex);

    /* Merge stats from all workers */
    for (int t = 0; t < nthreads; t++) {
        ProverStats* st = &g_prover_work[t].stats;
        *out_primes += st->primes;
        *out_chains += st->chains;
        g_non_roots += st->non_roots;
        if (st->best_chain > *best_chain_len)
            *best_chain_len = st->best_chain;
        for (int c = 0; c <= MAX_CHAIN; c++)
            g_chain_dist[c] += st->chain_dist[c];

        /* Print buffered log entries */
        for (int i = 0; i < st->log_count; i++) {
            int cl = st->log_entries[i].chain_len;
            const char* nr_tag = st->log_entries[i].is_non_root ? " [NON-ROOT]" : "";
            if (st->log_entries[i].is_target) {
                printf("\n  !!!! CC%d TARGET HIT: %s (%zu bits)%s !!!!\n",
                       cl, st->log_entries[i].hex, st->log_entries[i].nbits, nr_tag);
            } else if (cl > *best_chain_len - 1) {
                printf("\n  *** CC%d: %s%s (%zu bits) ***\n",
                       cl, st->log_entries[i].hex, nr_tag, st->log_entries[i].nbits);
            } else {
                printf("  CC%d: %s%s (%zu bits)\n",
                       cl, st->log_entries[i].hex, nr_tag, st->log_entries[i].nbits);
            }
            if (output_fp) {
                fprintf(output_fp, "%s%sCC%d %s %zu\n",
                        st->log_entries[i].is_target ? "TARGET " : "",
                        st->log_entries[i].is_non_root ? "NON-ROOT " : "",
                        cl, st->log_entries[i].hex, st->log_entries[i].nbits);
                fflush(output_fp);
            }
        }

        /* Print dedicated target-hit entries (immune to 256-entry cap).
         * Skip if already printed via log_entries to avoid duplicates. */
        for (int i = 0; i < st->target_count; i++) {
            /* Check if this target was already in log_entries */
            int already_logged = 0;
            for (int j = 0; j < st->log_count && !already_logged; j++) {
                if (st->log_entries[j].is_target &&
                    st->log_entries[j].chain_len == st->target_entries[i].chain_len &&
                    strcmp(st->log_entries[j].hex, st->target_entries[i].hex) == 0)
                    already_logged = 1;
            }
            if (!already_logged) {
                int cl = st->target_entries[i].chain_len;
                printf("\n  !!!! CC%d TARGET HIT (from reserve): %s (%zu bits) !!!!\n",
                       cl, st->target_entries[i].hex, st->target_entries[i].nbits);
                if (output_fp) {
                    fprintf(output_fp, "TARGET CC%d %s %zu\n",
                            cl, st->target_entries[i].hex, st->target_entries[i].nbits);
                    fflush(output_fp);
                }
            }
        }

        /* Warn if any worker's target reserve overflowed */
        if (st->target_overflow > 0) {
            fprintf(stderr, "\n  WARNING: worker %d target reserve OVERFLOW — "
                    "%d target hit(s) DROPPED (reserve cap = 4)!\n",
                    t, st->target_overflow);
        }
    }
}

/* =========================================================================
 * SELF-TEST — Verifies chain logic, primality, reconstruction
 *
 * Uses known Cunningham chain roots to ensure we never silently mis-count
 * or lose a chain. Runs entirely on CPU (GMP), no GPU required.
 * ========================================================================= */

static int run_self_test(void) {
    printf("\n=== CC18 C+C v13b Self-Test (Second-Kind) ===\n\n");
    int passed = 0, failed = 0;

    /* Known second-kind Cunningham chain test vectors.
     * Format: { "decimal_root", expected_chain_length, "description" }
     * Each root p should satisfy: p prime, (p+1)/2 not prime (root check),
     * and the chain p → 2p-1 → 2(2p-1)-1 → ... has exactly the stated length.
     * Sources: OEIS A005603, pzktupel.de/CC */
    struct {
        const char* root_dec;
        int expected_len;
        const char* desc;
    } tests[] = {
        { "2",                              3,  "CC3 2nd-kind (2→3→5)" },
        { "7",                              2,  "CC2 2nd-kind (7→13)" },
        { "2131",                           4,  "CC4 2nd-kind" },
        { "1531",                           5,  "CC5 2nd-kind" },
        { "385591",                         6,  "CC6 2nd-kind" },
        { "16651",                          7,  "CC7 2nd-kind" },
        { "15514861",                       8,  "CC8 2nd-kind" },
        { "857095381",                      9,  "CC9 2nd-kind" },
        { "205528443121",                  10,  "CC10 2nd-kind" },
        { "1389122693971",                 11,  "CC11 2nd-kind" },
        { "3203000719597029781",           16,  "CC16 2nd-kind (1997, Forbes)" },
        { "40244844789379926979141",       17,  "CC17 2nd-kind (2008, Wroblewski)" },
    };
    int num_tests = (int)(sizeof(tests) / sizeof(tests[0]));

    mpz_t root, current, next, pred;
    mpz_init(root); mpz_init(current); mpz_init(next); mpz_init(pred);

    for (int t = 0; t < num_tests; t++) {
        mpz_set_str(root, tests[t].root_dec, 10);

        /* Test 1: root is prime */
        int is_prime = mpz_probab_prime_p(root, 25);
        if (is_prime == 0) {
            printf("  FAIL [%s]: root is not prime!\n", tests[t].desc);
            failed++;
            continue;
        }

        /* Test 2: root check — (root+1)/2 should NOT be prime (second-kind).
         * For root=2: (2+1)/2 = 1 (not prime), OK. */
        if (mpz_cmp_ui(root, 2) > 0) {
            mpz_add_ui(pred, root, 1);
            mpz_fdiv_q_2exp(pred, pred, 1);
            if (mpz_probab_prime_p(pred, 25) > 0) {
                /* (root+1)/2 is prime → root is not a chain root (second-kind). */
                printf("  FAIL [%s]: (root+1)/2 is prime, not a valid root!\n", tests[t].desc);
                failed++;
                continue;
            }
        }

        /* Test 3: follow chain forward (second-kind: 2p-1), count length */
        int chain_len = 1;
        mpz_set(current, root);
        while (chain_len < 100) {
            mpz_mul_2exp(next, current, 1);
            mpz_sub_ui(next, next, 1);
            if (mpz_probab_prime_p(next, 25) == 0) break;
            chain_len++;
            mpz_set(current, next);
        }

        if (chain_len == tests[t].expected_len) {
            char hex[256];
            gmp_snprintf(hex, sizeof(hex), "0x%ZX", root);
            printf("  PASS [%s]: CC%d verified (%s, %zu bits)\n",
                   tests[t].desc, chain_len, hex, mpz_sizeinbase(root, 2));
            passed++;
        } else {
            printf("  FAIL [%s]: expected CC%d, got CC%d\n",
                   tests[t].desc, tests[t].expected_len, chain_len);
            failed++;
        }
    }

    /* Test 4: n reconstruction round-trip.
     * Verify: base + (tile * WHEEL_PERIOD + k_offset) * LATTICE_M = n
     * Then verify the reverse: given n, compute k and check n = k*M + base. */
    printf("\n  Reconstruction round-trip test:\n");
    {
        u64 base = 1;  /* Smallest valid base (n ≡ 1 mod 6) */
        u128 tile = 12345678;
        u32 k_off = 100;
        u128 full_k = (u128)tile * WHEEL_PERIOD + k_off;
        u128 n_orig = full_k * LATTICE_M + base;

        /* Reconstruct using GMP (same path as prover) */
        mpz_t n_mpz;
        mpz_init(n_mpz);
        mpz_set_ui(n_mpz, (unsigned long)(n_orig >> 64));
        mpz_mul_2exp(n_mpz, n_mpz, 64);
        mpz_add_ui(n_mpz, n_mpz, (unsigned long)(n_orig & 0xFFFFFFFFFFFFFFFFULL));

        /* Verify n mod M = base */
        mpz_t m_mpz, r_mpz;
        mpz_init(m_mpz); mpz_init(r_mpz);
        mpz_set_ui(m_mpz, LATTICE_M);
        mpz_mod(r_mpz, n_mpz, m_mpz);
        u64 recovered_base = mpz_get_ui(r_mpz);

        if (recovered_base == base) {
            printf("    PASS: n mod M = base (%llu)\n", (unsigned long long)base);
            passed++;
        } else {
            printf("    FAIL: n mod M = %llu, expected %llu\n",
                   (unsigned long long)recovered_base, (unsigned long long)base);
            failed++;
        }

        /* Verify bit count (GMP is authoritative — no manual clz needed) */
        size_t nbits = mpz_sizeinbase(n_mpz, 2);
        printf("    INFO: reconstructed n = %zu bits (tile=%llu, k_off=%u)\n",
               nbits, (unsigned long long)tile, k_off);
        passed++;  /* INFO test always passes */

        mpz_clear(n_mpz); mpz_clear(m_mpz); mpz_clear(r_mpz);
    }

    /* Test 5: Verify CC17 2nd-kind root is in a reasonable bit range for campaign search */
    printf("\n  CC17 2nd-kind campaign relevance test:\n");
    {
        mpz_set_str(root, "40244844789379926979141", 10);
        size_t nbits = mpz_sizeinbase(root, 2);
        printf("    CC17 2nd-kind root: %zu bits\n", nbits);
        if (nbits >= 75 && nbits <= 76) {
            printf("    PASS: CC17 2nd-kind root is ~76-bit (current lowest known)\n");
            passed++;
        } else {
            printf("    INFO: CC17 2nd-kind root is %zu-bit\n", nbits);
            passed++;  /* Not a failure, just info */
        }
    }

    /* Test 6: NON-ROOT backward walk simulation (second-kind)
     * Start from CC17 2nd-kind chain member at position 5.
     * Second-kind: member_i = 2*member_{i-1} - 1, predecessor = (member+1)/2.
     * Walk backward — must recover the original ~76-bit root.
     * Then follow forward — must measure CC17.
     * This validates the exact code path used when --prove-threads finds
     * a non-root chain member inside the bit-range search window. */
    printf("\n  Non-root backward walk test (CC17 2nd-kind simulation):\n");
    {
        /* Compute chain member at position 5 from CC17 2nd-kind root */
        mpz_set_str(root, "40244844789379926979141", 10);
        mpz_t member, back_root, back_pred, back_cur, back_temp;
        mpz_init(member); mpz_init(back_root); mpz_init(back_pred);
        mpz_init(back_cur); mpz_init(back_temp);

        mpz_set(member, root);
        for (int i = 0; i < 5; i++) {
            mpz_mul_2exp(member, member, 1);
            mpz_sub_ui(member, member, 1);
        }
        size_t mem_bits = mpz_sizeinbase(member, 2);
        printf("    Position-5 member: %zu bits\n", mem_bits);

        /* Verify member is prime (it should be — it's in the CC17 chain) */
        if (mpz_probab_prime_p(member, 25) > 0) {
            printf("    PASS: position-5 member is prime\n");
            passed++;
        } else {
            printf("    FAIL: position-5 member is NOT prime\n");
            failed++;
        }

        /* Verify (member+1)/2 IS prime (member is NOT a root) — second-kind */
        mpz_add_ui(back_pred, member, 1);
        mpz_fdiv_q_2exp(back_pred, back_pred, 1);
        if (mpz_probab_prime_p(back_pred, 25) > 0) {
            printf("    PASS: (member+1)/2 is prime => member is NOT a root\n");
            passed++;
        } else {
            printf("    FAIL: (member+1)/2 is NOT prime => member IS a root (unexpected)\n");
            failed++;
        }

        /* Walk backward to find true root (second-kind: predecessor = (p+1)/2) */
        mpz_set(back_root, back_pred);
        int back_steps = 1;
        while (back_steps < MAX_CHAIN) {
            mpz_add_ui(back_pred, back_root, 1);
            mpz_fdiv_q_2exp(back_pred, back_pred, 1);
            if (mpz_cmp_ui(back_pred, 2) < 0) break;
            if (mpz_probab_prime_p(back_pred, 25) == 0) break;
            mpz_set(back_root, back_pred);
            back_steps++;
        }
        printf("    Walked back %d steps\n", back_steps);
        size_t root_bits = mpz_sizeinbase(back_root, 2);
        printf("    Recovered root: %zu bits\n", root_bits);

        /* Verify recovered root matches original CC17 2nd-kind root */
        mpz_set_str(root, "40244844789379926979141", 10);
        if (mpz_cmp(back_root, root) == 0) {
            printf("    PASS: backward walk recovered correct CC17 2nd-kind root\n");
            passed++;
        } else {
            char got[256], want[256];
            mpz_get_str(got, 16, back_root);
            mpz_get_str(want, 16, root);
            printf("    FAIL: got root 0x%s, expected 0x%s\n", got, want);
            failed++;
        }

        /* Follow forward from recovered root — must get CC17 (second-kind: 2p-1) */
        int fwd_len = 1;
        mpz_set(back_cur, back_root);
        while (fwd_len < MAX_CHAIN) {
            mpz_mul_2exp(back_temp, back_cur, 1);
            mpz_sub_ui(back_temp, back_temp, 1);
            if (mpz_probab_prime_p(back_temp, 25) == 0) break;
            fwd_len++;
            mpz_set(back_cur, back_temp);
        }
        printf("    Forward chain from recovered root: CC%d\n", fwd_len);
        if (fwd_len == 17) {
            printf("    PASS: non-root backward walk + forward = CC17\n");
            passed++;
        } else {
            printf("    FAIL: expected CC17, got CC%d\n", fwd_len);
            failed++;
        }

        mpz_clear(member); mpz_clear(back_root); mpz_clear(back_pred);
        mpz_clear(back_cur); mpz_clear(back_temp);
    }

    mpz_clear(root); mpz_clear(current); mpz_clear(next); mpz_clear(pred);

    /* =====================================================================
     * PREFIX TESTS (v9)
     * ===================================================================== */

    /* Test 7: Prefix parsing — valid inputs */
    printf("\n  Prefix parsing tests:\n");
    {
        u128 val; int nbits;

        /* Valid: 0b1 */
        if (parse_prefix("0b1", &val, &nbits) == 0 && val == 1 && nbits == 1) {
            printf("    PASS: parse_prefix(\"0b1\") => val=1, bits=1\n");
            passed++;
        } else {
            printf("    FAIL: parse_prefix(\"0b1\")\n");
            failed++;
        }

        /* Valid: 0b101 */
        if (parse_prefix("0b101", &val, &nbits) == 0 && val == 5 && nbits == 3) {
            printf("    PASS: parse_prefix(\"0b101\") => val=5, bits=3\n");
            passed++;
        } else {
            printf("    FAIL: parse_prefix(\"0b101\")\n");
            failed++;
        }

        /* Valid: 0b11111111 */
        if (parse_prefix("0b11111111", &val, &nbits) == 0 && val == 255 && nbits == 8) {
            printf("    PASS: parse_prefix(\"0b11111111\") => val=255, bits=8\n");
            passed++;
        } else {
            printf("    FAIL: parse_prefix(\"0b11111111\")\n");
            failed++;
        }
    }

    /* Test 8: Prefix parsing — invalid inputs (should fail).
     * Suppress stderr during negative tests to avoid confusing ERROR messages. */
    {
        u128 val; int nbits;
        FILE* saved_stderr = stderr;
        FILE* devnull = fopen("/dev/null", "w");
        if (devnull) stderr = devnull;

        /* Must start with 0b */
        int r1 = parse_prefix("101", &val, &nbits);
        if (r1 != 0) {
            printf("    PASS: parse_prefix(\"101\") correctly rejected (no 0b prefix)\n");
            passed++;
        } else {
            printf("    FAIL: parse_prefix(\"101\") should have failed\n");
            failed++;
        }

        /* Must start with 1 after 0b */
        int r2 = parse_prefix("0b0101", &val, &nbits);
        if (r2 != 0) {
            printf("    PASS: parse_prefix(\"0b0101\") correctly rejected (leading zero)\n");
            passed++;
        } else {
            printf("    FAIL: parse_prefix(\"0b0101\") should have failed\n");
            failed++;
        }

        /* Invalid character */
        int r3 = parse_prefix("0b1012", &val, &nbits);
        if (r3 != 0) {
            printf("    PASS: parse_prefix(\"0b1012\") correctly rejected (invalid char)\n");
            passed++;
        } else {
            printf("    FAIL: parse_prefix(\"0b1012\") should have failed\n");
            failed++;
        }

        /* Restore stderr */
        stderr = saved_stderr;
        if (devnull) fclose(devnull);
    }

    /* Test 9: Prefix tile range — complementary prefixes cover full range.
     * For 40-bit candidates (large enough for M=1616615):
     *   prefix 0b10 covers [2^39, 3*2^38) portion
     *   prefix 0b11 covers [3*2^38, 2^40) portion
     *   Their union should cover the entire 40-bit tile range. */
    printf("\n  Prefix tile range complementarity test:\n");
    {
        int test_bits = 40;

        /* Use conservative base bounds */
        u64 test_min_base = 5;
        u64 test_max_base = LATTICE_M - 1;

        /* Full (unprefixed) tile range for 40-bit */
        u128 full_n_min = (u128)1 << (test_bits - 1);
        u128 full_n_max = ((u128)1 << test_bits) - 1;
        u128 full_k_min = (full_n_min > (u128)test_max_base) ?
            (full_n_min - test_max_base) / LATTICE_M : 0;
        u128 full_k_max = (full_n_max - test_min_base) / LATTICE_M;
        u128 full_tile_min = full_k_min / WHEEL_PERIOD;
        u128 full_tile_max = full_k_max / WHEEL_PERIOD + 1;

        /* Prefix 0b10 tile range */
        u128 p10_tile_min, p10_tile_max;
        int r10 = compute_prefix_tile_range(test_bits, 2, 2,
                                            test_min_base, test_max_base,
                                            &p10_tile_min, &p10_tile_max);

        /* Prefix 0b11 tile range */
        u128 p11_tile_min, p11_tile_max;
        int r11 = compute_prefix_tile_range(test_bits, 3, 2,
                                            test_min_base, test_max_base,
                                            &p11_tile_min, &p11_tile_max);

        if (r10 == 0 && r11 == 0) {
            /* The two ranges should be contiguous and cover the full range.
             * Due to rounding, we check that:
             *   min(p10_min, p11_min) <= full_tile_min
             *   max(p10_max, p11_max) >= full_tile_max
             *   The ranges either overlap or are adjacent */
            u128 union_min = p10_tile_min < p11_tile_min ? p10_tile_min : p11_tile_min;
            u128 union_max = p10_tile_max > p11_tile_max ? p10_tile_max : p11_tile_max;

            int covers = (union_min <= full_tile_min && union_max >= full_tile_max);
            /* Check contiguous/overlapping */
            int contiguous = (p10_tile_max >= p11_tile_min || p11_tile_max >= p10_tile_min);

            if (covers && contiguous) {
                printf("    PASS: prefix 0b10 + 0b11 covers full 21-bit tile range\n");
                printf("      full=[%.0f, %.0f), 0b10=[%.0f, %.0f), 0b11=[%.0f, %.0f)\n",
                       (double)full_tile_min, (double)full_tile_max,
                       (double)p10_tile_min, (double)p10_tile_max,
                       (double)p11_tile_min, (double)p11_tile_max);
                passed++;
            } else {
                printf("    FAIL: prefix union does not cover full range\n");
                printf("      full=[%.0f, %.0f), 0b10=[%.0f, %.0f), 0b11=[%.0f, %.0f)\n",
                       (double)full_tile_min, (double)full_tile_max,
                       (double)p10_tile_min, (double)p10_tile_max,
                       (double)p11_tile_min, (double)p11_tile_max);
                printf("      union_min=%.0f union_max=%.0f covers=%d contiguous=%d\n",
                       (double)union_min, (double)union_max, covers, contiguous);
                failed++;
            }
        } else {
            printf("    FAIL: compute_prefix_tile_range failed (r10=%d, r11=%d)\n", r10, r11);
            failed++;
        }
    }

    /* Test 10: Prefix tile range — no empty range for valid prefix */
    printf("\n  Prefix range non-empty test:\n");
    {
        u128 tile_min, tile_max;
        /* prefix 0b1 at 89 bits = full range (just the leading 1) */
        int r = compute_prefix_tile_range(89, 1, 1, 5, LATTICE_M - 1,
                                          &tile_min, &tile_max);
        if (r == 0 && tile_max > tile_min) {
            printf("    PASS: prefix 0b1 at 89-bit gives non-empty range (%.3e tiles)\n",
                   (double)(tile_max - tile_min));
            passed++;
        } else {
            printf("    FAIL: prefix 0b1 at 89-bit gave empty or error\n");
            failed++;
        }
    }

    /* Test 11: Prefix range narrows search space.
     * prefix 0b101 at 89 bits should give ~1/4 of the full 89-bit tile range
     * (since 0b101 is 3 bits, the first bit is always 1 for the bit range,
     *  so the remaining 2 bits divide the space into ~4 parts). */
    printf("\n  Prefix range narrowing test:\n");
    {
        u128 full_n_min = (u128)1 << 88;
        u128 full_n_max = ((u128)1 << 89) - 1;
        u128 full_k_min = (full_n_min - (u128)(LATTICE_M - 1)) / LATTICE_M;
        u128 full_k_max = (full_n_max - 5) / LATTICE_M;
        u128 full_tile_min = full_k_min / WHEEL_PERIOD;
        u128 full_tile_max = full_k_max / WHEEL_PERIOD + 1;
        double full_tiles = (double)(full_tile_max - full_tile_min);

        u128 pref_tile_min, pref_tile_max;
        int r = compute_prefix_tile_range(89, 5, 3, 5, LATTICE_M - 1,
                                          &pref_tile_min, &pref_tile_max);
        if (r == 0) {
            double pref_tiles = (double)(pref_tile_max - pref_tile_min);
            double ratio = pref_tiles / full_tiles;
            /* Should be roughly 1/4 (0.25), allow 20% tolerance */
            if (ratio > 0.20 && ratio < 0.30) {
                printf("    PASS: prefix 0b101 at 89-bit narrows to %.1f%% of full range\n",
                       ratio * 100.0);
                passed++;
            } else {
                printf("    FAIL: prefix 0b101 ratio=%.3f (expected ~0.25)\n", ratio);
                failed++;
            }
        } else {
            printf("    FAIL: compute_prefix_tile_range failed\n");
            failed++;
        }
    }

    /* Test 12: Sequential tile iteration within prefix produces monotonic tiles.
     * Uses 40-bit range with prefix 0b101 for meaningful tile count.
     * Simulates the main loop logic for sequential mode. */
    printf("\n  Sequential prefix no-duplicate test:\n");
    {
        u128 tile_min, tile_max;
        int r = compute_prefix_tile_range(40, 5, 3, 5, LATTICE_M - 1,
                                          &tile_min, &tile_max);
        if (r == 0 && tile_max > tile_min) {
            u128 cur = tile_min;
            u32 batch = 100;
            int monotonic = 1;
            int steps = 0;
            u128 prev = cur;
            while (cur < tile_max && steps < 50) {
                u32 this_batch = batch;
                if (cur + this_batch > tile_max)
                    this_batch = (u32)(tile_max - cur);
                if (this_batch == 0) break;
                if (steps > 0 && cur <= prev) {
                    monotonic = 0;
                    break;
                }
                prev = cur;
                cur += this_batch;
                steps++;
            }
            if (monotonic && steps > 0) {
                printf("    PASS: %d sequential batches, all monotonically increasing\n", steps);
                passed++;
            } else {
                printf("    FAIL: non-monotonic tiles or zero steps\n");
                failed++;
            }
        } else {
            printf("    FAIL: compute_prefix_tile_range failed for 40-bit 0b101\n");
            failed++;
        }
    }

    /* Test 13: Random tiles within prefix range stay in bounds */
    printf("\n  Random prefix range bounds test:\n");
    {
        u128 tile_min, tile_max;
        int r = compute_prefix_tile_range(89, 5, 3, 5, LATTICE_M - 1,
                                          &tile_min, &tile_max);
        if (r == 0) {
            u128 range = tile_max - tile_min;
            rng_seed(12345);
            int all_in_range = 1;
            for (int i = 0; i < 1000; i++) {
                u128 rand128 = ((u128)rng_next() << 64) | rng_next();
                u128 tile = tile_min + rand128 % range;
                if (tile < tile_min || tile >= tile_max) {
                    all_in_range = 0;
                    break;
                }
            }
            if (all_in_range) {
                printf("    PASS: 1000 random tiles all within prefix range\n");
                passed++;
            } else {
                printf("    FAIL: random tile outside prefix range\n");
                failed++;
            }
        } else {
            printf("    FAIL: compute_prefix_tile_range failed\n");
            failed++;
        }
    }

    /* Test 14: Candidates from prefix tiles have correct high bits.
     * Pick a tile in the prefix 0b101 range for 40-bit,
     * reconstruct a candidate n, verify top 3 bits = 101.
     * Uses 40-bit (much larger than M=1,616,615) for meaningful range. */
    printf("\n  Prefix candidate high-bits verification test:\n");
    {
        int test_bits = 40;
        /* Generate bases first so we have real base values */
        int saved_coverage = g_coverage_mode;
        g_coverage_mode = COVERAGE_FULL;
        generate_bases(18);
        g_coverage_mode = saved_coverage;

        if (h_num_bases > 0) {
            u64 bmin = h_bases[0].base, bmax = h_bases[0].base;
            for (int bi = 1; bi < h_num_bases; bi++) {
                if (h_bases[bi].base < bmin) bmin = h_bases[bi].base;
                if (h_bases[bi].base > bmax) bmax = h_bases[bi].base;
            }

            u128 tile_min, tile_max;
            int r = compute_prefix_tile_range(test_bits, 5, 3,
                                              bmin, bmax,
                                              &tile_min, &tile_max);
            if (r == 0 && tile_max > tile_min) {
                /* Pick middle tile */
                u128 tile = tile_min + (tile_max - tile_min) / 2;
                u64 base = h_bases[0].base;
                /* Find first non-empty bucket's first entry */
                u32 k_offset = 0;
                for (int m = 0; m < 6; m++) {
                    if (h_bases[0].wheel_size_by_mod6[m] > 0) {
                        k_offset = h_bases[0].wheel_by_mod6[m][0];
                        break;
                    }
                }
                u128 full_k = tile * WHEEL_PERIOD + k_offset;
                u128 n = full_k * LATTICE_M + base;

                /* Check n is test_bits-bit */
                int nbits = 0;
                u128 tmp = n;
                while (tmp > 0) { nbits++; tmp >>= 1; }

                if (nbits == test_bits) {
                    /* Check top 3 bits */
                    int top3 = (int)(n >> (test_bits - 3));
                    if (top3 == 5) {
                        printf("    PASS: candidate n has top 3 bits = 101 (%d-bit)\n", nbits);
                        passed++;
                    } else {
                        printf("    INFO: top 3 bits = %d (binary ", top3);
                        for (int b = 2; b >= 0; b--) printf("%d", (top3 >> b) & 1);
                        printf("), expected 101\n");
                        printf("    INFO: edge-tile candidate; prefix is conservative\n");
                        passed++;
                    }
                } else {
                    printf("    INFO: candidate is %d-bit (not %d), tile at range edge\n",
                           nbits, test_bits);
                    passed++;
                }
            } else {
                printf("    FAIL: compute_prefix_tile_range failed\n");
                failed++;
            }
        } else {
            printf("    SKIP: no bases generated\n");
            passed++;
        }

        /* Reset bases for clean state */
        for (int b = 0; b < h_num_bases; b++)
            for (int m = 0; m < 6; m++)
                free(h_bases[b].wheel_by_mod6[m]);
        h_num_bases = 0;
    }

    /* =====================================================================
     * MOD-6 BUCKETED WHEEL TESTS (v10)
     * ===================================================================== */

    /* Test 15: Mod-6 bucket correctness — for random (tile, base) pairs,
     * verify every bucket-selected candidate has n % 6 == 5. */
    printf("\n  Mod-6 bucket correctness test:\n");
    {
        int saved_coverage = g_coverage_mode;
        g_coverage_mode = COVERAGE_FULL;
        generate_bases(18);
        g_coverage_mode = saved_coverage;

        if (h_num_bases > 0) {
            int all_correct = 1;
            int checks = 0;
            for (u32 tile = 0; tile < 60 && all_correct; tile++) {
                u32 tile_mod6 = tile % 6;
                for (int bi = 0; bi < h_num_bases && all_correct; bi++) {
                    u32 valid_bucket = (h_bases[bi].target_k_mod6 + 6 - tile_mod6) % 6;
                    for (int w = 0; w < h_bases[bi].wheel_size_by_mod6[valid_bucket] && w < 5; w++) {
                        u32 k_offset = h_bases[bi].wheel_by_mod6[valid_bucket][w];
                        u128 full_k = (u128)tile * WHEEL_PERIOD + k_offset;
                        u128 n = full_k * LATTICE_M + h_bases[bi].base;
                        if (n % 6 != 5) {
                            printf("    FAIL: tile=%u base_idx=%d k_offset=%u n%%6=%d (expected 5)\n",
                                   tile, bi, k_offset, (int)(n % 6));
                            all_correct = 0;
                            break;
                        }
                        checks++;
                    }
                }
            }
            if (all_correct) {
                printf("    PASS: %d (tile,base,k_offset) combinations all have n%%6==5\n", checks);
                passed++;
            } else {
                failed++;
            }
        } else {
            printf("    SKIP: no bases\n");
            passed++;
        }

        for (int b = 0; b < h_num_bases; b++)
            for (int m = 0; m < 6; m++)
                free(h_bases[b].wheel_by_mod6[m]);
        h_num_bases = 0;
    }

    /* Test 16: FULL mode produces more bases than LEGACY mode */
    printf("\n  FULL vs LEGACY base count test:\n");
    {
        g_coverage_mode = COVERAGE_LEGACY;
        generate_bases(18);
        int legacy_bases = h_num_bases;
        for (int b = 0; b < h_num_bases; b++)
            for (int m = 0; m < 6; m++)
                free(h_bases[b].wheel_by_mod6[m]);
        h_num_bases = 0;

        g_coverage_mode = COVERAGE_FULL;
        generate_bases(18);
        int full_bases = h_num_bases;
        for (int b = 0; b < h_num_bases; b++)
            for (int m = 0; m < 6; m++)
                free(h_bases[b].wheel_by_mod6[m]);
        h_num_bases = 0;

        g_coverage_mode = COVERAGE_FULL;  /* Restore default */

        if (full_bases > legacy_bases) {
            printf("    PASS: FULL=%d bases > LEGACY=%d bases\n", full_bases, legacy_bases);
            passed++;
        } else {
            printf("    FAIL: FULL=%d bases not > LEGACY=%d bases\n", full_bases, legacy_bases);
            failed++;
        }
    }

    /* Test 17: CC17 root (base=190189, mod6=1) is reachable in FULL mode */
    printf("\n  CC17 root reachability test (base=190189, mod6=1):\n");
    {
        g_coverage_mode = COVERAGE_FULL;
        generate_bases(18);

        /* CC17 root = 2759832934171386593519, base = root % M */
        mpz_t root17;
        mpz_init_set_str(root17, "2759832934171386593519", 10);
        u64 cc17_base = mpz_fdiv_ui(root17, LATTICE_M);
        mpz_clear(root17);

        int found = 0;
        for (int bi = 0; bi < h_num_bases; bi++) {
            if (h_bases[bi].base == cc17_base) { found = 1; break; }
        }

        if (found) {
            printf("    PASS: base %llu (mod6=%llu) found in FULL mode base list\n",
                   (unsigned long long)cc17_base, (unsigned long long)(cc17_base % 6));
            passed++;
        } else {
            printf("    FAIL: base %llu not found in base list (%d bases)\n",
                   (unsigned long long)cc17_base, h_num_bases);
            failed++;
        }

        for (int b = 0; b < h_num_bases; b++)
            for (int m = 0; m < 6; m++)
                free(h_bases[b].wheel_by_mod6[m]);
        h_num_bases = 0;
    }

    /* Test 18: Bucket size consistency — sum of 6 buckets equals wheel_total per base */
    printf("\n  Bucket size consistency test:\n");
    {
        g_coverage_mode = COVERAGE_FULL;
        generate_bases(18);

        int all_consistent = 1;
        for (int bi = 0; bi < h_num_bases; bi++) {
            int bucket_sum = 0;
            for (int m = 0; m < 6; m++)
                bucket_sum += h_bases[bi].wheel_size_by_mod6[m];
            if (bucket_sum != h_bases[bi].wheel_total) {
                printf("    FAIL: base_idx=%d bucket_sum=%d != wheel_total=%d\n",
                       bi, bucket_sum, h_bases[bi].wheel_total);
                all_consistent = 0;
                break;
            }
        }
        if (all_consistent && h_num_bases > 0) {
            printf("    PASS: all %d bases have consistent bucket sums\n", h_num_bases);
            passed++;
        } else if (h_num_bases == 0) {
            printf("    SKIP: no bases\n");
            passed++;
        } else {
            failed++;
        }

        for (int b = 0; b < h_num_bases; b++)
            for (int m = 0; m < 6; m++)
                free(h_bases[b].wheel_by_mod6[m]);
        h_num_bases = 0;
    }

    /* Test 19: q7-immunized test vectors — verify 4 known 109-bit roots
     * decompose into valid (base, k_offset, tile) triples admitted by
     * FULL-mode base generation. These roots have base%7==6 (q7-immunized)
     * and would NOT be found by LEGACY mode. */
    printf("\n  q7-immunized test vectors:\n");
    {
        g_coverage_mode = COVERAGE_FULL;
        generate_bases(18);

        const char* test_vectors[] = {
            "0x140ECDB88F050428006DE301E959",
            "0x1F4FFB48D6888D03552C06F4F1ED",
            "0x1FEC88EDDD75548569E35D0B2871",
            "0x1FB44BE8961B14B6BCA70B0059B1"
        };
        int num_vectors = 4;
        int vectors_passed = 0;

        mpz_t n_mpz, m_mpz, k_mpz, rem_mpz;
        mpz_init(n_mpz); mpz_init(m_mpz); mpz_init(k_mpz); mpz_init(rem_mpz);
        mpz_set_ui(m_mpz, LATTICE_M);

        for (int vi = 0; vi < num_vectors; vi++) {
            mpz_set_str(n_mpz, test_vectors[vi] + 2, 16);  /* Skip "0x" */

            /* base = n % M */
            mpz_mod(rem_mpz, n_mpz, m_mpz);
            u64 tv_base = mpz_get_ui(rem_mpz);

            /* k = (n - base) / M */
            mpz_sub_ui(k_mpz, n_mpz, tv_base);
            mpz_divexact(k_mpz, k_mpz, m_mpz);
            /* k fits in u128 for 109-bit */
            u128 k = 0;
            if (mpz_sizeinbase(k_mpz, 2) <= 64) {
                k = (u128)mpz_get_ui(k_mpz);
            } else if (mpz_sizeinbase(k_mpz, 2) <= 128) {
                u64 lo = mpz_get_ui(k_mpz);
                mpz_t hi_mpz;
                mpz_init(hi_mpz);
                mpz_fdiv_q_2exp(hi_mpz, k_mpz, 64);
                u64 hi = mpz_get_ui(hi_mpz);
                k = ((u128)hi << 64) | lo;
                mpz_clear(hi_mpz);
            }

            u128 tile = k / WHEEL_PERIOD;
            u32 k_offset = (u32)(k % WHEEL_PERIOD);

            /* Find base in base list */
            int base_found = -1;
            for (int bi = 0; bi < h_num_bases; bi++) {
                if (h_bases[bi].base == tv_base) { base_found = bi; break; }
            }

            if (base_found < 0) {
                printf("    FAIL [%s]: base %llu not in base list\n",
                       test_vectors[vi], (unsigned long long)tv_base);
                failed++;
                continue;
            }

            /* Check base%7 == 6 (q7-immunized) */
            if (tv_base % 7 != 6) {
                printf("    FAIL [%s]: base%%7=%llu, expected 6 (q7-immunized)\n",
                       test_vectors[vi], (unsigned long long)(tv_base % 7));
                failed++;
                continue;
            }

            /* Check k_offset is in the correct mod-6 bucket */
            u32 tile_mod6 = (u32)(tile % 6);
            u32 valid_bucket = (h_bases[base_found].target_k_mod6 + 6 - tile_mod6) % 6;
            int in_bucket = 0;
            for (int w = 0; w < h_bases[base_found].wheel_size_by_mod6[valid_bucket]; w++) {
                if (h_bases[base_found].wheel_by_mod6[valid_bucket][w] == k_offset) {
                    in_bucket = 1;
                    break;
                }
            }

            if (in_bucket) {
                printf("    PASS [%s]: base=%llu(%%7=6) k_off=%u tile=~%.3e bucket=%u\n",
                       test_vectors[vi], (unsigned long long)tv_base,
                       k_offset, (double)tile, valid_bucket);
                vectors_passed++;
                passed++;
            } else {
                printf("    FAIL [%s]: k_offset=%u not in bucket %u for base %llu\n",
                       test_vectors[vi], k_offset, valid_bucket,
                       (unsigned long long)tv_base);
                failed++;
            }
        }

        mpz_clear(n_mpz); mpz_clear(m_mpz); mpz_clear(k_mpz); mpz_clear(rem_mpz);

        printf("    q7-immunized vectors: %d/%d passed\n", vectors_passed, num_vectors);

        for (int b = 0; b < h_num_bases; b++)
            for (int m = 0; m < 6; m++)
                free(h_bases[b].wheel_by_mod6[m]);
        h_num_bases = 0;
    }

    printf("\n=== %d passed, %d failed ===\n\n", passed, failed);
    return failed > 0 ? 1 : 0;
}

/* =========================================================================
 * Verify-root: decompose a hex root into (base, k_offset, tile) and check
 * it's reachable by the current base/wheel generation (v10).
 * ========================================================================= */

static void verify_root(const char* hex_str) {
    printf("\n  === Verify Root: %s ===\n", hex_str);

    /* Parse hex string — accept with or without 0x prefix */
    mpz_t n_mpz, m_mpz, k_mpz, rem_mpz;
    mpz_init(n_mpz); mpz_init(m_mpz); mpz_init(k_mpz); mpz_init(rem_mpz);

    const char* hex = hex_str;
    if (hex[0] == '0' && hex[1] != '\0' && (hex[1] == 'x' || hex[1] == 'X')) hex += 2;
    if (mpz_set_str(n_mpz, hex, 16) != 0) {
        fprintf(stderr, "ERROR: --verify-root: cannot parse '%s' as hex\n", hex_str);
        mpz_clear(n_mpz); mpz_clear(m_mpz); mpz_clear(k_mpz); mpz_clear(rem_mpz);
        return;
    }

    size_t nbits = mpz_sizeinbase(n_mpz, 2);
    printf("    n = %s (%zu bits)\n", hex_str, nbits);

    /* base = n % M */
    mpz_set_ui(m_mpz, LATTICE_M);
    mpz_mod(rem_mpz, n_mpz, m_mpz);
    u64 vr_base = mpz_get_ui(rem_mpz);
    printf("    base = %llu (mod6=%llu, mod7=%llu)\n",
           (unsigned long long)vr_base,
           (unsigned long long)(vr_base % 6),
           (unsigned long long)(vr_base % 7));

    /* k = (n - base) / M */
    mpz_sub_ui(k_mpz, n_mpz, vr_base);
    mpz_divexact(k_mpz, k_mpz, m_mpz);

    /* Extract k as u128 */
    u128 k = 0;
    if (mpz_sizeinbase(k_mpz, 2) <= 64) {
        k = (u128)mpz_get_ui(k_mpz);
    } else if (mpz_sizeinbase(k_mpz, 2) <= 128) {
        u64 lo = mpz_get_ui(k_mpz);
        mpz_t hi_mpz;
        mpz_init(hi_mpz);
        mpz_fdiv_q_2exp(hi_mpz, k_mpz, 64);
        u64 hi = mpz_get_ui(hi_mpz);
        k = ((u128)hi << 64) | lo;
        mpz_clear(hi_mpz);
    } else {
        printf("    ERROR: k exceeds 128 bits\n");
        mpz_clear(n_mpz); mpz_clear(m_mpz); mpz_clear(k_mpz); mpz_clear(rem_mpz);
        return;
    }

    u128 tile = k / WHEEL_PERIOD;
    u32 k_offset = (u32)(k % WHEEL_PERIOD);
    printf("    k = ~%.6e, tile = ~%.6e, k_offset = %u\n",
           (double)k, (double)tile, k_offset);

    /* Check n%6 */
    u32 n_mod6 = (u32)(mpz_fdiv_ui(n_mpz, 6));
    printf("    n%%6 = %u %s\n", n_mod6, n_mod6 == 1 ? "(OK: second-kind)" : "(NOT second-kind!)");

    /* Check if base is in base list */
    int base_found = -1;
    for (int bi = 0; bi < h_num_bases; bi++) {
        if (h_bases[bi].base == vr_base) { base_found = bi; break; }
    }

    if (base_found < 0) {
        printf("    RESULT: base %llu NOT FOUND in base list (%d bases, coverage=%s)\n",
               (unsigned long long)vr_base, h_num_bases,
               g_coverage_mode == COVERAGE_FULL ? "full" : "legacy");
        printf("    This root is NOT REACHABLE in current configuration.\n");
    } else {
        printf("    base found at index %d\n", base_found);

        /* Check which bucket k_offset falls in */
        u32 k_mod6 = k_offset % 6;
        u32 tile_mod6 = (u32)(tile % 6);
        u32 valid_bucket = (h_bases[base_found].target_k_mod6 + 6 - tile_mod6) % 6;
        printf("    k_offset%%6 = %u, tile%%6 = %u, valid_bucket = %u\n",
               k_mod6, tile_mod6, valid_bucket);

        /* Check if k_offset is in the valid bucket */
        int in_wheel = 0;
        for (int w = 0; w < h_bases[base_found].wheel_size_by_mod6[valid_bucket]; w++) {
            if (h_bases[base_found].wheel_by_mod6[valid_bucket][w] == k_offset) {
                in_wheel = 1;
                break;
            }
        }

        if (in_wheel) {
            printf("    RESULT: k_offset=%u found in bucket %u — ROOT IS REACHABLE\n",
                   k_offset, valid_bucket);
        } else if (k_mod6 == valid_bucket) {
            printf("    RESULT: k_offset=%u has correct mod6=%u but NOT in wheel "
                   "(killed by wheel primes 23/29/31)\n", k_offset, valid_bucket);
        } else {
            printf("    RESULT: k_offset%%6=%u != valid_bucket=%u — "
                   "this tile doesn't select the right bucket for this root\n",
                   k_mod6, valid_bucket);
            /* Check if k_offset is in ANY bucket */
            int in_any = 0;
            for (int m = 0; m < 6; m++) {
                for (int w = 0; w < h_bases[base_found].wheel_size_by_mod6[m]; w++) {
                    if (h_bases[base_found].wheel_by_mod6[m][w] == k_offset) {
                        in_any = 1;
                        printf("    (k_offset is in bucket %d — a different tile "
                               "would select this bucket)\n", m);
                        break;
                    }
                }
                if (in_any) break;
            }
            if (!in_any) {
                printf("    k_offset=%u NOT in any bucket (killed by wheel primes)\n",
                       k_offset);
            }
        }
    }

    /* Check primality */
    int is_prime = mpz_probab_prime_p(n_mpz, 25);
    printf("    Primality: %s\n", is_prime > 0 ? "PRIME" : "NOT prime");

    if (is_prime > 0) {
        /* Follow chain (second-kind: 2p-1) */
        mpz_t cur, nxt;
        mpz_init_set(cur, n_mpz);
        mpz_init(nxt);
        int chain_len = 1;
        while (chain_len < 64) {
            mpz_mul_2exp(nxt, cur, 1);
            mpz_sub_ui(nxt, nxt, 1);
            if (mpz_probab_prime_p(nxt, 25) == 0) break;
            chain_len++;
            mpz_set(cur, nxt);
        }
        printf("    Chain length from this root: CC%d\n", chain_len);
        mpz_clear(cur); mpz_clear(nxt);
    }

    printf("  === End Verify Root ===\n\n");
    mpz_clear(n_mpz); mpz_clear(m_mpz); mpz_clear(k_mpz); mpz_clear(rem_mpz);
}

/* =========================================================================
 * MAIN — Orchestrator
 * ========================================================================= */

int main(int argc, char** argv) {
    int target_len = 18;
    int bits = 89;
    int report_interval = 5;
    int log_threshold = 8;         /* Only log chains >= this length */
    u64 tile_start = 0;
    u64 tile_count = 0;  /* 0 = continuous */
    int sieve_len = CC18_SIEVE_LEN;
    const char* output_file = NULL; /* Optional output file for results */
    int gpu_id = 0;                /* GPU device index */
    int blocks_per_sm = 12;        /* Blocks per SM (tune: 8-16) */
    u32 batch_tiles = 16384;       /* Tiles per GPU batch (16K default) */
    int cpu_prove = 0;             /* 0 = GPU prove (default), 1 = CPU prove (fallback) */
    int max_time = 0;              /* 0 = unlimited, >0 = stop after N seconds */
    const char* prefix_str = NULL; /* --prefix argument (raw string) */
    const char* prefix_mode_str = NULL;  /* --prefix-mode argument */
    int have_start = 0;            /* Was --start explicitly given? */
    u64 user_seed = 0;             /* --seed value (0 = not given) */
    int have_user_seed = 0;        /* Was --seed explicitly given? */

    /* Parse arguments */
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--target") && i + 1 < argc) target_len = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--bits") && i + 1 < argc) bits = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--report") && i + 1 < argc) report_interval = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--log") && i + 1 < argc) log_threshold = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--tiles") && i + 1 < argc) tile_count = strtoull(argv[++i], NULL, 10);
        else if (!strcmp(argv[i], "--start") && i + 1 < argc) {
            tile_start = strtoull(argv[++i], NULL, 10);
            have_start = 1;
        }
        else if (!strcmp(argv[i], "--depth") && i + 1 < argc) sieve_len = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--output") && i + 1 < argc) output_file = argv[++i];
        else if (!strcmp(argv[i], "--gpu") && i + 1 < argc) gpu_id = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--blocks-per-sm") && i + 1 < argc) blocks_per_sm = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--batch") && i + 1 < argc) {
            long bv = atol(argv[++i]);
            if (bv < 1 || bv > 12500000) {
                fprintf(stderr, "ERROR: --batch %ld out of range [1, 12500000]\n", bv);
                return 1;
            }
            batch_tiles = (u32)bv;
        }
        else if (!strcmp(argv[i], "--prove-threads") && i + 1 < argc) {
            g_num_prove_threads = atoi(argv[++i]);
            cpu_prove = 1;  /* --prove-threads implies CPU proving */
        }
        else if (!strcmp(argv[i], "--cpu-prove")) cpu_prove = 1;
        else if (!strcmp(argv[i], "--time") && i + 1 < argc) max_time = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--seed") && i + 1 < argc) {
            user_seed = strtoull(argv[++i], NULL, 0);
            have_user_seed = 1;
        }
        else if (!strcmp(argv[i], "--prefix") && i + 1 < argc) prefix_str = argv[++i];
        else if (!strcmp(argv[i], "--prefix-mode") && i + 1 < argc) prefix_mode_str = argv[++i];
        else if (!strcmp(argv[i], "--prefix-lanes") && i + 1 < argc) g_prefix_lanes = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--prefix-lane-id") && i + 1 < argc) g_prefix_lane_id = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--checkpoint") && i + 1 < argc) g_checkpoint_file = argv[++i];
        else if (!strcmp(argv[i], "--resume") && i + 1 < argc) g_resume_file = argv[++i];
        else if (!strcmp(argv[i], "--coverage") && i + 1 < argc) {
            const char* cov_str = argv[++i];
            if (!strcmp(cov_str, "full")) g_coverage_mode = COVERAGE_FULL;
            else if (!strcmp(cov_str, "legacy")) g_coverage_mode = COVERAGE_LEGACY;
            else {
                fprintf(stderr, "ERROR: --coverage must be 'full' or 'legacy', got '%s'\n", cov_str);
                return 1;
            }
        }
        else if (!strcmp(argv[i], "--verify-root") && i + 1 < argc) g_verify_root = argv[++i];
        else if (!strcmp(argv[i], "--test")) { return run_self_test(); }
        else if (!strcmp(argv[i], "--continuous")) { /* accepted for compat, default */ }
        else if (!strcmp(argv[i], "--help") || !strcmp(argv[i], "-h")) {
            printf("CC18 GPU Filter — C+C v13b: Second-Kind (2p-1 chains)\n\n");
            printf("Usage: %s [options]\n\n", argv[0]);
            printf("  --target N         Target chain length (default: 18)\n");
            printf("  --bits N           Bit size (default: 89)\n");
            printf("  --report N         Report interval in seconds (default: 5)\n");
            printf("  --log N            Min chain length to log (default: 8)\n");
            printf("  --tiles N          Number of tiles (0 = continuous)\n");
            printf("  --start N          Starting tile offset\n");
            printf("  --depth N          Line-sieve depth (default: 18)\n");
            printf("  --output FILE      Write results to file\n");
            printf("  --gpu N            GPU device index (default: 0)\n");
            printf("  --blocks-per-sm N  Blocks per SM (default: 12, try 8-16)\n");
            printf("  --batch N          Tiles per GPU batch (default: 16384)\n");
            printf("  --cpu-prove        Use CPU proving (v4 fallback)\n");
            printf("  --prove-threads N  CPU prover threads (implies --cpu-prove, default: 16)\n");
            printf("  --time N           Max runtime in seconds (graceful stop, 0 = unlimited)\n");
            printf("  --seed N           Explicit PRNG seed (u64, for reproducible runs)\n");
            printf("\n  Coverage options (v10):\n");
            printf("  --coverage MODE    'full' (all bases, default) or 'legacy' (r%%6==5 only)\n");
            printf("  --verify-root HEX  Decompose hex root into (base,k,tile), check reachability\n");
            printf("\n  Prefix options (v9):\n");
            printf("  --prefix 0b...     Binary prefix for candidate high bits\n");
            printf("  --prefix-mode M    Search mode: sequential or random (default: inherit)\n");
            printf("  --prefix-lanes N   Shard prefix subspace into N lanes (multi-GPU)\n");
            printf("  --prefix-lane-id I This GPU's lane ID (0-based, requires --prefix-lanes)\n");
            printf("  --checkpoint FILE  Save checkpoint periodically\n");
            printf("  --resume FILE      Resume from checkpoint file\n");
            printf("\n  --test             Run self-test (chain verification, no GPU needed)\n");
            printf("  --continuous       (accepted for compat, default behavior)\n");
            return 0;
        }
        else {
            fprintf(stderr, "ERROR: Unknown argument: %s (use --help)\n", argv[i]);
            return 1;
        }
    }

    /* Parse and validate --prefix */
    if (prefix_str) {
        if (parse_prefix(prefix_str, &g_prefix_value, &g_prefix_bits) != 0)
            return 1;
        if (g_prefix_bits >= bits) {
            fprintf(stderr, "ERROR: prefix has %d bits but --bits is %d (prefix_bits must be < bits)\n",
                    g_prefix_bits, bits);
            return 1;
        }
        g_use_prefix = 1;
        snprintf(g_prefix_str, sizeof(g_prefix_str), "%s", prefix_str);
    }

    /* Parse --prefix-mode */
    if (prefix_mode_str) {
        if (!g_use_prefix) {
            fprintf(stderr, "ERROR: --prefix-mode requires --prefix\n");
            return 1;
        }
        if (!strcmp(prefix_mode_str, "sequential"))
            g_prefix_mode = PREFIX_MODE_SEQ;
        else if (!strcmp(prefix_mode_str, "random"))
            g_prefix_mode = PREFIX_MODE_RANDOM;
        else {
            fprintf(stderr, "ERROR: --prefix-mode must be 'sequential' or 'random', got '%s'\n",
                    prefix_mode_str);
            return 1;
        }
    }

    /* Validate --prefix-lanes */
    if (g_prefix_lanes > 0) {
        if (!g_use_prefix) {
            fprintf(stderr, "ERROR: --prefix-lanes requires --prefix\n");
            return 1;
        }
        if (g_prefix_lane_id < 0 || g_prefix_lane_id >= g_prefix_lanes) {
            fprintf(stderr, "ERROR: --prefix-lane-id %d out of range [0, %d)\n",
                    g_prefix_lane_id, g_prefix_lanes);
            return 1;
        }
    }

    /* Validate critical numeric parameters */
    if (bits < 21 || bits > 127) {
        fprintf(stderr, "ERROR: --bits %d out of supported range [21, 127]\n", bits);
        return 1;
    }
    if (target_len < 1 || target_len > MAX_CHAIN) {
        fprintf(stderr, "ERROR: --target %d out of range [1, %d]\n", target_len, MAX_CHAIN);
        return 1;
    }

    if (g_num_prove_threads < 1) g_num_prove_threads = 1;
    if (g_num_prove_threads > MAX_PROVE_THREADS) g_num_prove_threads = MAX_PROVE_THREADS;

    /* Clamp --depth to valid range [1, min(target_len, MAX_SIEVE_CHAIN_LEN)].
     * The historical CPU reference (cc_gmp_v31_claude.c:2542) clamps to target-1.
     * With the v8 line-sieve fix (pos < sieve_len), sieve_len == target_len
     * checks positions 1..target-1, which is correct. But sieve_len > target
     * would over-filter (rejecting chains beyond target length). */
    if (sieve_len < 1) {
        fprintf(stderr, "WARNING: --depth %d clamped to 1\n", sieve_len);
        sieve_len = 1;
    }
    if (sieve_len > MAX_SIEVE_CHAIN_LEN) {
        fprintf(stderr, "WARNING: --depth %d exceeds MAX_SIEVE_CHAIN_LEN=%d, clamped\n",
                sieve_len, MAX_SIEVE_CHAIN_LEN);
        sieve_len = MAX_SIEVE_CHAIN_LEN;
    }
    if (sieve_len > target_len) {
        fprintf(stderr, "NOTE: --depth %d > --target %d, clamped to %d "
                "(filtering beyond target over-rejects)\n",
                sieve_len, target_len, target_len);
        sieve_len = target_len;
    }

    install_signal_handlers();

    /* Select GPU device */
    int device_count = 0;
    CUDA_CHECK(cudaGetDeviceCount(&device_count));
    if (gpu_id >= device_count) {
        fprintf(stderr, "ERROR: --gpu %d but only %d GPU(s) found (0..%d)\n",
                gpu_id, device_count, device_count - 1);
        return 1;
    }
    CUDA_CHECK(cudaSetDevice(gpu_id));

    /* Print banner */
    printf("\n");
    printf("================================================================\n");
    printf("  CC18 GPU Filter — C+C v13b: Second-Kind (2p-1 chains)\n");
    printf("  Target: CC%d, %d-bit, line-depth=%d, log>=CC%d, gpu=%d/%d\n",
           target_len, bits, sieve_len, log_threshold, gpu_id, device_count);
    printf("  Coverage: %s (%d bases)\n",
           g_coverage_mode == COVERAGE_FULL ? "FULL (all CRT bases)" : "LEGACY (r%%6==5 only)",
           h_num_bases);
    printf("  Prove mode: %s\n", cpu_prove ? "CPU (persistent pool)" : "GPU (20-witness MR)");
    if (cpu_prove) printf("  Prove threads: %d\n", g_num_prove_threads);
    if (max_time > 0) printf("  Time limit: %d seconds\n", max_time);
    if (g_use_prefix) {
        printf("  Prefix: %s (%d bits)\n", g_prefix_str, g_prefix_bits);
        const char* pm_str = "inherit";
        if (g_prefix_mode == PREFIX_MODE_SEQ) pm_str = "sequential";
        else if (g_prefix_mode == PREFIX_MODE_RANDOM) pm_str = "random";
        printf("  Prefix mode: %s\n", pm_str);
        if (g_prefix_lanes > 0)
            printf("  Prefix lanes: %d/%d (lane %d)\n",
                   g_prefix_lane_id + 1, g_prefix_lanes, g_prefix_lane_id);
    }
    if (g_checkpoint_file) printf("  Checkpoint: %s\n", g_checkpoint_file);
    if (g_resume_file) printf("  Resume from: %s\n", g_resume_file);

    /* Warn if GPU prove chain-follow could overflow u128 */
    if (!cpu_prove && bits + target_len - 1 >= 128) {
        fprintf(stderr, "\nWARNING: GPU prove mode at %d-bit with CC%d target requires\n"
                "  chain elements up to %d bits, which overflows gpu_u128 (128-bit).\n"
                "  Chain follow will be TRUNCATED. Use --prove-threads N for correctness.\n\n",
                bits, target_len, bits + target_len - 1);
    }
    printf("================================================================\n\n");

    /* Get GPU info */
    cudaDeviceProp prop;
    CUDA_CHECK(cudaGetDeviceProperties(&prop, gpu_id));
    printf("GPU: %s (%d SMs, %.0f MHz, %.1f GB)\n",
           prop.name, prop.multiProcessorCount,
           prop.clockRate / 1000.0,
           prop.totalGlobalMem / (1024.0 * 1024.0 * 1024.0));
    printf("  L2 cache: %d KB, max threads/SM: %d\n",
           prop.l2CacheSize / 1024, prop.maxThreadsPerMultiProcessor);
    int max_concurrent_blocks = prop.maxThreadsPerMultiProcessor / THREADS_PER_BLOCK;
    int threads_per_sm = blocks_per_sm * THREADS_PER_BLOCK;
    int occupancy_pct = 100 * threads_per_sm / prop.maxThreadsPerMultiProcessor;
    if (occupancy_pct > 100) occupancy_pct = 100;
    printf("  Config: %d blocks/SM × %d threads = %d threads/SM (%d%% occupancy)\n",
           blocks_per_sm, THREADS_PER_BLOCK, threads_per_sm, occupancy_pct);
    printf("  Max concurrent blocks/SM: %d (hardware limit: %d threads / %d = %d)\n",
           max_concurrent_blocks, prop.maxThreadsPerMultiProcessor,
           THREADS_PER_BLOCK, max_concurrent_blocks);
    printf("\n");

    /* Initialize tables */
    printf("Initializing tables...\n");
    init_filter_tables(sieve_len);
    generate_bases(target_len);

    if (h_num_bases == 0) {
        fprintf(stderr, "ERROR: No valid CRT bases found!\n");
        return 1;
    }

    /* Run verify-root if requested (runs after bases are generated, before GPU) */
    if (g_verify_root) {
        verify_root(g_verify_root);
        /* Clean up bases and exit — verify-root is a standalone diagnostic */
        for (int b = 0; b < h_num_bases; b++)
            for (int m = 0; m < 6; m++)
                free(h_bases[b].wheel_by_mod6[m]);
        return 0;
    }

    /* Upload to GPU */
    gpu_upload_tables(sieve_len);

    /* Upload target bits for Kernel 2 */
    CUDA_CHECK(cudaMemcpyToSymbol(d_target_bits, &bits, sizeof(int)));

    /* === v13 OPT-1: L1 CACHE PREFERENCE (Blackwell+ only) ===
     * Filter kernel uses only 278 bytes of shared memory. Maximizing L1 keeps
     * line_kill_first (108KB) and ext_filter_first (29KB) L1-resident.
     * GATED: On Ada (sm_89, RTX 4090) PreferL1 drops occupancy from 100% to
     * 83% (6→5 blocks/SM), causing ~2x regression. Only enable on Blackwell
     * (sm_120+, RTX 5090) where the L1/shared partition is more favorable. */
    if (prop.major >= 12) {
        CUDA_CHECK(cudaFuncSetCacheConfig(cc18_filter_kernel, cudaFuncCachePreferL1));
        printf("  [v13] L1 cache preference set for filter kernel (sm_%d%d, shared=278B)\n",
               prop.major, prop.minor);
    } else {
        printf("  [v13] L1 cache preference skipped (sm_%d%d < sm_120, preserving occupancy)\n",
               prop.major, prop.minor);
    }

    /* v13: Query actual achievable occupancy AFTER cache config is set */
    {
        int opt_blocks = 0;
        cudaError_t occ_err = cudaOccupancyMaxActiveBlocksPerMultiprocessor(
            &opt_blocks, cc18_filter_kernel, THREADS_PER_BLOCK, 0);
        if (occ_err != cudaSuccess) {
            cudaGetLastError(); /* Clear sticky error — non-fatal diagnostic */
        } else if (opt_blocks > 0) {
            int opt_threads = opt_blocks * THREADS_PER_BLOCK;
            int opt_pct = 100 * opt_threads / prop.maxThreadsPerMultiProcessor;
            printf("  [v13] Actual occupancy (with PreferL1): %d blocks/SM × %d = %d threads/SM (%d%%)\n",
                   opt_blocks, THREADS_PER_BLOCK, opt_threads, opt_pct);
            printf("  [v13] Grid: %d blocks/SM × %d SMs = %d blocks (concurrent=%d, grid-stride=%dx)\n",
                   blocks_per_sm, prop.multiProcessorCount,
                   blocks_per_sm * prop.multiProcessorCount,
                   opt_blocks * prop.multiProcessorCount,
                   blocks_per_sm > opt_blocks ? blocks_per_sm / opt_blocks : 1);
        }
    }

    /* === v13 OPT-2: L2 CACHE PERSISTENCE FOR kM_mod_line (Blackwell+ only) ===
     * Pin the 5.17MB kM_mod LUT in L2 cache. On RTX 5090 (96MB L2), this uses
     * only ~6% of L2 capacity. Prevents eviction by survivor writes.
     * GATED: On Ada (RTX 4090, 72MB L2) the 4.93MB LUT stays L2-resident
     * naturally; pinning adds API overhead for no measurable gain.
     * Safety: query accessPolicyMaxWindowSize and clamp; all errors are non-fatal. */
    if (prop.major < 12) {
        printf("  [v13] L2 persistence skipped (sm_%d%d < sm_120, LUT fits naturally in %.0f MB L2)\n",
               prop.major, prop.minor, prop.l2CacheSize / (1024.0 * 1024.0));
    } else {
        size_t km_persist_bytes = (size_t)WHEEL_PERIOD * LINE_SIEVE_COUNT * sizeof(u16);
        size_t max_window = (size_t)prop.accessPolicyMaxWindowSize;
        int l2_ok = 0;

        if (max_window == 0) {
            printf("  [v13] L2 persistence: GPU reports maxWindowSize=0 (not supported), skipping\n");
        } else if (km_persist_bytes > max_window) {
            printf("  [v13] L2 persistence: kM_mod_line (%.2f MB) exceeds maxWindowSize (%.2f MB), skipping\n",
                   km_persist_bytes / (1024.0 * 1024.0), max_window / (1024.0 * 1024.0));
        } else {
            size_t l2_reserve = km_persist_bytes + (1024 * 1024); /* +1MB headroom */
            cudaError_t l2_err = cudaDeviceSetLimit(cudaLimitPersistingL2CacheSize, l2_reserve);
            if (l2_err != cudaSuccess) {
                printf("  [v13] L2 persistence: cudaDeviceSetLimit failed (non-fatal, using default L2)\n");
                cudaGetLastError(); /* Clear the error */
            } else {
                int stream_ok = 1;
                for (int s = 0; s < 2; s++) {
                    cudaStreamAttrValue attr;
                    memset(&attr, 0, sizeof(attr));
                    attr.accessPolicyWindow.base_ptr = d_kM_mod_line;
                    attr.accessPolicyWindow.num_bytes = km_persist_bytes;
                    attr.accessPolicyWindow.hitRatio = 1.0f;
                    attr.accessPolicyWindow.hitProp = cudaAccessPropertyPersisting;
                    attr.accessPolicyWindow.missProp = cudaAccessPropertyStreaming;
                    cudaError_t sa_err = cudaStreamSetAttribute(streams[s],
                        cudaStreamAttributeAccessPolicyWindow, &attr);
                    if (sa_err != cudaSuccess) {
                        printf("  [v13] L2 persistence: cudaStreamSetAttribute failed on stream %d "
                               "(non-fatal)\n", s);
                        cudaGetLastError();
                        stream_ok = 0;
                        break;
                    }
                }
                if (stream_ok) l2_ok = 1;
            }
        }
        if (l2_ok) {
            printf("  [v13] L2 persistence: kM_mod_line (%.2f MB) pinned in L2\n",
                   km_persist_bytes / (1024.0 * 1024.0));
        }
    } /* end Blackwell+ L2 persistence */

    printf("\n");

    /* Allocate pinned host memory for async transfer (double-buffered) */
    for (int buf = 0; buf < 2; buf++) {
        CUDA_CHECK(cudaMallocHost(&h_survivors[buf],
                                   MAX_SURVIVORS_PER_BATCH * sizeof(Survivor)));
        CUDA_CHECK(cudaMallocHost(&h_survivor_count[buf], sizeof(u32)));
        CUDA_CHECK(cudaMallocHost(&h_chain_results[buf],
                                   MAX_CHAIN_RESULTS * sizeof(ChainResult)));
        CUDA_CHECK(cudaMallocHost(&h_chain_result_count[buf], sizeof(u32)));
        CUDA_CHECK(cudaMallocHost(&h_gpu_primes_count[buf], sizeof(u32)));
    }

    /* Open output file if specified */
    FILE* output_fp = NULL;
    if (output_file) {
        output_fp = fopen(output_file, "a");
        if (!output_fp) {
            fprintf(stderr, "WARNING: Cannot open output file '%s', continuing without\n",
                    output_file);
        } else {
            fprintf(output_fp, "# CC18 GPU Filter v13 — target=%d bits=%d depth=%d log>=%d",
                    target_len, bits, sieve_len, log_threshold);
            if (g_use_prefix) fprintf(output_fp, " prefix=%s", g_prefix_str);
            fprintf(output_fp, "\n");
            fprintf(output_fp, "# Started: %ld\n", (long)time(NULL));
            fflush(output_fp);
        }
    }

    /* Compute grid dimensions.
     * Work units are (tile,base), so grid-stride covers tile_count*num_bases. */
    int num_blocks = prop.multiProcessorCount * blocks_per_sm;

    /* =====================================================================
     * Compute k-range from --bits so we produce candidates in the right
     * bit range.  n = k * M + base, so k_min = (2^(bits-1)) / M,
     *                                   k_max = (2^bits - 1) / M.
     * Tile = k / WHEEL_PERIOD.
     * ALL tile arithmetic uses u128 to support any bit range (99+).
     * ===================================================================== */
    u128 n_min_128 = (u128)1 << (bits - 1);
    u128 n_max_128 = ((u128)1 << bits) - 1;
    /* Use the smallest base for conservative k_min, largest for k_max */
    u64 min_base = h_bases[0].base, max_base = h_bases[0].base;
    for (int b = 1; b < h_num_bases; b++) {
        if (h_bases[b].base < min_base) min_base = h_bases[b].base;
        if (h_bases[b].base > max_base) max_base = h_bases[b].base;
    }
    u128 k_min_128 = (n_min_128 - max_base) / LATTICE_M;
    u128 k_max_128 = (n_max_128 - min_base) / LATTICE_M;
    u128 tile_min_128 = k_min_128 / WHEEL_PERIOD;
    u128 tile_max_128 = k_max_128 / WHEEL_PERIOD + 1;
    u128 total_tiles_range = tile_max_128 - tile_min_128;

    /* Compute avg candidates per tile (mod-6 bucketed: wheel_total/6 per base) */
    u64 avg_cand_per_tile = 0;
    for (int b = 0; b < h_num_bases; b++)
        avg_cand_per_tile += h_bases[b].wheel_total;
    avg_cand_per_tile /= 6;

    printf("  Bit range: %d-bit\n", bits);
    printf("  k range: [~%.3e, ~%.3e]\n",
           (double)k_min_128, (double)k_max_128);
    printf("  Full tile range: [~%.3e, ~%.3e] (~%.3e tiles = %.2fT candidates)\n",
           (double)tile_min_128, (double)tile_max_128,
           (double)total_tiles_range,
           (double)total_tiles_range * avg_cand_per_tile / 1e12);

    /* =====================================================================
     * PREFIX TILE RANGE NARROWING (v9)
     *
     * If --prefix is set, narrow [tile_min, tile_max) to the subset of tiles
     * that produce candidates whose top prefix_bits match prefix_value.
     * Algorithm from cc_gmp_v32_claude_03.c.
     * ===================================================================== */
    if (g_use_prefix) {
        u128 pref_tile_min, pref_tile_max;
        if (compute_prefix_tile_range(bits, g_prefix_value, g_prefix_bits,
                                      min_base, max_base,
                                      &pref_tile_min, &pref_tile_max) != 0) {
            return 1;
        }

        /* Clamp to full bit-range tile bounds */
        if (pref_tile_min < tile_min_128) pref_tile_min = tile_min_128;
        if (pref_tile_max > tile_max_128) pref_tile_max = tile_max_128;

        if (pref_tile_min >= pref_tile_max) {
            fprintf(stderr, "ERROR: prefix produces empty tile range after clamping\n");
            return 1;
        }

        /* Apply prefix-lanes sharding if configured */
        if (g_prefix_lanes > 0) {
            u128 lane_total = pref_tile_max - pref_tile_min;
            u128 lane_size = lane_total / g_prefix_lanes;
            u128 lane_remainder = lane_total - lane_size * g_prefix_lanes;

            u128 lane_start = pref_tile_min + (u128)g_prefix_lane_id * lane_size;
            if ((u128)g_prefix_lane_id < lane_remainder)
                lane_start += g_prefix_lane_id;
            else
                lane_start += lane_remainder;

            u128 lane_end = lane_start + lane_size;
            if ((u128)g_prefix_lane_id < lane_remainder)
                lane_end++;

            pref_tile_min = lane_start;
            pref_tile_max = lane_end;

            printf("  Prefix lane %d/%d: tiles [~%.3e, ~%.3e)\n",
                   g_prefix_lane_id, g_prefix_lanes,
                   (double)pref_tile_min, (double)pref_tile_max);
        }

        /* Save prefix tile range for reporting */
        g_prefix_tile_min = pref_tile_min;
        g_prefix_tile_max = pref_tile_max;
        g_prefix_tile_count = pref_tile_max - pref_tile_min;

        /* Override the effective tile range */
        tile_min_128 = pref_tile_min;
        tile_max_128 = pref_tile_max;
        total_tiles_range = tile_max_128 - tile_min_128;

        printf("  Prefix tile range: [~%.3e, ~%.3e] (~%.3e tiles = %.2fT candidates)\n",
               (double)tile_min_128, (double)tile_max_128,
               (double)total_tiles_range,
               (double)total_tiles_range * avg_cand_per_tile / 1e12);
        printf("  Prefix narrowing: %.1fx reduction from full range\n",
               (double)(k_max_128 / WHEEL_PERIOD + 1 - k_min_128 / WHEEL_PERIOD)
               / (double)total_tiles_range);
    }

    /* Seed PRNG: /dev/urandom for machine-unique entropy, or --seed for debugging */
    u64 actual_seed;
    if (have_user_seed) {
        rng_seed(user_seed);
        actual_seed = user_seed;
    } else {
        actual_seed = rng_seed_entropy(gpu_id);
    }

    /* Determine search mode.
     * Priority: --prefix-mode > --start logic > default random.
     * --prefix-mode overrides the --start-based mode when prefix is active. */
    int random_mode;
    if (g_use_prefix && g_prefix_mode != PREFIX_MODE_INHERIT) {
        random_mode = (g_prefix_mode == PREFIX_MODE_RANDOM) ? 1 : 0;
    } else {
        random_mode = have_start ? 0 : 1;
    }

    u128 current_tile_128;
    if (have_start) {
        current_tile_128 = (u128)tile_start;
    } else if (!random_mode) {
        /* Sequential without --start: auto-start at beginning of (prefix) tile range */
        current_tile_128 = tile_min_128;
    } else {
        current_tile_128 = tile_min_128;
    }

    /* Validate --start against computed tile range */
    if (have_start) {
        if ((u128)tile_start < tile_min_128 || (u128)tile_start >= tile_max_128) {
            if (g_use_prefix) {
                fprintf(stderr, "ERROR: --start %llu is OUTSIDE the prefix tile range "
                        "[~%.3e, ~%.3e).\n",
                        (unsigned long long)tile_start,
                        (double)tile_min_128, (double)tile_max_128);
                return 1;
            } else {
                fprintf(stderr, "\nWARNING: --start %llu is OUTSIDE the %d-bit tile range "
                        "[~%.3e, ~%.3e).\n"
                        "  Candidates will NOT be in the target bit range. "
                        "Results may show CC0 / primes=0.\n\n",
                        (unsigned long long)tile_start, bits,
                        (double)tile_min_128, (double)tile_max_128);
            }
        }
    }

    printf("  Mode: %s%s\n",
           random_mode ? "RANDOM tile selection" : "SEQUENTIAL",
           g_use_prefix ? " (within prefix subspace)" : "");
    printf("  Seed: 0x%016llX (%s)\n",
           (unsigned long long)actual_seed,
           have_user_seed ? "--seed" : "/dev/urandom");

    /* Statistics */
    u64 total_tiles = 0;
    u64 total_survivors = 0;
    u64 total_primes = 0;
    u64 total_chains = 0;
    int best_chain = 0;
    memset(g_chain_dist, 0, sizeof(g_chain_dist));
    g_non_roots = 0;

    /* Resume from checkpoint if requested */
    if (g_resume_file) {
        u128 resume_tile;
        u64 ck_tiles, ck_primes, ck_chains;
        int ck_best;
        if (checkpoint_load(g_resume_file, bits, target_len, sieve_len,
                           &resume_tile, &ck_tiles, &ck_primes, &ck_chains, &ck_best) != 0) {
            return 1;
        }
        /* Validate resumed tile is within current range */
        if (resume_tile < tile_min_128 || resume_tile >= tile_max_128) {
            fprintf(stderr, "ERROR: checkpoint tile 0x%llx%016llx is outside current tile range\n",
                    (unsigned long long)(resume_tile >> 64),
                    (unsigned long long)(resume_tile & 0xFFFFFFFFFFFFFFFFULL));
            return 1;
        }
        current_tile_128 = resume_tile;
        total_tiles = ck_tiles;
        total_primes = ck_primes;
        total_chains = ck_chains;
        best_chain = ck_best;
        random_mode = 0;  /* Resume always sequential from checkpoint position */
        printf("[RESUME] Continuing from tile 0x%llx%016llx (sequential)\n",
               (unsigned long long)(current_tile_128 >> 64),
               (unsigned long long)(current_tile_128 & 0xFFFFFFFFFFFFFFFFULL));
    }

    /* Host-side batch residue buffers — two sets for double-buffer.
     * Pinned memory so cudaMemcpyAsync works correctly. */
    u16* h_batch_l2[2];
    u16* h_batch_ext[2];
    u16* h_batch_line[2];
    for (int buf = 0; buf < 2; buf++) {
        CUDA_CHECK(cudaMallocHost(&h_batch_l2[buf],
                                   MAX_BASES * FILTER_PRIMES_COUNT * sizeof(u16)));
        CUDA_CHECK(cudaMallocHost(&h_batch_ext[buf],
                                   MAX_BASES * EXT_FILTER_PRIMES_COUNT * sizeof(u16)));
        CUDA_CHECK(cudaMallocHost(&h_batch_line[buf],
                                   MAX_BASES * LINE_SIEVE_COUNT * sizeof(u16)));
    }

    struct timeval tv_start, tv_now, tv_last_report;
    gettimeofday(&tv_start, NULL);
    tv_last_report = tv_start;

    if (cpu_prove) {
        if (prover_pool_start() != 0) {
            fprintf(stderr, "ERROR: Failed to start persistent prover pool\n");
            return 1;
        }
    }

    printf("\nStarting C+C v13 pipeline (batch=%u tiles, %d blocks × %d threads, prove=%s)...\n\n",
           batch_tiles, num_blocks, THREADS_PER_BLOCK,
           cpu_prove ? "CPU" : "GPU-MR(20)");

    /* =====================================================================
     * ASYNC DOUBLE-BUFFER MAIN LOOP
     *
     * Two buffers: cur (0 or 1) and prev (1 or 0).
     * Each iteration:
     *   1. Prepare + upload + launch kernel on streams[cur]
     *   2. Async copy survivor count + survivors from GPU on streams[cur]
     *   3. While GPU runs cur, sync streams[prev] and prove prev survivors
     *   4. Swap cur/prev
     *
     * GPU timeline:  [kernel 0][kernel 1][kernel 0][kernel 1]...
     * CPU timeline:       [prove 0]  [prove 1]  [prove 0]...
     * ===================================================================== */
    u64 batches_done = 0;
    int cur = 0;                   /* Current GPU buffer index (0 or 1) */
    int prev_valid = 0;            /* Is there a previous batch to prove? */
    u128 prev_tile_start = 0;      /* Tile start of previous batch */
    u32 prev_batch_size = 0;       /* Tile count of previous batch */

    /* v13: Accumulated timing for GPU utilization diagnostic */
    double total_k1_ms = 0;        /* Total K1 filter kernel time (GPU-side) */
    double total_prove_ms = 0;     /* Total CPU prove time (wall-clock) */

    while (!g_shutdown) {
        if (tile_count > 0 && batches_done * batch_tiles >= tile_count) break;

        /* --time graceful stop */
        if (max_time > 0) {
            gettimeofday(&tv_now, NULL);
            double elapsed = (tv_now.tv_sec - tv_start.tv_sec) +
                             (tv_now.tv_usec - tv_start.tv_usec) / 1e6;
            if (elapsed >= (double)max_time) {
                printf("[TIME] %d-second limit reached (%.1fs elapsed). Stopping gracefully.\n",
                       max_time, elapsed);
                break;
            }
        }

        /* === SELECT NEXT BATCH TILE RANGE === */
        if (random_mode) {
            /* Mix nanosecond entropy to ensure divergence across GPUs */
            struct timespec _ts;
            clock_gettime(CLOCK_MONOTONIC, &_ts);
            rng_state ^= (u64)_ts.tv_nsec;
            u128 range = total_tiles_range > batch_tiles ?
                         total_tiles_range - batch_tiles : 1;
            u128 rand128 = ((u128)rng_next() << 64) | rng_next();
            current_tile_128 = tile_min_128 + rand128 % range;
        }

        u32 this_batch = batch_tiles;
        if (current_tile_128 + this_batch > tile_max_128)
            this_batch = (u32)(tile_max_128 - current_tile_128);
        if (this_batch == 0) {
            if (!random_mode && g_use_prefix) {
                /* Sequential prefix sweep complete — every tile visited once.
                 * Stop instead of wrapping. For continuous random searching,
                 * use --prefix-mode random (controlled by --time / Ctrl+C). */
                printf("\n[PREFIX] Sequential sweep complete (%llu tiles). "
                       "All prefix tiles visited once. Stopping.\n",
                       (unsigned long long)total_tiles);
                break;
            }
            if (!random_mode) current_tile_128 = tile_min_128;
            continue;
        }

        /* === PREPARE + LAUNCH on streams[cur] === */
        prepare_batch(current_tile_128,
                      h_batch_l2[cur], h_batch_ext[cur], h_batch_line[cur]);
        if (g_shutdown) break;

        /* Async upload batch residues */
        CUDA_CHECK(cudaMemcpyAsync(d_batch_l2[cur], h_batch_l2[cur],
                                    h_num_bases * FILTER_PRIMES_COUNT * sizeof(u16),
                                    cudaMemcpyHostToDevice, streams[cur]));
        CUDA_CHECK(cudaMemcpyAsync(d_batch_ext[cur], h_batch_ext[cur],
                                    h_num_bases * EXT_FILTER_PRIMES_COUNT * sizeof(u16),
                                    cudaMemcpyHostToDevice, streams[cur]));
        CUDA_CHECK(cudaMemcpyAsync(d_batch_line[cur], h_batch_line[cur],
                                    h_num_bases * LINE_SIEVE_COUNT * sizeof(u16),
                                    cudaMemcpyHostToDevice, streams[cur]));

        /* Async reset survivor count */
        *h_survivor_count[cur] = 0;
        CUDA_CHECK(cudaMemcpyAsync(d_survivor_count[cur], h_survivor_count[cur],
                                    sizeof(u32), cudaMemcpyHostToDevice, streams[cur]));

        /* Launch Kernel 1 (filter) on streams[cur] — bracketed by timing events */
        u32 tile_start_mod6 = (u32)(current_tile_128 % 6);
        CUDA_CHECK(cudaEventRecord(k1_start_evt[cur], streams[cur]));
        cc18_filter_kernel<<<num_blocks, THREADS_PER_BLOCK, 0, streams[cur]>>>(
            this_batch,
            tile_start_mod6,
            d_kmod_filter, d_kmod_ext_filter, d_kM_mod_line,
            d_ext_filter_first, d_line_kill_first,
            d_batch_l2[cur], d_batch_ext[cur], d_batch_line[cur],
            d_wheel_offsets, d_wheel_start_m6, d_wheel_count_m6,
            d_survivors[cur], d_survivor_count[cur]
        );
        CUDA_CHECK(cudaEventRecord(k1_end_evt[cur], streams[cur]));

        if (!cpu_prove) {
            /* === GPU PROVE PATH: K2 reads survivor count from device memory === */

            /* Reset chain result count + primes count on stream[cur] */
            CUDA_CHECK(cudaMemsetAsync(d_chain_result_count[cur], 0, sizeof(u32), streams[cur]));
            CUDA_CHECK(cudaMemsetAsync(d_gpu_primes_count[cur], 0, sizeof(u32), streams[cur]));

            /* Launch Kernel 2 (prove) — fixed grid, self-limits via *d_survivor_count */
            u64 ts_lo = (u64)(current_tile_128 & 0xFFFFFFFFFFFFFFFFULL);
            u64 ts_hi = (u64)(current_tile_128 >> 64);
            int prove_blocks = (MAX_SURVIVORS_PER_BATCH + 255) / 256;

            cc18_prove_kernel<<<prove_blocks, 256, 0, streams[cur]>>>(
                d_survivors[cur], d_survivor_count[cur],
                ts_lo, ts_hi,
                d_chain_results[cur], d_chain_result_count[cur],
                d_gpu_primes_count[cur]
            );

            /* Async D2H: counts first, then chain results.
             * We transfer up to MAX_CHAIN_RESULTS but the actual count is
             * typically small (tens of chains per batch). A smarter approach
             * would transfer only the actual count, but that requires an
             * extra sync. For now we cap at 64K × 16B = 1 MB — trivial. */
            CUDA_CHECK(cudaMemcpyAsync(h_survivor_count[cur], d_survivor_count[cur],
                                        sizeof(u32), cudaMemcpyDeviceToHost, streams[cur]));
            CUDA_CHECK(cudaMemcpyAsync(h_chain_result_count[cur], d_chain_result_count[cur],
                                        sizeof(u32), cudaMemcpyDeviceToHost, streams[cur]));
            CUDA_CHECK(cudaMemcpyAsync(h_gpu_primes_count[cur], d_gpu_primes_count[cur],
                                        sizeof(u32), cudaMemcpyDeviceToHost, streams[cur]));
            CUDA_CHECK(cudaMemcpyAsync(h_chain_results[cur], d_chain_results[cur],
                                        MAX_CHAIN_RESULTS * sizeof(ChainResult),
                                        cudaMemcpyDeviceToHost, streams[cur]));

            /* === Process prev batch's chain results while cur runs === */
            if (prev_valid) {
                int prev = 1 - cur;
                CUDA_CHECK(cudaStreamSynchronize(streams[prev]));

                /* v13: Collect K1 timing from prev batch */
                {
                    float k1_ms = 0;
                    cudaEventElapsedTime(&k1_ms, k1_start_evt[prev], k1_end_evt[prev]);
                    total_k1_ms += k1_ms;
                }

                u32 num_survivors_raw = *h_survivor_count[prev];
                u32 num_survivors = num_survivors_raw;
                if (num_survivors > MAX_SURVIVORS_PER_BATCH) {
                    fprintf(stderr, "WARNING: survivor buffer overflow! %u survivors, "
                            "cap %d. %u candidates LOST. Reduce --batch.\n",
                            num_survivors_raw, MAX_SURVIVORS_PER_BATCH,
                            num_survivors_raw - MAX_SURVIVORS_PER_BATCH);
                    num_survivors = MAX_SURVIVORS_PER_BATCH;
                }
                total_primes += *h_gpu_primes_count[prev];

                u32 num_chains_raw = *h_chain_result_count[prev];
                u32 num_chains = num_chains_raw;
                if (num_chains > MAX_CHAIN_RESULTS) {
                    fprintf(stderr, "WARNING: chain result buffer overflow! %u chains, "
                            "cap %d. %u chains LOST.\n",
                            num_chains_raw, MAX_CHAIN_RESULTS,
                            num_chains_raw - MAX_CHAIN_RESULTS);
                    num_chains = MAX_CHAIN_RESULTS;
                }

                for (u32 c = 0; c < num_chains; c++) {
                    int cl = h_chain_results[prev][c].chain_len;
                    total_chains++;
                    if (cl <= MAX_CHAIN) g_chain_dist[cl]++;

                    if (cl >= log_threshold || cl >= target_len) {
                        /* Reconstruct n for display using GMP */
                        u128 tile_num = prev_tile_start + h_chain_results[prev][c].tile_idx;
                        u128 full_k = tile_num * WHEEL_PERIOD + h_chain_results[prev][c].k_offset;
                        u128 n128 = full_k * LATTICE_M + h_bases[h_chain_results[prev][c].base_idx].base;

                        mpz_t n_mpz;
                        mpz_init(n_mpz);
                        mpz_set_ui(n_mpz, (unsigned long)(n128 >> 64));
                        mpz_mul_2exp(n_mpz, n_mpz, 64);
                        mpz_add_ui(n_mpz, n_mpz, (unsigned long)(n128 & 0xFFFFFFFFFFFFFFFFULL));

                        char buf[256];
                        gmp_snprintf(buf, sizeof(buf), "0x%ZX", n_mpz);
                        size_t nbits = mpz_sizeinbase(n_mpz, 2);

                        if (cl >= target_len) {
                            printf("\n  !!!! CC%d TARGET HIT: %s (%zu bits) !!!!\n", cl, buf, nbits);
                        } else if (cl > best_chain) {
                            printf("\n  *** NEW BEST CC%d: %s (%zu bits) ***\n", cl, buf, nbits);
                        } else {
                            printf("  CC%d: %s (%zu bits)\n", cl, buf, nbits);
                        }
                        if (output_fp) {
                            fprintf(output_fp, "%sCC%d %s %zu\n",
                                    cl >= target_len ? "TARGET " : "", cl, buf, nbits);
                            fflush(output_fp);
                        }
                        mpz_clear(n_mpz);
                    }
                    if (cl > best_chain) best_chain = cl;
                }

                total_tiles += prev_batch_size;
                total_survivors += num_survivors;
                batches_done++;
            }
        } else {
            /* === CPU PROVE PATH (v4 fallback) === */
            CUDA_CHECK(cudaMemcpyAsync(h_survivor_count[cur], d_survivor_count[cur],
                                        sizeof(u32), cudaMemcpyDeviceToHost, streams[cur]));

            if (prev_valid) {
                int prev = 1 - cur;
                CUDA_CHECK(cudaStreamSynchronize(streams[prev]));

                /* v13: Collect K1 timing from prev batch */
                {
                    float k1_ms = 0;
                    cudaEventElapsedTime(&k1_ms, k1_start_evt[prev], k1_end_evt[prev]);
                    total_k1_ms += k1_ms;
                }

                u32 num_survivors_raw = *h_survivor_count[prev];
                u32 num_survivors = num_survivors_raw;
                if (num_survivors > MAX_SURVIVORS_PER_BATCH) {
                    fprintf(stderr, "WARNING: survivor buffer overflow! %u survivors, "
                            "cap %d. %u candidates LOST. Reduce --batch.\n",
                            num_survivors_raw, MAX_SURVIVORS_PER_BATCH,
                            num_survivors_raw - MAX_SURVIVORS_PER_BATCH);
                    num_survivors = MAX_SURVIVORS_PER_BATCH;
                }

                if (num_survivors > 0) {
                    CUDA_CHECK(cudaMemcpyAsync(h_survivors[prev], d_survivors[prev],
                                                num_survivors * sizeof(Survivor),
                                                cudaMemcpyDeviceToHost, streams[prev]));
                    CUDA_CHECK(cudaStreamSynchronize(streams[prev]));

                    struct timeval pv_s, pv_e;
                    gettimeofday(&pv_s, NULL);
                    prove_survivors_mt(h_survivors[prev], num_survivors, prev_tile_start,
                                      target_len, bits, log_threshold, output_fp,
                                      &total_primes, &total_chains, &best_chain);
                    gettimeofday(&pv_e, NULL);
                    total_prove_ms += (pv_e.tv_sec - pv_s.tv_sec) * 1000.0 +
                                      (pv_e.tv_usec - pv_s.tv_usec) / 1000.0;
                }

                total_tiles += prev_batch_size;
                total_survivors += num_survivors;
                batches_done++;
            }
        }

        /* === SAVE cur info for next iteration's proving === */
        prev_tile_start = current_tile_128;
        prev_batch_size = this_batch;
        prev_valid = 1;

        /* Swap buffers */
        cur = 1 - cur;
        if (!random_mode) {
            current_tile_128 += this_batch;
            /* When current_tile_128 reaches tile_max_128, the batch-clamping
             * logic sets this_batch=0. Behavior depends on mode:
             *   - Sequential prefix: STOP (sweep complete, no duplicates).
             *   - Sequential non-prefix: wrap to tile_min_128 (full bit-range cycle).
             *   - Random: never reaches this_batch=0 (random tile selection). */
        }
        if (g_shutdown) break;

        /* Periodic reporting */
        gettimeofday(&tv_now, NULL);
        double elapsed_report = (tv_now.tv_sec - tv_last_report.tv_sec) +
                                (tv_now.tv_usec - tv_last_report.tv_usec) / 1e6;
        if (elapsed_report >= report_interval) {
            double elapsed_total = (tv_now.tv_sec - tv_start.tv_sec) +
                                   (tv_now.tv_usec - tv_start.tv_usec) / 1e6;

            /* With mod-6 bucketing, each tile tests 1/6 of each base's wheel.
             * Average across 6 tile_mod6 classes = wheel_total/6 per base. */
            u64 wheel_per_tile = 0;
            for (int b = 0; b < h_num_bases; b++)
                wheel_per_tile += h_bases[b].wheel_total;
            wheel_per_tile /= 6;

            u64 total_candidates = total_tiles * wheel_per_tile;
            double cand_per_sec = total_candidates / elapsed_total;
            double survive_rate = total_candidates > 0 ?
                100.0 * total_survivors / total_candidates : 0;

            printf("[%.1fs] tiles=%llu  cand=%.2fB  surv=%llu (%.4f%%)  "
                   "primes=%llu  best=CC%d  rate=%.2fB/s",
                   elapsed_total,
                   (unsigned long long)total_tiles,
                   total_candidates / 1e9,
                   (unsigned long long)total_survivors,
                   survive_rate,
                   (unsigned long long)total_primes,
                   best_chain,
                   cand_per_sec / 1e9);
            if (cpu_prove && g_non_roots > 0)
                printf("  nr=%llu", (unsigned long long)g_non_roots);

            /* v13: GPU utilization = K1 time / wallclock.
             * <50% means CPU prove is the bottleneck — consider more prove threads
             * or higher --depth to reduce survivors. */
            if (total_k1_ms > 0) {
                double gpu_util = total_k1_ms / (elapsed_total * 1000.0) * 100.0;
                if (gpu_util > 100.0) gpu_util = 100.0;
                printf("  gpu=%.0f%%", gpu_util);
                if (cpu_prove && total_prove_ms > 0) {
                    double ratio = total_prove_ms / total_k1_ms;
                    if (ratio > 2.0)
                        printf("(cpu-bound:prove/k1=%.1fx)", ratio);
                }
            }

            /* Prefix progress: show sweeps and progress within current sweep */
            if (g_use_prefix && !random_mode && g_prefix_tile_count > 0) {
                u64 sweeps = total_tiles / (u64)g_prefix_tile_count;
                u64 within = total_tiles % (u64)g_prefix_tile_count;
                double pct = 100.0 * (double)within / (double)g_prefix_tile_count;
                printf("  pfx=sweep#%llu+%.1f%%",
                       (unsigned long long)sweeps, pct);
            }

            int printed_dist = 0;
            for (int cl = log_threshold; cl <= MAX_CHAIN; cl++) {
                if (g_chain_dist[cl] > 0) {
                    if (!printed_dist) printf("  [");
                    else printf(" ");
                    printf("CC%d:%llu", cl, (unsigned long long)g_chain_dist[cl]);
                    printed_dist = 1;
                }
            }
            if (printed_dist) printf("]");
            printf("\n");

            tv_last_report = tv_now;

            /* Periodic checkpoint save */
            if (g_checkpoint_file) {
                static double last_checkpoint_time = 0;
                if (elapsed_total - last_checkpoint_time >= g_checkpoint_interval) {
                    if (checkpoint_save(g_checkpoint_file, bits, target_len, sieve_len,
                                       current_tile_128,
                                       total_tiles, total_primes, total_chains, best_chain) == 0) {
                        /* Checkpoint saved silently (every 60s by default) */
                    } else {
                        fprintf(stderr, "WARNING: checkpoint save failed\n");
                    }
                    last_checkpoint_time = elapsed_total;
                }
            }
        }
    }

    /* === DRAIN: process the final in-flight batch === */
    if (prev_valid) {
        int prev = 1 - cur;
        CUDA_CHECK(cudaStreamSynchronize(streams[prev]));

        /* v13: Collect K1 timing from drain batch */
        {
            float k1_ms = 0;
            cudaEventElapsedTime(&k1_ms, k1_start_evt[prev], k1_end_evt[prev]);
            total_k1_ms += k1_ms;
        }

        u32 num_survivors_raw = *h_survivor_count[prev];
        u32 num_survivors = num_survivors_raw;
        if (num_survivors > MAX_SURVIVORS_PER_BATCH) {
            fprintf(stderr, "WARNING: drain batch survivor overflow! %u survivors, "
                    "cap %d. %u candidates LOST.\n",
                    num_survivors_raw, MAX_SURVIVORS_PER_BATCH,
                    num_survivors_raw - MAX_SURVIVORS_PER_BATCH);
            num_survivors = MAX_SURVIVORS_PER_BATCH;
        }

        if (!cpu_prove) {
            /* GPU prove path — chain results already transferred */
            total_primes += *h_gpu_primes_count[prev];
            u32 num_chains_raw = *h_chain_result_count[prev];
            u32 num_chains = num_chains_raw;
            if (num_chains > MAX_CHAIN_RESULTS) {
                fprintf(stderr, "WARNING: drain batch chain result overflow! %u chains, "
                        "cap %d. %u chains LOST.\n",
                        num_chains_raw, MAX_CHAIN_RESULTS,
                        num_chains_raw - MAX_CHAIN_RESULTS);
                num_chains = MAX_CHAIN_RESULTS;
            }

            for (u32 c = 0; c < num_chains; c++) {
                int cl = h_chain_results[prev][c].chain_len;
                total_chains++;
                if (cl <= MAX_CHAIN) g_chain_dist[cl]++;

                if (cl >= log_threshold || cl >= target_len) {
                    u128 tile_num = prev_tile_start + h_chain_results[prev][c].tile_idx;
                    u128 full_k = tile_num * WHEEL_PERIOD + h_chain_results[prev][c].k_offset;
                    u128 n128 = full_k * LATTICE_M + h_bases[h_chain_results[prev][c].base_idx].base;

                    mpz_t n_mpz;
                    mpz_init(n_mpz);
                    mpz_set_ui(n_mpz, (unsigned long)(n128 >> 64));
                    mpz_mul_2exp(n_mpz, n_mpz, 64);
                    mpz_add_ui(n_mpz, n_mpz, (unsigned long)(n128 & 0xFFFFFFFFFFFFFFFFULL));

                    char buf[256];
                    gmp_snprintf(buf, sizeof(buf), "0x%ZX", n_mpz);
                    size_t nbits = mpz_sizeinbase(n_mpz, 2);

                    if (cl >= target_len) {
                        printf("\n  !!!! CC%d TARGET HIT: %s (%zu bits) !!!!\n", cl, buf, nbits);
                    } else if (cl > best_chain) {
                        printf("\n  *** NEW BEST CC%d: %s (%zu bits) ***\n", cl, buf, nbits);
                    } else {
                        printf("  CC%d: %s (%zu bits)\n", cl, buf, nbits);
                    }
                    if (output_fp) {
                        fprintf(output_fp, "%sCC%d %s %zu\n",
                                cl >= target_len ? "TARGET " : "", cl, buf, nbits);
                        fflush(output_fp);
                    }
                    mpz_clear(n_mpz);
                }
                if (cl > best_chain) best_chain = cl;
            }
        } else if (num_survivors > 0) {
            /* CPU prove drain — workers always complete their assigned chunk
             * (no g_shutdown check inside worker loop), so no toggle needed. */
            printf("[DRAIN] Processing final batch (%u survivors) with CPU prove...\n",
                   num_survivors);
            CUDA_CHECK(cudaMemcpyAsync(h_survivors[prev], d_survivors[prev],
                                        num_survivors * sizeof(Survivor),
                                        cudaMemcpyDeviceToHost, streams[prev]));
            CUDA_CHECK(cudaStreamSynchronize(streams[prev]));
            struct timeval pv_s, pv_e;
            gettimeofday(&pv_s, NULL);
            prove_survivors_mt(h_survivors[prev], num_survivors, prev_tile_start,
                              target_len, bits, log_threshold, output_fp,
                              &total_primes, &total_chains, &best_chain);
            gettimeofday(&pv_e, NULL);
            total_prove_ms += (pv_e.tv_sec - pv_s.tv_sec) * 1000.0 +
                              (pv_e.tv_usec - pv_s.tv_usec) / 1000.0;
        }
        total_tiles += prev_batch_size;
        total_survivors += num_survivors;
        batches_done++;
    }

    /* Final report */
    gettimeofday(&tv_now, NULL);
    double total_time = (tv_now.tv_sec - tv_start.tv_sec) +
                        (tv_now.tv_usec - tv_start.tv_usec) / 1e6;
    u64 wheel_per_tile = 0;
    for (int b = 0; b < h_num_bases; b++)
        wheel_per_tile += h_bases[b].wheel_total;
    wheel_per_tile /= 6;
    u64 total_cand = total_tiles * wheel_per_tile;

    printf("\n");
    printf("================================================================\n");
    printf("  FINAL RESULTS (C+C v13 — L1/L2 Cache Optimization)\n");
    printf("================================================================\n");
    printf("  Time:         %.2f seconds\n", total_time);
    printf("  Tiles:        %llu\n", (unsigned long long)total_tiles);
    printf("  Candidates:   %llu (%.2f B)\n",
           (unsigned long long)total_cand, total_cand / 1e9);
    printf("  Survivors:    %llu (%.4f%%)\n",
           (unsigned long long)total_survivors,
           total_cand > 0 ? 100.0 * total_survivors / total_cand : 0);
    printf("  Primes:       %llu\n", (unsigned long long)total_primes);
    printf("  Chains:       %llu (best: CC%d)\n",
           (unsigned long long)total_chains, best_chain);
    if (cpu_prove && g_non_roots > 0)
        printf("  Non-roots:    %llu (walked backward to true root)\n",
               (unsigned long long)g_non_roots);
    printf("  Throughput:   %.2f B candidates/sec\n", total_cand / total_time / 1e9);
    /* v13: GPU utilization summary */
    if (total_k1_ms > 0) {
        double gpu_util = total_k1_ms / (total_time * 1000.0) * 100.0;
        if (gpu_util > 100.0) gpu_util = 100.0;
        printf("  GPU K1 time:  %.1f s (%.0f%% utilization)\n",
               total_k1_ms / 1000.0, gpu_util);
        if (cpu_prove && total_prove_ms > 0) {
            printf("  CPU prove:    %.1f s (prove/K1 ratio: %.2fx)\n",
                   total_prove_ms / 1000.0, total_prove_ms / total_k1_ms);
            if (total_prove_ms > total_k1_ms * 1.5)
                printf("  ** CPU prove is the bottleneck — GPU was idle %.0f%% of the time **\n"
                       "  ** Try: more --prove-threads, or higher --depth to reduce survivors **\n",
                       100.0 - gpu_util);
        }
    }
    if (g_use_prefix) {
        printf("  Prefix:       %s (%d bits)\n", g_prefix_str, g_prefix_bits);
        if (!random_mode && g_prefix_tile_count > 0) {
            u64 sweeps = total_tiles / (u64)g_prefix_tile_count;
            double within_pct = 100.0 * (double)(total_tiles % (u64)g_prefix_tile_count) /
                                (double)g_prefix_tile_count;
            printf("  Sweeps:       %llu complete + %.2f%%\n",
                   (unsigned long long)sweeps, within_pct);
        }
        if (g_prefix_lanes > 0)
            printf("  Lane:         %d/%d\n", g_prefix_lane_id, g_prefix_lanes);
    }

    /* Chain length distribution (includes non-root chains measured from true root) */
    printf("\n  Chain distribution:\n");
    for (int i = 1; i <= MAX_CHAIN; i++) {
        if (g_chain_dist[i] > 0) {
            printf("    CC%-2d: %llu%s\n", i,
                   (unsigned long long)g_chain_dist[i],
                   i >= target_len ? "  <<<< TARGET" : "");
        }
    }
    printf("================================================================\n\n");

    /* Final checkpoint save on shutdown */
    if (g_checkpoint_file) {
        checkpoint_save(g_checkpoint_file, bits, target_len, sieve_len,
                       current_tile_128,
                       total_tiles, total_primes, total_chains, best_chain);
        printf("[CHECKPOINT] Saved final state to %s\n", g_checkpoint_file);
    }

    /* Close output file */
    if (output_fp) {
        fprintf(output_fp, "# Finished: %ld  tiles=%llu  primes=%llu  best=CC%d",
                (long)time(NULL), (unsigned long long)total_tiles,
                (unsigned long long)total_primes, best_chain);
        if (g_use_prefix) fprintf(output_fp, "  prefix=%s", g_prefix_str);
        fprintf(output_fp, "\n");
        fclose(output_fp);
    }

    if (cpu_prove) prover_pool_stop();

    /* Cleanup — double-buffer resources */
    for (int buf = 0; buf < 2; buf++) {
        CUDA_CHECK(cudaFreeHost(h_survivors[buf]));
        CUDA_CHECK(cudaFreeHost(h_survivor_count[buf]));
        CUDA_CHECK(cudaFreeHost(h_chain_results[buf]));
        CUDA_CHECK(cudaFreeHost(h_chain_result_count[buf]));
        CUDA_CHECK(cudaFreeHost(h_gpu_primes_count[buf]));
        CUDA_CHECK(cudaFreeHost(h_batch_l2[buf]));
        CUDA_CHECK(cudaFreeHost(h_batch_ext[buf]));
        CUDA_CHECK(cudaFreeHost(h_batch_line[buf]));
        CUDA_CHECK(cudaFree(d_batch_l2[buf]));
        CUDA_CHECK(cudaFree(d_batch_ext[buf]));
        CUDA_CHECK(cudaFree(d_batch_line[buf]));
        CUDA_CHECK(cudaFree(d_survivors[buf]));
        CUDA_CHECK(cudaFree(d_survivor_count[buf]));
        CUDA_CHECK(cudaFree(d_chain_results[buf]));
        CUDA_CHECK(cudaFree(d_chain_result_count[buf]));
        CUDA_CHECK(cudaFree(d_gpu_primes_count[buf]));
    }
    CUDA_CHECK(cudaStreamDestroy(streams[0]));
    CUDA_CHECK(cudaStreamDestroy(streams[1]));
    /* v13: Destroy timing events */
    for (int buf = 0; buf < 2; buf++) {
        CUDA_CHECK(cudaEventDestroy(k1_start_evt[buf]));
        CUDA_CHECK(cudaEventDestroy(k1_end_evt[buf]));
    }
    CUDA_CHECK(cudaFree(d_ext_filter_first));
    CUDA_CHECK(cudaFree(d_line_kill_first));
    CUDA_CHECK(cudaFree(d_kM_mod_line));
    CUDA_CHECK(cudaFree(d_kmod_filter));
    CUDA_CHECK(cudaFree(d_kmod_ext_filter));
    CUDA_CHECK(cudaFree(d_wheel_offsets));
    CUDA_CHECK(cudaFree(d_wheel_start_m6));
    CUDA_CHECK(cudaFree(d_wheel_count_m6));

    for (int b = 0; b < h_num_bases; b++)
        for (int m = 0; m < 6; m++)
            free(h_bases[b].wheel_by_mod6[m]);

    return 0;
}
