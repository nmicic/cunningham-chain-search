# Master Index — Cunningham Chain Program Suite

Comprehensive audit and evaluation of all programs, focused on algorithms, features, and techniques relevant for building an optimized GMP/CUDA hybrid chain searcher.

---

## Status Key
- **Core CC**: Implements Cunningham-targeted construction/search logic (CC1, CC2, BiTwin, or CC-specific sieve/design).
- **CC-related visualization**: Not a pure CC search engine, but includes CC2 recurrence behavior in its tracked paths.
- **Non-CC**: Prime/2-adic analytics without enforcing Cunningham recurrence checks.

## Families & Overlaps
- **Constructor family**: `cc_constructor*`, `cc_cuda`, `cc_lattice_search`, `cc_line_sieve`, `cc_bitwin_prototype`, `cc_hybrid_lattice_design`.
- **Search family**: `cc_search*`, `cc_search_unified*`, `cc_search_256bit`.
- **Algebraic family**: `cunningham` (primorial-based k*P#*2^m-1 form).
- **CC visualization family**: `prime_rails`, `prime_tree_chains` (right-edge CC2 behavior).
- **Non-CC analytics family**: `Non_CC/*` folders (`cuda_prime*`, `prime_analysis`, `prime_lines`, `twoadic_cuda`).

---

## Grand Feature Comparison Matrix

### Core CC Programs — Constructor Family

| Feature | cc_constructor v20 | cc_alg_hybrid v05b | cc_hybrid v07 | cc_cuda v10 | cc_lattice_search | cc_line_sieve | cc_bitwin_proto |
|---|---|---|---|---|---|---|---|
| **LOC** | ~1481 | ~963 | ~889 | ~1183 | ~951 | ~304 | ~725 |
| **Arithmetic** | 128-bit + GMP | 128-bit | 128-bit | 128-bit | 128-bit | 64-bit only | 128-bit |
| **MR Witnesses (64b)** | 7 (det.) | 7 (det.) | 7 (det.) | 7 (det.) | 7 (det.) | stub | 7 (det.) |
| **MR Witnesses (128b)** | 12 (prob.) | 4 (prob.) | 4 + trial div | 4 + trial div | 7 (prob.) | N/A | 16 (prob.) |
| **GMP Deep Verify** | Yes | No | No | No | No | No | No |
| **CRT Lattice** | 25 primes (3-101) | 6 primes | 6 primes | 6 primes | Grid-based | No | 6 primes |
| **Wheel Opt** | No | No | No | 3 primes (23,29,31) | No | No | No |
| **Bitmask Filter** | Periodic tables (25 primes) | 18 primes | 10 primes (<=61) | 7 primes | Bitmap prefilter | Line sieve (14 primes) | Quick sieve |
| **Chain Type** | CC1 only | CC1 only | CC1 only | CC1 only | CC1 only | CC1 only | **BiTwin (CC1+CC2)** |
| **k-value Tracking** | **Yes** | No | No | No | No | No | No |
| **Backward Walk** | No | No | **Yes** | **Yes** | No | No | No |
| **Checkpoint/Resume** | No | No | No | **Yes** | No | No | No |
| **File Output** | promising_k1.txt | No | Optional | Optional | No | No | No |
| **Search Mode** | Lattice+2-adic | CRT+primorial+2-adic | Random Monte Carlo | Sequential+Random | Grid bitmap | Brute-force | CRT lattice |
| **Completeness** | Complete | Complete | Complete | **Production** | Complete | **Prototype** | Complete |

### Core CC Programs — Search Family

| Feature | cc_search v7 | cc_search_256bit v3 | cc_search_cuda v1 | cc_search_merged v6 | cc_search_unified v10 | cc_unified_c7 v25/v27 |
|---|---|---|---|---|---|---|
| **LOC** | ~365 | ~1060 | ~490 | ~502 | ~2475 | **~4217 / ~4598** |
| **Arithmetic** | __int128 | Custom 256-bit + Montgomery | __int128 | __int128 | __int128 + binary mulmod | + fast-path mulmod |
| **Max bits** | 64 | **256** | 64 | 64 | 128 | 128 |
| **MR Witnesses (64b)** | 12 (det.) | 12 (det.) | 12 (det.) | 12 (det.) | **7 (proven det.)** | **7 (proven det.)** |
| **MR Witnesses (128b)** | N/A | 12 + Montgomery | N/A | N/A | **BPSW** | **BPSW** |
| **Sieve** | Death tables | Bitwise div + pow tables | Death tables | Death tables | Death residues | **Bitmask sieve + sieve-only line filter (v27)** |
| **Chain Types** | CC1 only | CC1 only | CC1 only | CC1 only | **CC1 + CC2** | **CC1 + CC2** |
| **Warp Reduction** | No | **Yes** | No (shared mem) | No | **Yes** | **Yes** |
| **Batch Processing** | No | No | No | No | **Yes (2^22)** | **Yes (2^22)** |
| **Batch Ordering** | N/A | N/A | N/A | N/A | Seq/Random | **Seq/Random/2-adic/Deep** |
| **CSV Logging** | No | No | No | No | No | **Yes (heatmap)** |
| **Checkpoint** | No | No | No | No | **Yes** | **Yes (timestamped)** |
| **Self-Tests** | No | Basic | Verify known | No | **5 levels** | **5 levels + GPU sim** |
| **CPU Verify** | No | No | Yes | No | **BPSW re-verify** | **BPSW re-verify** |
| **Signal Handling** | No | Yes | No | No | **Yes** | **Yes** |
| **Recommend Engine** | No | No | No | No | No | **Yes** |
| **Completeness** | Complete | Complete | Fixed-range | Bridge | **Production** | **Most mature (v27: +lattice constructor)** |

### Algebraic + Visualization + Non-CC

| Feature | cunningham_alg v3 | prime_rails v1-2 | prime_tree v1 | cuda_prime_256 | cuda_prime_analysis | prime_analysis | prime_lines v14 | twoadic_cuda |
|---|---|---|---|---|---|---|---|---|
| **LOC** | 727 | 774 | 409 | 620 | 620 | 1118 | 845 | 592 |
| **Arithmetic** | 64+128 | 64 | 64 trial div | **256-bit Montgomery** | 64 | 64+256 | 64 | 32-bit |
| **Chain Logic** | CC1 + 128b follow | CC1 + factor analysis | CC1 L/R/MIXED tree | None | None | None | None | None |
| **Sieve** | Death residues, mod6 | Bitwise div 3/5/7/11/13 | Eratosthenes | Trial div | None | None | Per-number MR | None |
| **CSV Output** | No | 2 files | 2 files | Yes | 22-column | Dual format | 3 files | No |
| **Checkpoint** | No | No | No | No | No | No | **Binary** | No |
| **Status** | Core CC | CC-related | CC-viz | Non-CC | Non-CC | Non-CC | Non-CC | Non-CC |

---

## Detailed Per-Program Evaluations

---

### 1. cc_constructor v20 — Score: 9.2/10

**Path**: `cc_constructor/00_cc_constructor_v20__20260122-035258.cu`

**Algorithm**: CRT lattice construction using 25 primes (3-101). Computes forbidden residues for each chain depth: `P mod q == (2^{-k}-1) mod q` means k-th element is divisible by q. Single-residue primes combined via CRT into base/modulus pair. Constructs candidates of form `(prefix << k) | ((1 << k) - 1)` — terminals with k trailing 1-bits. Three search modes: standard CRT, 2-adic subtree, 2-adic bit-pattern with CRT intersection.

**Key Strengths**:
- **k-value tracking**: Distinguishes true CC roots (k=1) from fragments (k>1) — critical for CC20 hunting since only k=1 chains can be true CC20s
- **Periodic forbidden residue tables**: Precomputed in `__constant__` memory using multiplicative orders of 2 mod q, eliminating per-step computation. 25 arrays for primes 7-107
- **GMP deep analysis**: Extends chains beyond 128-bit using GMP with 25-50 rounds, configurable threshold
- **2-adic pattern intersection with CRT**: Computes CRT lattice points within specified binary prefix range

**Weaknesses**: CC1-only, no checkpoint/resume, no CSV output.

**Unique technique for GMP version**: The periodic forbidden residue table approach (precomputing `(2^{-i}-1) mod p` for each prime using multiplicative orders) is highly efficient and directly portable to GMP.

---

### 2. cc_constructor_algebraic_hybrid v05b — Score: 8.4/10

**Path**: `cc_constructor_algebraic_hybrid/00_cc_constructor_algebraic_hybrid_v05b.cu`

**Algorithm**: Three-strategy hybrid: (1) CRT line-lattice with 6 primes producing valid bases mod M=1,616,615; (2) Primorial construction `k = a * P# * 2^s` where P# is configurable up to 31#; (3) 2-adic exploration of the `a` parameter space using random anchor + neighborhood. Candidate formula: `base + k * M`.

**Key Strengths**:
- **Primorial-structured k values**: Directly mimics historical record-finding methodology on GPU
- **Odd-only a-value fix** (v05b): Eliminates redundant search caused by even `a` values — since (a, s) and (a/2, s+1) give identical k
- **Three-way strategy fusion**: CRT + primorial + 2-adic is unique

**Weaknesses**: Only 4 MR witnesses for 128-bit, no file output, no checkpoint.

**Unique technique for GMP version**: The algebraic form `k * P# * 2^m - 1` with parameterized primorial is the classic record-finding approach. The odd-only deduplication fix is important.

---

### 3. cc_constructor_hybrid v07 — Score: 8.1/10

**Path**: `cc_constructor_hybrid/00_cc_constructor_hybrid_v07_1.cu`

**Algorithm**: CRT lattice + random Monte Carlo sampling. Generates random k values and base indices, forming candidates `n = base + k * M`. Optional prefix steering constrains k range for bit-pattern targeting.

**Key Strengths**:
- **Backward chain walking**: Finds true root of any prime found, even non-roots, and tracks full chain lengths
- **Non-root statistics**: Counts non-roots and their full chain lengths — useful for search coverage analysis
- **Conservative bitmask safety**: Limits bitmask primes to <=61 to prevent undefined behavior on 64-bit shifts
- **Extended trial division**: 62 primes 67-409 (v07 fix)

**Weaknesses**: Random-only (no exhaustive mode), limited bitmask primes.

**Unique technique for GMP version**: Backward walking to find true roots + non-root statistics give valuable data about search completeness.

---

### 4. cc_cuda v10 tuned — Score: 7.8/10

**Path**: `cc_cuda/00_cc_cuda_v10_tuned.cu`

**Algorithm**: CRT lattice + per-base wheel optimization. Innovation: 3 additional "wheel primes" {23, 29, 31} with period W=20,677. For each CRT base, precomputes which k offsets within each wheel period are valid. Supports sequential (exhaustive) and random modes.

**Key Strengths**:
- **Per-base wheel optimization**: Precomputes valid k offsets per CRT base, integrating CRT + wheel + mod6 into one step — most sophisticated filtering in constructor family
- **Three-tier filtering** (L1 wheel / L2 bitmask / L3 MR): Achieves 5-6x improvement over v08
- **Integrated mod6 filtering**: `target_k_mod6` eliminates runtime mod6 check
- **Tile-aligned sequential search**: Cache-friendly access patterns
- **Checkpoint/resume**: The only constructor with this feature

**Weaknesses**: CC1-only, no GMP integration.

**Unique technique for GMP version**: The three-tier filtering architecture (CRT → wheel → bitmask → MR) is the most throughput-optimized approach and should inform the GMP candidate generation pipeline.

---

### 5. cc_hybrid_lattice_design — Score: 4.0/10 (design doc)

**Path**: `cc_hybrid_lattice_design/00_cc_hybrid_lattice_design__20260122-124303.cu`

**Algorithm**: Design specification, not compilable. Describes root-based lattice vs terminal-based lattice.

**Key Insight**: `cc_constructor`'s form `a * 2^k - 1` optimizes for terminals with k trailing 1-bits, but CC17's root has only 4 trailing 1s while the terminal has 20. A root-based lattice would catch all chains regardless of trailing 1-bit structure.

**Unique technique for GMP version**: The root vs terminal insight is foundational — GMP version should search for roots, not terminals.

---

### 6. cc_lattice_search — Score: 7.4/10

**Path**: `cc_lattice_search/00_cc_lattice_search_1.cu`

**Algorithm**: 2-adic grid sieve with precomputed bitmap. Constructs a 2D grid of size D x D (default 1024x1024). For each small prime q and chain position k, marks forbidden residues. Builds "CC-possible" bitmap by walking each L-step path. At search time, a single bit lookup determines if a candidate is worth testing.

**Key Strengths**:
- **Grid-based chain path feasibility**: Models entire chain path on finite grid, precomputing whether a CC of length L can exist at each position
- **Structural prime exclusion**: Small primes eliminated by grid structure itself — no runtime divisibility check
- **Configurable grid resolution**: `d_bits` trades memory for filter precision (1M to 16M cells)

**Weaknesses**: High memory usage for large grids, sequential candidate enumeration.

**Unique technique for GMP version**: The positional sieve / bitmap prefilter concept — precompute which low-bit patterns can start a CC of target length, then filter candidates by low bits before expensive GMP primality testing.

---

### 7. cc_line_sieve — Score: 3.5/10 (incomplete prototype)

**Path**: `cc_line_sieve/00_cc_line_sieve.cu`

**Algorithm**: Brute-force line extension sieve. For each candidate, sequentially checks chain elements for divisibility by small primes 3-47. Simplest possible approach.

**Key Strengths**:
- **Progressive sieve with bucketing**: Counts how many chain steps survive and buckets candidates — useful for analyzing "almost good enough" candidates

**Weaknesses**: Miller-Rabin is a stub (`return false`), 64-bit only, tiny test range, no CLI.

**Unique technique for GMP version**: The progressive bucketing idea (tracking how far a candidate gets before failing) is useful for diagnostics.

---

### 8. cc_bitwin_prototype — Score: 5.9/10

**Path**: `cc_bitwin_prototype/00_cc_bitwin_prototype__20260121-023053.cu`

**Algorithm**: CRT lattice for BiTwin chains. For even center n, both n-1 (CC1) and n+1 (CC2) must start chains. BiTwin length = min(first_kind_len, second_kind_len). Combined constraint system for both chain kinds with evenness constraint.

**Key Strengths**:
- **Only BiTwin searcher** in the entire family
- **Combined constraint system**: Forbidden masks combine CC1, CC2, and evenness — mathematically correct formulation
- **Dual chain following**: Simultaneously follows 2p+1 and 2p-1 chains
- **2-adic spine reference**: Computes spine numbers `(4^d-1)/3`
- **16 MR witnesses for 128-bit**: Most thorough probabilistic testing

**Weaknesses**: Basic filtering (no bitmask line filter), potential intermediate overflow in CRT combination.

**Unique technique for GMP version**: If BiTwin search is desired, this is the only reference implementation. The combined constraint system formulation is directly portable.

---

### 9. cc_search v7 — Score: 6.6/10

**Path**: `cc_search/00_cc_search_v7__20260118-100858.cu`

**Algorithm**: Brute-force sieve-then-test over single shell `[2^(N-1), 2^N)`. Death-table lookahead prefilter. Dual-target build (CUDA + CPU).

**Key Strengths**:
- **Death tables in constant memory**: Precomputed steps-before-death for each prime/residue pair
- **Dual CUDA/CPU build**: Single `.cu` source compiles for either target via preprocessor macros

**Weaknesses**: No batching, no checkpoint, CC1 only, per-candidate atomics (contention).

**Unique technique for GMP version**: Death table concept — precompute for each small prime p and residue r, how many chain steps survive before hitting a multiple of p.

---

### 10. cc_search_256bit v3 — Score: 8.6/10

**Path**: `cc_search_256bit/00_cc_search_256bit_v3.cu`

**Algorithm**: Sieve-then-test for candidates up to 256 bits using custom 256-bit integer type with Montgomery multiplication.

**Key Strengths**:
- **Division-free 256-bit composite filter**: popcount tricks for mod 3, byte-sum for mod 5, chunk-sum for mod 7/11/13, precomputed power tables for 17-47
- **Full Montgomery pipeline**: R^2 precomputation, stays in Montgomery form during all MR squarings
- **`__umul64hi` PTX intrinsic**: Avoids `__int128` dependency
- **Warp-level reduction**: Dramatically reduces atomic contention
- **Bidirectional chain trace**: Walks both backward (find root) and forward (find tip)

**Weaknesses**: R^2 computation (512 doublings per candidate) is a known bottleneck, CC1 only, no checkpoint.

**Unique technique for GMP version**: The division-free composite filter (popcount/chunk-sum tricks) is valuable even in GMP context for quick pre-screening. The death residue survival check for 256-bit numbers shows how to scale the technique.

---

### 11. cc_search_cuda v1 — Score: 6.2/10

**Path**: `cc_search_cuda/00_cc_search_cuda_1__20260117-230736.cu`

**Algorithm**: CUDA-only search over `[2^31, 2^32)` with shared-memory result buffering.

**Key Strengths**:
- **Shared-memory result buffering**: Block-local results collected before single atomic to global — reduces contention
- **CPU verification pass**: GPU results independently re-verified on CPU

**Weaknesses**: Hardcoded to 32-bit range, no CLI range selection.

---

### 12. cc_search_merged v6 — Score: 6.5/10

**Path**: `cc_search_merged/00_cc_search_merged_v6__20260118-013418.cu`

**Algorithm**: Bridge version between cc_search and cc_search_unified. Functionally equivalent to cc_search v7 with slightly different constant-memory management.

**Weaknesses**: No novel features beyond cc_search v7.

---

### 13. cc_search_unified v10 — Score: 8.8/10

**Path**: `cc_search_unified/00_cc_search_unified_v10.cu`

**Algorithm**: Most architecturally mature 64/128-bit search engine. Batch-based (2^22 candidates/batch) with configurable ordering.

**Key Strengths**:
- **BPSW for 128-bit**: Full Baillie-PSW (base-2 MR + Lucas probable prime test with Selfridge Method A). Gold standard — no known counterexamples
- **Proven 7-witness deterministic MR** for 64-bit: `{2, 325, 9375, 28178, 450775, 9780504, 1795265022}` — proven deterministic for all n < 2^64
- **Both chain kinds** via `--kind` parameter (CC1: 6k+5, CC2: 6k+1)
- **Packed best-result tracking**: Single `atomicMax` on packed u64 (8-bit length + 56-bit inverted start)
- **Extensive test suite**: 5 levels (SANITY through STRESS), CC5-CC19 test vectors, shell boundary tests
- **Checkpoint with signal-safe shutdown**
- **Warp-level reduction**: Two-phase (warp → shared → block → global)

**Weaknesses**: No CSV logging, no 2-adic ordering, 128-bit max.

**Unique technique for GMP version**: BPSW test implementation is the gold standard and should be the primality test used in GMP version. The 7-witness deterministic MR for 64-bit pre-screening is also optimal.

---

### 14. cc_search_unified_c7 v25/v27 — Score: 9.5/10 (MOST ADVANCED)

**Baseline**: `cc_search_unified_c7/00_cc_search_unified_c7_v25.cu` (~4217 lines)
**Latest production**: `cc_search_unified_c7/01_cc_search_unified_c7_v27__20260122-121913.cu` (~4598 lines)

**Algorithm**: Fork of cc_search_unified, specialized for CC7+ hunting. All features of v10 plus:

**Key Strengths (v25 baseline)**:
- **Pre-computed forbidden residue bitmasks** (v25): For each small prime q and chain length L, precomputes which residues `(n mod q)` kill any chain element. Stored as 64-bit bitmasks. **One mod + one bit-test per prime replaces L mod operations per prime.** Expected 3-8x faster sieving.
- **Hardcoded sieve checks**: Individual `if (n % 5 == 0)` etc. allows compiler magic-number multiplication (~10x faster than loop-based modulo)
- **Five batch-ordering modes**: Sequential, random, 2-adic (Roots First), 2-adic-reverse (Spines First), 2-adic-deep (random jumping + neighborhood)
- **CSV heatmap logging**: Per-batch chain statistics with 2-adic metadata for visualization
- **Timestamped output files**: Prevents overwrites across runs
- **CC7+ independent tracking**: Always tracks best long chain regardless of `--min`
- **GPU batch simulation testing** (v23): End-to-end kernel correctness testing
- **Recommendation engine**: Outputs optimal `--exp` for target chain lengths
- **Three search range modes**: Default (chains end < 2^exp), shell, full
- **Optimized 128-bit mulmod**: Fast-path when both operands fit in 64 bits

**Key Strengths (v27 additions)**:
- **Sieve-only line filter** (`METHOD-L17`): Tests chain p → 2p+1 → 4p+3 → ... using ONLY sieve checks — no MR until candidate survives all filter steps. 10-50x faster than v26's MR-based quick_prime filter.
- **Root-based lattice constructor** (`METHOD-L18`, `--lattice`): CRT enumeration of valid root residues across sieve primes. Generates candidates by construction instead of testing all p ≡ 5 (mod 6). Key insight: the line filter doubles as a generator.
- **`--line-filter N`**: Configurable sieve-only line filter depth (default 5).

**Weaknesses**: 128-bit max (no 256-bit), CC1+CC2 but no BiTwin.

**Unique technique for GMP version**: The forbidden residue bitmask approach (v25) is THE key optimization to port. v27's sieve-only line filter and root-based lattice constructor add two more high-value techniques — the lattice constructor is especially relevant since it implements root-based (not terminal-based) candidate generation via CRT. Also the 2-adic batch ordering strategies are worth exploring for targeted GMP searches.

---

### 15. cunningham_algebraic v3 — Score: 6.9/10

**Path**: `cunningham/00_cunningham_algebraic_v3.cu`

**Algorithm**: Searches for CC1 starting primes of form `k * P# * 2^m - 1`. Iterates over (m, k) pairs.

**Key Strengths**:
- **Death residue multi-step lookahead**: For 14 small primes, checks up to 10 future chain steps for guaranteed compositeness
- **Algebraic form**: `k * P# * 2^m - 1` guarantees small-factor properties from primorial
- **128-bit chain following**: Dispatches to fast 64-bit when value fits

**Weaknesses**: No CRT lattice, per-m kernel launch overhead, no checkpoint.

**Unique technique for GMP version**: The algebraic form with primorial is the classic record-finding strategy and should be a candidate generation option.

---

### 16. prime_rails v1-2 — Score: CC visualization tool

**Path**: `prime_rails/00_prime_rails_v1-2__20260117-175813.cu`

**Algorithm**: For every prime in a shell, follows CC1 chain, records chain length, smallest factor of interrupting composite, verifies mod-4 assertion.

**Unique value**: Interrupt factor analysis (which small primes break chains at which lengths) provides statistical insight into chain-ending patterns.

---

### 17. prime_tree_chains — Score: CC visualization tool

**Path**: `prime_tree_chains/00_prime_tree_chains__20260117-181212.cu`

**Algorithm**: Models numbers as binary tree (left=2n, right=2n+1). Traces prime chains with L/R/MIXED classification. RIGHT-only = CC1.

**Unique value**: Binary tree perspective clarifies CC structure. Sieve-based (Eratosthenes), limited to sieve range.

---

### 18-22. Non-CC Analytics Tools

| Program | Purpose | Key Value |
|---|---|---|
| **cuda_prime_256** | 256-bit Montgomery prime segment analyzer | Production-quality 256-bit Montgomery implementation on GPU |
| **cuda_prime_analysis_64** | 64-bit range prime statistics exporter | CAS-based atomic double add, 22-column rich CSV |
| **prime_analysis** | Unified 64/256-bit prime analyzer | Merged multi-precision tool with auto-detection |
| **prime_lines v14** | Dual-mode prime line tracker | Octagonal coordinate + angular ray analysis, binary checkpoint |
| **twoadic_cuda** | 2-adic structural analysis toolkit | GPU hash table for odd-core clustering, anomaly scoring |

---

## Maturity Ranking (re-evaluated)

### Tier 1 — Production Quality
| Rank | Program | Score | Why |
|---|---|---|---|
| 1 | **cc_search_unified_c7 v25/v27** | 9.5 | Most features, bitmask sieve, 5 batch modes, BPSW, CSV, checkpoint, tests, recommendations. v27 adds sieve-only line filter + lattice constructor. |
| 2 | **cc_search_unified v10** | 8.8 | BPSW, proven 7-witness MR, dual chain kinds, extensive tests, checkpoint |
| 3 | **cc_search_256bit v3** | 8.6 | 256-bit Montgomery, division-free filters, warp reduction, bidirectional trace |

### Tier 2 — Complete & Innovative
| Rank | Program | Score | Why |
|---|---|---|---|
| 4 | **cc_constructor v20** | 9.2 | k-value tracking, periodic tables, GMP deep analysis, 2-adic pattern intersection |
| 5 | **cc_constructor_algebraic_hybrid v05b** | 8.4 | Three-strategy fusion, primorial-structured k, odd dedup |
| 6 | **cc_constructor_hybrid v07** | 8.1 | Backward walking, non-root statistics, random Monte Carlo |
| 7 | **cc_cuda v10** | 7.8 | Three-tier filtering, wheel optimization, checkpoint |
| 8 | **cc_lattice_search** | 7.4 | Grid bitmap prefilter, structural prime exclusion |

### Tier 3 — Specialized / Earlier Versions
| Rank | Program | Score | Why |
|---|---|---|---|
| 9 | **cunningham_algebraic v3** | 6.9 | Classic k*P#*2^m-1 form, death residue lookahead |
| 10 | **cc_search v7** | 6.6 | Foundation death tables, dual CUDA/CPU build |
| 11 | **cc_search_merged v6** | 6.5 | Bridge version |
| 12 | **cc_search_cuda v1** | 6.2 | Shared-memory result buffering |
| 13 | **cc_bitwin_prototype** | 5.9 | Only BiTwin searcher, dual chain constraint system |

### Tier 4 — Prototypes / Design
| Rank | Program | Score | Why |
|---|---|---|---|
| 14 | **cc_hybrid_lattice_design** | 4.0 | Design doc only, but root-vs-terminal insight is valuable |
| 15 | **cc_line_sieve** | 3.5 | Incomplete prototype, MR stub |

---

## Technique Inventory — What to Port to GMP Version

### MUST-HAVE (present in best programs)

1. **Forbidden residue bitmask sieve** (cc_unified_c7 v25): For each small prime q and target chain length L, precompute 64-bit bitmask of residues `(n mod q)` that kill any chain element. One mod + one bit-test per prime. This is the single most impactful optimization.

1b. **Sieve-only line filter** (cc_unified_c7 v27, `METHOD-L17`): Tests chain elements p → 2p+1 → 4p+3 → ... using ONLY sieve checks — no MR until candidate survives all filter steps. 10-50x faster than MR-based filtering.

1c. **Root-based lattice constructor** (cc_unified_c7 v27, `METHOD-L18`): CRT enumeration of valid root residues across sieve primes. Generates candidates by construction instead of testing all p ≡ 5 (mod 6). The line filter doubles as a generator — root-based (not terminal-based) lattice search.

2. **BPSW primality test** (cc_unified v10, c7 v25): Base-2 Miller-Rabin + Lucas test with Selfridge Method A. No known counterexamples. In GMP: use `mpz_probab_prime_p` or implement BPSW directly.

3. **7-witness deterministic MR for 64-bit** (cc_unified v10+): Witnesses `{2, 325, 9375, 28178, 450775, 9780504, 1795265022}` are proven deterministic for all n < 2^64. Use for any 64-bit pre-screening.

4. **k-value tracking** (cc_constructor v20): Distinguish true chain roots (k=1) from fragments. Only k=1 chains can be full-length Cunningham chains.

5. **Death residue lookahead** (cc_search family, cunningham_algebraic): For each small prime p, precompute which residues at each chain step guarantee compositeness. Reject candidates that cannot sustain target length.

### SHOULD-HAVE (significant speedup)

6. **CRT lattice construction** (constructor family): Combine forbidden residue constraints via CRT. Walk only lattice points guaranteed to survive initial small-prime checks.

7. **Per-base wheel optimization** (cc_cuda v10): For additional primes beyond CRT, precompute valid k offsets within each wheel period. Integrated mod6 filtering.

8. **Grid bitmap prefilter** (cc_lattice_search): Precompute bitmap of which low-bit patterns can start a CC of target length. Single bit lookup before expensive tests.

9. **Periodic forbidden residue tables** (cc_constructor v20): Using multiplicative orders of 2 mod q, precompute periodic tables in constant memory. Eliminates per-step computation.

10. **Algebraic form candidates** (cunningham_algebraic, cc_alg_hybrid): Generate candidates as `k * P# * 2^m - 1` to leverage primorial divisibility properties.

### NICE-TO-HAVE (diagnostics / special cases)

11. **Backward chain walking** (cc_hybrid v07, cc_cuda v10): Find true root of any discovered prime. Track full chain lengths even from non-roots.

12. **BiTwin constraint system** (cc_bitwin_prototype): Combined CC1+CC2 forbidden residues for BiTwin search.

13. **2-adic batch ordering** (cc_unified_c7 v25): Search by 2-adic valuation order for structured exploration.

14. **Division-free composite filters** (cc_search_256bit): Popcount tricks for mod 3, byte-sum for mod 5, chunk-sum for 7/11/13.

15. **Progressive sieve bucketing** (cc_line_sieve): Track how far candidates get before failing — useful for diagnostics.

### TECHNIQUES PRESENT IN GMP BUT NOT IN ANY CUDA PROGRAM

Check your GMP version for these — if missing, they are unique GMP advantages:
- Arbitrary precision (no 128/256-bit ceiling)
- Exact primality via mpz_probab_prime_p with high rounds
- GMP's highly optimized mulmod/powmod
- Potential for primorial sieve of Eratosthenes at GMP scale

---

## Architectural Lineage

```
cc_line_sieve (prototype, brute-force line sieve)
    │
    ├── cc_hybrid_lattice_design (design doc: root vs terminal insight)
    │       │
    │       ├── cc_constructor_hybrid v07 (root-based CRT + random)
    │       │       │
    │       │       ├── cc_constructor_algebraic_hybrid v05b (+ primorial k)
    │       │       │
    │       │       └── cc_cuda v10 (+ wheel optimization)
    │       │
    │       └── cc_constructor v20 (prefix-based + periodic tables + GMP)
    │
    ├── cc_lattice_search (grid bitmap prefilter)
    │
    └── cc_bitwin_prototype (BiTwin extension)

cc_search_cuda v1 (pioneer CUDA search)
    │
    ├── cc_search v7 (+ death tables, dual CUDA/CPU)
    │       │
    │       └── cc_search_merged v6 (bridge)
    │
    └── cc_search_unified v10 (batch processing, BPSW, dual kinds)
            │
            └── cc_search_unified_c7 v25 (bitmask sieve, 2-adic ordering, CSV)
                    │
                    └── cc_search_unified_c7 v27 (+ sieve-only line filter, lattice constructor)

cc_search_256bit v3 (256-bit Montgomery, independent branch)

cunningham_algebraic v3 (algebraic form k*P#*2^m-1, independent)
```

---

## High-Similarity Clusters
- `cc_constructor`, `cc_constructor_hybrid`, `cc_constructor_algebraic_hybrid`, `cc_cuda`, `cc_lattice_search`, `cc_line_sieve` are the same CC-construction family with different heuristics.
- `cc_search`, `cc_search_merged`, `cc_search_unified`, `cc_search_unified_c7`, `cc_search_256bit` are the same CC-search family across maturity/bit-width specialization.
- `prime_rails` and `prime_tree_chains` are visualization-heavy tools that still expose CC2 behavior.
- `Non_CC/*` folders should stay separate from CC search/constructor workflows.

---

## Folder-Level README Links
- cc_bitwin_prototype: `cc_bitwin_prototype/README.md`
- cc_constructor: `cc_constructor/README.md`
- cc_constructor_algebraic_hybrid: `cc_constructor_algebraic_hybrid/README.md`
- cc_constructor_hybrid: `cc_constructor_hybrid/README.md`
- cc_cuda: `cc_cuda/README.md`
- cc_hybrid_lattice_design: `cc_hybrid_lattice_design/README.md`
- cc_lattice_search: `cc_lattice_search/README.md`
- cc_line_sieve: `cc_line_sieve/README.md`
- cc_search: `cc_search/README.md`
- cc_search_256bit: `cc_search_256bit/README.md`
- cc_search_cuda: `cc_search_cuda/README.md`
- cc_search_merged: `cc_search_merged/README.md`
- cc_search_unified: `cc_search_unified/README.md`
- cc_search_unified_c7: `cc_search_unified_c7/README.md`
- cunningham: `cunningham/README.md`
- prime_rails: `prime_rails/README.md`
- prime_tree_chains: `prime_tree_chains/README.md`
- Non_CC/cuda_prime: `Non_CC/cuda_prime/README.md`
- Non_CC/cuda_prime_analysis: `Non_CC/cuda_prime_analysis/README.md`
- Non_CC/prime_analysis: `Non_CC/prime_analysis/README.md`
- Non_CC/prime_lines: `Non_CC/prime_lines/README.md`
- Non_CC/twoadic_cuda: `Non_CC/twoadic_cuda/README.md`
