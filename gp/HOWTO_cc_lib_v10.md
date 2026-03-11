# HOWTO: cc_lib_v10.gp — Interactive Cunningham Chain Toolkit

A practical guide to the PARI/GP library for Cunningham Chain exploration,
analysis, search, and safe prime generation. Runs anywhere GP runs: Linux,
macOS, iSH (iPhone).

---

## Quick Start

```bash
# Launch GP and load the library
gp -q
\r cc_lib_v10.gp
```

On load, 37 self-tests run automatically. You should see:

```
  All self-tests passed (37 tests). Ready.
```

**Tip for phone/iSH:** The library fits in ~64 MB stack. If you hit stack
errors on big primes, increase with `default(parisizemax, 512000000)` (done
automatically on load).

**Minimal smoke check:**

```gp
cc_chain_test(1122659)      \\ should return 7
cc_is_root(1122659)         \\ should return 1
cc_sp_verify(11)            \\ should return [11, 23, 5]
```

---

## Table of Contents

1. [Analyzing a Known Chain](#1-analyzing-a-known-chain)
2. [Understanding Chain Breakers](#2-understanding-chain-breakers)
3. [Shadow Analysis — Why Chains Break](#3-shadow-analysis--why-chains-break)
4. [Residue Forensics — Kill Positions](#4-residue-forensics--kill-positions)
5. [Searching for Chains (Layer 7 — CRT)](#5-searching-for-chains-layer-7--crt)
6. [Constructing Chains (Layer 8 — Primorial)](#6-constructing-chains-layer-8--primorial)
7. [Comparing CRT vs Constructor Search](#7-comparing-crt-vs-constructor-search)
8. [Planning a GPU Campaign with cc_con_info / cc_search_info](#8-planning-a-gpu-campaign)
9. [Analyzing GPU Hits — Post-Discovery Workflow](#9-analyzing-gpu-hits--post-discovery-workflow)
10. [Target Parameter Tuning](#10-target-parameter-tuning)
11. [Safe Prime Search (Layer 9)](#11-safe-prime-search-layer-9)
12. [Function Reference](#12-function-reference)
12b. [OpenSSL Comparison (Layer 9b)](#12b-openssl-comparison-layer-9b)
13. [Alternative Algorithms — Layer 10 (cc_x_ prefix)](#13-alternative-algorithms--layer-10)
14. [Combining Layer 10 with Other Layers](#14-combining-layer-10-with-other-layers)

---

## 1. Analyzing a Known Chain

The CC7 root `1122659` is a classic test case (7 consecutive primes in the
2p+1 chain).

```gp
\\ Walk the chain — shows every member, residues, primality
cc_walk(1122659)
```

Output shows:
- Root, chain length, back-steps from input to root
- Each member with bit-length, popcount, NAF weight, trailing ones
- Residues mod {3,5,7,11,13,17,19,23} — bracketed `[r]` means immune
- Post-break lookahead (6 composites after the chain ends)

**Start from any member — the library finds the root automatically:**

```gp
\\ Give it the 4th member, it walks back to find root 1122659
cc_walk(2*2*2*(1122659+1)-1)
```

**Silent walk (returns data, no display):**

```gp
v = cc_walk_quiet(1122659)
\\ return shape: [chain_vec, chain_len, root, breaker]
v[2]                 \\ chain length
v[3]                 \\ root
v[4]                 \\ first composite after the chain
```

**Exact-start forward walk (do not backtrack to root):**

```gp
m = 2*(1122659 + 1) - 1;    \\ member at position 1 of the chain
w = cc_walk_from_quiet(m)
\\ return shape: [chain_vec, chain_len, breaker]
w[2]                        \\ run from this exact member onward
```

### Exact Input vs Normalized Subject

This distinction matters when your input is not the true root.

```gp
\\ A non-root chain member
m = 2*(1122659 + 1) - 1;

\\ Exact-input analysis: treat m as the subject
cc_candidate_profile(m, 1, 0)

\\ Normalized analysis: walk back to the true root first
cc_candidate_profile(m, 1, 1)
```

---

## 2. Understanding Chain Breakers

When a chain breaks, you want to know *why*. `cc_autopsy` dissects the
breaker composite.

```gp
cc_autopsy(1122659)
```

---

## 3. Shadow Analysis — Why Chains Break

```gp
cc_shadow(1122659)
```

Output table shows for each of the 24 analysis primes {3..97}:
- Residue: `p mod q`
- Immune: YES if `p mod q == q-1` (first kind)
- Kill position: the chain position where `q` first divides a member

---

## 4. Residue Forensics — Kill Positions

```gp
cc_kill_pos(3, 7)           \\ where does 7 kill the chain?
cc_forbidden_res(7, 18)     \\ residues mod 7 that kill before pos 18
cc_valid_res(7, 18)         \\ residues that survive to pos 18
cc_immune_residue(7)        \\ = 6  (first kind)
```

---

## 5. Searching for Chains (Layer 7 — CRT)

```gp
cc_search(40, 5)                    \\ CC5+ at 40 bits
cc_search(60, 7, 1, 5, 3)          \\ CC7+ at 60 bits, prefix 0b101
cc_qsearch(35, 4)                  \\ quick search (50K candidates)
cc_search_info(89, 18)             \\ search space analysis
```

---

## 6. Constructing Chains (Layer 8 — Primorial)

```gp
cc_construct(35, 4, 13)            \\ CC4+ at 35 bits
cc_construct(50, 5, 19)            \\ CC5+ at 50 bits
cc_construct(89, 10, 31)           \\ CC10+ at 89 bits
cc_con_info(89, 18)                \\ constructor parameters
```

---

## 7. Comparing CRT vs Constructor Search

| Feature | Layer 7 (CRT) | Layer 8 (Constructor) |
|---------|---------------|----------------------|
| Candidate form | `n = k*M + base` | `p = S*R - 1` |
| Immunization | CRT bases avoid forbidden residues | S = primorial gives automatic immunity |
| Best for | Matching GPU pipeline | Algebraic exploration, phone |

---

## 8. Planning a GPU Campaign

```gp
cc_search_info(89, 18)
cc_con_info(89, 18)
```

---

## 9. Analyzing GPU Hits — Post-Discovery Workflow

```gp
p = 0x4300E6A5560F71F54B7E593;    \\ paste hex root from GPU log
cc_walk(p)
cc_shadow(p)
cc_autopsy(p)
cc_chain_verify(p)
```

---

## 10. Target Parameter Tuning

```gp
for(t = 18, 23, \
  my(BM = cc_build_bases(t, 1)); \
  printf("target=%d: %d bases\n", t, #BM[1]); \
)
```

---

## 11. Safe Prime Search (Layer 9)

Layer 9 is new in v9: a complete safe prime search engine that implements
the same combined q+2q+1 trial division algorithm as `cc_gmp_v35_03.c`.

### Why this matters

At 4096 bits, BPSW primality testing costs ~5-20 ms per candidate.
Without trial division, you test BPSW on every random `q` — most of which
have a small factor in `q` or `2q+1`. Combined trial division against
thousands of small primes eliminates ~99.4% of candidates cheaply, giving
a **~160x speedup** before any expensive BPSW call.

This is exactly what OpenSSL does for `dhparam` generation, and it's what
v35-03 adds to the C codebase.

### Quick safe prime generation

```gp
\\ Find safe primes at various bit sizes
cc_sp_search(512)           \\ 512-bit, ~instant
cc_sp_search(1024)          \\ 1024-bit, ~1 second
cc_sp_search(2048)          \\ 2048-bit, ~seconds
cc_sp_search(4096)          \\ 4096-bit, ~30-120s in GP
```

### Understanding the search space

```gp
\\ Before searching: understand the math
cc_sp_info(4096)
```

Output shows:
- Trial prime count and limit (auto-scaled: bits*25)
- Estimated trial division survival rate
- Estimated PRP and safe prime rates
- Cost comparison: BPSW tests per safe prime with/without trial division

### How it works: combined q+2q+1 trial division

The core optimization, reimplemented here as `cc_sp_trial_check()`:

```gp
\\ For each trial prime p, compute:
\\   rp = q mod p          — if 0, q is composite, skip
\\   rp2 = (2*rp + 1) mod p — if 0, 2q+1 is composite, skip
\\ This is O(1) per prime — no bignum arithmetic needed.

\\ Example: q=11, primes=[3,5,7]
cc_sp_trial_check(11, [3, 5, 7])    \\ returns 1 (both survive)

\\ q=13, 2q+1=27=3*9
cc_sp_trial_check(13, [3, 5, 7])    \\ returns 0 (2q+1 div by 3)
```

### Verified safe primes

```gp
\\ Prove both q and 2q+1 with APRCL/ECPP
cc_sp_verify(11)                    \\ returns [11, 23, 5]

\\ Or search with proven=1 (slower but proven)
cc_sp_search(512, 0, 10000000, 1, 1)
```

### Benchmarking trial vs no-trial

```gp
\\ Compare search with and without trial division
cc_sp_bench(64)                     \\ 64-bit: too fast to measure
cc_sp_bench(256)                    \\ 256-bit: visible speedup
```

At large bit sizes, the trial division speedup is dramatic because BPSW
cost grows as O(n^3) while trial division stays O(n * num_primes).

### Tuning the trial limit

```gp
\\ Default: auto-scaled to max(1000, bits*25)
cc_sp_search(4096)                  \\ uses trial_limit = 102400

\\ Override with a larger sieve
cc_sp_search(4096, 200000)          \\ primes up to 200000

\\ Smaller sieve (faster per candidate, but more reach BPSW)
cc_sp_search(4096, 10000)
```

### Safe prime search pipeline

Each candidate goes through:
1. Random odd `q` in `[2^(bits-2), 2^(bits-1)-1]`
2. Force `q ≡ 5 (mod 6)` — guarantees `q` odd, not div by 3, `2q+1` not div by 3
3. Combined trial division against primes 3 to trial_limit
4. BPSW test on `q` (rejects ~99.96% at 4096 bits)
5. BPSW test on `2q+1`
6. Optional proven primality via `isprime()`

### Correspondence with C code (v35-03)

| GP (cc_lib_v10) | C (cc_gmp_v35_03) | Purpose |
|---|---|---|
| `cc_sp_trial_primes(3, N)` | `init_safeprime_trial_primes(N)` | Sieve of Eratosthenes |
| `cc_sp_trial_check(q, primes)` | Inner loop in worker | Combined q+2q+1 mod check |
| `cc_sp_search(bits)` | `--safeprime --bits N` | Full search pipeline |
| `cc_sp_info(bits)` | Startup banner output | Search space analysis |
| Auto-scale: `bits*25` | `target_bits * 25` | Trial limit formula |
| `--trial-limit N` | `--trial-limit N` | CLI override |

This GP implementation serves as both POC and reference for verifying
the C code produces correct results.

---

## 12. Function Reference

### Layer 1 — Primitive Helpers

| Function | Description |
|----------|-------------|
| `cc_next(p, {kind=1})` | Next chain member: 2p+1 or 2p-1 |
| `cc_prev(p, {kind=1})` | Previous member: (p-1)/2 or (p+1)/2 |
| `cc_immune_residue(q, {kind})` | Immunizing residue: q-1 (first) or 1 (second) |
| `cc_is_immune(p, q, {kind})` | Check if p is immune to prime q |
| `cc_kill_pos(r, q, {kind}, {mx=32})` | Kill position for residue r mod q |
| `cc_root(p, {kind})` | Walk back to find chain root |
| `cc_is_root(p, {kind})` | Check if p is a chain root |

### Layer 2 — Chain Walking

| Function | Description |
|----------|-------------|
| `cc_walk(p, {kind}, {lookahead=6})` | Full chain walk with display |
| `cc_walk_quiet(p, {kind})` | Silent walk, returns `[chain_vec, len, root, breaker]` |
| `cc_walk_from_quiet(p, {kind})` | Exact-start walk, returns `[chain_vec, len, breaker]` |

### Layer 3 — Analysis

| Function | Description |
|----------|-------------|
| `cc_shadow(p, {kind}, {primes})` | Shadow spiral: immune/kill per prime |
| `cc_autopsy(p, {kind}, {depth})` | Deep analysis of chain breaker |
| `cc_deep_factor(n, {depth})` | Recursive p-1/p+1 factor tree |
| `cc_full(p, {kind})` | Everything at once |
| `cc_compare(p1, p2, {kind})` | Side-by-side comparison |

### Layer 7 — CRT Search Engine

| Function | Description |
|----------|-------------|
| `cc_search(bits, target, ...)` | Main CRT search |
| `cc_qsearch(bits, target, ...)` | Quick wrapper (50K candidates) |
| `cc_search_info(bits, target, ...)` | Search configuration |
| `cc_chain_test(n, {kind})` | Quick PRP chain length |
| `cc_chain_verify(n, {kind})` | Proven primality chain length |

### Layer 8 — Constructor Search Engine

| Function | Description |
|----------|-------------|
| `cc_construct(bits, target, ...)` | Main constructor search |
| `cc_qconstruct(bits, target, ...)` | Quick wrapper (50K candidates) |
| `cc_con_info(bits, target, {S_limit})` | Constructor parameters |
| `cc_primorial(n)` | Product of primes <= n |

### Layer 9 — Safe Prime Search Engine

| Function | Description |
|----------|-------------|
| `cc_sp_search(bits, {trial_limit}, {max_cand}, {verbose}, {proven})` | Main safe prime search |
| `cc_sp_qsearch(bits, {trial_limit})` | Quick wrapper (1M candidates) |
| `cc_sp_info(bits, {trial_limit})` | Search space analysis |
| `cc_sp_trial_primes(lo, hi)` | Generate trial prime list |
| `cc_sp_trial_check(q, trial_primes)` | Combined q+2q+1 trial division |
| `cc_sp_verify(q)` | Prove q and 2q+1 prime |
| `cc_sp_bench(bits, {max_cand})` | Benchmark trial vs no-trial |

### Layer 9b — OpenSSL Comparison

| Function | Description |
|----------|-------------|
| `cc_sp_openssl_trial_count(bits)` | OpenSSL's trial prime count (from calc_trial_divisions) |
| `cc_sp_delta_sieve(bits, {n_trial}, {max_cand}, {max_delta}, {verbose})` | OpenSSL-style delta-sieve search |
| `cc_sp_compare(bits, {runs}, {max_cand})` | Head-to-head comparison: OpenSSL vs more primes |
| `cc_sp_trial_analysis(bits)` | Theoretical analysis of trial prime count impact |

---

## 12b. OpenSSL Comparison (Layer 9b)

Layer 9b reimplements OpenSSL's safe prime algorithm (Zimmermann quick sieve)
in GP/PARI to enable head-to-head comparison with our approach.

### OpenSSL's algorithm

1. Generate random `p₀ ≡ 3 (mod 4)` at target bit size
2. Compute `mods[i] = p₀ mod prime[i]` for all trial primes (ONCE)
3. Try offsets `p = p₀ + delta` (step +4, maintaining `p ≡ 3 mod 4`):
   - Check `(mods[i] + delta) % prime[i] <= 1` for all trial primes
   - The `<= 1` trick: residue 0 means `p` divisible, residue 1 means `q=(p-1)/2` divisible
4. Survivor: BPSW test on `p`, then BPSW on `q = (p-1)/2`

Key difference from our v35-03: OpenSSL computes residues ONCE per random start
and adjusts cheaply per delta, while we do fresh `mpz_fdiv_ui` per candidate.
But OpenSSL uses only 1024 trial primes at 4096-bit while we use ~9805.

### Analyzing the trial prime impact

```gp
\\ Show how survival rate changes with trial prime count
cc_sp_trial_analysis(4096)
```

Output shows that going from OpenSSL's 1024 primes to our 9805 primes
eliminates ~39% more candidates before expensive BPSW testing.

### Running the comparison

```gp
\\ Head-to-head at 512 bits (fast, good for demos)
cc_sp_compare(512, 5)

\\ Full comparison at larger bit sizes
cc_sp_compare(1024, 3)
```

### Using the delta-sieve directly

```gp
\\ OpenSSL-style with OpenSSL's prime count
cc_sp_delta_sieve(512)

\\ Delta-sieve with more trial primes (our recommendation)
cc_sp_delta_sieve(512, 1000)

\\ Silent mode
result = cc_sp_delta_sieve(512, 0, 10000000, 0, 0);
```

### Key finding

OpenSSL issue #17956 externally validated that using 20M trial primes
was 2x faster than OpenSSL's default. Our analysis shows even a modest
increase from 1024 to ~10000 primes yields ~39% fewer BPSW tests at 4096 bits.
The extra sieve cost per prime is ~2 ns (negligible vs ~5 ms BPSW cost).

### C implementation: v35-04 (validated)

The delta-sieve algorithm from this GP layer was ported to C as
`cc_gmp_v35_04_safeprime.c` — a dedicated multi-threaded safe prime generator.
Initial benchmarks confirm it **beats OpenSSL's `dhparam`** at 4096-bit
single-threaded, with further gains from multi-threading and thread pinning.
The GP `cc_sp_delta_sieve()` function implements the same algorithm and can be
used to prototype, verify, and explain the approach before running the C code.

---

## 13. Alternative Algorithms — Layer 10

Layer 10 contains 7 functions mined from the 22 CUDA programs in `cuda_legacy/`.
All use the `cc_x_` prefix (experimental/alternative).  These implement techniques
explored during GPU development that are mathematically interesting and complement
the production pipeline.

### Bitmask sieve — O(1) candidate check per prime

`cc_forbidden_res()` returns a vector.  `cc_x_forbidden_mask()` returns an
integer bitmask where bit `r` is set if residue `r` is forbidden.  This enables
O(1) lookup via `bittest(mask, p % q)`.

```gp
\\ Build mask for prime 5, chain depth 3
mask = cc_x_forbidden_mask(5, 3);      \\ = 13 (bits 0, 2, 3)

\\ Test a candidate — 0 means survives, 1 means killed
bittest(mask, 1122659 % 5)             \\ 0 — survives (CC7a root)
bittest(mask, 25 % 5)                  \\ 1 — killed (25 mod 5 = 0)

\\ Batch-screen many candidates against multiple primes
sp = [5, 7, 11, 13, 17, 19, 23];
masks = vector(#sp, i, cc_x_forbidden_mask(sp[i], 18));

screen(p) = {
  for(i = 1, #sp,
    if(bittest(masks[i], p % sp[i]), return(0))
  ); 1;
};

screen(1122659)    \\ 1 — survives all 7 primes for CC18
```

### Periodic forbidden tables — exploit multiplicative order

The forbidden residues repeat with period `ord_q(2)` (multiplicative order
of 2 mod q).  For primes with short periods, this saves precomputation.

```gp
cc_x_periodic_table(7)     \\ [[0, 3, 1], 3]  — period 3
cc_x_periodic_table(31)    \\ [table, 5]       — period only 5!
cc_x_periodic_table(89)    \\ [table, 11]      — period 11

\\ Check forbidden residue at chain position 100 for prime 31
pt = cc_x_periodic_table(31);
table = pt[1]; period = pt[2];
table[(100 % period) + 1]      \\ forbidden residue at position 100
```

### Sieve-only line filter — fast pre-screening without isprime()

`cc_x_line_filter()` follows a chain using only trial division by small primes.
No BPSW, no `isprime()`.  Returns how many steps survive the sieve.

```gp
\\ CC7a root: all 7 members survive small-prime sieve
cc_x_line_filter(1122659, 7)           \\ = 7

\\ Random composite: fails quickly
cc_x_line_filter(1122660, 7)           \\ = 0 (even, fails immediately)

\\ Use as pre-filter before expensive isprime()
candidate_worth_testing(p, depth) = cc_x_line_filter(p, depth) >= depth;
```

### BiTwin chains — simultaneous first + second kind

BiTwin chains start from an even center `n` with a first-kind chain from `n-1`
and a second-kind chain from `n+1` running simultaneously.

```gp
\\ Walk BiTwin from center 6
cc_x_bitwin_walk(6)        \\ [4, 2, 2] — CC4a from 5, CC2b from 7

\\ Find interesting BiTwin centers
forstep(n = 2, 200, 2,
  my(bt = cc_x_bitwin_walk(n));
  if(bt[3] >= 2,
    printf("  center=%d: CC%da from %d, CC%db from %d, BiTwin=%d\n",
           n, bt[1], n-1, bt[2], n+1, bt[3])
  )
)

\\ BiTwin forbidden mask: combined constraints for both kinds
cc_x_bitwin_mask(5, 2)     \\ = 15 (only residue 4 out of 5 survives)
cc_x_bitwin_mask(7, 3)     \\ how constrained is depth 3?
```

### Primorial algebraic scan — how historical records were found

Historical CC world records used the form `p = k * P# * 2^m - 1` where `P#`
is a primorial.  This form guarantees divisibility by the first N primes is
already handled.

```gp
\\ Scan 20-bit range with 7# = 210, find CC3+
cc_x_primorial_scan(20, 4, 3)

\\ Scan 30-bit range with 11# = 2310, find CC4+
cc_x_primorial_scan(30, 5, 4)

\\ Larger scan: 40-bit, 13# = 30030, CC5+
cc_x_primorial_scan(40, 6, 5, 50000)
```

### Root depth — k-value tracking

How deep in the chain is a given prime?  `cc_x_root_depth()` counts backward
steps to the true root.

```gp
\\ Chain: 2 → 5 → 11 → 23 → 47
cc_x_root_depth(2)     \\ 0 — IS the root
cc_x_root_depth(5)     \\ 1
cc_x_root_depth(11)    \\ 2
cc_x_root_depth(23)    \\ 3
cc_x_root_depth(47)    \\ 4

\\ Verify consistency with cc_is_root
cc_is_root(2)          \\ 1 (depth 0 = root)
cc_is_root(47)         \\ 0 (depth 4 = not root)
```

### Layer 10 function reference

| Function | Description |
|----------|-------------|
| `cc_x_forbidden_mask(q, depth, {kind=1})` | Bitmask of forbidden residues mod q (O(1) lookup via bittest) |
| `cc_x_periodic_table(q, {kind=1})` | Periodic forbidden table `[table_vec, period]` via `znorder(Mod(2,q))` |
| `cc_x_line_filter(p, depth, {kind=1}, {max_prime=107})` | Sieve-only chain pre-screen (no isprime), returns surviving steps |
| `cc_x_bitwin_mask(q, depth)` | Combined forbidden mask for BiTwin (both kinds + evenness) |
| `cc_x_bitwin_walk(n, {max_depth=50})` | BiTwin chain walk, returns `[CC_a_len, CC_b_len, min]` |
| `cc_x_primorial_scan(bits, prim_n, {target=5}, {k_max=10000})` | Algebraic form `k*P#*2^m-1` search, returns `[[p, len], ...]` |
| `cc_x_root_depth(p, {kind=1})` | Backward steps to chain root (0 = is root, -1 = not prime) |

---

## 14. Combining Layer 10 with Other Layers

The `cc_x_` functions are most powerful when combined with existing layers.

### Pre-filter + full verification pipeline

```gp
\\ Use line_filter (Layer 10) as fast pre-screen,
\\ then cc_chain_verify (Layer 7) for proven chains

p = 1122659;
depth = 18;

\\ Step 1: cheap sieve-only check (microseconds)
if(cc_x_line_filter(p, depth) >= depth,
  \\ Step 2: only if sieve passes, do expensive proven check
  printf("Sieve survived %d steps, verifying...\n", depth);
  cc_chain_verify(p),

  printf("Sieve killed at step %d, skipping\n", cc_x_line_filter(p, depth))
)
```

### Bitmask sieve + CRT search (Layer 7)

```gp
\\ Build bitmask table once, use in custom search loop
sp = [5, 7, 11, 13, 17, 19, 23, 29, 31];
masks = vector(#sp, i, cc_x_forbidden_mask(sp[i], 10));

\\ Screen candidates: O(1) per prime via bittest
bitmask_screen(p) = {
  for(i = 1, #sp, if(bittest(masks[i], p % sp[i]), return(0))); 1;
};

\\ Combine with standard search for CC5+ at 30 bits
hits = cc_search(30, 5, 1, 0, 0, 50000, 0);
printf("Standard search found %d hits\n", #hits);

\\ Cross-validate: every hit should pass bitmask screen for depth 5
ok = 1;
for(i = 1, #hits, if(!bitmask_screen(hits[i][1]), ok = 0));
printf("All hits pass bitmask screen: %s\n", if(ok, "YES", "NO"));
```

### BiTwin discovery + full analysis (Layers 1-6)

```gp
\\ Find BiTwin chains, then analyze the best one with full toolkit
best_center = 0; best_bt = 0;
forstep(n = 2, 10000, 2,
  my(bt = cc_x_bitwin_walk(n));
  if(bt[3] > best_bt, best_bt = bt[3]; best_center = n)
);
printf("Best BiTwin center in [2,10000]: %d (BT%d)\n", best_center, best_bt);

\\ Deep analysis of both chains
printf("\nFirst-kind chain from %d:\n", best_center - 1);
cc_walk(best_center - 1, 1);
printf("\nSecond-kind chain from %d:\n", best_center + 1);
cc_walk(best_center + 1, 2);
```

### Periodic table for long chain analysis

```gp
\\ For CC18 analysis: which primes have short periods?
\\ Short period = forbidden residues repeat quickly = fewer unique constraints
printf("Prime  Period  Unique_forbidden\n");
forprime(q = 3, 50,
  my(pt = cc_x_periodic_table(q));
  printf("  %3d  %5d  %d residues\n", q, pt[2], min(pt[2], 18))
)
```

### Root depth census — fragment vs true root statistics

```gp
\\ Find CC5+ in a range and classify: true roots (depth 0) vs fragments
hits = cc_search(35, 5, 1, 0, 0, 100000, 0);
my(roots = 0, frags = 0);
for(i = 1, #hits,
  my(d = cc_x_root_depth(hits[i][1]));
  if(d == 0, roots++, frags++)
);
printf("CC5+ hits: %d roots, %d fragments (%.0f%% are true roots)\n",
       roots, frags, 100.0 * roots / (roots + frags));
```

---

## Cheat Sheet

```gp
\r cc_lib_v10.gp                          \\ load library (37 tests)

\\ ---- Analyze a discovery ----
cc_walk(p)                               \\ chain members + residues
cc_shadow(p)                             \\ immunization map
cc_autopsy(p)                            \\ breaker forensics
cc_full(p)                               \\ everything

\\ ---- Plan a search ----
cc_search_info(89, 18)                   \\ CRT search space
cc_con_info(89, 18)                      \\ constructor search space
cc_sp_info(4096)                         \\ safe prime search space

\\ ---- Run a search ----
cc_qsearch(40, 5)                        \\ CRT quick search
cc_construct(50, 5, 19)                  \\ constructor search
cc_sp_search(512)                        \\ safe prime search

\\ ---- Safe prime workflow ----
cc_sp_search(4096)                       \\ find 4096-bit safe prime
cc_sp_search(4096, 200000)               \\ with larger trial sieve
result = cc_sp_search(2048, , , 0);      \\ silent search
cc_sp_verify(result[1])                  \\ prove it

\\ ---- OpenSSL comparison ----
cc_sp_trial_analysis(4096)               \\ trial prime impact table
cc_sp_compare(512, 5)                    \\ head-to-head comparison
cc_sp_delta_sieve(512)                   \\ OpenSSL-style search

\\ ---- Alternative algorithms (Layer 10, cc_x_) ----
cc_x_forbidden_mask(5, 18)               \\ bitmask sieve for CC18
cc_x_line_filter(p, 10)                  \\ sieve-only pre-screen
cc_x_bitwin_walk(6)                      \\ BiTwin from even center
cc_x_primorial_scan(20, 4, 3)            \\ algebraic k*P#*2^m-1
cc_x_root_depth(p)                       \\ backward steps to root
cc_x_periodic_table(31)                  \\ periodic forbidden table

\\ ---- Verify a hit ----
cc_chain_verify(p)                       \\ proven chain primality
cc_sp_verify(q)                          \\ proven safe prime
```
