# cc20 first-kind immune (v16) — design

CUDA Cunningham Chain search engine. Successor to
`cc18_filter_cuda_CpC_v15.cu`. Target: CC20/CC21-scale first-kind
Cunningham chains.

> **Lineage in one line:** v13 and v15 are CRT-lattice search engines that
> enumerate every candidate at a single sieve shape `D = 5·7·11·13·17·19`.
> v16 (this engine) is a **fingerprint-pool Q-iteration** engine that
> draws a randomly-sampled `Q` from each of ~4,200 empirically-derived
> sieve shapes `D_V` and tests `p = Q · D_V − 1`.

---

## 1. Architectural lineage

| version | search architecture | candidate model | sieve shapes | record era |
|---|---|---|---|---|
| `cc18_filter_cuda_CpC_v13.cu` | CRT lattice + wheel | `n = (tile·WP + k)·M + base` | 1 (`D = 5·7·11·13·17·19`) | CC17–CC18 |
| `cc18_filter_cuda_CpC_v15.cu` | CRT lattice + wheel + full mod-6 coverage | as v13, `~36` CRT bases | 1 | CC17–CC18 |
| `cc20_first_kind_immune_v16.cu` (this engine) | **Fingerprint-pool Q-iteration** | `p = Q · D_V − 1` | **~4,200** (empirical pool) | CC18–CC20 |

v13 and v15 sieve a single arithmetic progression and hope the wheel hits
chain roots. v16 inverts the question: *what fingerprints have produced
chain roots in the historical data, and how do we enumerate primes that
share those fingerprints?*

---

## 2. The Q-iteration construction

Given a target Cunningham chain root `p_0`, we have
`p_0 + 1 = ∏ q_i^{e_i}`  for small primes `q_i`. We call this product
the chain's **fingerprint** `D_V`. By construction, for any prime `q | D_V`
we have `p_0 ≡ −1 (mod q)`, which means small primes never divide `p_0`
itself. The fingerprint **immunizes** `p_0` against small-prime
compositeness.

Empirical observation (from `data/cc10plus_roots_snapshot_2026-03-19.csv`,
1.55 M CC10+ root samples): the set of distinct `D_V` fingerprints in real
chains is small (~4,200 with frequency ≥ 10) and dominated by a long tail
of `2·3·5·11·13·19 · …` shapes. We bake these into a **pool file**:

```
pools/known_fingerprints_exp_d.txt
  # via scripts/build_known_fingerprints_exp.py
  # variant D: freq >= 10  (4248 entries, 96.78% coverage of CC10+)
  2*3*5*11*13*19                  # freq=271525 bits=16.31
  2*3^2*5*11*13*19                # freq=89509  bits=17.90
  2*3*5*7*11*13*19                # freq=74063  bits=19.12
  ...
```

Each line is a `D_V` written as `q1^e1 * q2^e2 * …`.

For each variant `V` (one line), the engine **iterates Q** over the bit
band: `Q_min(V) = ⌈2^{N1−1} / D_V⌉`, `Q_max(V) = ⌊(2^{N2}−1) / D_V⌋`,
and tests `p = Q · D_V − 1` for primality and chain extension. `Q` is
sampled either sequentially or via xor-shift-64* random sampling
(`--q-order random`).

The throughput unit is **G(Q,V)/s** — billions of (Q, variant) pairs
evaluated per second.

### 2.1. Per-V Q-walker mode (`--q-band-mode`)

The Q-range bit-width `bits(Q_max - Q_min + 1)` varies wildly across a
pool: tiny for high-`bits(D_V)` primorial-quotient variants (every
integer in the range can be enumerated in seconds), astronomical for the
empirical exp-d shapes (random sampling is the only practical option).

A single global `--q-order` flag forces a bad trade-off either way:
random on the narrow variants produces birthday-style repeat samples
within hours, while sequential on the wide ones never reaches the bits
where chains live. `--q-band-mode exhaustive` resolves this by
classifying each variant individually at startup:

- `bits(Q_max - Q_min + 1) ≤ --exhaustive-max-q-bits` ⇒ sequential walk
  (deterministic coverage, no duplicates).
- Otherwise ⇒ random sampling.

`--q-band-mode fixed` (default) preserves the legacy behaviour: all
variants follow `--q-order`. The chosen split is reported in the banner
and persisted into the checkpoint so a resume cannot silently flip
between modes.

### 2.2. Emit dedup (`--dedup-p0`)

Two independent mechanisms can surface the same `p_0` more than once:

1. **Random with-replacement sampling.** xor-shift-64* has no
   visited-Q memory; narrow Q-ranges produce birthday collisions.
2. **Structural multiplicity in primorial-quotient pools.** A given
   `p_0 + 1` factors as `Q · D_V` for several `(Q, V)` pairs when the
   pool contains overlapping primorial quotients — each framing
   "discovers" the same chain.

The engine fingerprints `p_0` (xor of GMP limbs + Murmur mixer) and
emits the first sighting verbose (full fp/pool buffer); subsequent
sightings collapse to a compact `HIT-DUP` line so the hit log still
records the duplicate without re-printing the multi-V framing. Disable
with `--no-dedup-p0` if downstream tooling expects every emit verbose.

---

## 3. Why v13/v15's CRT model couldn't get here

v15's CRT lattice picks a *single* sieve shape and uses ~36 bases to
cover its residue classes. To get the same coverage v16 has, v15 would
need ~4,200 separate runs (one per fingerprint). The Q-iteration model
interleaves them naturally on the GPU — every kernel launch touches all
active variants. The per-variant work is amortized via a forbid-mask
table that fits comfortably in L2.

A second consequence: v13/v15 ignore the empirical fingerprint
distribution. v16 weights effort toward shapes that have produced known
chains. This is also why the recent **"pool 7-pin" bug** mattered so
much: a script-level bug forced every `D_V` to include the prime 7, but
empirically only 24.16 % of CC10+ chains have `7 | (p+1)`. The fix —
augmenting only with the universal floor `{2, 3, 5}` — increased
structural coverage from 24 % to 97 % of the historical chain space.

---

## 4. Pipeline

```
GPU                                           CPU
─────────────────────────────────────────     ─────────────────────────────
K1   Filter kernel                           ProvePool (pthread workers,
     - For each (block, V), iterate Q          async queue, qd=8..16)
       over its window.                       - Pop survivor (Q, V)
     - Reject Q whose residues mod the        - Reconstruct p_0 = Q·D_V − 1
       sieve primes (7..503) would make       - Top-MR p_(target-1)
       any of p_0..p_{depth-1} divisible        (skipped if K1.5 ran)
       by a small prime.                      - Root MR p_0
     - Active-prime pruning per variant:      - Birth-certificate root
       skip primes q where q | D_V.             check (v34 OPT-I)
                                              - Chain follow → emit
K1.5 GPU top-MR pre-filter
     - Miller-Rabin witnesses {2,3,(5)} on
       p_{target-1}; reject composites
       before CPU sees them.

Stream pool: N=2..16 GPU streams in flight, each owns a pinned host
survivor buffer; CPU prove drains them asynchronously.
```

### Filtering layers (GPU)

| layer | primes | notes |
|---|---|---|
| variant immunity (D_V) | per pool entry | `p_0 ≡ −1` mod every prime in `D_V` |
| forbid-mask sieve | 7..503 (93) | bitmask of bad Q residues per (V, q) |
| active-prime pruning | per V | skip masks where q ∈ D_V |
| K1.5 GPU top-MR | — | 2- or 3-witness MR on `p_{target-1}` |

`V16_Q_MASK_WORDS = 8` (i.e. 512 bits) caps individual sieve primes at
≤ 503. Going wider needs either a bigger mask (4× memory) or a separate
out-of-shared-mem secondary sieve.

### Prove layers (CPU)

| stage | cost | notes |
|---|---|---|
| top MR (reverse-prove) | ~30 µs | skipped when K1.5 already ran |
| root MR (`mpz_probab_prime_p(p, 5)`) | ~30 µs | drops 90 % of survivors |
| root check `(p−1)/2` | ~30 µs | **birth-certificate fast path** (v34 OPT-I): `p ≡ 1 (mod q)` for any small `q` ⇒ `q | (p−1)/2` ⇒ `p` is the chain root, skip MR. ~70 % coverage in our prime range. |
| chain follow | ~30 µs × `chain_len` | forward MR walk |
| re-verification | ~100 µs × `chain_len` | high-rep MR (`reps=25`) when `chain_len ≥ --min-report-len`. |

`--prove-order reverse` is the default and trims the dominant case
(`p_{target-1}` composite, ~98.6 % of survivors at CC20 / 90+ bits).

**Why reverse here when v13 used forward.** v13 ran forward to harvest
the 1.5 M CC10+ chain corpus we now use as ground truth: its CRT lattice
produced relatively few survivors, so the CPU could absorb a full
chain-walk per candidate and emit every short-chain side-find. v16's
fingerprint-pool model produces 10²–10³× more survivors per second
(4,248 variants in flight vs. 1), so forward would crater the CPU prove
queue. Reverse pre-gates on `p_{target-1}` and discards 98.6 % of
survivors before any chain walk — a CPU-protection move forced by the
new architecture, not a primality-testing preference. `--prove-order
forward` is still available if a future workload (small pool, sparse
band) wants more side-emits.

---

## 5. CLI surface (delta vs v15)

| flag | meaning | v15 analog |
|---|---|---|
| `--q-iter` | enable q-iteration mode (default for v16 wrappers) | — |
| `--seed-pool-file FILE` | path to `known_fingerprints_exp_d.txt` | — |
| `--immune-prime N` | sieve floor; pool `D_V` must be a superset | similar |
| `--target N` | required chain length | same |
| `--depth N` | sieve depth (positions 0..N−1 must avoid small primes) | `--depth` |
| `--bits-min N1` / `--bits-max N2` | `p_0` bit band | `--bits` (single value) |
| `--exp-start E` | shift `Q` window by `2^E` | — |
| `--q-order sequential\|random` | global Q-walker mode (active in `fixed` band-mode) | similar |
| `--q-band-mode fixed\|exhaustive` | per-V auto-pick: sequential when Q-range narrow, random otherwise (default `fixed`) | — |
| `--exhaustive-max-q-bits N` | Q-range bit threshold for the auto-pick (default 40) | — |
| `--dedup-p0` / `--no-dedup-p0` | collapse repeat-`p_0` emits to `HIT-DUP` (default on) | — |
| `--streams N` | GPU stream pool depth (default 4) | `--streams 2` (fixed) |
| `--prove-threads N` | CPU pool size | same |
| `--prove-queue-depth N` | async queue depth (default 8) | — |
| `--test` (default) / `--no-test` | startup self-test gate (CC5..CC17 + 5 q-iter round-trip vectors) | `--test` (opt-in) |
| `--min-report-len N` | emit HITs ≥ this length | `--log` |
| `--immune-factor-set 2,3,5,11,…` | explicit sieve floor | — |

---

## 6. Self-test

`run_v16_self_test()` runs at startup (CPU/GMP only, no GPU init) and
asserts the engine's chain-follow machinery against literature:

1. **Classic vectors** (CC5..CC14, CC17 Wroblewski 2008) — root primality,
   chain-root check, forward chain length.
2. **Q-iteration round-trip vectors** — five known CC17/CC18 roots,
   including both world-record CC18s, with their (Q, D_V) decomposition.
   Asserts `Q · D_V − 1 = p_0`, that `p_0` is prime, the chain extends
   to the declared length, and the `D_V` shape is present in the loaded
   pool.

If any vector fails the engine returns 2 from `main()` *before* any GPU
init or pool dispatch. This catches the class of structural bug
exemplified by the May 2026 pool 7-pin incident.

Run with `--no-test` to skip (saves ~30 s on hot-restart loops).

---

## 7. Pool curation

`scripts/build_known_fingerprints_exp.py` reads the CC10+ snapshot
CSV (1.55 M chain roots) and emits the pool with these rules:

1. For each root, factor `p_0 + 1` over the palette
   `{3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61}` (exponent-aware).
2. Augment with the universal floor `{2, 3, 5}` (every CC10+ chain has
   2, 3, 5 dividing `p_0+1`). `--pin-7` reproduces the historical
   `{2,3,5,7}` augmentation for legacy 7-immune campaigns.
3. Drop shapes with frequency < `--min-freq` (default 10).
4. Emit in `q^e * q^e * …` notation; the engine's parser accepts both
   `q*q` (repeated) and `q^e` forms.

The structural blind-spot bug (2026-05-14): forcing `7 ∈ D_V` for every
variant silently rejected ~76 % of historical CC10+ chains, because q=7
is immune in only 24.16 % of them. Verified against the snapshot via
`tools/snapshot_replay.py` (the *external prefix test*).

---

## 8. Performance characteristics

On two RTX 5090 test systems (10-core and 30-core CPUs, May 2026):

| config | 10c CPU | 30c CPU | notes |
|---|---|---|---|
| v15 baseline | — | ~33.8 G/s | single CRT shape, mid-2024 baseline |
| v16 depth=20 target=16 | ~140 G/s | ~140 G/s | but silently filtered CC10..15 |
| v16 depth=16 target=16 + birth-cert + qd=16 | 65.9 G/s | 88.7 G/s | structurally correct, prove-bound on the 10-core host |
| v16 depth=20 target=20 (historical CC20 config) | **147 G/s** | **141 G/s** | reverse-prove gates on `p_19` |

Per-(Q,V) survivor density `s/c`:

- depth=16: 1.91e-5 (1 in 52K Qs survives sieve)
- depth=20: 1.54e-6 (1 in 650K)

depth=20 cuts survivors 12× — frees CPU prove to keep pace with GPU.

---

## 9. Known limitations / open work

Two pieces particularly worth knowing:

1. **`V16_Q_MASK_WORDS = 8` caps sieve primes at ≤ 503.** Going wider
   needs a 4× larger forbid-mask or a separate L3-style filter for cold
   primes 503..few-K. Not yet implemented.
2. **CPU trial-division pre-MR (primes 503..10 K)** is *net-negative*
   without bit-vector batching like `cc_gmp_v34_bit-vector_10.c`'s
   John Armitage L2 filter — Fermat is too cheap to pre-filter for one
   candidate at a time. The v34 trick batches 64 candidates per pass;
   porting that batching to v16's prove worker is the natural next
   speedup but a substantial refactor.

---

## 10. File map (key)

```
cc20_first_kind_immune_v16.cu              # this engine
pools/known_fingerprints_exp_d.txt       # pool (4248 entries)
scripts/build_known_fingerprints_exp.py    # pool generator
scripts/build_known_fingerprints_all.py    # alt (squarefree-shape) generator
tools/snapshot_replay.py                   # external prefix test / structural validator
scripts/run_cc20_hunter.sh                 # production wrapper
examples/                                 # one-command tmux launchers by pool
```

---

## 11. References

- `src/cuda/cc18_filter_cuda_CpC_v13.cu` —
  CRT baseline, mod-6 buckets, single sieve shape.
- `src/cuda/cc18_filter_cuda_CpC_v15.cu` —
  CRT + full mod-7 coverage, GPU K1+K2.
- `src/cpu/cc_gmp_v34_bit-vector_10.c` —
  John Armitage's bit-vector L2 filter, OPT-I birth certificate
  (ported here), super-ext-L2 and L1-resident kill tables (relevant
  for future v16 prove-worker batching).
