# Bugs and lessons — postmortems from development

Selected debugging stories worth reading if you're building or modifying
a CC search engine in this style. Each bug here was found empirically
and cost us hours; documenting them is the most useful thing we can pay
forward.

## 1. The pool 7-pin bug — 76 % of CC10+ structures silently rejected

**Symptom.** v16 was at ~140 G(Q,V)/s but found 0 hits in 4+ hours at
the [90, 127] bit band. v15 on the same hardware found chains within
minutes. We assumed it was sparsity at high bits, but the discrepancy
was structural, not statistical.

**Root cause.** `scripts/build_known_fingerprints_exp.py` had a
`PINNED_SIEVE_PRIMES = (2, 3, 5, 7)` constant. The script
**force-augmented every fingerprint** with `{2, 3, 5, 7}` to satisfy a
mistaken belief that the v16 sieve required `7 ∈ D_V` for every
variant. So *every* pool entry contained 7 as a factor of `D_V`.

The math problem: in the 1.55 M CC10+ snapshot, **only 24.16 % of
chains have `7 | (p+1)`**. The other 75.84 % were silently impossible
to find with the broken pool. We were searching one-fourth of the
known chain space and assuming we were searching all of it.

**Diagnosis.** A direct snapshot analysis showed q=7 immunization at
24.16 %, vs q=3, q=5, q=11, and q=13 at 100 %. The "universal floor"
for CC10+ is `{2, 3, 5, 11, 13}` (with q=19 nearly universal at
98.90 %). Pinning 7 was wrong.

**Fix.** Changed `PINNED_SIEVE_PRIMES_DEFAULT = (2, 3, 5)`. Added
`--pin-7` flag for explicit legacy behavior. Regenerated pools:
`known_fingerprints_exp_d.txt` went from 3,433 to **4,248 entries**;
top entry changed from `2*3*5*7*11*13*19` (freq 74,063) to
`2*3*5*11*13*19` (freq 271,525, the new #1).

**Validation.** Built `tools/snapshot_replay.py` to confirm. Coverage
went from 24.16 % to 96.78 % of historical CC10+ chains. Both known
world-record CC18s became structurally reachable.

**Lesson.** When building a search pool from historical data, **verify
coverage against the data**. The corpus is the source of truth, not
intuition. An "external prefix test" tool that round-trips known
solutions through the search structure is the only thing that catches
this class of bug reliably.

**Indirect lesson.** v15 (the predecessor) avoided this because it
didn't use a pool — its CRT lattice was an arithmetic progression that
naturally covered all residue classes. Pool-based search is more
powerful but introduces a new failure mode: pool construction.

## 2. The `depth > target` trap — silent CC10..(target−1) filter

**Symptom.** Switched from `target=10 depth=10` (validation config that
produced 75 hits/hr) to `target=16 depth=20 [90,127]` (production CC16
config). Throughput was great at 140 G(Q,V)/s, but 0 hits emerged in
4 hours.

**Root cause.** The GPU forbid-mask sieve encodes "no small prime
divides any of `p_0..p_{depth−1}`". With `depth=20` and `target=16`:

- The engine *reports* hits of length ≥ 10 (per `--min-report-len`)
- But the sieve *rejects* candidates where `p_10..p_19` would be
  divisible by a small prime
- A real CC10 chain where `p_10..p_15` happen to fail at small primes
  gets dropped at sieve time, never reaching prove

The depth filter doesn't know that the user wanted to *count* chains of
length 10. It pre-filters as if you wanted chains of length 20.

**Fix.** Set `depth = target` (or `depth < target` if you genuinely
want only the longest chains and accept losing all shorter ones).

**Lesson.** "Depth" and "target" sound interchangeable but they're not.
Depth is the sieve's structural constraint; target is the prover's
emission threshold. They should be **equal** unless you have a specific
reason and have done the math.

## 3. Reverse vs. forward prove — CPU protection in the pool architecture

**Symptom.** Tried `--prove-order forward` on the CC20 hunter, hoping
to see more CC10..CC19 side-hits. CPU prove queue immediately backed up
and throughput collapsed to ~25 G(Q,V)/s.

**Root cause.** Forward prove walks `p_0, p_1, p_2, ...` testing
primality at each step. The fingerprint-pool architecture produces
~10²–10³× more survivors per second than v15's CRT lattice (because
4,248 variants are in flight). Forward-walking each survivor through
~20 Miller-Rabin tests **before any rejection** saturates the CPU.

Reverse prove (`--prove-order reverse`, default) tests `p_{target−1}`
first. If composite, the whole candidate is dropped after one MR.
~98.6 % of survivors are rejected by this single check at CC20-level
target depth.

**Lesson.** Reverse prove is not a primality-testing preference — it's
a CPU-protection move forced by the fingerprint-pool's high survivor
volume. If you make changes that increase survivor rate (richer pool,
wider band, deeper bit range), reverse stays mandatory.

## 4. GMP's hidden trial division — wasted ~15 µs/survivor

**Symptom.** Even after the CPU extended-sieve patch reduced prove
load by ~5×, prove timings showed ~30 µs of Fermat work per survivor,
which seemed high for 90-bit numbers.

**Root cause.** `mpz_probab_prime_p(n, reps)` does an internal trial
division by primes up to ~5000 before BPSW. For survivors that have
already passed our GPU sieve (primes 7..503) and CPU ext-sieve
(509..1009), the trial-div pass redoes work we already did.

**Status.** Not yet fixed — when CPU is heavily idle (current
post-ext-sieve state), saving 15 µs/survivor doesn't move throughput.
Custom-MR replacement is a future cleanup.

**Lesson.** "Black-box" primality libraries duplicate work you've
already done structurally. In tight CPU regimes, replace
`mpz_probab_prime_p` with direct `mpz_powm_ui`-based MR. In our regime
this is currently cosmetic, but it'd be the highest-leverage CPU-side
change if the GPU ceiling lifts.

## 5. `V16_Q_MASK_WORDS = 8` — the hard architectural cap

**Symptom.** Tried `--sieve-max 1009` to catch more composites at the
GPU stage. Engine errored at startup: `sieve prime 521 exceeds
V16_Q_MASK_WORDS*64 = 512`.

**Root cause.** Each (variant, prime) pair stores forbidden-Q residues
as a bitmask of size `q` bits. The mask is allocated as 8 × u64 = 512
bits per pair. Primes above 512 don't fit. Raising the mask doubles
forbid-mask memory: 24 MiB → 78 MiB (still within RTX 5090's 96 MiB L2,
but margin shrinks).

**Workaround we shipped.** CPU-side "extended sieve" — primes 509..1009
checked in the prove worker before MR. Same math (forbid-mask bitmask),
just stored in pinned host RAM instead of L2.

**Lesson.** Constant-time array bounds in CUDA kernels are real
performance contracts. Bumping them isn't free, and a CPU-side mirror
of the same filter can be cheaper than enlarging the GPU footprint.

## 6. Per-variant vs aggregate throughput — the metric pitfall

**Symptom.** v16 reports `G (Q,V)/s` (e.g., 150 G/s). The natural
comparison is v15's "G/s" (e.g., 33.8 G/s) — looks like a 4.5× speedup.

**The correction.** v16's metric counts **(Q × variant) pairs**, not
raw Q values. At 4,248 variants in the active pool, per-variant
Q-rate is 150 G / 4248 = **36 MQ/s per variant**, *~1,000× slower
than v15 per variant*.

**The actual architectural bet.** v16 sacrifices per-variant speed for
breadth (4,248 productive fingerprints vs v15's single CRT shape).
Pool-size sweep showed aggregate throughput is flat across pool sizes
10..4248; smaller pools give proportionally higher per-variant rate
(pool=10 → 15.7 GQ/s per variant, ~50 % of v15's single-D rate).

**Lesson.** Throughput metrics in pool-architecture search must be
reported per-variant to compare against single-shape search. Aggregate
G(Q,V)/s flatters the engine when it's actually trading speed for
breadth.
