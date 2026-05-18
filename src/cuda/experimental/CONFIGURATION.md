# CONFIGURATION — flag reference & recommended configs

## Recommended configurations by intent

### "I want to run a CC21-scale long search"

```bash
--q-iter
--seed-pool-file pools/known_fingerprints_exp_d.txt   # 4248 entries
--immune-prime 5
--target 21 --depth 21            # match: depth=target, reverse-prove gate on p_20
--sieve-max 503
--bits-min 90 --bits-max 107       # u128-safe at target=21 (chain top p_20 ≤ 2^127)
--q-order random
--cpu-prove --prove-threads <cores - 4>
--prove-queue-depth 16
--streams 4
--min-report-len 10                # emit CC10+ side-hits too
--time 0                           # run forever
```

Throughput on a single RTX 5090: ~140–160 G(Q,V)/s.

### "I want exhaustive coverage of the high-V primorial slice"

For pools that mix narrow- and wide-Q-range variants (e.g. a primorial-
quotient pool spanning bits(V) ∈ [33, 85]), `--q-band-mode exhaustive`
auto-picks per V: variants with `bit-width(Q_max - Q_min + 1) ≤ 40` walk
every Q deterministically (no birthday collisions, completes in
seconds-to-hours); the rest stay on random sampling.

```bash
--q-iter
--seed-pool-file pools/known_fingerprints_primorial.txt   # 3000 V's
--immune-prime 5
--target 21 --depth 21
--sieve-max 503
--bits-min 90 --bits-max 107
--q-band-mode exhaustive --exhaustive-max-q-bits 40
--dedup-p0                           # default ON — collapses repeat-p_0 emits
--cpu-prove --prove-threads <cores - 4>
--prove-queue-depth 16
--streams 4
--min-report-len 10
--time 0
```

The banner reports the split, e.g. `q-band-mode: exhaustive
threshold=2^40 Q-range (sequential V's=1250 random V's=1750)`. Same-`p_0`
emits from multi-V framings collapse to `HIT-DUP` lines.

### "I want to focus on the most-productive fingerprints" (focused hunter)

```bash
--seed-pool-file pools/known_fingerprints_top10.txt   # 10 entries
# ...rest same as the long-run search preset
```

Same aggregate throughput; **per-variant Q-rate is 400× higher** —
each top fingerprint receives nearly v15-level GPU bandwidth. Trades
coverage for focus. World-record CC18 #1's fingerprint
(`2·3·5·11·13·19`) is the top entry, so this configuration is well-
suited if you believe the next CC20 will share that small-prime
signature.

### "I want fast hit-flow for validation/debugging"

```bash
--target 10 --depth 10
--sieve-max 211
--bits-min 86 --bits-max 90
--cpu-prove --prove-threads 4
--time 300
```

This is a diagnostic: low bits, low target, narrow band — produces
CC10+ hits within seconds, useful for sanity-checking the toolchain
before launching a multi-day campaign.

### "I want exhaustive coverage of all CC10+ shapes" (broad explorer)

The default `--seed-pool-file pools/known_fingerprints_exp_d.txt`
already does this: 4,248 fingerprints covering 96.78 % of historical
CC10+ structures. Useful if you suspect CC20 may live in a rare
fingerprint we haven't catalogued yet.

## CLI flag reference

### Core search

| flag | default | meaning |
|---|---|---|
| `--q-iter` | off | enable q-iteration mode (REQUIRED for v16 production) |
| `--target N` | 10 | target chain length |
| `--depth N` | = target | sieve depth. Set `depth == target` for sane reverse-prove behavior |
| `--sieve-max N` | 211 | sieve prime ceiling (cap: 503 due to V16_Q_MASK_WORDS=8) |
| `--bits-min N1` | 90 | p_0 bit-band lower edge |
| `--bits-max N2` | 127 | p_0 bit-band upper edge. Use 107 at target=21 to keep the chain top inside u128 (`p_20 ≤ 2^127`) with no fixture-validated >128-bit code paths |
| `--exp-start E` | 0 | shift Q-window by 2^E (for chains rooted at higher positions) |

### Q sampling

| flag | default | meaning |
|---|---|---|
| `--q-order random\|sequential` | sequential | global Q sampling strategy in fixed mode. `random` spreads bits across band uniformly; `sequential` walks every Q in `[Q_min(V), Q_max(V)]` in order |
| `--q-band-mode fixed\|exhaustive` | fixed | `fixed`: all V's use `--q-order`. `exhaustive`: per-V auto-pick — V's whose Q-range bit-width ≤ `--exhaustive-max-q-bits` get sequential (deterministic coverage), the rest get random (no realistic chance of exhausting). Recommended for primorial-style pools that mix wide and narrow Q-ranges |
| `--exhaustive-max-q-bits N` | 40 | threshold (in bits) for `--q-band-mode exhaustive` per-V auto-pick |
| `--dedup-p0\|--no-dedup-p0` | on | collapse repeat-`p_0` emits to a compact `HIT-DUP` line. Same `p_0` shows up via random with-replacement sampling **and** via different `(Q,V)` framings of the same chain in primorial-quotient pools |
| `--seed-pool-file FILE` | — | path to fingerprint pool (see `pools/`) |
| `--seed-name NAME` | unnamed | run identifier (appears in HIT lines) |

### Sieve floor

| flag | default | meaning |
|---|---|---|
| `--immune-prime P` | 41 | sieve floor: D contains all primes ≤ P. Use 5 for the v16-default pool |
| `--immune-factor-set 2,3,5,…` | — | explicit prime set (alternative to --immune-prime) |

### CPU prove

| flag | default | meaning |
|---|---|---|
| `--cpu-prove` | on | enable CPU prove pool (required when no GPU prove kernel) |
| `--prove-threads N` | 8 | pthread worker count. Sweet spot: `cores − 4`. Going over the core count oversubscribes |
| `--prove-queue-depth N` | 8 | survivor queue depth per stream. Raise to 16 for depth=16 configs |
| `--prove-order forward\|reverse` | reverse | reverse: test p_{target−1} first, walk down. Strongly recommended |

### GPU streams

| flag | default | meaning |
|---|---|---|
| `--streams N` | 4 | N-deep stream pool. Range [2, 16]. 4 is the empirical sweet spot — 8 doesn't help when GPU is compute-bound |
| `--gpu N` | 0 | CUDA device id |

### Self-test (default-ON)

| flag | default | meaning |
|---|---|---|
| `--test` | on | run startup self-test (~30 s) |
| `--no-test` | — | skip self-test (for hot restarts after a verified build) |

### Output

| flag | default | meaning |
|---|---|---|
| `--report SEC` | 60 | periodic report interval |
| `--time SEC` | 0 | wallclock budget (0 = run forever) |
| `--min-report-len N` | target | minimum chain length to emit to log |
| `--log FILE` | — | hit log path |
| `--checkpoint FILE` | — | resume-state path |
| `--resume FILE` | — | resume from a checkpoint |

## Pool files in `pools/`

| file | entries | coverage of CC10+ | use when |
|---|---|---|---|
| `known_fingerprints_exp_d.txt` | 4,248 | 96.78 % | default broad campaigns; CC10+ exhaustive coverage |
| `known_fingerprints_top10.txt` | 10 | high-frequency shapes only | focused per-variant high-rate hunting |
| `known_fingerprints_all.txt` | 1,372 | squarefree-shape pool | when you want one variant per distinct prime-support (no exponent variants) |

Pool format (every line):

```
<prime>[^exponent]*<prime>[^exponent]*...   # freq=N bits=B.B
```

Examples:
```
2*3*5*11*13*19          # freq=271525  bits=16.31
2*3^2*5*11*13*19        # freq=89509   bits=17.90
2*3*5*7^2*11*13*19      # freq=12345   bits=20.71
```

You can hand-edit these or generate variants via the build scripts
(`scripts/build_known_fingerprints_*.py`). Anything that satisfies the
format and contains `{2, 3, 5}` as a subset (the v16 sieve floor for
`--immune-prime 5`) is a valid pool entry.

## Hardware tuning

### RTX 5090 + ~10-core CPU

```
--prove-threads 8
--prove-queue-depth 16
--streams 4
```

Expected: ~150 G(Q,V)/s aggregate. Will be CPU-prove-bound at this
core count.

### RTX 5090 + ~30-core CPU

```
--prove-threads 24
--prove-queue-depth 16
--streams 4
```

Expected: ~140–150 G(Q,V)/s aggregate, GPU at 85–89 % saturation,
CPU at low utilization (extended sieve filters most survivors).

### RTX 4090 (Ada, sm_89)

Build with `-arch=sm_89` in Makefile. Throughput likely 60–80 % of
RTX 5090, otherwise identical config.

## Common parameter pitfalls

- **`--depth > --target`**: the sieve silently filters chains of length
  `< depth` whose `p_target..p_{depth−1}` are small-prime-composite.
  You'll see lower hit rates without obvious errors. Always set
  `depth ≤ target` unless you know exactly what you're doing.
- **`--prove-threads > cores`**: oversubscription. With cores=30,
  threads=28 is fine; threads=32 is worse than threads=24. Leave a few
  cores for K1 dispatch and IO.
- **`--immune-prime 7` with the v16 default pool**: the pool was built
  for `--immune-prime 5`. Using 7 will reject every D_V that doesn't
  contain 7 — see `BUGS_AND_LESSONS.md` (this was a real bug we hit).
- **`--sieve-max > 503`**: the GPU forbid_mask is hard-capped at 512
  bits (`V16_Q_MASK_WORDS=8`). Higher sieve-max values will error at
  startup.
