# HOWTO — build, run, observe, report

End-to-end walkthrough. Assumes you have a CUDA-capable GPU and a Linux
host. RTX 5090 (Blackwell, sm_120) is the reference; older cards work
with a sm-target change.

## 1. System requirements

| | |
|---|---|
| **GPU** | CUDA-capable, ≥ 8 GiB VRAM. RTX 5090 / Blackwell ideal; RTX 4090 / Ada works (change `-arch=sm_120` → `sm_89`). |
| **CPU** | 8+ cores. The CPU prove path scales linearly up to ~24 threads on a 30-core box. |
| **RAM** | 4 GiB free. Engine uses ~150 MiB host + ~80 MiB ext-sieve forbid table at the 4248-pool default. |
| **CUDA** | 13.x or newer (older may work with `-arch` tweaks). |
| **GMP** | `libgmp-dev`. Required for the CPU prove path. |

## 2. Build

```bash
# Edit Makefile if your GPU isn't Blackwell (change sm_120 → sm_89 for Ada,
# sm_86 for Ampere, etc.)
make -j

# Produces ./bin/cc20_first_kind_immune_v16

./bin/cc20_first_kind_immune_v16 --help | less
```

## 3. Verify your build

The engine has a startup self-test that runs by default. It asserts the
chain-follow machinery is correct against 11 known Cunningham chains
(CC5..CC14, CC17 Wroblewski 2008) and 5 q-iter round-trip vectors
including both known CC18 world records. Run with a minimal config:

```bash
./bin/cc20_first_kind_immune_v16 \
    --q-iter --immune-prime 5 \
    --seed-pool-file pools/known_fingerprints_exp_d.txt \
    --target 10 --depth 10 --sieve-max 211 \
    --bits-min 86 --bits-max 90 \
    --cpu-prove --prove-threads 2 \
    --time 1 --gpu 0
```

You should see:

```
=== v16 startup self-test ===
  PASS [CC5  classic (2->5->11->23->47)]: CC5 verified (2 bits)
  ...
  PASS [q-iter v0]: CC18 round-trip (88 bits, D_V in pool)
  ...
INFO: self-test PASSED (16/16 vectors)
```

If any vector FAILS, **do not run a real campaign**. Open an issue; the
chain-follow path is broken on your toolchain.

To skip the ~30 s self-test in hot-restart loops: `--no-test`.

## 4. Run the long-run search preset

The default-recommended configuration for CC21-scale searching (current
target; the engine supports any `--target N` ≤ 24):

```bash
./bin/cc20_first_kind_immune_v16 \
    --q-iter \
    --seed-pool-file pools/known_fingerprints_exp_d.txt \
    --immune-prime 5 \
    --target 21 --depth 21 \
    --sieve-max 503 \
    --bits-min 90 --bits-max 107 \
    --q-order random \
    --cpu-prove --prove-threads $(($(nproc) - 4)) \
    --prove-queue-depth 16 \
    --streams 4 \
    --min-report-len 10 \
    --time 0 \
    --gpu 0 \
    --log run/hit.log \
    --checkpoint run/checkpoint.txt
```

Or just use the wrapper:

```bash
./scripts/run_cc20_hunter.sh
```

The wrapper picks sensible defaults, creates a timestamped `run/`
directory, and starts. Edit it to suit your hardware.

> **Why `--bits-max 107` at target 21?** For a chain of length 21 the top
> prime is `p_20 ≈ 2^20 · p_0`. With `bits-max 107` the top stays
> `≤ 2^127` — at the tight edge of `u128` but inside it (and the engine
> validates this at startup). The wider `[90, 127]` band is mathematically valid but
> exercises code paths near the `u128` boundary that have no synthetic
> CC20 fixture covering them yet. Narrow to 107 unless you've added that
> fixture.

For the fastest hit-flow validation (chains in seconds, not hours), drop
the band to `--bits-min 86 --bits-max 90`. You'll find CC10+ in a
minute or two on most GPUs — useful to confirm the prove path is wired
up correctly. Chains found there are not records, just sanity checks.

### Exhaustive-mode (primorial-pool) preset

Primorial-quotient pools (e.g. `pools/known_fingerprints_primorial.txt`)
mix variants whose A-range can be exhausted in seconds (high-V, narrow
Q-range) with variants whose A-range is astronomical. Use exhaustive
band-mode so the engine picks per V:

```bash
./bin/cc20_first_kind_immune_v16 \
    --q-iter \
    --seed-pool-file pools/known_fingerprints_primorial.txt \
    --immune-prime 5 \
    --target 21 --depth 21 \
    --sieve-max 503 \
    --bits-min 90 --bits-max 107 \
    --q-band-mode exhaustive --exhaustive-max-q-bits 40 \
    --cpu-prove --prove-threads $(($(nproc) - 4)) \
    --prove-queue-depth 16 \
    --streams 4 \
    --min-report-len 10 \
    --time 0 --gpu 0 \
    --log run/hit.log \
    --checkpoint run/checkpoint.txt
```

Banner reports the per-V split, e.g.
`q-band-mode: exhaustive threshold=2^40 Q-range (sequential V's=1250
random V's=1750)`. The sequential V's complete coverage of their entire
A-range; the random V's continue probabilistic sampling. `--dedup-p0`
(default ON) collapses repeat-`p_0` emits — both random with-replacement
collisions and same-`p_0` multi-V framings — to compact `HIT-DUP` lines.

## 5. Reading the periodic report

Every 30 seconds the engine prints a line like:

```
[t=300.00s] launches=1585882 Q_eval=26606684864512 (88.69G (Q,V)/s)
  survivors=508155542 gpu=68% prove/k1=21.7x s/c=1.91e-05
  hits>=10=0 emitted=0 dups=0 active_V=4248/4248
```

What it means:

| field | meaning |
|---|---|
| `t=` | wall-clock seconds since start |
| `launches=` | total K1 kernel launches |
| `Q_eval=` | total (Q, variant) pairs evaluated |
| `(N G (Q,V)/s)` | aggregate throughput (Q × variant per second). **Not raw Q rate** — divide by active variants for per-variant Q/s. |
| `survivors=` | Q-candidates that survived the GPU sieve |
| `gpu=` | GPU utilization. 80–90 % is healthy and means GPU is the bottleneck |
| `prove/k1=` | CPU prove total time / K1 wall time. <1× means prove is keeping up; >4× means prove is the bottleneck |
| `s/c=` | survivors per K1 candidate, the sieve survival rate |
| `hits>=N=` | chain prove successes ≥ N (every chain prove counts, including dedup-collapsed repeats) |
| `emitted=` | distinct chains by `p_0` written verbose to `hit.log` (the `HIT seed=...` lines) |
| `dups=` | repeat-`p_0` sightings written as compact `HIT-DUP` lines. Total log lines = `emitted + dups`. Always 0 if `--no-dedup-p0`. |
| `active_V=N/M` | active variants of total. Variants exhaust when their Q-window is consumed |

## 6. Reading a HIT

`hit.log` lines look like:

```
HIT seed=q_iter_v18_top10_s1 len=12 pool=[2*3*5*11*13*19, 2*3^2*5*11*13*19, ...]
    A=0x... p_0=0x1bf67f52a82539b77a5e641
```

| field | meaning |
|---|---|
| `seed=` | run-tag from `--seed-name` |
| `len=` | confirmed chain length (Miller-Rabin verified) |
| `pool=[...]` | the variants whose `D_V` divides `p_0+1` — the chain's full fingerprint |
| `A=` | the `Q` value (in hex) that produced this hit |
| `p_0=` | the chain root in hex. The chain is `p_0, 2·p_0+1, 4·p_0+3, ...` |

A `HIT-DUP` line is the dedup-p0 short form for a chain whose `p_0` was
already reported verbose earlier in this run:

```
HIT-DUP seed=q_iter_v18_top10_s1 len=12 A=0x... p_0=0x1bf67f52a82539b77a5e641
```

Same `p_0` arrives twice either because random Q-order resampled it, or
(more commonly with primorial-quotient pools) because multiple V's frame
the same chain. Disable with `--no-dedup-p0` if you want the full verbose
emit every time.

## 7. Verify a chain found by the engine

```bash
python3 << EOF
p = int("0x1bf67f52a82539b77a5e641", 16)  # paste p_0 from hit.log
from sympy import isprime
k = 0
while isprime(p):
    print(f"p_{k} = 0x{p:x}  ({p.bit_length()} bits)")
    p = 2*p + 1
    k += 1
print(f"chain length = {k}")
EOF
```

## 8. Verify your build's structural coverage

The `tools/snapshot_replay.py` script does an **external prefix test**:
for every chain root in `data/cc10plus_roots_snapshot_2026-03-19.csv`,
it computes the (Q, D_V) decomposition and checks whether the loaded
pool covers that fingerprint.

The CSV is not bundled — fetch it from the data release first:

```bash
bash data/fetch_data.sh        # ~5 MB gzipped, ~52 MB uncompressed
python3 tools/snapshot_replay.py
```

You should see ~96.78 % coverage on the default pool. Significantly
lower coverage means your pool was built incorrectly — see
`BUGS_AND_LESSONS.md` for the canonical "we silently lost 76 % of
chains" pool-builder bug.

## 9. Report a chain you found

1. **Verify it externally.** Run the python snippet above to confirm
   chain length.
2. **Submit to the Cunningham Chain table:** Norman Luhn's CC table is at
   <https://www.pzktupel.de/CC/cc.php>. Contact:
   `pzktupel [at] pzktupel [dot] de`.
3. **Cite engine + pool:**
   > Found with `cc20_first_kind_immune_v16` (commit `<git-sha>`),
   > fingerprint pool `known_fingerprints_exp_d.txt` SHA-256
   > `<pool-sha>`, q-order random, seed `<value-from-banner>`.

4. **Keep the verification details.** Include the chain root, length,
   independent primality check, and any run parameters needed to reproduce
   the search context.

## 10. Troubleshooting

| symptom | check |
|---|---|
| `wheel size N exceeds cap` | usually `--immune-prime` is too high for current `D` config. Try `--immune-prime 5` |
| `gpu=` stuck at <30 % | CPU prove is starved — increase `--prove-threads` |
| `prove/k1=` > 10× | CPU prove is the bottleneck — reduce `--prove-threads` (avoid oversubscription), or raise `--prove-queue-depth` |
| 0 hits after 1 h at narrow bit band | check band density — `[90,107]` is sparse, expect tens of minutes per CC10+ hit; try `--bits-min 86 --bits-max 90` for fast hit-flow validation |
| Self-test fails | DO NOT run. Filing an issue with the failing vector helps everyone |
