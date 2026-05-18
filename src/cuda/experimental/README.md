# cc20-first-kind-immune — experimental Cunningham Chain search engine

A CUDA engine for hunting **Cunningham chains of the first kind** at
CC20/CC21 scale. Successor in spirit to
`cc18_filter_cuda_CpC_v15.cu`, but with a fundamentally different
architecture: instead of one sieve shape, it iterates Q-values across
~4,200 fingerprint shapes drawn from the **empirical distribution of
real CC10+ chain roots**.

## What it does, in one paragraph

For each variant `V` (a small-prime fingerprint `D_V` taken from the
historical CC10+ catalog), it samples a value `Q` in the configured
bit band and tests whether `p = Q · D_V − 1` starts a Cunningham chain
of the first kind. The construction guarantees small primes never
divide `p` (because `D_V | (p+1)`), so the candidates pass a strong
structural filter before the engine spends a cycle on Miller-Rabin.
Throughput at the current optimum: **~150 G(Q,V)/s per RTX 5090**.

## Why this exists

Open Cunningham-chain search code is conspicuously absent. Existing
record announcements publish numbers, not implementations. The
methodology here — fingerprint pool derived from a real chain corpus,
external-prefix structural validator, GPU sieve + CPU prove pipeline —
is the part worth contributing.

If you use this engine to find a new chain, please:

1. Submit it to Norman Luhn's Cunningham Chain table:
   <https://www.pzktupel.de/CC/cc.php>. Contact:
   `pzktupel [at] pzktupel [dot] de`.
2. Cite the engine version + fingerprint-pool snapshot SHA in the
   submission.
3. Include an independent verification of the chain root and length.

## Quick start

```bash
# 1. Build (requires CUDA 13.x, GCC supporting sm_120 for Blackwell/RTX 5090)
make -j

# 2. Run the long-run search preset on GPU 0, two prove threads per core,
#    output to ./run/ — runs forever until Ctrl+C
./scripts/run_cc20_hunter.sh

# 3. Watch hits roll in
tail -F run/*/hit.log
```

## At a glance

| | |
|---|---|
| **Target** | Cunningham chains of the first kind, length ≥ 10 |
| **Default bit band** | `p_0 ∈ [2^89, 2^107 − 1]` (90–107 bit roots) — keeps the target-21 chain top `< 2^127`, inside `u128` |
| **Default depth** | 21 (filters small-prime composites at chain positions 0..20) |
| **Pool size** | 4,248 fingerprints (96.78 % coverage of historical CC10+); plus a primorial-quotient pool (3,000 V's) for exhaustive-mode runs |
| **Hardware** | One CUDA GPU (RTX 5090 reference); 8+ CPU cores for prove |
| **License** | Apache-2.0 |

## Files in this folder

| | |
|---|---|
| `README.md` (this file) | What it is, why, quick start |
| `HOWTO.md` | Step-by-step: build, run, interpret, report findings |
| `CONFIGURATION.md` | CLI reference, recommended configs for different hardware |
| `DESIGN.md` | Architecture deep-dive: q-iteration construction, GPU pipeline, comparison vs v15 |
| `BUGS_AND_LESSONS.md` | Postmortems on bugs we hit so you can avoid them |
| `LICENSE` | Apache-2.0 |
| `Makefile` | `make -j` produces `bin/cc20_first_kind_immune_v16` |
| `cc20_first_kind_immune_v16.cu` | The engine. ~5300 lines, CUDA + GMP |
| `pools/` | Fingerprint pool files (text, human-readable) |
| `scripts/` | Pool generators + the production wrapper `run_cc20_hunter.sh` |
| `tools/snapshot_replay.py` | External-prefix structural validator |
| `data/` | Fetch script + README for the CC10+ root snapshot. The CSV itself lives in the [`cunningham-chain-data`](https://github.com/nmicic/cunningham-chain-data) release; only needed if you want to regenerate pools or run `snapshot_replay.py` |
| `examples/` | One-command tmux launchers for each pool size — read this first |

## Status

**Experimental.** This is research code under active iteration. The
default configuration finds chains; specific parameters and even the
file/binary basename may change between commits.
