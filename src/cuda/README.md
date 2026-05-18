# `src/cuda/` — CUDA search engines

Three production engines and one experimental branch live here. They
share the high-level idea — a GPU sieve feeds candidates into a CPU
Miller-Rabin prove path — but differ in what gets sieved and how the
hot kernel is structured.

## The variants

### `cc18_filter_cuda_CpC_v13.cu` — campaign engine, first-kind
Production code for the 2026 first-kind campaign. Produced the
**1,551,489 CC10+ roots** published as the
[`v2026-03-19-snapshot`](https://github.com/nmicic/cunningham-chain-data/releases/tag/v2026-03-19-snapshot)
in the [`cunningham-chain-data`](https://github.com/nmicic/cunningham-chain-data)
repo, including the **two CC18s that are currently the largest
first-kind entries listed on the public
[Cunningham chain tables](https://www.pzktupel.de/CC/cc.php)**
(88-bit root `214325014495971624590189129`, 16 Mar 2026, and
87-bit root `106103983461039119546815109`, 14 Mar 2026).
Three-stage modular-depth filter:
CRT base residues mod {5,7,11,13,17,19} → packed kmod sieves for
inner primes → wheel-lattice candidate enumeration. Steady ~57–65 G
candidates/sec on RTX 4090 / 5090 at depth 20. See
[`docs/SEARCH_PIPELINE.md`](../../docs/SEARCH_PIPELINE.md) for the
full pipeline description.

### `cc18_filter_cuda_CpC_v13b.cu` — second-kind specialization
Same architecture as v13, retargeted to **second-kind** chains
(`p_{k+1} = 2·p_k − 1`). The candidate generation and CRT layer carry
over almost unchanged; the prove path flips the recursion direction.
Found a handful of small CC16 second-kind entries. The bit-window
sweep runner is [`tools/run_cc17b_sequential.sh`](../../tools/run_cc17b_sequential.sh).

### `cc18_filter_cuda_CpC_v15.cu` — throughput-optimized first-kind
Same three-stage filter shape as v13, but the hot path is rebuilt for
Blackwell-era cache hierarchies: interleaved `kmod` tables for the L2
and ext-L2 stages, constant-memory ext-L2 bitmasks, packed line-sieve
kill masks, and cooperative per-prime residue setup inside the filter
kernel. On RTX 5090 at depth 20 in CC19-style configs it has crossed
**~100 B candidates/sec** (96–98 G/s steady). Treat it as an
optimization branch — the published campaign baseline is still v13.

### `experimental/` — fingerprint-pool Q-iteration (v16)
A different sampling strategy. Instead of one CRT shape per run, v16
iterates Q-values across ~4,200 small-prime fingerprint shapes drawn
from the **empirical distribution of real CC10+ chain roots** in
[`cunningham-chain-data`](https://github.com/nmicic/cunningham-chain-data).
Targets CC20/CC21 scale. Methodology is the contribution; records, if
any, would be a bonus. At depth 21 with reverse-prove and ~32 CPU
prove threads it has crossed **~145 G (Q,V)/s** in favorable configs
(note: `(Q, variant)` pairs/sec, not raw Q rate — divide by active
variants for per-Q throughput). See
[`experimental/README.md`](experimental/README.md) for the elevator
pitch and [`experimental/DESIGN.md`](experimental/DESIGN.md) for the
architecture deep-dive vs v15.

## Build

Each engine ships with a recommended `-arch` line and external deps
(`libgmp`, `pthread`). Quick reference for Blackwell / RTX 5090:

```bash
# v13 (Ada / sm_89 — adjust for Blackwell)
nvcc -O3 -arch=sm_89 src/cuda/cc18_filter_cuda_CpC_v13.cu \
     -o cc18_filter -lgmp -lpthread

# v13b (second-kind)
nvcc -O3 -arch=sm_89 src/cuda/cc18_filter_cuda_CpC_v13b.cu \
     -o cc18_filterb -lgmp -lpthread

# v15 (Blackwell)
nvcc -O3 -arch=sm_120 src/cuda/cc18_filter_cuda_CpC_v15.cu \
     -o cc18_filter_v15 -lgmp -lpthread

# experimental v16 — uses its own Makefile
cd src/cuda/experimental && make -j
```
