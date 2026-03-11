# Search Pipeline

## Overview

The published search engine is a staged heterogeneous pipeline:

1. precompute forbidden residues for the first `d` chain positions
2. build a CRT-and-wheel search lattice
3. reject almost all candidates on the GPU with a 3-stage modular sieve
4. hand the tiny survivor set to the proving path
5. record confirmed roots and derive analysis datasets

The speed comes from avoiding expensive primality work almost all the time.

## Phase 1: Precomputed Forbidden Residues

For a first-kind Cunningham chain, a root `p` generates:

`p, 2p+1, 4p+3, 8p+7, ...`

For each small prime `q`, the host precomputes the root residues that would force one of the first `d` projected chain positions to be divisible by `q`. These forbidden-residue tables are built once at startup and then uploaded to the GPU.

In the published v13 engine, the filtering primes are organized in three layers:

- CRT primes: `5, 7, 11, 13, 17, 19`
- wheel primes: `23, 29, 31`
- chain-screen primes for filtering and proving support: `37..863`

This is the logic behind `--depth`: the code is not asking only whether `p` is prime, but whether any of the first `d` chain positions are already doomed modulo tracked small primes.

**Implementation:** `src/cuda/cc18_filter_cuda_CpC_v13.cu`

## Phase 2: CRT Bases and Wheel Lattice

The engine does not scan all odd integers.

First, it enumerates CRT-surviving base residues modulo:

`LATTICE_M = 5 * 7 * 11 * 13 * 17 * 19 = 1,616,615`

Then, for each surviving base, it generates a wheel over:

`WHEEL_PERIOD = 23 * 29 * 31 = 20,677`

So the concrete candidate parameterization is:

`n = base + (tile * WHEEL_PERIOD + k_offset) * LATTICE_M`

where:

- `base` is one of the CRT-surviving residues
- `tile` is the coarse batch/sweep coordinate
- `k_offset` is a wheel position surviving the wheel primes

The code also partitions wheel positions into `mod 6` buckets so each `(tile, base)` pair only scans the bucket that can produce candidates in the desired `n ≡ 5 (mod 6)` class.

Two coverage modes exist:

- `full`: all CRT-surviving bases
- `legacy`: only the older `r % 6 == 5` subset

The published release defaults to `full`.

## Phase 3: Search Partitioning

Before each run, the lattice can be narrowed by:

- target bit size
- binary prefix (`--prefix`)
- prefix mode (`sequential` or `random`)
- prefix lanes (`--prefix-lanes`, `--prefix-lane-id`)

This is how the same engine is used for:

- full-range sequential sweeps
- prefix-constrained experiments
- multi-GPU non-overlapping lane scans

Sequential prefix mode stops after one complete sweep of its assigned prefix range. Random mode keeps sampling until a time limit or manual stop.

## Phase 4: GPU Filter Kernel

The hot path is Kernel 1 in `src/cuda/cc18_filter_cuda_CpC_v13.cu`.

Each block works on one `(tile, base)` work unit. Candidate wheel positions then pass through three rejection stages:

1. **L2 bitmask filter**
   Small-prime bitmask tests for `37..61`
2. **ext-L2 byte filter**
   Byte-table tests for `67..97`
3. **line sieve**
   A deeper per-candidate sieve for `101..863`

Only candidates that pass all three stages are emitted as survivors.

Important implementation details:

- per-tile residues are computed once in shared memory
- wheel positions are already bucketed by `mod 6`
- the large `kM_mod` lookup table is prefetched only after the first two filters pass
- survivor emission is warp-aggregated to reduce atomic traffic

This is the stage that delivers the campaign-scale throughput.

## Phase 5: Double-Buffered Handoff

The engine uses two CUDA streams and two copies of the batch/survivor buffers.

This lets it overlap:

- batch setup and upload
- GPU filtering
- survivor copy-back
- proving of the previous batch

In practice, the pipeline runs as:

- buffer A filtering while buffer B is being proved
- then swap

That overlap is essential for keeping the GPU saturated.

## Phase 6: Proving and Chain Confirmation

Two proving paths exist in v13:

- **GPU prove kernel**
  Fast, but incomplete for non-root handling and subject to correctness limits at larger bit/target combinations
- **CPU prove path**
  The recommended path for correctness and the path used for the published snapshot/analysis workflow

The CPU path does the following for each survivor:

1. reconstruct `n` from `(base_idx, tile_idx, k_offset)`
2. run trial division with the chain-screen primes
3. run GMP probable-prime testing on the candidate
4. walk backward if `(n-1)/2` is also prime, to find the true root
5. follow the chain forward from the true root to measure the full chain length

This is also where non-roots are identified and logged.

The published campaign path favors CPU confirmation because it preserves:

- true-root recovery
- non-root logging
- better post-campaign analysis visibility

**Implementations:**

- GPU prove path: `src/cuda/cc18_filter_cuda_CpC_v13.cu`
- CPU reference engine: `src/cpu/cc_gmp_v33_03.c`

## Phase 7: Output and Analysis

Confirmed roots are written out and later turned into the released analysis material, including:

- gap statistics
- kill-position summaries
- immunization / residue summaries
- `p+1` breaker analysis
- ghost-chain census
- twins / clusters

This is why the published code should be read as a campaign-and-analysis version, not only as a raw speed benchmark.

## Published vs Alternate Variants

The released code reflects the campaign version used to generate the published snapshot and derived analysis.

Some later or alternate optimizations were explored separately. One example is a reverse-check proving path (`cc18_filter_cuda_CpC_v14.cu`), which can reduce CPU prover work by checking from higher chain positions inward, leveraging lower prime density in the outer shells. That path was not used for the published snapshot because the forward path was better suited for logging, inspection, and post-campaign analysis.

## Performance

| Metric | Value |
|--------|-------|
| GPU throughput | 57-65 B candidates/sec |
| GPU survivor rate | 0.0012% |
| GPU rejection rate | 99.9988% |
| Primary search band | 89-91 bit |
| Published data snapshot | 929,574 roots |
