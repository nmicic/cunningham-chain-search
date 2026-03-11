# cc_cuda index

**Latest file**: `00_cc_cuda_v10_tuned.cu`

## Functional spec
CUDA v10 applies per-base wheels, mod6 filtering, and tile-aligned scanning for dramatic performance gains over v08.

## Problem / capability
Speed up core constructor search on modern GPUs while preserving correctness for sequential/prefix-driven runs.

## Key methods / implementation notes
- Per-base wheels for primes 23, 29, 31
- Mod6 filtering baked into the k offsets
- Tile-aligned sequential search aligned to WHEEL_PERIOD boundaries
- Supports random and sequential modes

## Build / run
- **Compile**: `nvcc -O3 -arch=sm_86 cc_cuda_v10.cu -o cc_cuda_v10`
- **Usage tip**: Example: `./cc_cuda_v10 --target 20 --bits 107 --continuous` or `--prefix`/`--sequential` for controlled sweeps.

## Notes
GPU-oriented sibling to cc_constructor_hybrid; same constructors but tuned for high-throughput scanning.

Historical versions (2 files) in `archive/` subfolder.

## Cunningham relation
- **Status**: Core CC program
- **Why**: Implements Cunningham construction/search logic (CC1/CC2/BiTwin or direct CC-targeted sieving).
