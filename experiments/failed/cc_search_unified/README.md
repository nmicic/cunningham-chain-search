# cc_search_unified index

**Latest file**: `00_cc_search_unified_v10.cu`

## Functional spec
Latest unified search harness (CUDA / CPU agnostic) emphasizing determinism and the exact Miller-Rabin witness set for robust chain discovery.

## Problem / capability
Keep search code maintainable while supporting CPU fallback, 4-bit grid logging, and deterministic MR coverage.

## Key methods / implementation notes
- Constant memory for sieve primes and death residues
- Deterministic witness set for 64-bit numbers
- Data structures for result tracking and atomic best-chain selection

## Build / run
- **Compile**: `nvcc -O3 -arch=sm_86 cc_search_unified_v10.cu -o cc_search_unified_v10` (same as other CUDA builds).
- **Usage tip**: Focus on structured stats and GPU reliability; refer to CLI for specific flags.

## Notes
Shares DNA with cc_search_unified_c7 but without the dedicated CC7 logging pipeline; part of the same search family.

Historical versions (14 files) in `archive/` subfolder.

## Cunningham relation
- **Status**: Core CC program
- **Why**: Implements Cunningham construction/search logic (CC1/CC2/BiTwin or direct CC-targeted sieving).
