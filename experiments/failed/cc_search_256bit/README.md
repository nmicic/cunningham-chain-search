# cc_search_256bit index

**Latest file**: `00_cc_search_256bit_v3.cu`

## Functional spec
256-bit search handles giant targets without relying on __uint128_t by using Montgomery, bitwise divisibility, and explicit carry management.

## Problem / capability
Search for chains that exceed 64 bits without dangerous intrinsics, while still keeping modular arithmetic efficient.

## Key methods / implementation notes
- Multiple modulus elimination tricks (mod 3, 5, 7, 11, 13) on bytes
- Warp-level reductions via `__shfl_down_sync`
- Montgomery arithmetic inside Miller-Rabin
- Precomputed pow tables and CUDA error checking

## Build / run
- **Compile**: `nvcc -O3 -arch=sm_86 cc_search_256bit_v3.cu -o cc_search_256`
- **Usage tip**: `./cc_search_256 --test` validates the pipeline before production runs.

## Notes
Search-focused cousin to constructor efforts; useful once targets grow beyond 64 bits.

Historical versions (13 files) in `archive/` subfolder.

## Cunningham relation
- **Status**: Core CC program
- **Why**: Implements Cunningham construction/search logic (CC1/CC2/BiTwin or direct CC-targeted sieving).
