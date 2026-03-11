# cc_search_merged index

**Latest file**: `00_cc_search_merged_v6__20260118-013418.cu`

## Functional spec
Search v6 with all CUDA fixes plus CPU fallback and unified kernels.

## Problem / capability
Eliminate CUDA errors from previous merged attempts while keeping CPU portability.

## Key methods / implementation notes
- Conditional CPU_ONLY sections for CPU builds
- Unified search logic for shared macros
- Robust error handling in host/device code

## Build / run
- **Compile**: `nvcc -O3 -arch=sm_86 cc_search_merged_v6.cu -o cc_search_v6 (or CPU path similar to cc_search).`
- **Usage tip**: Standard run once built; identical CLI to cc_search suite.

## Notes
Bridge between cc_search and cc_search_unified; should behave identically when configured as CUDA or CPU build.

Historical versions (2 files) in `archive/` subfolder.

## Cunningham relation
- **Status**: Core CC program
- **Why**: Implements Cunningham construction/search logic (CC1/CC2/BiTwin or direct CC-targeted sieving).
