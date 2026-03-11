# cc_search_cuda index

**Latest file**: `00_cc_search_cuda_1__20260117-230736.cu`

## Functional spec
Early CUDA search applying z-slice parallelism, death tables, and coalesced memory access for n ≡ 5 (mod 6) candidates.

## Problem / capability
Bring the search to GPU slices while minimizing garbage candidates via modular pruning.

## Key methods / implementation notes
- Restricted to n≡5 (mod 6)
- Pre-sieved death tables
- Parallel search across Z-slices and coalesced CUDA access patterns

## Build / run
- **Compile**: `nvcc -O3 -arch=sm_70 cc_search_cuda.cu -o cc_search_cuda`
- **Usage tip**: Straightforward GPU run once compiled; targeted at HPC verification.

## Notes
Pioneer for cc_search_unified* derivatives, but retains simpler structure for slicing.

Historical versions (2 files) in `archive/` subfolder.

## Cunningham relation
- **Status**: Core CC program
- **Why**: Implements Cunningham construction/search logic (CC1/CC2/BiTwin or direct CC-targeted sieving).
