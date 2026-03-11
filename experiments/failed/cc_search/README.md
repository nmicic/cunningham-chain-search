# cc_search index

**Latest file**: `00_cc_search_v7__20260118-100858.cu`

## Functional spec
Unified CUDA+CPU v7 search covering the essential death residue filtering, Miller-Rabin filters, and CPU fallbacks.

## Problem / capability
Locate long Cunningham chains while offering both GPU acceleration and CPU-only fallback for portability.

## Key methods / implementation notes
- Set of 15 sieve primes with death residues
- Deterministic Miller-Rabin (12 witnesses) after pre-sieve
- Optional CPU-only build with threading placeholders

## Build / run
- **Compile**: `CUDA: `nvcc -O3 -arch=sm_86 cc_search_v7.cu -o cc_search`; CPU: `g++ -O3 -DCPU_ONLY -x c++ cc_search_v7.cu -o cc_search -lpthread`.`
- **Usage tip**: Runs generically once built; see command-line help inside the binary.

## Notes
Baseline search engine; many later searches (merged/unified) derive from this core pipeline.

Historical versions (2 files) in `archive/` subfolder.

## Cunningham relation
- **Status**: Core CC program
- **Why**: Implements Cunningham construction/search logic (CC1/CC2/BiTwin or direct CC-targeted sieving).
