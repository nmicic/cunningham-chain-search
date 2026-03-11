# cc_lattice_search index

**Latest file**: `00_cc_lattice_search_1.cu`

## Functional spec
Implements the 2-adic grid sieve from CC20_LATTICE_SPEC to treat chain positions as bit-walks with precomputed exclusion masks.

## Problem / capability
Avoid Miller-Rabin on the many excluded candidates by modeling the chain path on a 2-adic lattice instead of terminal 1s.

## Key methods / implementation notes
- 2-adic grid where bits represent positions; chains follow bit shifts and ORs
- Precomputed composite patterns for small primes
- Mask operations to keep the search within feasible grid cells

## Build / run
- **Compile**: `nvcc -O3 -arch=sm_86 cc_lattice_search.cu -o cc_lattice_search`
- **Usage tip**: Example run: `./cc_lattice_search --target 20 --bits 100 --prime-bound 10000 --continuous`.

## Notes
Relates to cc_hybrid_lattice_design (spec) and complements cc_constructor_hybrid by supplying lattice-based candidate generation.

Historical versions (3 files) in `archive/` subfolder.

## Cunningham relation
- **Status**: Core CC program
- **Why**: Implements Cunningham construction/search logic (CC1/CC2/BiTwin or direct CC-targeted sieving).
