# cc_constructor_hybrid index

**Latest file**: `00_cc_constructor_hybrid_v07_1.cu`

## Functional spec
Hybrid v07 refines earlier lattice+2-adic prototypes with safe bitmask ranges, extended trial division up to prime 409, and removal of failed experiments.

## Problem / capability
Deliver robust hybrid constructor behavior without the undefined shifts and lattice noise that plagued v06.

## Key methods / implementation notes
- Bitmask filter capped at primes ≤61 to avoid undefined behavior
- Trial division for primes 67‑409 inside is_prime_128()
- Lattice search executed alongside 2-adic guidance

## Build / run
- **Compile**: `nvcc -O3 -arch=sm_86 cc_constructor_hybrid_v07.cu -o cc_hybrid_v07`
- **Usage tip**: Use like other constructor builds; additional diagnostics around lattice scalars were added in this revision.

## Notes
Part of the broader constructor family; tracks the same goal as cc_constructor_algebraic_hybrid but with different heuristics.

**Also kept**: `02_cc_constructor_hybrid_v06_lattice.cu` — lattice-based filtering (~30,000x fewer primality tests), removed in v07.

Historical versions (34 files) in `archive/` subfolder.

## Cunningham relation
- **Status**: Core CC program
- **Why**: Implements Cunningham construction/search logic (CC1/CC2/BiTwin or direct CC-targeted sieving).
