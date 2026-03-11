# cc_constructor_algebraic_hybrid index

**Latest file**: `00_cc_constructor_algebraic_hybrid_v05b.cu`

## Functional spec
Algebraic hybrid v05 fuses CRT line-lattices, primorial-based k forms, and 2-adic patterns to attack chains that historically came from k × P# × 2^m - 1.

## Problem / capability
Revisit record-setting k values by combining algebraic (primorial) structure with lattice sieving to reduce wasted primality tests.

## Key methods / implementation notes
- Primorial scaffolding (k = a × P# × 2^s) to reuse known fertile residues
- Bitmask-based CRT lattice filters up to 200+ primes
- 2-adic hints and low-overhead trial division

## Build / run
- **Compile**: `nvcc -O3 -arch=sm_86 cc_constructor_algebraic_hybrid_v05.cu -o cc_alg_hybrid_v05`
- **Usage tip**: Uses the same option vocabulary as other constructors (see cc_constructor). Ideal for k deep dives with primorial heuristics.

## Notes
Variant of the constructor family focused on algebraic+hybrid heuristics; sits next to cc_constructor_hybrid and cc_cuda branches.

Historical versions (1 file) in `archive/` subfolder.

## Cunningham relation
- **Status**: Core CC program
- **Why**: Implements Cunningham construction/search logic (CC1/CC2/BiTwin or direct CC-targeted sieving).
