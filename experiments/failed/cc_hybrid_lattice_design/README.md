# cc_hybrid_lattice_design index

**Latest file**: `00_cc_hybrid_lattice_design__20260122-124303.cu`

## Functional spec
Design document that explains why previous terminal-based searches miss chains like CC17 and describes the hybrid lattice fix.

## Problem / capability
Sequential tests for n ≡ 5 (mod 6) are too slow and fail to hit low-1-root chains whose lines still pass the sieve.

## Key methods / implementation notes
- Define lattice based on small prime set P = {5,7,11,13,17,19}
- Work with root lines rather than terminals
- Formulate constraints so the root lattice aligns with feasible 2-adic lines

## Build / run
- **Compile**: `Documentation-only; no build command inside the spec.`
- **Usage tip**: Use as HLD/architecture note when developing cc_constructor and lattice search variants.

## Notes
Essential background for understanding cc_constructor_hybrid and cc_lattice_search. Treat as high-level spec.

## Cunningham relation
- **Status**: Core CC program
- **Why**: Implements Cunningham construction/search logic (CC1/CC2/BiTwin or direct CC-targeted sieving).
