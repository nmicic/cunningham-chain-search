# cunningham index

**Latest file**: `00_cunningham_algebraic_v3.cu`

## Functional spec
128-bit fixed chain follow-up that allows searches of the form k × P# × 2^m - 1 beyond 64-bit limits.

## Problem / capability
Previous builds capped at 2^63; need ability to follow chains requiring 128-bit arithmetic.

## Key methods / implementation notes
- Uses custom uint128 type with hi/lo helpers
- Performs Miller-Rabin over 128-bit numbers
- Supports prod-grade runs with --target-len and --m-end flags

## Build / run
- **Compile**: `nvcc -O3 -arch=sm_86 cunningham_algebraic_v3.cu -o cunningham_alg3`
- **Usage tip**: `./cunningham_alg3 --test` for smoke, `--target-len 12 --m-end 50` for production.

## Notes
Provides reference for algebraic chain forms; complements cc_constructor_algebraic_hybrid.

Historical versions (11 files) in `archive/` subfolder.

## Cunningham relation
- **Status**: Core CC program
- **Why**: Implements Cunningham construction/search logic (CC1/CC2/BiTwin or direct CC-targeted sieving).
