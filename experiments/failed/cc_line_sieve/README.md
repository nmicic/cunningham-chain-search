# cc_line_sieve index

**Latest file**: `00_cc_line_sieve.cu`

## Functional spec
Line extension sieve shifts focus from terminals to roots by ensuring the line of successive doublings stays prime.

## Problem / capability
Classic constructor searches miss roots with few trailing 1s but long line lengths; this solver targets those roots directly.

## Key methods / implementation notes
- Sieve eliminates most candidates (~1 in 12 million left) using ~170 modular checks
- Only remaining seeds hit expensive primality tests, accelerating chain discovery by ~256x

## Build / run
- **Compile**: `Assumed `nvcc -O3 -arch=sm_86 cc_line_sieve.cu -o cc_line_sieve` (not explicit).`
- **Usage tip**: Designed for investigation; follow the sieve-first approach rather than terminal enumeration.

## Notes
Still part of the constructor/search family; emphasizes line-based heuristics instead of bit-heavy lattice work.

## Cunningham relation
- **Status**: Core CC program
- **Why**: Implements Cunningham construction/search logic (CC1/CC2/BiTwin or direct CC-targeted sieving).
