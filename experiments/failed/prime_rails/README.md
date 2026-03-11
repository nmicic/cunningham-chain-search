# prime_rails index

**Latest file**: `00_prime_rails_v1-2__20260117-175813.cu`

## Functional spec
Prime Rails v1 follows 2n+1 rails (Cunningham chain second kind) with examples and length breakdowns.

## Problem / capability
Investigate sequences p → 2p+1 that stay prime and highlight lengths/interruptions.

## Key methods / implementation notes
- Logs lengths for rails, e.g. 2→5→11→23→47
- Notes relation to CC2 and safe primes
- Emphasizes interrupts when composites appear

## Build / run
- **Compile**: `nvcc -O3 -arch=sm_86 prime_rails_v1.cu -o prime_rails`
- **Usage tip**: Run once compiled to log rails and their lengths; cross-check with CC2 heuristics.

## Notes
Focuses on second-kind chain behavior; mathematical cousin to main Cunningham chain efforts.

Historical versions (1 file) in `archive/` subfolder.

## Cunningham relation
- **Status**: CC-related (Cunningham second kind, CC2)
- **Why**: Uses the explicit right-edge recurrence `p -> 2p+1`, which is CC2.
