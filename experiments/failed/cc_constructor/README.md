# cc_constructor index

**Latest file**: `00_cc_constructor_v20__20260122-035258.cu`

## Functional spec
Cunningham Chain Constructor v20 tracks the k-value for every chain and applies GMP deep analysis for promising k=1 roots, with the ability to split fragments (k>1) from true starter chains.

## Problem / capability
Generate and validate Cunningham chain starters with insight into which roots can theoretically support CC20s.

## Key methods / implementation notes
- GPU/kernel search plus optional GMP verification for long k=1 chains
- statistics split between k=1 and k>1 fragments
- options such as --k1-only, --2adic-pattern, --ext-bits, and periodic tables carry over from earlier versions

## Build / run
- **Compile**: `nvcc -O3 -arch=sm_86 -DUSE_GMP cc_constructor_v20.cu -o cc_constructor -lgmp`
- **Usage tip**: Use flags like --k1-only to focus on true CC20 starters, or rely on default search for mixed fragments.

## Notes
Shares lineage with cc_constructor_hybrid, cc_constructor_algebraic_hybrid, and the CUDA-specific/bit-lattice variants; keeps the core constructor workflow.

**Also kept**: `12_cc_constructor_v13__20260121-023739.cu` — BiTwin + second-kind chain support not present in v20.

Historical versions (36 files) in `archive/` subfolder.

## Cunningham relation
- **Status**: Core CC program
- **Why**: Implements Cunningham construction/search logic (CC1/CC2/BiTwin or direct CC-targeted sieving).
