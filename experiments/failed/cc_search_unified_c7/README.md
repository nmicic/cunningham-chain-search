# cc_search_unified_c7 index

**Latest file**: `00_cc_search_unified_c7_v25.cu`

## Functional spec
CC7-focused unified search that logs CSV heatmaps, tracks every CC7+ chain, and reports best verified CC7s.

## Problem / capability
Observe the 2-adic heatmap of CC7+ activity while guaranteeing best-chain tracking and checkpoint serialisation.

## Key methods / implementation notes
- CSVLogger for per-batch v2/odd-core/chain counts
- Always tracks best CC7+ regardless of --min length
- Checkpointing + final verification of best chain

## Build / run
- **Compile**: `nvcc -O3 -arch=sm_86 cc_search_unified_c7_v25.cu -o cc_search_v25`
- **Usage tip**: Useful for detailed CC7 research; check CSV output for heatmap analysis.

## Notes
Direct companion to cc_search_unified; highlight: CC7+ focus, per-batch CSV logging, small carbohydrates.

**Also kept**: `01_cc_search_unified_c7_v27__20260122-121913.cu` — latest production (sieve-only line filter + lattice constructor).

Historical versions (31 files) in `archive/` subfolder.

## Cunningham relation
- **Status**: Core CC program
- **Why**: Implements Cunningham construction/search logic (CC1/CC2/BiTwin or direct CC-targeted sieving).
