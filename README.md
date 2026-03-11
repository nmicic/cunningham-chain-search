# CC18: Cunningham Chain Search Engine

A visualization-driven project that grew into a high-throughput GPU/CPU pipeline for searching long Cunningham chains (first-kind). More about HPC optimization and AI-assisted iteration than the math itself.

## Results

As of March 2026, this project has contributed:
- New CC16 and CC17 entries to the public [Cunningham chain tables](https://www.pzktupel.de/CC/cc.php)
- The current smallest known second-kind CC17
- The current largest known first-kind CC17

The original target (first-kind CC18 or CC19) has not yet been reached. This was compute-limited.

The campaign dataset is documented in [`data/`](data/). Summary statistics and analysis CSVs are tracked in the repo; the full raw dataset (~30 MB, 1.77M roots including 44 CC16, a CC17, and ~929K CC10+) is available via GitHub Releases.

## Quick Links

| What | Where |
|------|-------|
| GPU filter engine (CUDA) | [`src/cuda/`](src/cuda/) |
| CPU search engine (GMP) | [`src/cpu/`](src/cpu/) |
| GP/PARI library + MCP server | [`gp/`](gp/) |
| Interactive visualizations | [`visualizations/`](visualizations/) |
| Analysis notes | [`analysis/notes/`](analysis/notes/) |
| Analysis scripts | [`analysis/scripts/`](analysis/scripts/) |
| Campaign dataset | [`data/`](data/) |
| Failed experiments (17 approaches) | [`experiments/failed/`](experiments/failed/) |
| Repo map | [`docs/REPO_MAP.md`](docs/REPO_MAP.md) |
| Search pipeline overview | [`docs/SEARCH_PIPELINE.md`](docs/SEARCH_PIPELINE.md) |

## Search Pipeline

Two-phase GPU+CPU pipeline:
1. **GPU sieve** (CUDA): modular-residue filtering at 57-65 billion candidates/sec (RTX 4090 / RTX 5090). Rejects 99.9988% of candidates.
2. **CPU prove** (GMP): BPSW primality proving on the 0.0012% that survive.

See [`docs/SEARCH_PIPELINE.md`](docs/SEARCH_PIPELINE.md) for details.

## Interactive Visualizations

Standalone HTML/JS tools that shaped the search design:

- [Cunningham Chain Mesh](visualizations/chain-mesh/) - 2-adic square-perimeter map with CC1/CC2 edges
- [2-Adic Tree Explorer](visualizations/2adic-tree/) - inspectable tree with local chain structure
- [3D Fold](visualizations/3d-fold/) - shell structure across levels
- [p+1 Analysis](visualizations/p1-analysis/) - factor coloring and mod-p views
- [Chain Analyzer](visualizations/chain-analyzer/) - single-chain analysis and breaker autopsy
- [Campaign Dashboard](visualizations/campaign-dashboard/) - campaign tracking
- [Immunization Dashboard](visualizations/immunization-dashboard/) - residue immunity analysis

## Build

```bash
# GPU engine (requires CUDA toolkit + GMP)
nvcc -O3 -arch=sm_89 src/cuda/cc18_filter_cuda_CpC_v13.cu -o cc18_filter -lgmp -lpthread

# CPU engine (requires GMP)
gcc -O3 -march=native -flto src/cpu/cc_gmp_v33_03.c -o cc_search -lgmp -lpthread -lm

# Run tests
./cc_search --test
```

## GP/PARI Library

```bash
gp -q
\r gp/cc_lib_v10.gp
```

Self-tests (37) run automatically on load.

See [`gp/HOWTO_cc_lib_v10.md`](gp/HOWTO_cc_lib_v10.md) for usage.

## Project Stack

- **Search:** CUDA filtering, GMP proving, prefix sharding, checkpointing
- **Visualization:** standalone HTML/JS tools
- **Analysis:** Python + HTML for autopsy and fingerprinting
- **AI-assisted:** LLMs used for coding iteration; math direction and search decisions are human-driven
- **Extra:** PARI/GP library and MCP server for interactive analysis

## External Links

- [Cunningham chain tables (pzktupel.de)](https://www.pzktupel.de/CC/cc.php)
- [CC16 history](https://www.pzktupel.de/CC/HCC16.php)
- [CC17 history](https://www.pzktupel.de/CC/HCC17.php)

## Author

Nenad Micic, Belgium — [LinkedIn](https://be.linkedin.com/in/nenadmicic)

## License

See [LICENSE](LICENSE).
