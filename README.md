# CC18: Cunningham Chain Search Engine

> **Update (March 14–17, 2026):** Just days after publishing this repository, continued runs found two first-kind CC18s, starting with `106103983461039119546815109` and `214325014495971624590189129`. A later prior-art review showed that first-kind CC18 had already been documented in John Armitage's 2021 Oxford thesis via the smallest known example, so these are not the first known CC18s. They do, however, remain the largest listed first-kind CC18 entries on the current public Cunningham tables. The campaign has now finished because my available GPU compute time ran out.
>
> The text below is the original public release README from March 11, kept largely intact as the baseline snapshot. The March 11 dataset and release counts are unchanged here and do not include the later CC18 findings.

A visualization-driven project that grew into a high-throughput GPU/CPU pipeline for searching long Cunningham chains (first-kind). More about HPC optimization and AI-assisted iteration than the math itself.

## Results

As of March 2026, this project has contributed:
- New CC16 and CC17 entries to the public [Cunningham chain tables](https://www.pzktupel.de/CC/cc.php)
- The current largest listed first-kind CC16 and CC17 entries on the public tables

The campaign ultimately did reach its CC18 target, but not the more ambitious CC19 goal. It remained compute-limited.

For future searches, this code is more efficient at higher target/depth settings (for example CC19-oriented runs) than in the CC18-oriented configuration used here, because stricter filtering reduces survivor pressure on the downstream confirmation path.

The campaign dataset lives in a separate repo: [cunningham-chain-data](https://github.com/nmicic/cunningham-chain-data). Summary statistics and analysis CSVs are in [`data/`](data/) here; the full raw dataset (~30 MB, 929K CC10+ roots including 44 CC16 and a CC17) is available as a [release download](https://github.com/nmicic/cunningham-chain-data/releases/tag/v2026-03-12-snapshot).

## Quick Links

| What | Where |
|------|-------|
| GPU filter engine (CUDA) | [`src/cuda/`](src/cuda/) |
| CPU search engine (GMP) | [`src/cpu/`](src/cpu/) |
| GP/PARI library + MCP server | [`gp/`](gp/) |
| Interactive visualizations | [`visualizations/`](visualizations/) |
| Analysis notes | [`analysis/notes/`](analysis/notes/) |
| Analysis scripts | [`analysis/scripts/`](analysis/scripts/) |
| Campaign dataset | [cunningham-chain-data](https://github.com/nmicic/cunningham-chain-data) (separate repo) |
| Analysis CSVs | [`data/`](data/) |
| Failed experiments (17 approaches) | [`experiments/failed/`](experiments/failed/) |
| Repo map | [`docs/REPO_MAP.md`](docs/REPO_MAP.md) |
| Search pipeline overview | [`docs/SEARCH_PIPELINE.md`](docs/SEARCH_PIPELINE.md) |

## Search Pipeline

High-level staged pipeline:
1. **Search lattice**: enumerate CRT-surviving bases and wheel positions on the published first-kind lattice.
2. **GPU filter** (CUDA): three-stage modular depth filtering at 57-65 billion candidates/sec (public v13 baseline on RTX 4090 / RTX 5090). Rejects 99.9988% of candidates. Experimental CUDA branch [`src/cuda/cc18_filter_cuda_CpC_v15.cu`](src/cuda/cc18_filter_cuda_CpC_v15.cu) has reached roughly 96-98B candidates/sec on RTX 5090 in some CC19-style runs, but it still needs more validation before replacing the public baseline.
3. **CPU confirmation** (GMP): probable-prime testing, true-root recovery for non-roots, and full chain-length confirmation on the 0.0012% that survive.

See [`docs/SEARCH_PIPELINE.md`](docs/SEARCH_PIPELINE.md) for details.

## Interactive Visualizations

Standalone HTML/JS tools that shaped the search design. Try them live on GitHub Pages:

- [Cunningham Chain Mesh](https://nmicic.github.io/cunningham-chain-search/visualizations/chain-mesh/) — 2-adic square-perimeter map with CC1/CC2 edges
- [2-Adic Tree Explorer](https://nmicic.github.io/cunningham-chain-search/visualizations/2adic-tree/) — inspectable tree with local chain structure
- [3D Fold](https://nmicic.github.io/cunningham-chain-search/visualizations/3d-fold/) — shell structure across levels
- [p+1 Analysis](https://nmicic.github.io/cunningham-chain-search/visualizations/p1-analysis/) — factor coloring and mod-p views
- [Chain Analyzer](https://nmicic.github.io/cunningham-chain-search/visualizations/chain-analyzer/) — single-chain analysis and breaker autopsy
- [Campaign Dashboard](https://nmicic.github.io/cunningham-chain-search/visualizations/campaign-dashboard/) — campaign tracking
- [Immunization Dashboard](https://nmicic.github.io/cunningham-chain-search/visualizations/immunization-dashboard/) — residue immunity analysis

## Build

```bash
# GPU engine (requires CUDA toolkit + GMP)
# Adjust -arch=sm_XX for your GPU / CUDA toolchain
# Example: Ada may use sm_89; some setups may still require sm_86
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
- [CC18 history](https://www.pzktupel.de/CC/HCC18.php)

## Author

Nenad Micic, Belgium — [LinkedIn](https://be.linkedin.com/in/nenadmicic)

## License

See [LICENSE](LICENSE).
