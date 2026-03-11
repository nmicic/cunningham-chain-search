# Repository Map

```
CC18/
├── README.md                  Top-level overview
├── LICENSE
│
├── src/
│   ├── cuda/                  GPU filter engine (CUDA)
│   │   ├── cc18_filter_cuda_CpC_v13.cu    Production GPU sieve + prover
│   │   └── cc18_filter_cuda_CpC_v13b.cu   Variant
│   └── cpu/                   CPU search engine (C + GMP)
│       └── cc_gmp_v33_03.c   Production CPU sieve + prover (wide-mode, >127-bit)
│
├── gp/                        GP/PARI library and MCP server
│   ├── cc_lib_v10.gp          Library (37 self-tests, Layers 1-10)
│   ├── HOWTO_cc_lib_v10.md    Usage guide
│   ├── HOWTO_cc_lib_v10_SESSIONS.md  Worked examples
│   └── mcp/                   Model Context Protocol server (Node.js)
│       ├── cc_math_server.js  Server implementation
│       ├── SETUP.md           Setup instructions
│       └── ...                Docker config, knowledge base
│
├── analysis/
│   ├── notes/                 Post-campaign analysis documents
│   │   ├── CC_GAP_ANALYSIS.md          Root spacing and uniformity
│   │   ├── CC_IMMUNIZATION_ANALYSIS.md   Immunization analysis overview
│   │   ├── CC_IMMUNIZATION_STATISTICS.md  Residue immunity patterns
│   │   ├── CC_IMMUNE_COMBINATIONS.md   Immune fingerprint catalog
│   │   ├── CC_TWINS.md                 Closest CC root pairs
│   │   ├── CC_CLUSTERS.md             Triplets, quadruplets, ghost chains
│   │   └── ghost_chains_top.txt       Ghost chain census (depth 20)
│   └── scripts/               Analysis tools (Python + C)
│       ├── cc_gap_analysis.py          Gap distribution analysis
│       ├── analyze_cc_immunization.py  Immunization analysis
│       ├── cc_immunization_analysis.c  Immunization analysis → JSON (C)
│       ├── cc_immunization_heatmap.c   Immunization heatmap CSV (C)
│       ├── cc_chain_autopsy_combo.py   Chain breaker autopsy
│       ├── cc_fingerprint.py           Immune fingerprint extraction
│       ├── ghost_chains.py             Ghost chain scanner
│       └── cc_dashboard_data.py        Dashboard data generator
│
├── visualizations/            Interactive HTML tools (GitHub Pages ready)
│   ├── chain-analyzer/        Single-chain analysis, breaker autopsy
│   ├── 2adic-tree/            2-adic tree explorer + generator script
│   ├── 3d-fold/               3D shell structure visualization
│   ├── p1-analysis/           p+1 factor coloring + generator script
│   ├── chain-mesh/            2-adic square-perimeter CC mesh
│   ├── campaign-dashboard/    Campaign tracking dashboard
│   └── immunization-dashboard/  Residue immunity dashboard
│
├── docs/
│   ├── REPO_MAP.md            This file
│   └── SEARCH_PIPELINE.md     How the GPU+CPU pipeline works
│
├── tools/                     Utilities
│   ├── run_cc18_campaign.sh   Campaign runner script (CPU)
│   ├── run_cc18_random_bits.sh  GPU runner, random bit sizes
│   ├── run_cc17b_sequential.sh  GPU runner, second-kind sequential
│   ├── validate_chains.c      Chain validator (C)
│   ├── validate_chains.py     Chain validator (Python)
│   └── validate_chains.gp     Chain validator (GP/PARI)
│
├── data/                      Campaign dataset
│   ├── README.md              What's included and excluded
│   ├── sample/                Small sample (1000 roots)
│   ├── *.csv                  Gap statistics, twins, heatmaps
│   └── (full dataset via GitHub Releases)
│
└── experiments/
    └── failed/                17 failed/exploratory CUDA approaches
        ├── MASTER_INDEX.md    Index with one-line summaries
        └── <name>/README.md   Per-experiment documentation
```
