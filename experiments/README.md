# Experiments

## Failed / Exploratory Approaches

The `failed/` directory contains 17 CUDA-based approaches that were explored and ultimately set aside. Each has a README with context. See `failed/MASTER_INDEX.md` for the full evaluation.

| Experiment | One-line summary |
|-----------|-----------------|
| cc_constructor | Baseline CUDA constructor — superseded by CpC filter approach |
| cc_constructor_algebraic_hybrid | Algebraic pre-screening — overhead exceeded filtering gains |
| cc_constructor_hybrid | CPU/GPU hybrid split — synchronization bottleneck |
| cc_cuda | Early CUDA port — naive parallelism, poor occupancy |
| cc_hybrid_lattice_design | Lattice-based candidate generation — too sparse |
| cc_lattice_search | Lattice enumeration — combinatorial explosion at depth |
| cc_line_sieve | Line sieve on GPU — memory bandwidth limited |
| cc_search | General search framework — too generic, no depth filtering |
| cc_search_256bit | 256-bit variant — arithmetic too slow on GPU |
| cc_search_cuda | Refactored CUDA search — still lacked depth filtering |
| cc_search_merged | Merged CPU+GPU kernel — divergence killed throughput |
| cc_search_unified | Unified memory approach — page faults dominated |
| cc_search_unified_c7 | C7 primorial variant — wrong granularity |
| cunningham_algebraic | Algebraic chain construction — mathematically interesting but slow |
| feature_map | Feature extraction for ML-guided search — insufficient signal |
| prime_rails | Rail-based prime enumeration — novel but not competitive |
| prime_tree_chains | Tree-structured chain search — branching overhead |

The approach that worked (`cc18_filter_cuda_CpC_v13.cu`) uses modular depth filtering: test all chain positions against small primes before any expensive primality test. This achieves 99.9988% rejection on GPU.
