# cc_gmp_v34_bit-vector_07: CPU Cunningham Chain Search (Bit-Vector L2)

Experimental CPU search engine for long first-kind Cunningham chains.
Built on `cc_gmp_v33_03.c` with John Armitage's bit-vector L2 filter and a
small set of additional optimizations. Intended use is high-bit CPU searches
(notably 153-bit CC15+) where GPU wide-mode buffering makes GPU runs slower,
so the CPU path can dominate.

## Build

```bash
gcc -O3 -march=native -flto -o cc_v34_bit-vector_07 \
    cc_gmp_v34_bit-vector_07_public.c -lgmp -lpthread -lm
```

Requires: GCC with `__builtin_ctzll`, GMP, pthreads.

## Quick Start

```bash
# Run self-tests
./cc_v34_bit-vector_07 --test

# 153-bit CC15 search (intended use case)
./cc_v34_bit-vector_07 --prefix 0b1 --sequential \
    --bits 153 --target 15 --log 15 --threads 16

# Sieve-only throughput (no proving)
./cc_v34_bit-vector_07 --prefix 0b1 --sequential \
    --bits 89 --target 13 --threads 4 --log 50000 --sieve-only

# Baseline comparison (disable bit-vector filter)
./cc_v34_bit-vector_07 --no-bitvec --prefix 0b1 --sequential \
    --bits 89 --target 13 --threads 4 --log 50000
```

## Key Flags

| Flag | Purpose |
|------|---------|
| `--bits N` | Bit-size of candidate roots |
| `--target N` | Minimum chain length to report |
| `--log N` | Minimum chain length to log (set = target for reverse-prove gain) |
| `--threads N` | Worker threads |
| `--prefix 0bXXX` | Binary prefix for sequential search |
| `--sequential` | Sequential prefix scan (vs random chunks) |
| `--no-bitvec` | Disable bit-vector L2 filter (fallback to v33 behavior) |
| `--no-ext2` | Disable OPT-G super-ext-L2 filter |
| `--ext2` | Force-enable OPT-G (overrides auto-tune at low bits) |
| `--sieve-only` | Run sieve filters only (no proving), for throughput checks |
| `--line-depth N` | Line-sieve depth (default 17; use 14 for CC15) |
| `--screen-primes N` | Number of chain screening primes (default 132) |

## What’s New vs v33_03

Short summary of the effective changes:

- Bit-vector L2 filter (John Armitage) for 64 candidates at once.
- OPT-A: incremental `tile_M_line` update.
- OPT-B: `ctzll` survivor iteration.
- OPT-E: reverse-prove at top chain position.
- OPT-I: birth-certificate root detection when `n == 1 (mod q)`.
- OPT-G: super-ext-L2 bitmask for 101–127.
- OPT-C: packed L1 line-sieve kill table.

The primary benchmark summary is included in the header of
`cc_gmp_v34_bit-vector_07_public.c`.

## Correctness

Compatibility checks against `cc_gmp_v33_03.c` are included in the source
header and in the built-in `--test` suite.

## Credits

- Nenad Micic — base engine `cc_gmp_v33_03.c` and infrastructure
- John Armitage — bit-vector L2 filter algorithm (Oxford, 2021)

Apache 2.0. See LICENSE in the main repo.
