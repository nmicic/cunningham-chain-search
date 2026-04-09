# CPU Search Engines

This folder contains the CPU/GMP search engines.

## Files

- `cc_gmp_v33_03.c`
  Legacy baseline engine from the March 2026 public release. Use this as the
  reference implementation for correctness and comparisons.

- `cc_gmp_v34_bit-vector_07_public.c`
  Experimental CPU engine based on John Armitage's bit-vector L2 filter plus
  additional optimizations. Tuned for high-bit searches (notably 153-bit CC15+)
  where GPU wide-mode buffering can dominate.

- `README_v34_bit-vector.md`
  Focused documentation for the v34 bit-vector engine.

