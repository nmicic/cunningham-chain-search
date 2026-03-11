# Search Pipeline

## Overview

Two-phase pipeline: GPU rejects candidates cheaply, CPU proves the survivors.

```
Candidate space (~2^74 grid positions per bit size)
    │
    ▼
┌─────────────────────────────┐
│  GPU Sieve (CUDA)           │  57-65 B candidates/sec
│  Modular residue filtering  │  Rejects 99.9988%
│  across d chain positions   │
└─────────────┬───────────────┘
              │ 0.0012% survive
              ▼
┌─────────────────────────────┐
│  CPU Prover (GMP)           │  BPSW primality proving
│  Full chain verification    │  on each survivor
└─────────────┬───────────────┘
              │
              ▼
         CC10+ roots
```

## GPU Phase: Modular Depth Filtering

For a target of CC18, the GPU tests each candidate root against small primes across all 18 chain positions. If any chain position is divisible by a small prime, the candidate is rejected without expensive primality testing.

The key insight: a candidate that fails at chain position 11 due to a small factor can never be CC18, regardless of whether it's prime at positions 1-10. This depth filtering is what achieves the 99.9988% rejection rate.

**Implementation:** `src/cuda/cc18_filter_cuda_CpC_v13.cu`
- Prefix-sharded search: binary prefix selects the search region
- Sequential sweep within each prefix
- Two-kernel architecture: filter kernel + prove kernel
- Tested on RTX 4090 and RTX 5090

## CPU Phase: Primality Proving

Survivors from the GPU sieve undergo full BPSW primality testing using GMP. The chain is extended position by position until a composite is found.

**Implementation:** `src/cpu/cc_gmp_v33_03.c`
- GMP-based BPSW proving
- Thread-pinned for V-Cache locality (AMD 7950X3D)
- Prefix sharding matches GPU regions
- Supports arbitrary bit widths (tested 59-153 bit)

## The Mandatory Grid

All standard CC10+ first-kind roots satisfy `p ≡ -1 (mod 40755)` where 40755 = lcm(3, 5, 11, 13, 19). This creates a mandatory grid with spacing ~2^15.3. Only grid positions need testing, which the search exploits.

## Note on Published Code

The published search code reflects the campaign version used to generate the released snapshot and analysis outputs. Some later or alternate optimizations were explored separately. One example is a reverse-check proving path (see `cc18_filter_cuda_CpC_v14.cu`), which can reduce CPU prover work by checking from higher chain positions inward, leveraging lower prime density in the outer shells. This path was not used for the published snapshot because the forward path was better suited for logging, inspection, and post-campaign analysis.

## Performance

| Metric | Value |
|--------|-------|
| GPU throughput | 57-65 B candidates/sec |
| GPU survivor rate | 0.0012% |
| CPU wheel rate | ~2.65 B/sec (30 threads, 89-bit) |
| Primary search band | 89-91 bit |
| Dataset produced | 1.77M roots, 929K CC10+ |
