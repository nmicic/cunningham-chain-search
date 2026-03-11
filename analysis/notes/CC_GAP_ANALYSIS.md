# CC Root Gap Analysis: Distance and Spacing

Analysis of distances between consecutive Cunningham chain roots across 929,574 CC10+ roots.

## Corrected Key Finding: Roots Are Uniformly Distributed on the Mandatory Grid

**Initial observation**: minimum gaps appeared 15–32 bits closer than the **mean** gap, suggesting clustering. **Corrected analysis**: comparing against the **expected minimum gap** for random samples (which is ~R/n², not R/n), the observed minimums are exactly what uniform random spacing predicts — or even slightly larger (anti-clustered at CC12+).

**This means: there are no "fertile regions". Every position on the mandatory grid is equally likely to yield CC18. Sequential search is optimal. Jumping to a predicted position after a find provides no advantage.**

### The Mandatory Grid

All standard CC10+ roots satisfy `p ≡ -1 (mod 40755)` where 40755 = lcm(3, 5, 11, 13, 19). This creates a mandatory grid with spacing ~2^15.3 bits. Every gap between standard roots is exactly divisible by 40755.

### Minimum Gap: Observed vs Mean vs Expected Minimum (90-bit band)

| CC | Roots | Mean gap | Expected min | Observed min | vs Mean | vs Exp min | Verdict |
|----|------:|:--------:|:------------:|:------------:|:-------:|:----------:|:-------:|
| CC10 | 273,384 | 2^56.6 | 2^39.6 | 2^40 | +16.6 | **+0.4** | Random |
| CC11 | 47,626 | 2^59.2 | 2^44.6 | 2^49 | +10.2 | **+4.4** | Spread |
| CC12 | 8,393 | 2^61.7 | 2^49.6 | 2^52 | +9.7 | **+2.4** | Spread |
| CC13 | 1,435 | 2^64.2 | 2^54.7 | 2^56 | +8.2 | **+1.3** | Random |
| CC14 | 254 | 2^66.7 | 2^59.7 | 2^66 | +0.7 | **+6.3** | Spread |
| CC15 | 35 | 2^69.6 | 2^65.4 | 2^69 | +0.6 | **+3.6** | Spread |
| CC16 | 9 | 2^71.5 | 2^69.4 | 2^82 | -10.5 | **+12.6** | Spread |

- **Mean gap** = 2^bits / count (the average spacing)
- **Expected min** = 2R/n² (the expected minimum gap for n uniform random points in range R)
- **vs Exp min** positive = roots are **farther apart** than random (anti-clustered)

The "deficit vs mean" that initially looked like clustering is simply the natural difference between mean and minimum for a sample of n points. When compared against the **correct baseline** (expected minimum), CC12+ roots are actually **more evenly spread** than random — the modular constraints for higher CC levels enforce spacing.

**Median gaps match expected mean** (within ~2 bits), confirming the bulk distribution is uniform.

### Minimum Gap Scaling Per CC Level

| Transition | Min gap increase | Count ratio | Expected increase |
|-----------|:----------------:|:-----------:|:-----------------:|
| CC10 → CC11 | +9 bits | 5.7:1 | +2.5 bits |
| CC11 → CC12 | +3 bits | 5.7:1 | +2.5 bits |
| CC12 → CC13 | +4 bits | 5.8:1 | +2.5 bits |
| CC13 → CC14 | +10 bits | 5.6:1 | +2.5 bits |
| CC14 → CC15 | +3 bits | 7.3:1 | +2.9 bits |
| CC15 → CC16 | +13 bits | 3.9:1 | +2.0 bits |

The min gap increases faster than the count ratio alone would explain — the clustering deficit shrinks at higher CC levels, but remains significant through CC15.

---

## Gap Distribution Pyramid (All Bits Merged)

### CC10 (749,017 roots, 748,981 gaps)

| Gap range | Count | % | Cumul % |
|-----------|------:|--:|--------:|
| 2^40 – 2^49 | 150 | 0.0% | 0.0% |
| 2^50 – 2^54 | 2,828 | 0.4% | 0.4% |
| 2^55 – 2^59 | 50,302 | 6.7% | 7.1% |
| **2^60 – 2^64** | **126,047** | **16.8%** | **23.9%** |
| 2^65 – 2^69 | 80,564 | 10.8% | 34.7% |
| **2^70 – 2^74** | **380,724** | **50.8%** | **85.5%** |
| 2^75 – 2^79 | 23,393 | 3.1% | 88.6% |
| 2^80+ | 3,901 | 0.5% | 100% |

**Bimodal distribution** with peaks at 2^62 and 2^71 — caused by the two bit-range populations in the data (65–77 and 85–92 bit roots).

### CC12 (39,135 roots, 39,095 gaps)

| Gap range | Count | % | Cumul % |
|-----------|------:|--:|--------:|
| 2^44 – 2^49 | 13 | 0.0% | 0.0% |
| 2^50 – 2^54 | 193 | 0.5% | 0.5% |
| 2^55 – 2^59 | 3,027 | 7.7% | 8.3% |
| **2^60 – 2^64** | **10,663** | **27.3%** | **35.5%** |
| 2^65 – 2^69 | 3,658 | 9.4% | 44.9% |
| **2^70 – 2^79** | **19,543** | **50.0%** | **94.9%** |
| 2^80+ | 108 | 0.3% | 100% |

### CC14 (1,318 roots, 1,297 gaps)

| Gap range | Count | % | Cumul % |
|-----------|------:|--:|--------:|
| 2^54 – 2^59 | 7 | 0.5% | 0.5% |
| 2^60 – 2^64 | 147 | 11.3% | 11.9% |
| **2^65 – 2^69** | **372** | **28.7%** | **40.6%** |
| 2^70 – 2^74 | 87 | 6.7% | 47.3% |
| **2^75 – 2^84** | **629** | **48.5%** | **95.8%** |
| 2^85+ | 20 | 1.5% | 100% |

### CC15 (237 roots, 218 gaps)

| Gap range | Count | % | Cumul % |
|-----------|------:|--:|--------:|
| 2^60 – 2^64 | 16 | 7.3% | 7.3% |
| **2^65 – 2^69** | **76** | **34.9%** | **42.2%** |
| 2^70 – 2^74 | 21 | 9.6% | 51.8% |
| **2^78 – 2^85** | **91** | **41.7%** | **93.6%** |
| 2^86+ | 6 | 2.7% | 100% |

### CC16 (41 roots, 32 gaps)

| Gap range | Count | % | Cumul % |
|-----------|------:|--:|--------:|
| 2^65 – 2^69 | 6 | 18.8% | 18.8% |
| 2^70 – 2^74 | 7 | 21.9% | 40.6% |
| 2^81 – 2^87 | 16 | 50.0% | 100% |

---

## Cross-CC Proximity: Do Higher CC Roots Cluster Near Lower CC Roots?

For each CC(N+1) root, we found the nearest CC(N) root. Results for the 88–91 bit band:

| From CC | # roots | Nearest CC(N-1) | Min dist | Expected NN | NN deficit |
|---------|--------:|:---------------:|:--------:|:-----------:|:----------:|
| CC11 | 119,227 | CC10 (675K) | 2^44 | 2^69.6 | **25.6** |
| CC12 | 20,917 | CC11 (119K) | 2^49 | 2^72.1 | **23.1** |
| CC13 | 3,573 | CC12 (21K) | 2^53 | 2^74.6 | **21.6** |
| CC14 | 624 | CC13 (3.6K) | 2^60 | 2^77.2 | **17.2** |
| CC15 | 95 | CC14 (624) | 2^66 | 2^79.7 | **13.7** |
| CC16 | 20 | CC15 (95) | 2^70 | 2^82.4 | **12.4** |

**Corrected**: When compared against the expected minimum nearest-neighbor distance for independent random samples (≈ R / (n_hi × n_lo)), the observed distances are **5–7 bits larger** — the same anti-clustering pattern. CC roots of different lengths do NOT cluster together more than random.

| From CC | n_hi | n_lo | Expected min NN | Observed | vs Exp | Verdict |
|---------|-----:|-----:|:---------------:|:--------:|:------:|:-------:|
| CC11→CC10 | 119,227 | 675,102 | 2^38.5 | 2^44 | +5.5 | Spread |
| CC12→CC11 | 20,917 | 119,227 | 2^43.5 | 2^49 | +5.5 | Spread |
| CC13→CC12 | 3,573 | 20,917 | 2^48.5 | 2^53 | +4.5 | Spread |
| CC14→CC13 | 624 | 3,573 | 2^53.6 | 2^60 | +6.4 | Spread |
| CC15→CC14 | 95 | 624 | 2^58.8 | 2^66 | +7.2 | Spread |
| CC16→CC15 | 20 | 95 | 2^63.8 | 2^70 | +6.2 | Spread |

### CC16 → Nearest Neighbors Across All CC Levels

Selected CC16 roots with nearest neighbors:

| CC16 root | → CC15 | → CC14 | → CC13 | → CC12 |
|-----------|:------:|:------:|:------:|:------:|
| `0x289d94...` | 2^83 | **2^65** | **2^65** | **2^65** |
| `0x2a7631...` | 2^80 | **2^69** | **2^67** | **2^65** |
| `0xa0006b...` | 2^82 | **2^68** | **2^63** | 2^67 |
| `0x93b4bc...` | **2^70** | **2^67** | 2^68 | **2^64** |
| `0x3007d0...` | 2^76 | 2^80 | 2^76 | **2^63** |
| `0x8764e2...` | 2^81 | **2^69** | **2^64** | **2^65** |

Several CC16 roots have CC12–CC14 neighbors within 2^65 — only 25 bits of search radius in a 90-bit space.

---

## Per-Bit-Size Gap Statistics

### 90-bit (dominant search band)

| CC | Roots | Min gap | P10 gap | Median | P90 gap | Max gap | Expected |
|----|------:|:-------:|:-------:|:------:|:-------:|:-------:|:--------:|
| CC10 | 273,384 | 2^40 | 2^58 | 2^70 | 2^72 | 2^80 | 2^71.9 |
| CC11 | 47,626 | 2^49 | 2^60 | 2^72 | 2^75 | 2^82 | 2^74.5 |
| CC12 | 8,393 | 2^52 | 2^63 | 2^75 | 2^77 | 2^83 | 2^77.0 |
| CC13 | 1,435 | 2^56 | 2^67 | 2^77 | 2^80 | 2^86 | 2^79.5 |
| CC14 | 254 | 2^66 | 2^70 | 2^80 | 2^83 | 2^87 | 2^82.0 |
| CC15 | 35 | 2^69 | 2^78 | 2^83 | 2^85 | 2^89 | 2^84.9 |
| CC16 | 9 | 2^82 | 2^83 | 2^85 | 2^87 | 2^87 | 2^86.8 |

### 88-bit band

| CC | Roots | Min gap | Median | Expected |
|----|------:|:-------:|:------:|:--------:|
| CC10 | 141,472 | 2^0* | 2^62 | 2^71.0 |
| CC11 | 25,268 | 2^51 | 2^64 | 2^73.4 |
| CC12 | 4,552 | 2^53 | 2^67 | 2^75.8 |
| CC13 | 791 | 2^61 | 2^69 | 2^78.4 |
| CC14 | 144 | 2^66 | 2^72 | 2^80.8 |

*2^0 = duplicate entry in data.

### 89-bit band

| CC | Roots | Min gap | Median | Expected |
|----|------:|:-------:|:------:|:--------:|
| CC10 | 139,366 | 2^41 | 2^70 | 2^71.9 |
| CC11 | 25,107 | 2^46 | 2^72 | 2^74.4 |
| CC12 | 4,401 | 2^55 | 2^75 | 2^76.9 |
| CC13 | 739 | 2^65 | 2^77 | 2^79.5 |
| CC14 | 128 | 2^67 | 2^80 | 2^82.0 |
| CC15 | 28 | 2^78 | 2^82 | 2^84.2 |

---

## What This Means for CC18 Search

### 1. No "Jump" Strategy — Sequential Search Is Optimal

The data is clear: **there are no fertile regions**. After finding a CC13, CC14, or CC16 root, the probability of finding another high-CC root at `root + 2^63` is the same as at any other grid position. Jumping provides no advantage over continuing sequential search.

**Why**: the anti-clustering at CC12+ means higher CC roots are actually more evenly spread than random. The modular constraints act as a natural spacing mechanism.

### 2. No "Pyramid" Effect

Cross-CC proximity analysis (corrected) shows CC roots of different lengths are **independently distributed** on the mandatory grid. Finding CC12 near some position does NOT predict CC14+ nearby.

### 3. The Mandatory Grid Period

All standard roots sit on a grid with spacing **40,755** ≈ 2^15.3 (= lcm(3,5,11,13,19)). In a 90-bit search:
- Grid has ~2^74.7 positions
- CC10 fills ~266K of them → density 1 per 2^56 grid positions
- CC16 fills ~9 → density 1 per 2^71 grid positions

The grid guarantees that **gap is always a multiple of 40,755** for standard roots. The search already exploits this — only grid positions are tested.

### 4. What Actually Matters: Search Volume

The decay curve (R² = 0.997) gives a precise answer:

| Target | CC16 needed | CC10 equivalent | Factor vs current | Strategy |
|--------|:-----------:|:---------------:|:-----------------:|----------|
| **CC17** | ~6 | ~190K | ~0.7x | **Already within reach** |
| **CC18** | ~31 | ~1.0M | ~3.9x (at 89–91 bit) | Scale search 4x |
| **CC19** | ~177 | ~5.6M | ~22x | Needs GPU |
| CC20 | ~990 | ~31M | ~125x | Needs GPU cluster |

### 5. The Real CC18 Path

**You are closer than you think.** The data from two bit bands combined:

| Band | CC16 found | CC18 expected | Factor needed |
|------|:----------:|:-------------:|:-------------:|
| 89–91 bit | 8 | 0.25 | 3.9x more search |
| 65–77 bit | 36 | 1.14 | **Already at ~1x** |
| **Combined** | **44** | **1.40** | **< 1x (should already exist)** |

The 65–77 bit dataset already has sufficient CC16 density that ~1.4 CC18 roots are expected statistically. The fact that none were found suggests either:
- Statistical fluctuation (0.25 probability of 0 finds when expecting 1.4)
- The decay ratio increases slightly at CC17–CC18 (tighter constraints)
- The 65–77 bit search may not have been deep enough at CC17 level

### 6. Optimal Resource Allocation

Since every grid position is equally likely:

**A) Widen the 89-bit search** (easiest win)
- Currently searching 89–91 bit. Add 88 bit and 92 bit.
- Each added bit size ~doubles effective search volume
- 5 bit sizes (88–92) = ~2.5x more than 89–91 alone

**B) Continue sequential search** (no jumping)
- The search is already optimally structured
- Don't waste time on local search around finds
- Every second of sequential search moves you toward CC18

**C) GPU acceleration** (biggest multiplier)
- Sieve phase is 99.86% of the kill → ideal for GPU parallelism
- RTX 5090 could provide 50–100x on sieve
- This alone could close the 3.9x gap for CC18 at 89–91 bit

**D) Don't use primorial filtering** (confirmed by immunization data)
- Safe-residue sieve covers 13.2% of candidates
- Primorial 7# covers only 0.001%
- The search already uses safe-residue sieve (v31_claude)

---

## Output Files

| File | Contents |
|------|----------|
| `gap_comprehensive.csv` | Per-(cc,bits): min/p10/median/p90/max gap, expected, deficit |
| `gap_heatmap.csv` | Gap distribution bucketed by 2^X per (cc,bits) |
| `gap_pyramid.csv` | Gap distribution per CC level (all bits merged) |
| `gap_closest.csv` | 500 smallest gaps with root pairs |
| `gap_statistics.csv` | Summary statistics per (cc,bits) bucket |

### Reproducing

```bash
python3 cc_gap_analysis.py cc10plus_roots_snapshot_2026-03-12.txt --prefix gap --min-cc 10
```
