# CC-Clusters: Triplets, Quadruplets, Quintuplets & Rainbows

When Cunningham chain roots of different lengths cluster together in a tiny region of number space, we call them **CC-clusters**. These are beautiful structures — three, four, or five independent chains starting almost at the same number.

Dataset: 929,574 CC10+ roots from the NEW4 campaign (59–153 bits).

---

## Record CC-Triplets

### Overall Champion — "The Ghost CC16 Triplet"

Three CC10 roots packed into a span of just **2^52.65** in 88-bit space:

```
CC10  0x9627b724253d4c1f827cc5
CC10  0x9627b72437bf339da0039f
CC10  0x9627b7243e60393fed82dd
```

All three share the prefix `0x9627b724` — the top **44 bits are identical**. Three independent length-10 Cunningham chains starting within a window of ~7 trillion.

**Root #1 has an interesting ghost chain** (see Ghost Chains section below). Tested to depth 20, it has **14 out of 20 links prime** — the chain breaks at link 11, but links 12–13 and 15–16 are prime:

```
CC10 ■■■■■■■■■■ □ ■■ □ ■■ □ □ □ □  (■=prime, □=composite, depth 20)
     1-------10 11 1213 14 1516 17-20
```

Link 14 was killed by 1579 × 2287 (small factors). Root #2 has link 11 killed by 40213, root #3 by 1181 × 2707.

### Champion Cross-CC Triplet

Three different chain lengths in a span of **2^53.11**:

```
CC12  0x799112DAC1B4F335A0C091B
CC11  0x799112DAC3ACE1947B4D219
CC10  0x799112DAC3DDB3C88C5D923
```

A CC12, CC11, and CC10 sharing the prefix `0x799112DAC` — a staircase of chain lengths in adjacent numbers. Full 23-hex-digit roots, differing only in the last 14 hex digits.

---

## Top CC-Triplets by Category

### Any Mix — Top 10

| Rank | Span | CC levels | Bit size | Notable |
|:----:|:----:|:---------:|:--------:|---------|
| 1 | 2^52.65 | 10+10+10 | 88b | Champion same-CC |
| 2 | 2^52.73 | 10+10+10 | 68b | |
| 3 | 2^52.88 | 10+10+12 | 68b | Mixed with CC12 |
| 4 | 2^53.11 | 12+11+10 | 91b | Champion cross-CC |
| 5 | 2^53.14 | 10+10+10 | 88b | |
| 6 | 2^53.42 | 10+10+10 | 91b | |
| 7 | 2^53.63 | 10+10+10 | 91b | |
| 8 | 2^53.76 | 10+10+10 | 91b | |
| 9 | 2^53.82 | 10+10+10 | 68b | |
| 10 | 2^53.86 | 10+10+10 | 68b | |

### Cross-CC Triplets (All 3 Different Levels) — Top 10

| Rank | Span | CC levels | Hex prefix | Notable |
|:----:|:----:|:---------:|:----------:|---------|
| 1 | 2^53.11 | 12+11+10 | `0x799112dac` | Staircase: 12→11→10 |
| 2 | 2^54.87 | 12+10+11 | `0x80` | |
| 3 | 2^55.34 | **12+15+13** | `0x1e3d9` | High-value! CC12+CC15+CC13 |
| 4 | 2^55.40 | 11+12+10 | `0x26c4f82a` | |
| 5 | 2^55.41 | 10+11+12 | `0x2a762ff7` | Ascending staircase |
| 6 | 2^55.58 | 12+11+10 | `0x1351` | |
| 7 | 2^56.00 | 10+12+11 | `0x1350` | |
| 8 | 2^56.06 | 12+10+11 | `0x809` | |
| 9 | 2^56.08 | 10+12+11 | `0x1296` | |
| 10 | 2^56.11 | **12+10+13** | `0x1157` | CC13 in the mix |

### "Staircase" Triplets

Some cross-CC triplets form perfect ascending or descending sequences:

- **CC10→CC11→CC12** at `0x2a762ff7`: span 2^55.41 — ascending chain lengths
- **CC12→CC11→CC10** at `0x799112dac`: span 2^53.11 — descending chain lengths
- **CC13→CC12→CC12** at `0xe6`: span 2^54.10 — high-value descending

---

## High-Value Triplets (All Roots CC12+)

The top triplets where every root is CC12 or higher:

| Rank | Span | CC levels | Root region |
|:----:|:----:|:---------:|:-----------:|
| 1 | 2^53.88 | 12+12+**13** | `0x48b0` |
| 2 | 2^53.94 | 12+12+12 | `0x1884c` |
| 3 | 2^54.10 | **13**+12+12 | `0xe6–0xea` |
| 4 | 2^54.18 | 12+12+12 | `0x13d7` |
| 5 | 2^54.47 | **13**+12+12 | `0x377f` |

### CC13+ Triplets (All Roots CC13+)

The rarest and most valuable — every root is CC13 or higher:

| Rank | Span | CC levels | Roots |
|:----:|:----:|:---------:|-------|
| **1** | **2^57.58** | **15+14+13** | `0x2c248e798ab8d24e8d` (CC15), `0x2c24e437661e4c625f` (CC14), `0x2c273ff1191c249261` (CC13) |
| 2 | 2^57.82 | 13+13+13 | `0x120394d4...` cluster |
| 3 | 2^57.96 | 13+13+13 | `0x130e06...` cluster |
| 4 | 2^58.75 | 13+14+13 | `0x7c8–0x7cd` |
| 5 | 2^58.86 | 13+13+13 | `0x17b0–0x17b5` |

**The CC15+CC14+CC13 triplet** at `0x2c24` is remarkable — three high-level chains of different lengths, all within 2^57.58 of each other in 70-bit space. This region also appears in the CC-Twins records.

---

## CC-Quadruplets

### Overall Champion

Four roots in a span of just **2^53.97**:

```
CC10  0x1003b93a4510814fd1
CC10  0x1003da31fbb82064f1
CC10  0x1003f1c82e201994e9
CC12  0x1003f7bb1fa2473b41    ← CC12 sneaks in!
```

Three CC10s and a CC12, all sharing the prefix `0x1003`. The CC12 is only 2^50.57 away from the third CC10.

### Same-CC Champion

Four CC10 roots in **2^54.13**:

```
CC10  0x100488814affbdd7b5
CC10  0x1004b44acbcce976a7
CC10  0x1004c982cde458b943
CC10  0x1004cec02b971f5cd1
```

### High-Value Quadruplets (All CC11+)

| Rank | Span | CC levels | Notable |
|:----:|:----:|:---------:|---------|
| 1 | 2^55.55 | **13**+12+12+12 | CC13 leads |
| 2 | 2^56.00 | 12+**13**+**13**+12 | Two CC13s! |
| 3 | 2^56.47 | 12+12+12+**13** | |
| 4 | 2^56.49 | 12+12+12+12 | All CC12 |
| 5 | 2^56.59 | 12+**13**+12+12 | |

---

## CC-Quintuplets

Five roots in a tight window — the largest confirmed clusters:

| Rank | Span | CC levels | Region |
|:----:|:----:|:---------:|--------|
| 1 | 2^54.87 | 11+10+10+10+12 | `0x1003` |
| 2 | 2^55.40 | 10+10+10+10+10 | `0x1004` |
| 3 | 2^55.49 | 11+11+10+10+10 | `0x1003` |
| 4 | 2^55.52 | 10+10+10+10+10 | `0x101b` |
| 5 | 2^55.72 | 10+10+10+12+10 | `0x1003` |

The **champion quintuplet** spans 2^54.87 and contains a CC11, a CC12, and three CC10s — five independent chains of length 10+ starting in a window of ~32 quadrillion.

The `0x1003–0x1005` region appears in multiple quintuplets — the densest concentration of CC10+ roots in this dataset.

---

## Rainbow Clusters (Most Distinct CC Levels)

The maximum number of distinct CC levels found in a consecutive window:

| Window size | Max distinct | Levels | Span | Region |
|:-----------:|:------------:|:------:|:----:|--------|
| 10 roots | **6** | CC10–CC15 | 2^62.61 | `0x1042–0x10a4` |
| 20 roots | **6** | CC10–CC15 | 2^62.92 | `0x101b–0x1095` |
| 50 roots | **6** | CC10–CC15 | 2^64.09 | `0x6f23–0x8026` |

**Six levels** (CC10 through CC15) is the maximum — CC16 and CC17 are too rare to appear in tight windows. The densest rainbow (window=10) packs all 6 levels into just 2^62.61:

```
Window of 10 consecutive roots near 0x1042–0x10a4:
  CC15: 1 root    ← crown jewel
  CC14: 1 root
  CC13: 1 root
  CC12: 5 roots
  CC11: 1 root
  CC10: 1 root
```

---

## Multi-Record Region at 0x2c24 (70-bit)

This region appears in several record tables across both CC_TWINS.md and this document:

| Record | What's there |
|--------|-------------|
| CC-Twin (cross) | CC14+CC15 at gap 2^54 |
| CC13+ Triplet #1 | CC15+CC14+CC13 at span 2^57.58 |
| Cross-CC Triplet #13 | CC12+CC15+CC14 at span 2^56.43 |

```
0x2c238c26bb18ae89a5  CC12
0x2c248e798ab8d24e8d  CC15
0x2c24e437661e4c625f  CC14
0x2c273ff1191c249261  CC13
```

Four different chain lengths (CC12–CC15) cluster within 2^57.66 of each other. This is the highest-value cluster in the dataset. Consistent with random placement — the gap analysis confirms uniform distribution — but interesting as a descriptive feature.

---

## The "Dense Block" at 0x1003–0x1005 (68-bit)

This region hosts the overall champion quadruplet and quintuplet:

```
0x100382a1a95fde5adf  CC11
0x1003b93a4510814fd1  CC10
0x1003da31fbb82064f1  CC10
0x1003f1c82e201994e9  CC10
0x1003f7bb1fa2473b41  CC12  ← CC12 joins the party
0x100488814affbdd7b5  CC10
0x1004b44acbcce976a7  CC10
0x1004c982cde458b943  CC10
0x1004cec02b971f5cd1  CC10
0x1005310879ce023b85  CC10
```

Ten roots spanning 2^56.88 — containing 7 CC10, 1 CC11, 1 CC12, with the tightest internal gaps at 2^50.

---

## Cluster Hierarchy

| Name | N | Record span | Champion levels |
|------|:-:|:-----------:|:---------------:|
| CC-Twin | 2 | 2^40 | CC10+CC10 |
| CC-Triplet | 3 | 2^52.65 | CC10×3 |
| CC-Quadruplet | 4 | 2^53.97 | CC10×3+CC12 |
| CC-Quintuplet | 5 | 2^54.87 | CC11+CC10×3+CC12 |
| Rainbow-6 | 10 | 2^62.61 | CC10–CC15 |
| Dense Block | 10 | 2^56.88 | CC10×7+CC11+CC12 |

### Record Spans

- Twin: 2^40, Triplet: 2^52.65, Quad: 2^53.97, Quint: 2^54.87

The triplet→quad→quint spans grow slowly (+1.32, +0.90 bits per added root). With 929K roots this is likely close to the tightest possible for this dataset size, but no theoretical bound is claimed.

---

## Ghost Chains: The Hidden Structure

A **ghost chain** is a root whose Cunningham chain has prime links beyond the official break point. The official CC length counts only consecutive primes from the root, but some roots have additional primes further down — separated by just one or two composite "gaps."

Full census: 929,574 roots scanned to depth 20 (see `ghost_chains_top.txt`).

### Ghost Chain Distribution

| Primes/20 | Count | % | Significance |
|:---------:|------:|:---:|-------------|
| **18/20** | **7** | 0.0008% | The "Magnificent Seven" |
| **17/20** | **121** | 0.013% | Elite ghost chains |
| 16/20 | 1,290 | 0.14% | Very strong |
| 15/20 | 8,525 | 0.92% | Strong |
| 14/20 | 38,744 | 4.2% | Above average |
| 13/20 | 121,273 | 13.0% | Average |
| 12/20 | 255,491 | 27.5% | Typical |
| 11/20 | 318,389 | 34.3% | Bulk |
| 10/20 | 181,420 | 19.5% | Minimum (CC10 only) |

Post-snapshot discovery: an 88-bit CC16 root whose first 21 links contain 20 primes, with only link 17 composite and links 18-21 prime again.

### The Magnificent Seven (18/20 Primes)

Only 7 roots in the initial dataset snapshot achieve 18 out of 20 prime links:

```
                                                              Pattern (P=prime .=composite)
#1  CC17  0xFE5F018B29D68951358F47     88b  PPPPPPPPPPPPPPPPP.P.  [L18: 36473]
#2  CC16  0x21807352AB07B8FAA5BF2C9    90b  PPPPPPPPPPPPPPPP.P.P  [L17: 1163, L19: 53]
#3  CC13  0x5EC5E93EFE5CC3B4DE6DBF     87b  PPPPPPPPPPPPP..PPPPP  [L14: 4159]
#4  CC13  0xBD13430898FADB3825B137     88b  PPPPPPPPPPPPP.PPP.PP  [L14: 1193, L18: 3557]
#5  CC12  0x215694578FD6F308168EFD     86b  PPPPPPPPPPPP.PPPPP.P  [L13: 21187, L19: 2129]
#6  CC11  0x42797BD6F7DE776520B2CA9    91b  PPPPPPPPPPP.P.PPPPPP  [L12: 8419]
#7  CC11  0xA9DAF6920434C87630276D     88b  PPPPPPPPPPP.PPP.PPPP
```

**Root #1** is the CC17 — the rarest chain in the dataset. It has 17 consecutive primes, then link 18 is killed by just **36473** (a single small factor), link 19 is prime again. One tiny factor away from CC19.

**Root #2** is a CC16 killed at link 17 by **1163** and link 19 by **53**. Two small-factor assassinations prevent CC20.

**Root #3** is a CC13 where after breaking at link 14 (killed by 4159), links 16–20 are all prime — a 5-prime "second wind."

### Hidden Gems: CC10 Roots That Are Secretly Elite

The best CC10 ghost chains have **17/20 primes** — only 3 composites in 20 links despite being officially "just CC10":

```
#1  CC10  0x11D53A292A58460DDB5F9A9  89b  PPPPPPPPPP..P.PPPPPP  [L11: 5189]
#2  CC10  0x13E8373FAC2155DBA361925  89b  PPPPPPPPPP.PP..PPPPP  [L15: 3821]
#3  CC10  0x140009267466F75A5535919  89b  PPPPPPPPPP.PP.PPPP.P  [L19: 71]
```

Root #1 has a "second wind" of 6 consecutive primes (links 15–20) — a hidden CC6 beyond the break. Link 11 was killed by just **5189**.

### "Second Wind" Champions

Roots with the longest consecutive prime run AFTER their first break:

| 2nd wind | Root | Official | Pattern |
|:--------:|------|:--------:|---------|
| **7** | `0x6D1DCDEB6AFBB0477F58B1` | CC10 | `PPPPPPPPPP..PPPPPPP.` |
| **7** | `0x8001B2248CE829CA637717` | CC10 | `PPPPPPPPPP.PPPPPPP..` |
| **7** | `0xAEBC33B59A175E1F41A7B` | CC10 | `PPPPPPPPPP..PPPPPPP.` |
| 6 | `0x42797BD6F7DE776520B2CA9` | CC11 | `PPPPPPPPPPP.P.PPPPPP` |

Three CC10 roots contain a hidden CC7 after their first break — 7 consecutive primes in a row beyond link 10.

### Ghost Chain Notation

```
P = prime link, . = composite link (all patterns measured to depth 20)

Champion triplet root:   PPPPPPPPPP.PP.PP.PP.  14/20 primes
                         1--------10 11 1213 14 1516 1718 1920
```

### Descriptive Summary

| Root | Official | Primes/20 | Ghost pattern | Killer factors |
|------|:--------:|:---------:|---------------|----------------|
| `0xFE5F018B29D68951358F47` | **CC17** | **18/20** | `PPPPPPPPPPPPPPPPP.P.` | L18: 36473 |
| `0x21807352AB07B8FAA5BF2C9` | **CC16** | **18/20** | `PPPPPPPPPPPPPPPP.P.P` | L17: 1163, L19: 53 |
| `0x5EC5E93EFE5CC3B4DE6DBF` | CC13 | **18/20** | `PPPPPPPPPPPPP..PPPPP` | L14: 4159 |
| `0x9627b724253d4c1f827cc5` | CC10 | 14/20 | `PPPPPPPPPP.PP.PP.PP.` | L14: 1579×2287 |

Ghost chains are descriptive — they show which roots had many prime links beyond the official break. Whether ghost-chain density correlates with anything useful (e.g., neighborhood of higher CC roots) is an open question not yet tested.

---

## Open Questions

Clusters and ghost chains are descriptive features of this dataset. The gap analysis (see CC_GAP_ANALYSIS.md) confirms roots are uniformly distributed on the mandatory grid — there is **no validated search heuristic** based on clusters or ghost chains.

Questions that could be explored post-campaign:

- **Ghost chain census done** (depth 20, see `ghost_chains_top.txt`). Open: does ghost-chain density in a region correlate with the presence of higher CC roots? Not yet tested.
- **CC16 neighborhood analysis**: Do the 50 nearest roots to each CC16 have different CC-level or ghost-chain distributions than a random sample?
- **Residue fingerprinting**: Do the multi-record regions (0x2c24, 0x1003) share any residue patterns mod higher primes, or are they simply random coincidences?

---

## Related

- [CC_TWINS.md](CC_TWINS.md) — CC-twin records (pairs)
- `cc_twins.csv` — CC-twin pairs (98 rows, in `data/`)
- [CC_GAP_ANALYSIS.md](CC_GAP_ANALYSIS.md) — Gap distribution and uniformity proof
