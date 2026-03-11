# CC-Twins: Closest Cunningham Chain Pairs

A **CC-twin** is a pair of Cunningham chain roots that are unusually close in value. Two flavors:

- **Same-CC twin**: two roots of the same chain length (e.g., two CC14 roots near each other)
- **Cross-CC twin**: roots of different lengths near each other (e.g., a CC12 next to a CC15)

Dataset: 929,574 CC10+ roots from the NEW4 campaign (59–153 bits).

---

## Record CC-Twins (Same-CC)

The tightest same-CC pair at each chain length:

| CC | Gap | Bits | Root 1 | Root 2 |
|----|:---:|:----:|--------|--------|
| **CC10** | **2^40** | 90b | `0x22f0b53edd7f07eb0953771` | `0x22f0b53edd7f2370ca67953` |
| **CC11** | **2^46** | 89b | `0x11851fa21d07278b912782d` | `0x11851fa21d0bcda3d675a6b` |
| **CC12** | **2^44** | 73b | `0x1486cb8f4411cacae3f` | `0x1486cb90bc68166ec7f` |
| **CC13** | **2^50** | 72b | `0x83e20e6fed56062713` | `0x83e216231498a70cd1` |
| **CC14** | **2^54** | 75b | `0x62a278ff41e70d1f725` | `0x62a27fbe19809f0f889` |
| **CC15** | **2^60** | 66b | `0x2c45205cf640715fb` | `0x2dc0bdebf384b0685` |
| **CC16** | **2^65** | 72b | `0xe2b9f8a983549c7d1b` | `0xe50be879c713d72a65` |

### What the Gap Means

The CC10 twin pair shares **12 hex digits** of prefix — they differ only in the last 10 hex digits. In a 90-bit number, the top 48 bits are identical. The gap of 2^40 means these two independent CC10 chains start only ~1.1 trillion apart.

The CC16 twin at 72-bit: two CC16 roots with a gap of only 2^65 in a 72-bit space. These are just 7 bits apart relative to the number size.

### Same-CC Twins — Top 3 Per Level

**CC10** (749,024 roots)

| Rank | Gap | Bits | Shared prefix | Gap hex |
|:----:|:---:|:----:|:-------------:|---------|
| 1 | 2^40 | 90b | 12 hex digits | `0x1b85c1141e2` |
| 2 | 2^41 | 89b | 12 hex digits | `0x2ff46d9d2ec` |
| 3 | 2^43 | 91b | 10 hex digits | `0xb5b40402008` |

**CC11** (132,580 roots)

| Rank | Gap | Bits | Shared prefix |
|:----:|:---:|:----:|:-------------:|
| 1 | 2^46 | 89b | 11 hex digits |
| 2 | 2^49 | 90b | 10 hex digits |
| 3 | 2^50 | 90b | 10 hex digits |

**CC12** (39,138 roots)

| Rank | Gap | Bits | Shared prefix |
|:----:|:---:|:----:|:-------------:|
| 1 | 2^44 | 73b | 6 hex digits |
| 2 | 2^46 | 74b | 7 hex digits |
| 3 | 2^46 | 73b | 6 hex digits |

**CC13** (7,226 roots)

| Rank | Gap | Bits | Shared prefix |
|:----:|:---:|:----:|:-------------:|
| 1 | 2^50 | 72b | 4 hex digits |
| 2 | 2^52 | 70b | 4 hex digits |
| 3 | 2^53 | 65b | 3 hex digits |

**CC14** (1,323 roots)

| Rank | Gap | Bits | Shared prefix |
|:----:|:---:|:----:|:-------------:|
| 1 | 2^54 | 75b | 5 hex digits |
| 2 | 2^57 | 64b | 1 hex digit |
| 3 | 2^57 | 73b | 3 hex digits |

**CC15** (238 roots)

| Rank | Gap | Bits | Shared prefix |
|:----:|:---:|:----:|:-------------:|
| 1 | 2^60 | 66b | 1 hex digit |
| 2 | 2^61 | 76b | 3 hex digits |
| 3 | 2^62 | 69b | 1 hex digit |

**CC16** (44 roots)

| Rank | Gap | Bits | Shared prefix |
|:----:|:---:|:----:|:-------------:|
| 1 | 2^65 | 72b | 1 hex digit |
| 2 | 2^65 | 73b | 2 hex digits |
| 3 | 2^66 | 73b | 1 hex digit |

---

## Record CC-Twins (Cross-CC)

The closest pair between different CC levels:

| Pair | Gap | Bits | Lower root | Higher root |
|------|:---:|:----:|------------|-------------|
| **CC10 + CC11** | **2^44** | 90b | `0x326bd358a6f48d989771745` (CC10) | `0x326bd358a6f36cf7153f87f` (CC11) |
| **CC10 + CC12** | **2^47** | 90b | `0x2042d722bc15714311989c3` (CC10) | `0x2042d722bc1f6651abe5497` (CC12) |
| **CC10 + CC13** | **2^46** | 90b | `0x29d8a26f0a3db164b9eebf7` (CC10) | `0x29d8a26f0a4373f466a014b` (CC13) |
| **CC11 + CC12** | **2^49** | 91b | `0x6bc0f9b5633016963f36bad` (CC11) | `0x6bc0f9b5630c8171c0b59f5` (CC12) |
| **CC11 + CC13** | **2^52** | 88b | `0x9d8b691e136a5dbff0b2b5` (CC11) | `0x9d8b691e0351899eb9ee29` (CC13) |
| **CC11 + CC14** | **2^54** | 88b | `0xa4ed2fa97129f610aceca3` (CC11) | `0xa4ed2fa9db0c1d569074c9` (CC14) |
| **CC12 + CC13** | **2^45** | 62b | `0x3e776b256e09e1bb` (CC12) | `0x3e7732f4e08941a3` (CC13) |
| **CC12 + CC14** | **2^50** | 68b | `0x88732cac87eb7dc2f` (CC12) | `0x8872bec362b82ee39` (CC14) |
| **CC12 + CC15** | **2^49** | 68b | `0xe1aef0e8eda975357` (CC12) | `0xe1af247aaf841e8fb` (CC15) |
| **CC13 + CC14** | **2^54** | 69b | `0x11932525db8fb7c9f7` (CC13) | `0x1192db45ae2bcfd8b7` (CC14) |
| **CC13 + CC15** | **2^50** | 73b | `0x1e3d99df12741e11c1d` (CC13) | `0x1e3d9993b12c06c0ec5` (CC15) |
| **CC13 + CC16** | **2^56** | 74b | `0x3256094058a8221e019` (CC13) | `0x3256258e3e9ef6d0b8d` (CC16) |
| **CC14 + CC15** | **2^54** | 70b | `0x2c24e437661e4c625f` (CC14) | `0x2c248e798ab8d24e8d` (CC15) |
| **CC14 + CC16** | **2^59** | 72b | `0xb3d106dba2693a091d` (CC14) | `0xb3da5a03eda307f5ef` (CC16) |
| **CC15 + CC16** | **2^60** | 72b | `0xe2cd9bb8317bfd2fc5` (CC15) | `0xe2b9f8a983549c7d1b` (CC16) |
| **CC14 + CC17** | **2^75** | 88b | `0xfe4f342126cd30076ef96d` (CC14) | `0xfe5f018b29d68951358f47` (CC17) |
| **CC16 + CC17** | **2^85** | 88b | `0xdd253a1b54851aba0b6b81` (CC16) | `0xfe5f018b29d68951358f47` (CC17) |

### Highlights

**CC12 + CC15 twin at 2^49**: A CC12 root and a CC15 root just 2^49 apart in 68-bit space. Only 19 bits of gap relative to number size.

**CC13 + CC15 twin at 2^50**: A CC13 and CC15 sharing the same 73-bit neighborhood, gap only 2^50.

**CC13 + CC16 twin at 2^56**: A CC13 next to a CC16 in 74-bit space. Gap of 2^56 = only 18 bits relative to number size.

**CC10 + CC13 twin at 2^46**: In 90-bit space, a CC10 and CC13 just 2^46 apart — 44 bits of shared prefix.

---

## Notable Cross-CC Pairs

Some of the tightest cross-CC twins also appear near other high-CC roots. For triplets and larger clusters, see [CC_CLUSTERS.md](CC_CLUSTERS.md).

| Pair | Gap | Bit size |
|------|:---:|:--------:|
| CC10+CC11 at `0x326bd358` | 2^44 | 90b |
| CC10+CC12 at `0x2042d722` | 2^47 | 90b |
| CC12+CC13 at `0x3e77` | 2^45 | 62b |

---

## Gap Progression: How Tight Can CC-Twins Get?

| CC twin | Record gap | Bit size | Gap / bit size ratio |
|---------|:----------:|:--------:|:--------------------:|
| CC10+CC10 | 2^40 | 90 | **44%** of bits shared |
| CC11+CC11 | 2^46 | 89 | **48%** of bits shared |
| CC12+CC12 | 2^44 | 73 | **40%** of bits shared |
| CC13+CC13 | 2^50 | 72 | **31%** of bits shared |
| CC14+CC14 | 2^54 | 75 | **28%** of bits shared |
| CC15+CC15 | 2^60 | 66 | **9%** of bits shared |
| CC16+CC16 | 2^65 | 72 | **10%** of bits shared |
| CC10+CC11 | 2^44 | 90 | **51%** of bits shared |
| CC12+CC15 | 2^49 | 68 | **28%** of bits shared |
| CC13+CC16 | 2^56 | 74 | **24%** of bits shared |

"% of bits shared" = (bit_size - gap_log2) / bit_size. Higher = tighter twin.

---

## Data Files

| File | Contents |
|------|----------|
| `cc_twins.csv` | All top-3 same-CC and cross-CC twins with full hex roots |
| `gap_closest.csv` | 500 tightest same-CC pairs (all levels combined) |

### Re-running

```bash
python3 analysis/scripts/cc_gap_analysis.py cc10plus_roots_snapshot_2026-03-12.txt --prefix gap --min-cc 10
# Then extract twins from the gap_closest.csv
```
