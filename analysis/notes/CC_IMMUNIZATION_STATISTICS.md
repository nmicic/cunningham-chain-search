# CC Immunization Statistics

Immunization analysis of 929,574 roots from the NEW4 dataset (CC10–CC17, 59–153 bits).

## How Roots Survive Small Primes

For a Cunningham chain of length L starting at root p, every chain member `p, 2p+1, 4p+3, ...` must be prime. A small prime q can kill the chain if it divides any member. There are three survival mechanisms:

1. **Immune**: `(p+1) ≡ 0 (mod q)` — q never divides any chain member
2. **Coset safe**: `(p+1) mod q` is not in the multiplicative subgroup `⟨2⟩ mod q` — q is structurally unable to kill
3. **Kill beyond chain**: the kill position `j ≥ L` — q would kill at position j, but the chain only has L members

## Group Theory Constants

| Prime q | ord(2,q) | \|⟨2⟩\| | Safe residues | Safe set | P(immune\|safe) |
|--------:|:--------:|:--------:|:-------------:|----------|:----------------:|
| 3 | 2 | 2 | 1 | {0} | 100% |
| 5 | 4 | 4 | 1 | {0} | 100% |
| 7 | 3 | 3 | 4 | {0,3,5,6} | 25.0% |
| 11 | 10 | 10 | 1 | {0} | 100% |
| 13 | 12 | 12 | 1 | {0} | 100% |
| 17 | 8 | 8 | 9 | {0,3,5,6,7,10,11,12,14} | 11.1% |
| 19 | 18 | 18 | 1 | {0} | 100% |
| 23 | 11 | 11 | 12 | {0,5,7,10,11,14,15,17,19,20,21,22} | 8.3% |
| 29 | 28 | 28 | 1 | {0} | 100% |
| 31 | 5 | 5 | 26 | {0,3,5,6,7,9,10,11,12,13,...} | 3.8% |

**Key insight:** For q = 3, 5, 11, 13, 19 — the subgroup `⟨2⟩` covers all of `(Z/qZ)*`, so the only safe residue is 0 (immune). These primes are **mandatory immune** for CC10+.

For q = 7, 17, 23, 29, 31 — the subgroup is smaller than `(Z/qZ)*`, leaving non-zero safe residues (coset-safe). These primes are **optional** — roots can survive without immunity.

## Overall Immunization Rates

### Aggregated across all CC lengths (929,574 roots)

| Prime q | Immune % | Safe % | Random 1/q | ord(2,q) |
|--------:|:--------:|:------:|:----------:|:--------:|
| 3 | **99.54%** | 99.54% | 33.3% | 2 |
| 5 | **99.54%** | 99.54% | 20.0% | 4 |
| 7 | 23.56% | **99.78%** | 14.3% | 3 |
| 11 | **99.54%** | 99.54% | 9.1% | 10 |
| 13 | **99.54%** | 99.54% | 7.7% | 12 |
| 17 | 12.24% | **99.74%** | 5.9% | 8 |
| 19 | **97.93%** | 100.00% | 5.3% | 18 |
| 23 | 8.30% | **99.77%** | 4.3% | 11 |
| 29 | 9.93% | **99.84%** | 3.4% | 28 |
| 31 | 3.83% | **99.93%** | 3.2% | 5 |

**Why not 100% immune for q=3,5,11,13?** The 0.46% non-immune roots are CC12+ "sub-chain" roots that survive through coset safety and kill-beyond-chain (see CC12+ anomaly below).

**Why not 100% immune for q=19?** With ord(2,19)=18, CC10–CC17 roots can survive without immunity if the kill position is ≥ chain length.

## Per-CC Immunization Breakdown

### Mandatory primes: q = 3, 5, 11, 13

| CC | Roots | Immune to {3,5,11,13} | Not immune % | Notes |
|----|------:|:----------------------:|:------------:|-------|
| CC10 | 749,024 | 100.00% | 0.00% | All immune as expected |
| CC11 | 132,580 | 100.00% | 0.00% | All immune as expected |
| CC12 | 39,138 | 91.14% | **8.86%** | Sub-chains appear |
| CC13 | 7,226 | 90.63% | **9.37%** | Growing fraction |
| CC14 | 1,323 | 90.17% | **9.83%** | Consistent |
| CC15 | 238 | 86.55% | **13.45%** | More surprising roots |
| CC16 | 44 | 86.36% | **13.64%** | 6 of 44 lack full base |
| CC17 | 1 | 100.00% | 0.00% | Only root has full base |

### q = 19 (near-mandatory, ord = 18)

| CC | Roots | Immune to 19 | Not immune % | Why they survive |
|----|------:|:------------:|:------------:|-----------------|
| CC10 | 749,024 | 99.04% | **0.96%** | Kill at pos 10–17 |
| CC11 | 132,580 | 99.05% | **0.95%** | Kill at pos 11–17 |
| CC12 | 39,138 | 77.95% | **22.05%** | Kill at pos 12–17 |
| CC13 | 7,226 | 76.27% | **23.73%** | Kill at pos 13–17 |
| CC14 | 1,323 | 76.04% | **23.96%** | Kill at pos 14–17 |
| CC15 | 238 | 69.75% | **30.25%** | Kill at pos 15–17 |
| CC16 | 44 | 68.18% | **31.82%** | Kill at pos 16–17 |
| CC17 | 1 | 100.00% | 0.00% | Must be immune (only pos 17 left, but ord=18) |

**The q=19 anomaly**: at CC12+, roughly 22–32% of roots survive without immunity to 19, because their kill position falls beyond the actual chain length. This fraction grows with chain length because fewer kill positions remain "in range."

### q = 7 (optional, ord = 3)

| CC | Immune to 7 | Safe (coset) | Notes |
|----|:-----------:|:------------:|-------|
| CC10 | 23.61% | 100.00% | 76.4% survive via coset {3,5,6} |
| CC11 | 23.74% | 100.00% | Same pattern |
| CC12 | 22.16% | 95.81% | Sub-chains slightly shift |
| CC15 | 20.59% | 91.60% | |
| CC16 | 25.00% | 95.45% | Small sample |

### q = 17 (optional, ord = 8)

| CC | Immune to 17 | Safe (coset) |
|----|:------------:|:------------:|
| CC10 | 12.28% | 100.00% |
| CC11 | 12.43% | 100.00% |
| CC12 | 11.13% | 95.00% |
| CC15 | 10.50% | 92.02% |
| CC16 | 13.64% | 97.73% |

### q = 31 (optional, ord = 5)

| CC | Immune to 31 | Safe (coset) |
|----|:------------:|:------------:|
| CC10 | 3.87% | 100.00% |
| CC12 | 3.40% | 98.63% |
| CC15 | 3.78% | 95.80% |
| CC16 | **0.00%** | 97.73% |
| CC17 | **0.00%** | 100.00% |

The CC16–CC17 roots survive q=31 entirely through coset safety — zero immunity, yet 100% safe.

---

## Residue Distribution: (p+1) mod q

The actual residue distribution confirms the group-theoretic predictions.

### q = 7 (ord = 3, subgroup {1,2,4})

| Residue | Meaning | Count | % |
|:-------:|---------|------:|--:|
| 0 | **Immune** | 218,996 | 23.56% |
| 3 | Coset safe | 236,622 | 25.45% |
| 5 | Coset safe | 235,300 | 25.31% |
| 6 | Coset safe | 236,621 | 25.45% |
| 1 | Killed (in ⟨2⟩) | 0 | 0% |
| 2 | Killed (in ⟨2⟩) | 985 | 0.11% |
| 4 | Killed (in ⟨2⟩) | 1,050 | 0.11% |

The 4 safe residues {0,3,5,6} capture 99.78% of all roots. The ~0.22% in killed residues {2,4} are CC12+ sub-chains surviving via kill-beyond-chain.

### q = 17 (ord = 8, subgroup {1,2,4,8,16,15,13,9})

| Category | Residues | Count | % |
|----------|----------|------:|--:|
| Immune | {0} | 113,740 | 12.24% |
| Coset safe | {3,5,6,7,10,11,12,14} | 813,411 | 87.51% |
| Killed (in ⟨2⟩) | {2,4,8,9,13,15,16} | 2,423 | 0.26% |

### q = 23 (ord = 11)

| Category | Count | % |
|----------|------:|--:|
| Immune | 77,160 | 8.30% |
| Coset safe | 849,413 | 91.37% |
| In subgroup (killed) | 3,001 | 0.32% |

### q = 31 (ord = 5, subgroup {1,2,4,8,16})

| Category | Count | % |
|----------|------:|--:|
| Immune | 35,627 | 3.83% |
| Coset safe (26 residues) | 893,277 | 96.10% |
| In subgroup (killed) | 670 | 0.07% |

**Pattern:** for optional primes, the fraction of roots that are immune closely matches `1/|safe_set|`:
- q=7: 23.56% ≈ 1/4 = 25%
- q=17: 12.24% ≈ 1/9 = 11.1%
- q=23: 8.30% ≈ 1/12 = 8.3%
- q=31: 3.83% ≈ 1/26 = 3.8%

This confirms roots distribute uniformly across safe residues — the sieve selects for safety, not specifically for immunity.

---

## Immune Combination Fingerprints

Each root's immune fingerprint is the set of primes q where `(p+1) ≡ 0 (mod q)`.

### Global Distribution (929,574 roots, 58 unique fingerprints)

| Rank | Fingerprint | # primes | Count | % |
|-----:|------------|:--------:|------:|--:|
| 1 | {3, 5, 11, 13, 19} | 5 | 480,840 | **51.73%** |
| 2 | {3, 5, **7**, 11, 13, 19} | 6 | 151,399 | 16.29% |
| 3 | {3, 5, 11, 13, **17**, 19} | 6 | 69,995 | 7.53% |
| 4 | {3, 5, 11, 13, 19, **29**} | 6 | 52,967 | 5.70% |
| 5 | {3, 5, 11, 13, 19, **23**} | 6 | 43,477 | 4.68% |
| 6 | {3, 5, 11, 13, 19, **31**} | 6 | 19,110 | 2.06% |
| 7 | {3, 5, **7**, 11, 13, **17**, 19} | 7 | 18,490 | 1.99% |
| 8 | {3, 5, **7**, 11, 13, 19, **29**} | 7 | 16,976 | 1.83% |
| 9 | {3, 5, **7**, 11, 13, 19, **23**} | 7 | 13,888 | 1.49% |
| 10 | {3, 5, 11, 13, **17**, 19, **29**} | 7 | 6,378 | 0.69% |
| 11 | {3, 5, 11, 13, **17**, 19, **23**} | 7 | 5,172 | 0.56% |
| 12 | {3, 5, **7**, 11, 13, 19, **31**} | 7 | 4,975 | 0.54% |
| 13 | {3, 5, 11, 13, 19, **23**, **29**} | 7 | 4,009 | 0.43% |
| 14 | {3, 5, 11, 13} | **4** | 8,852 | **0.95%** |
| 15 | {3, 5, 11, 13, **17**, 19, **31**} | 7 | 2,340 | 0.25% |

### Rarest Fingerprints

| Rank | Fingerprint | # primes | Count |
|-----:|------------|:--------:|------:|
| 54 | {3, 5, 11, 13, 23, 31} | 6 | 24 |
| 55 | {3, 5, 7, 11, 13, 17, 31} | 7 | 11 |
| 56 | {3, 5, 7, 11, 13, 23, 31} | 7 | 3 |
| 57 | {3, 5, 7, 11, 13, 23, 29, 31} | 8 | 2 |
| 58 | {3, 5, 7, 11, 13, 29, 31} | 7 | 1 |

**Zero roots** were observed with immunity to all 10 primes simultaneously.

### "Surprising" Fingerprints — Roots Missing Mandatory Primes

| CC | Total | Missing ≥1 of {3,5,11,13} | {} (none immune) | {19} only |
|----|------:|:--------------------------:|:------------------:|:---------:|
| CC10 | 749,024 | 0 | 0 | 0 |
| CC11 | 132,580 | 0 | 0 | 0 |
| CC12 | 39,138 | 10,324 (26.4%) | 1,577 (4.0%) | 1,585 (4.1%) |
| CC13 | 7,226 | 2,059 (28.5%) | 301 (4.2%) | 324 (4.5%) |
| CC14 | 1,323 | 388 (29.3%) | 53 (4.0%) | 62 (4.7%) |
| CC15 | 238 | 89 (37.4%) | 13 (5.5%) | 14 (5.9%) |
| CC16 | 44 | 16 (36.4%) | 1 (2.3%) | 2 (4.5%) |

At CC12+, these roots survive purely through coset safety and kill-beyond-chain mechanics.

### Fingerprint Diversity by Chain Length

| CC | Roots | Unique combos | Top fingerprint (%) | Missing base (%) |
|----|------:|--------------:|:-------------------:|:----------------:|
| CC10 | 749,024 | 61 | {3,5,11,13,19} (52.5%) | 0.96% |
| CC11 | 132,580 | 58 | {3,5,11,13,19} (52.2%) | 0.95% |
| CC12 | 39,138 | 60 | {3,5,11,13,19} (38.8%) | **26.4%** |
| CC13 | 7,226 | 51 | {3,5,11,13,19} (37.7%) | **28.5%** |
| CC14 | 1,323 | 34 | {3,5,11,13,19} (38.5%) | **29.3%** |
| CC15 | 238 | 28 | {3,5,11,13,19} (30.3%) | **37.4%** |
| CC16 | 44 | 10 | {3,5,11,13,19} (36.4%) | **36.4%** |
| CC17 | 1 | 1 | {3,5,11,13,19} (100%) | 0% |

**The longer the chain, the more likely it is to have a "surprising" fingerprint.**

---

## CC12+ Combination Breakdown

### CC14 (1,323 roots, 34 unique combos)

| Fingerprint | Primes | Count | % |
|------------|:------:|------:|--:|
| {3,5,11,13,19} | 5 | 510 | 38.55% |
| {3,5,7,11,13,19} | 6 | 162 | 12.24% |
| {3,5,11,13} | 4 | 132 | 9.98% |
| {3,5,11,13,17,19} | 6 | 70 | 5.29% |
| {19} | 1 | 62 | 4.69% |
| {} | 0 | 53 | 4.01% |
| {3,5,7,11,13} | 5 | 51 | 3.85% |
| {3,5,11,13,19,29} | 6 | 47 | 3.55% |
| {3,5,11,13,19,23} | 6 | 37 | 2.80% |

### CC15 (238 roots, 28 unique combos)

| Fingerprint | Primes | Count | % |
|------------|:------:|------:|--:|
| {3,5,11,13,19} | 5 | 72 | 30.25% |
| {3,5,11,13} | 4 | 29 | 12.18% |
| {3,5,7,11,13,19} | 6 | 23 | 9.66% |
| {19} | 1 | 14 | 5.88% |
| {} | 0 | 13 | 5.46% |
| {3,5,11,13,17,19} | 6 | 13 | 5.46% |
| {3,5,11,13,19,29} | 6 | 13 | 5.46% |
| {3,5,7,11,13} | 5 | 10 | 4.20% |
| {3,5,11,13,19,23} | 6 | 8 | 3.36% |

### CC16 (44 roots, 10 unique combos)

| Fingerprint | Primes | Count | % |
|------------|:------:|------:|--:|
| {3,5,11,13,19} | 5 | 16 | 36.36% |
| {3,5,11,13} | 4 | 6 | 13.64% |
| {3,5,7,11,13,19} | 6 | 6 | 13.64% |
| {3,5,11,13,17,19} | 6 | 5 | 11.36% |
| {3,5,7,11,13} | 5 | 3 | 6.82% |
| {29} | 1 | 3 | 6.82% |
| {19} | 1 | 2 | 4.55% |
| {} | 0 | 1 | 2.27% |
| {3,5,7,11,13,17,19} | 7 | 1 | 2.27% |
| {3,5,7,11,13,23} | 6 | 1 | 2.27% |

**Notable:** 1 CC16 root has fingerprint `{}` — immune to NONE of the 10 primes, surviving entirely through coset safety and kill-beyond-chain.

---

## Primorial Search Miss Rate

Primorial-based searches force `p+1 ≡ 0 (mod primorial)`, which excludes roots that survive via coset safety.

| Strategy | Forces immunity to | Miss % (all) | Miss % (CC14+) |
|----------|-------------------|:------------:|:---------------:|
| 5# = 30 | {3, 5} | 0.5% | 10.5% |
| **7# = 210** | **{3, 5, 7}** | **76.4%** | **78.8%** |
| 11# = 2,310 | {3, 5, 7, 11} | 76.4% | 78.8% |
| 13# = 30,030 | {3, 5, 7, 11, 13} | 76.4% | 78.8% |
| **17# = 510,510** | **{3, 5, 7, 11, 13, 17}** | **97.4%** | **98.1%** |
| 19# | {3, 5, 7, 11, 13, 17, 19} | 97.5% | 98.5% |
| 23# | all through 23 | 99.8% | 100.0% |
| All 10 | All 10 primes | 100.0% | 100.0% |

**The moment you include 7 in the primorial, you lose 76% of all valid chains.**

### Why 7# is catastrophic

For q=7, the safe residues are {0, 3, 5, 6}. A primorial search forces residue 0 (immune), discarding the 3 coset-safe residues. Since 3/4 of safe roots use coset safety for q=7, you lose 75% of the search space at this prime alone.

### Optimal CC18 search

The data-driven approach for CC18:
1. **Force** immunity to {3, 5, 11, 13, 19} — group theory mandates these
2. **Allow** any safe residue for {7, 17, 23, 29, 31} — coset safe is sufficient
3. This preserves `4/7 × 9/17 × 12/23 × 1 × 26/31 ≈ 13.6%` of search space that a primorial approach would discard

---

## The CC17 Root

The single CC17 root in the dataset:

| Property | Value |
|----------|-------|
| Fingerprint | {3, 5, 11, 13, 19} |
| Immune to 7? | No (coset safe) |
| Immune to 17? | No (coset safe) |
| Immune to 23? | No (coset safe) |
| Immune to 29? | No (coset safe) |
| Immune to 31? | No (coset safe) |

This root has the **minimum mandatory fingerprint** — it would be completely invisible to any search using `7#` or higher primorial multiplier.
