# Immune Combination Fingerprints

## What is a fingerprint?

Every CC root p has an **immune fingerprint**: the set of small primes q where `(p+1) ≡ 0 (mod q)`.

For example, if p+1 is divisible by 3, 5, 11, 13, and 19, the fingerprint is `{3,5,11,13,19}`.

## Mandatory base

For CC10+ first-kind chains, group theory forces certain primes to always appear in the fingerprint:

| Prime q | ord(2,q) | 2 is primitive root? | Mandatory for CC10+? |
|---------|----------|---------------------|---------------------|
| 3 | 2 | Yes | **Yes** — must be immune |
| 5 | 4 | Yes | **Yes** — must be immune |
| 11 | 10 | Yes | **Yes** — kill within 10 positions |
| 13 | 12 | Yes | **Yes** — kill within 12 positions |
| 19 | 18 | Yes | **Almost** — kill within 18, but pos ≥ 10 OK |

So the **mandatory base** for most CC10+ roots is `{3, 5, 11, 13, 19}`.

Primes 7, 17, 23, 29, 31 are optional — roots survive these via coset safety or kill-beyond-chain.

## Observed distribution (NEW4 dataset, 929K roots, 59–153 bits)

**58 unique fingerprints** observed across 929,574 roots.

### Top 10 (most common)

| Rank | Fingerprint | # primes | Count | % |
|------|------------|----------|-------|---|
| 1 | {3, 5, 11, 13, 19} | 5 | 480,840 | **51.73%** |
| 2 | {3, 5, **7**, 11, 13, 19} | 6 | 151,399 | 16.29% |
| 3 | {3, 5, 11, 13, **17**, 19} | 6 | 69,995 | 7.53% |
| 4 | {3, 5, 11, 13, 19, **29**} | 6 | 52,967 | 5.70% |
| 5 | {3, 5, 11, 13, 19, **23**} | 6 | 43,477 | 4.68% |
| 6 | {3, 5, 11, 13, 19, **31**} | 6 | 19,110 | 2.06% |
| 7 | {3, 5, **7**, 11, 13, **17**, 19} | 7 | 18,490 | 1.99% |
| 8 | {3, 5, **7**, 11, 13, 19, **29**} | 7 | 16,976 | 1.83% |
| 9 | {3, 5, **7**, 11, 13, 19, **23**} | 7 | 13,888 | 1.49% |
| 10 | {3, 5, 11, 13} | **4** | 8,852 | 0.95% |

### Bottom 5 (rarest observed)

| Rank | Fingerprint | # primes | Count |
|------|------------|----------|-------|
| 54 | {3, 5, 11, 13, 23, 31} | 6 | 24 |
| 55 | {3, 5, 7, 11, 13, 17, 31} | 7 | 11 |
| 56 | {3, 5, 7, 11, 13, 23, 31} | 7 | 3 |
| 57 | {3, 5, 7, 11, 13, 23, 29, 31} | 8 | 2 |
| 58 | {3, 5, 7, 11, 13, 29, 31} | 7 | 1 |

**Zero roots** were observed with full immunity to all 10 primes.

### Why these percentages?

Among roots with the mandatory base `{3,5,11,13,19}`, the probability of also being immune to prime q is `1/|safe_set(q)|`:

| Extra prime | P(immune \| safe) | Expected add-on % | Observed |
|-------------|-------------------|-------------------|----------|
| 7 | 1/4 = 25.0% | ~25% of base | ~31% (16.29/51.73) |
| 17 | 1/9 = 11.1% | ~11% of base | ~15% (7.53/51.73) |
| 23 | 1/12 = 8.3% | ~8% of base | ~9% (4.68/51.73) |
| 29 | 1/1 = 100% | Special (see below) | ~11% |
| 31 | 1/26 = 3.8% | ~4% of base | ~4% (2.06/51.73) |

**q=29**: ord(2,29) = 28, so `⟨2⟩` covers all of (Z/29Z)*. For CC10, kill position can be 0-27,
but only positions 0-9 matter. So 18/28 = 64% of non-immune residues have kill position ≥ 10.
Many roots survive q=29 without immunity.

---

## The critical finding: roots that BREAK the "mandatory base" pattern

**2.29% of all roots** (21,322) are missing at least one prime from {3,5,11,13,19}.

This seems impossible — group theory says these primes should be mandatory. What's happening?

### The q=19 anomaly

Prime 19 has `ord(2,19) = 18`. For a CC10 chain, the kill position only matters if it's < 10.
Since ord=18, the kill position for non-immune residues cycles through 0-17. Positions 10-17
are **beyond the chain length**, so the root survives.

This creates a loophole: **CC10-CC17 chains can exist without being immune to 19.**

| CC length | Not immune to 19 | % of that CC | Why they survive |
|-----------|------------------|-------------|-----------------|
| CC10 | 7,180 | 0.96% | Kill at pos 10-17 (beyond CC10) |
| CC11 | 1,254 | 0.95% | Kill at pos 11-17 |
| CC12 | 8,617 | **22.0%** | Kill at pos 12-17 |
| CC13 | 1,715 | **23.7%** | Kill at pos 13-17 |
| CC14 | 317 | **24.0%** | Kill at pos 14-17 |
| CC15 | 72 | **30.3%** | Kill at pos 15-17 |
| CC16 | 14 | **31.8%** | Kill at pos 16-17 |
| CC17 | 0 | 0% | Must be immune (ord=18, only pos 17 left) |

**Key insight**: for CC12-CC16, roughly **25-32%** of roots survive without immunity to 19.
These roots have `(p+1) mod 19` in a position that kills at chain position 12-17 — beyond
the actual chain length but before position 18.

### The CC12+ "sub-chain" roots

At CC12+, we also see roots with **completely empty fingerprints** `{}` — immune to NONE of
the 10 small primes. These are CC12 sub-chains within the data that survive purely through
coset safety and kill-beyond-chain for every prime:

| CC | Fingerprint `{}` (no immunity) | Fingerprint `{19}` only | Combined "surprising" |
|----|-------------------------------|------------------------|----------------------|
| CC12 | 1,577 (4.0%) | 1,585 (4.1%) | 3,162 (8.1%) |
| CC13 | 301 (4.2%) | 324 (4.5%) | 625 (8.7%) |
| CC14 | 53 (4.0%) | 62 (4.7%) | 115 (8.7%) |
| CC15 | 13 (5.5%) | 14 (5.9%) | 27 (11.3%) |
| CC16 | 1 (2.3%) | 2 (4.5%) | 3 (6.8%) |

These roots would be **completely invisible** to any search that forces even basic immunity.

---

## Why primorial searches miss most chains

A common approach in CC searching uses a **primorial multiplier**: construct candidates as
`p = k × primorial# × 2 - 1` to force `(p+1) ≡ 0` for all primes in the primorial.

This is a catastrophic constraint. The data proves it:

### Miss rate by primorial level

| Search strategy | Forces immunity to | Miss % (all) | Miss % (CC14+) |
|----------------|-------------------|-------------|----------------|
| 5# = 30 | {3, 5} | 0.5% | 10.5% |
| **7# = 210** | **{3, 5, 7}** | **76.4%** | **78.8%** |
| 11# = 2310 | {3, 5, 7, 11} | 76.4% | 78.8% |
| 13# = 30030 | {3, 5, 7, 11, 13} | 76.4% | 78.8% |
| **17# = 510510** | **{3, 5, 7, 11, 13, 17}** | **97.4%** | **98.1%** |
| 19# | {3, 5, 7, 11, 13, 17, 19} | 97.5% | 98.5% |
| 23# | {3, 5, 7, 11, 13, 17, 19, 23} | 99.8% | 100.0% |
| All 10 primes | All of them | 100.0% | 100.0% |

**The moment you include 7 in the primorial, you lose 76% of all valid chains.**

This is because only 1/4 of safe residues for q=7 are the immune residue 0.
The other 3/4 survive via coset safety: `(p+1) mod 7 ∈ {3, 5, 6}`.

Adding 17 pushes the miss rate to 97% — only 1/9 of safe residues for q=17 are immune.

### Why this matters for CC18/CC19

For CC18, the constraints tighten: `ord(2,19) = 18 = chain_length`, so 19 becomes
mandatory. But primes 7, 17, 23, 31 remain optional (survived via coset).

A primorial search forcing `19# = {3,5,7,11,13,17,19}` would:
- Correctly force immunity to 3, 5, 11, 13, 19 (mandatory for CC18)
- **Unnecessarily** force immunity to 7 and 17
- Miss **97.5%** of CC18-capable residue classes

The optimal CC18 search should:
1. **Force** immunity to {3, 5, 11, 13, 19} (group theory mandates these)
2. **Allow** any safe residue for {7, 17, 23, 29, 31} (coset safe is fine)
3. This preserves `4/7 × 9/17 × 12/23 × 19/29 × 26/31 ≈ 13.6%` of the search space
   that a primorial search would discard

For CC19, q=19 is still mandatory (ord=18 < 19), plus q=23 becomes tighter
(ord(2,23)=11, so kill positions cycle faster). But 7, 17, 31 remain coset-safe.

---

## Fingerprint diversity by chain length

| CC | Roots | Unique combos | Top fingerprint (%) | Missing base (%) | Trend |
|----|-------|--------------|--------------------|-----------------|----|
| CC10 | 748,989 | 50 | {3,5,11,13,19} (52.5%) | 0.96% | Base dominates |
| CC11 | 132,575 | 50 | {3,5,11,13,19} (52.2%) | 0.95% | Same |
| CC12 | 39,117 | 50 | {3,5,11,13,19} (38.8%) | **26.4%** | CC12 sub-chains appear |
| CC13 | 7,225 | 50 | {3,5,11,13,19} (37.7%) | **28.5%** | Sub-chains grow |
| CC14 | 1,323 | 34 | {3,5,11,13,19} (38.5%) | **29.3%** | Consistent |
| CC15 | 238 | 28 | {3,5,11,13,19} (30.3%) | **37.4%** | Minority have full base |
| CC16 | 44 | 10 | {3,5,11,13,19} (36.4%) | **36.4%** | A third lack full base |
| CC17 | 1 | 1 | {3,5,11,13,19} (100%) | 0% | Sample size = 1 |

**The longer the chain, the MORE likely it is to have a "surprising" fingerprint.**

At CC15-CC16, over a third of roots are missing at least one "mandatory" prime — they
survive through coset safety and kill-beyond-chain. This fraction would only grow at CC18+
as the search space gets sparser and unusual residue patterns become more common
among survivors.

---

## Re-running the analysis

```bash
# When new data arrives:
./cc_immunization_analysis new_dataset.txt -o immunization_new.json

# Compare 80+ bit vs all:
./cc_immunization_analysis cc10plus_roots_snapshot_2026-03-12.txt -o imm_all.json
./cc_immunization_analysis cc10plus_roots_snapshot_2026-03-12.txt -o imm_80plus.json --min-bits 80

# Open dashboard, load each JSON to compare
open cc_immunization_dashboard.html
```
