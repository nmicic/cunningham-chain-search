# Cunningham Chain Immunization Analysis

## What is immunization?

For a **first-kind Cunningham Chain** starting at prime p, the chain members are:

```
member 0 = p
member 1 = 2p + 1
member 2 = 4p + 3
   ...
member j = 2^j * (p + 1) - 1
```

A small prime **q** divides chain member j when:

```
2^j * (p + 1) ≡ 1  (mod q)
```

If this happens at position j < chain_length, then member j is composite and the chain breaks.

A root p is **immune** to prime q when `(p + 1) ≡ 0 (mod q)`. In that case, `2^j * (p + 1) ≡ 0 (mod q)` for all j, so `2^j * (p + 1) - 1 ≡ -1 (mod q)`, which is never zero. The chain can never be broken by q.

## Three survival mechanisms

For a valid CC chain, every member must be prime. A root survives prime q in one of three ways:

| Mechanism | Condition | Effect |
|-----------|-----------|--------|
| **Immune** | `(p+1) ≡ 0 (mod q)` | q never divides any member — permanent protection |
| **Coset safe** | `(p+1) mod q ∉ ⟨2⟩ mod q` | 2^j * (p+1) never reaches 1 mod q — q never kills |
| **Beyond chain** | kill position j ≥ chain length | q would kill at position j, but chain ends before j |

## The role of ord(2, q)

The **multiplicative order** of 2 modulo q determines the cycle length of `2^j mod q`.
The subgroup `⟨2⟩ mod q` contains exactly `ord(2,q)` elements.

- If `(p+1) mod q` is in `⟨2⟩`, then some j < ord(2,q) kills the chain
- If `(p+1) mod q` is NOT in `⟨2⟩`, then q never kills (coset safe)
- If `(p+1) mod q = 0`, the root is immune (strongest protection)

### Key primes for CC10+

| q | ord(2,q) | \|⟨2⟩\| | Safe residues | CC10+ constraint |
|---|----------|---------|---------------|------------------|
| 3 | 2 | 2/2 | {0} | Must be immune |
| 5 | 4 | 4/4 | {0} | Must be immune |
| 7 | 3 | 3/6 | {0,3,5,6} | Immune or coset safe |
| 11 | 10 | 10/10 | {0} | Must be immune |
| 13 | 12 | 12/12 | {0} | Must be immune |
| 17 | 8 | 8/16 | {0,3,5,6,7,10,11,12,14} | Immune or coset safe |
| 19 | 18 | 18/18 | {0} | Must be immune (but kill at j≥10 OK) |
| 23 | 11 | 11/22 | {0,2,4,6,7,9,14,17,18,20,21,22} | Immune or coset safe |
| 29 | 28 | 28/28 | {0} | Kill pos often > 10, so many survive |
| 31 | 5 | 5/30 | {0} + 25 coset residues | Most roots coset safe |

When `⟨2⟩ = (Z/qZ)*` (2 is a primitive root mod q), the only safe residue is 0 (immune).
This happens for q = 3, 5, 11, 13, 19, 29 — for these primes, CC10+ roots are forced to be immune.

When `⟨2⟩` is a proper subgroup, there are additional safe residues in other cosets.
This happens for q = 7, 17, 23, 31 — roots can survive without being immune.

## Tools

### C analyzer (no GMP needed)

```bash
gcc -O2 -o cc_immunization_analysis cc_immunization_analysis.c -lm

# Full analysis → JSON
./cc_immunization_analysis cc10plus_roots_snapshot_2026-03-12.txt -o immunization.json

# Filter to 80+ bit roots
./cc_immunization_analysis cc10plus_roots_snapshot_2026-03-12.txt -o imm_80plus.json --min-bits 80

# Filter to CC14+
./cc_immunization_analysis cc10plus_roots_snapshot_2026-03-12.txt -o imm_cc14plus.json --min-cc 14
```

### Interactive dashboard

```bash
# Open dashboard, load JSON via file picker
open cc_immunization_dashboard.html

# Or bake JSON into standalone HTML
./embed_immunization.sh immunization.json cc_immunization_dashboard.html > standalone.html
```

### Input format

The tool reads the combined CC file format:
```
CC10 0x1A3F5B7 25
CC12 0xABCDEF123 36
CC14 0x1F 5
```

Each line: `CC<length> 0x<hex_root> <bit_size>`

### Output (JSON)

The JSON contains:
- `meta` — source file, filters, prime list
- `group_theory` — ord(2,q), subgroup, safe set per prime
- `per_cc` — immune% and safe% aggregated per chain length
- `per_bucket` — full detail per (cc, bits) pair
- `kill_positions` — which position each prime kills at
- `combinations` — immune fingerprint distribution (most common sets)
- `residue_distribution` — (p+1) mod q counts per prime

## Performance

The tool computes `p mod q` from hex digits directly (no big number library needed).
Processing 930K roots takes ~1.5 seconds.
