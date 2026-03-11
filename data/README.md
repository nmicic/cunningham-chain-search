# Campaign Dataset

## What's Included

### In this repository

| File | Description |
|------|-------------|
| `sample/cc10plus_roots_sample_1000.txt` | First 1000 roots from the full dataset (for testing scripts) |
| `gap_closest.csv` | 500 tightest root pairs with hex values |
| `gap_comprehensive.csv` | Per-(CC, bits) gap statistics |
| `gap_heatmap.csv` | Gap distribution bucketed by 2^X |
| `gap_pyramid.csv` | Gap distribution per CC level |
| `gap_statistics.csv` | Summary statistics per bucket |
| `cc_twins.csv` | Top CC-twin pairs (same-CC and cross-CC) |
| `kill_positions.csv` | Chain-breaker positions by small prime, per (CC, bits) bucket |
| `immunization_heatmap.csv` | Residue immunity and safe-residue rates per (CC, bits) bucket |

### Via GitHub Releases

The full raw snapshot (`cc10plus_roots_snapshot_2026-03-12.txt.gz`, ~30 MB uncompressed) is available as a GitHub Release asset. It contains 1,773,994 recorded roots in the format:

```
CC10 0xHEXVALUE DIGITS
CC11 0xHEXVALUE DIGITS
...
```

The published CC10+ analysis subset derived from that snapshot contains 929,574 roots.

**CC10+ subset summary:**

| CC level | Count |
|----------|------:|
| CC10 | 749,024 |
| CC11 | 132,580 |
| CC12 | 39,138 |
| CC13 | 7,226 |
| CC14 | 1,323 |
| CC15 | 238 |
| CC16 | 44 |
| CC17 | 1 |

Primary search band: 89-91 bits. Additional data from 59-153 bit ranges.

### CSV column notes

**`kill_positions.csv`** — Where chains break, by small prime.

| Column | Meaning |
|--------|---------|
| `cc` | Official chain length of the root |
| `bits` | Bit size of the root |
| `prime` | Small prime that divides the composite link |
| `position` | Chain position of the break (1 = first link after the root, 2 = second, etc.) |
| `count` | Number of roots in this bucket killed at this position by this prime |
| `pct_of_bucket` | Percentage of the (cc, bits) bucket |

**`immunization_heatmap.csv`** — Residue immunity summary per (cc, bits) bucket.

- `imm_<p>` / `imm_<p>_pct`: count and percentage of roots immune to prime `p` (i.e., `p` cannot kill any chain position).
- `safe_<p>` / `safe_<p>_pct`: count and percentage of roots where the residue mod `p` is "safe" — the root avoids the dangerous residue classes for that prime.
- `safe_all` / `safe_all_pct`: roots that satisfy all tracked safe-residue conditions simultaneously.

Both CSVs are derived analysis outputs generated from the released dataset snapshot, not hand-maintained files.

## What's Excluded

- Raw operational logs and checkpoint files
- Intermediate sieve outputs
- Per-prefix lane reports
- Server configuration and deployment details

## Reproducing Analysis

```bash
# Gap analysis (requires dataset file)
python3 analysis/scripts/cc_gap_analysis.py cc10plus_roots_snapshot_2026-03-12.txt --prefix gap --min-cc 10

# Ghost chain scan (~12 min on 8 cores)
python3 analysis/scripts/ghost_chains.py cc10plus_roots_snapshot_2026-03-12.txt

# Immunization analysis
python3 analysis/scripts/analyze_cc_immunization.py cc10plus_roots_snapshot_2026-03-12.txt
```
