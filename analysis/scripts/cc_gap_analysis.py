#!/usr/bin/env python3
# Copyright (c) 2026 Nenad Mićić <nenad@micic.be>
# SPDX-License-Identifier: Apache-2.0
"""CC Gap Analysis — Distance between consecutive roots per (CC, bit_size).

Reads combined CC file (format: "CC10 0xHEX <digits>").
For each (cc, bits) bucket, sorts roots and computes all consecutive gaps.

Outputs:
  1. gap_heatmap.csv     — gap distribution bucketed by 2^X
  2. gap_statistics.csv   — per-(cc,bits) summary: min/median/mean/max gap
  3. gap_raw_top.csv      — smallest gaps (closest pairs) across all CC10+
  4. gap_pyramid.csv      — aggregated gap distribution per CC level (all bits merged)

Usage:
    python3 cc_gap_analysis.py /path/to/CC_x.txt [--prefix gap] [--min-cc 10]
"""

import argparse
import math
import os
import re
import sys
from collections import defaultdict


def read_combined_file(path, min_cc=10):
    """Parse CC file. Returns dict: (cc, bits) -> sorted list of roots."""
    buckets = defaultdict(list)
    skipped = 0

    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 2:
                skipped += 1
                continue

            m = re.match(r"CC(\d+)$", parts[0])
            if not m:
                skipped += 1
                continue
            cc = int(m.group(1))
            if cc < min_cc:
                skipped += 1
                continue

            hex_str = parts[1]
            if hex_str.startswith(("0x", "0X")):
                hex_str = hex_str[2:]
            try:
                root = int(hex_str, 16)
            except ValueError:
                skipped += 1
                continue

            bits = root.bit_length()
            buckets[(cc, bits)].append(root)

    # Sort each bucket
    for key in buckets:
        buckets[key].sort()

    if skipped > 0:
        print(f"  (skipped {skipped} lines below CC{min_cc})", file=sys.stderr)
    return buckets


def gap_bucket(gap):
    """Return floor(log2(gap)) as the bucket index."""
    if gap <= 0:
        return 0
    return gap.bit_length() - 1


def analyze_gaps(buckets):
    """Compute gap statistics and distributions."""
    # Per-bucket statistics
    stats = {}
    # Per-bucket gap distribution (bucketed by 2^X)
    heatmap = {}
    # All gaps with metadata (for finding smallest)
    all_gaps = []

    for (cc, bits), roots in sorted(buckets.items()):
        n = len(roots)
        if n < 2:
            stats[(cc, bits)] = {
                "count": n, "pairs": 0,
                "min_gap": None, "max_gap": None,
                "median_gap": None, "mean_gap": None,
                "min_gap_log2": None, "median_gap_log2": None,
            }
            continue

        gaps = []
        for i in range(n - 1):
            g = roots[i + 1] - roots[i]
            gaps.append(g)
            # Track smallest gaps with context
            if g > 0:
                all_gaps.append((g, cc, bits, roots[i], roots[i + 1]))

        gaps.sort()
        total = sum(gaps)
        mid = len(gaps) // 2
        median = gaps[mid] if len(gaps) % 2 == 1 else (gaps[mid - 1] + gaps[mid]) // 2

        stats[(cc, bits)] = {
            "count": n,
            "pairs": len(gaps),
            "min_gap": gaps[0],
            "max_gap": gaps[-1],
            "median_gap": median,
            "mean_gap": total // len(gaps),
            "min_gap_log2": gap_bucket(gaps[0]),
            "median_gap_log2": gap_bucket(median),
            "mean_gap_log2": gap_bucket(total // len(gaps)),
        }

        # Build gap distribution for heatmap
        dist = defaultdict(int)
        for g in gaps:
            b = gap_bucket(g)
            dist[b] += 1
        heatmap[(cc, bits)] = dict(dist)

    return stats, heatmap, all_gaps


def write_gap_statistics(stats, path):
    """Write per-(cc,bits) gap summary."""
    with open(path, "w") as f:
        f.write("cc,bits,count,pairs,min_gap_log2,median_gap_log2,mean_gap_log2,"
                "min_gap_hex,median_gap_hex,mean_gap_hex,max_gap_hex\n")
        for (cc, bits) in sorted(stats.keys()):
            s = stats[(cc, bits)]
            if s["pairs"] == 0:
                f.write(f"{cc},{bits},{s['count']},0,,,,,,\n")
                continue
            f.write(f"{cc},{bits},{s['count']},{s['pairs']},"
                    f"{s['min_gap_log2']},{s['median_gap_log2']},{s['mean_gap_log2']},"
                    f"{hex(s['min_gap'])},{hex(s['median_gap'])},{hex(s['mean_gap'])},"
                    f"{hex(s['max_gap'])}\n")
    print(f"Wrote {path}", file=sys.stderr)


def write_gap_heatmap(heatmap, path):
    """Write gap distribution heatmap: cc,bits,gap_log2,count,pct."""
    with open(path, "w") as f:
        f.write("cc,bits,gap_log2,count,pct\n")
        for (cc, bits) in sorted(heatmap.keys()):
            dist = heatmap[(cc, bits)]
            total = sum(dist.values())
            for log2_bucket in sorted(dist.keys()):
                cnt = dist[log2_bucket]
                pct = 100.0 * cnt / total if total > 0 else 0
                f.write(f"{cc},{bits},{log2_bucket},{cnt},{pct:.2f}\n")
    print(f"Wrote {path}", file=sys.stderr)


def write_gap_pyramid(heatmap, stats, path):
    """Write aggregated gap distribution per CC level (all bits merged)."""
    # Merge across bit sizes for each CC
    cc_dist = defaultdict(lambda: defaultdict(int))
    cc_total_pairs = defaultdict(int)
    cc_total_roots = defaultdict(int)

    for (cc, bits), dist in heatmap.items():
        for log2_bucket, cnt in dist.items():
            cc_dist[cc][log2_bucket] += cnt
            cc_total_pairs[cc] += cnt
        cc_total_roots[cc] += stats[(cc, bits)]["count"]

    with open(path, "w") as f:
        f.write("cc,roots,pairs,gap_log2,count,pct,cumulative_pct\n")
        for cc in sorted(cc_dist.keys()):
            dist = cc_dist[cc]
            total = cc_total_pairs[cc]
            roots = cc_total_roots[cc]
            cum = 0.0
            for log2_bucket in sorted(dist.keys()):
                cnt = dist[log2_bucket]
                pct = 100.0 * cnt / total if total > 0 else 0
                cum += pct
                f.write(f"{cc},{roots},{total},{log2_bucket},{cnt},{pct:.2f},{cum:.2f}\n")
    print(f"Wrote {path}", file=sys.stderr)


def write_top_gaps(all_gaps, path, top_n=500):
    """Write the smallest gaps (closest pairs)."""
    all_gaps.sort(key=lambda x: x[0])
    with open(path, "w") as f:
        f.write("gap_log2,gap_hex,gap_decimal_approx,cc,bits,root1_hex,root2_hex\n")
        for i, (g, cc, bits, r1, r2) in enumerate(all_gaps[:top_n]):
            gl2 = gap_bucket(g)
            # Decimal approximation for readability
            if g.bit_length() > 64:
                dec_approx = f"~2^{g.bit_length()-1}"
            else:
                dec_approx = str(g)
            f.write(f"{gl2},{hex(g)},{dec_approx},{cc},{bits},{hex(r1)},{hex(r2)}\n")
    print(f"Wrote {path} ({min(top_n, len(all_gaps))} closest pairs)", file=sys.stderr)


def print_summary(stats, heatmap, all_gaps):
    """Print summary to stderr."""
    print("\n=== GAP ANALYSIS SUMMARY ===", file=sys.stderr)

    # Per-CC summary
    cc_data = defaultdict(lambda: {"roots": 0, "pairs": 0, "min_gaps": [], "median_gaps": []})
    for (cc, bits), s in stats.items():
        d = cc_data[cc]
        d["roots"] += s["count"]
        d["pairs"] += s["pairs"]
        if s["min_gap"] is not None:
            d["min_gaps"].append(s["min_gap"])
            d["median_gaps"].append(s["median_gap"])

    print(f"\n{'CC':>4} {'Roots':>8} {'Pairs':>8} {'Min gap':>12} {'Min log2':>9} "
          f"{'Median gap':>12} {'Med log2':>9}", file=sys.stderr)
    print(f"{'----':>4} {'--------':>8} {'--------':>8} {'------------':>12} {'---------':>9} "
          f"{'------------':>12} {'---------':>9}", file=sys.stderr)

    for cc in sorted(cc_data.keys()):
        d = cc_data[cc]
        if not d["min_gaps"]:
            print(f"CC{cc:>2} {d['roots']:>8,} {d['pairs']:>8} {'N/A':>12} {'':>9} {'N/A':>12}", file=sys.stderr)
            continue
        overall_min = min(d["min_gaps"])
        # Weighted median approximation: use min of per-bucket medians
        overall_med = sorted(d["median_gaps"])[len(d["median_gaps"]) // 2]
        print(f"CC{cc:>2} {d['roots']:>8,} {d['pairs']:>8,} "
              f"{'2^' + str(gap_bucket(overall_min)):>12} {gap_bucket(overall_min):>9} "
              f"{'2^' + str(gap_bucket(overall_med)):>12} {gap_bucket(overall_med):>9}", file=sys.stderr)

    # Top 10 closest pairs
    all_gaps.sort(key=lambda x: x[0])
    print(f"\nTop 20 closest pairs:", file=sys.stderr)
    print(f"{'Gap':>14} {'log2':>5} {'CC':>4} {'Bits':>5}  Root1 → Root2", file=sys.stderr)
    for g, cc, bits, r1, r2 in all_gaps[:20]:
        gl2 = gap_bucket(g)
        print(f"{'2^'+str(gl2):>14} {gl2:>5} CC{cc:>2} {bits:>5}  "
              f"{hex(r1)[:24]}... → {hex(r2)[:24]}...", file=sys.stderr)

    # Expected vs observed gap analysis
    print(f"\n--- Expected gap (uniform random) vs observed ---", file=sys.stderr)
    for cc in sorted(cc_data.keys()):
        d = cc_data[cc]
        if d["pairs"] == 0:
            continue
        # For uniform random in [2^(b-1), 2^b), expected gap ≈ 2^b / count
        # Group by dominant bit size
        bit_counts = defaultdict(int)
        for (c, b), s in stats.items():
            if c == cc and s["count"] > 10:
                bit_counts[b] += s["count"]
        if not bit_counts:
            continue
        dom_bits = max(bit_counts, key=bit_counts.get)
        dom_count = bit_counts[dom_bits]
        expected_gap_log2 = dom_bits - 1 - (dom_count.bit_length() - 1)  # 2^bits / count
        actual_min = min(s["min_gap_log2"] for (c, b), s in stats.items()
                        if c == cc and s["min_gap_log2"] is not None and b == dom_bits)
        actual_med = None
        for (c, b), s in stats.items():
            if c == cc and b == dom_bits and s["median_gap_log2"] is not None:
                actual_med = s["median_gap_log2"]
                break
        print(f"  CC{cc} ({dom_bits}-bit, {dom_count:,} roots): "
              f"expected gap ~2^{expected_gap_log2}, "
              f"min observed 2^{actual_min}, "
              f"median 2^{actual_med}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description="CC Gap Analysis")
    parser.add_argument("input_file", help="Combined CC file (e.g., CC_x.txt)")
    parser.add_argument("--prefix", "-p", default="gap",
                        help="Output file prefix (default: gap)")
    parser.add_argument("--min-cc", type=int, default=10,
                        help="Minimum CC level to analyze (default: 10)")
    parser.add_argument("--top-gaps", type=int, default=500,
                        help="Number of smallest gaps to export (default: 500)")
    args = parser.parse_args()

    if not os.path.isfile(args.input_file):
        print(f"ERROR: {args.input_file} not found", file=sys.stderr)
        sys.exit(1)

    print(f"Reading {args.input_file} (CC{args.min_cc}+)...", file=sys.stderr)
    buckets = read_combined_file(args.input_file, min_cc=args.min_cc)

    total_roots = sum(len(v) for v in buckets.values())
    total_buckets = len(buckets)
    print(f"Loaded {total_roots:,} roots in {total_buckets} (cc,bits) buckets", file=sys.stderr)

    print("Computing gaps...", file=sys.stderr)
    stats, heatmap, all_gaps = analyze_gaps(buckets)
    print(f"Computed {len(all_gaps):,} gaps", file=sys.stderr)

    # Write outputs
    write_gap_statistics(stats, f"{args.prefix}_statistics.csv")
    write_gap_heatmap(heatmap, f"{args.prefix}_heatmap.csv")
    write_gap_pyramid(heatmap, stats, f"{args.prefix}_pyramid.csv")
    write_top_gaps(all_gaps, f"{args.prefix}_closest.csv", top_n=args.top_gaps)

    print_summary(stats, heatmap, all_gaps)


if __name__ == "__main__":
    main()
