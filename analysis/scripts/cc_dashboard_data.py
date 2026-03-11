#!/usr/bin/env python3
# Copyright (c) 2026 Nenad Mićić <nenad@micic.be>
# SPDX-License-Identifier: Apache-2.0
"""CC Campaign Dashboard — Data Pipeline.

Reads all CC*_*.txt from a data directory, produces cc_dashboard_data.json
for the companion cc_dashboard.html dashboard.

Usage:
    python3 cc_dashboard_data.py /path/to/NEW2 [--output cc_dashboard_data.json]
"""

import argparse
import glob
import json
import math
import os
import re
import sys
from collections import Counter, defaultdict


# ---------------------------------------------------------------------------
# File parsing
# ---------------------------------------------------------------------------

def parse_filename(path):
    """Extract (cc_length, bits) from a filename like CC10_89.txt."""
    base = os.path.basename(path)
    m = re.match(r"CC(\d+)_(\d+)\.txt$", base)
    if not m:
        return None, None
    return int(m.group(1)), int(m.group(2))


def read_roots(path):
    """Read hex roots from a file. Handles 0x prefix, bare hex, mixed case."""
    roots = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            # Take first column if space-separated
            token = line.split()[0]
            # Strip 0x/0X prefix
            if token.startswith(("0x", "0X")):
                token = token[2:]
            try:
                val = int(token, 16)
                roots.append(val)
            except ValueError:
                continue
    return roots


# ---------------------------------------------------------------------------
# Immunization depth (pure Python, no GP)
# ---------------------------------------------------------------------------

CRT_PRIMES = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31]


def immunization_depth(root, primes=CRT_PRIMES):
    """For a first-kind chain root p, compute how many chain elements
    survive each small prime q.  Chain element j = 2^j * p + (2^j - 1).
    Element j ≡ 0 mod q when 2^j * (p+1) ≡ 1 mod q.
    Returns dict {q: first_j_killed} or {q: None} if immune through 18."""
    result = {}
    for q in primes:
        p_mod = root % q
        # chain element j: 2^j * p + 2^j - 1 = 2^j*(p+1) - 1
        # killed at j if (2^j*(p+1) - 1) % q == 0
        # i.e. 2^j*(p+1) ≡ 1 (mod q)
        pp1 = (p_mod + 1) % q
        if pp1 == 0:
            # p+1 ≡ 0 mod q → 2^j*(p+1) ≡ 0 mod q → never ≡ 1
            # but element 0 = p itself: if p%q==q-1 then p+1≡0, element 0 = p ≡ q-1 ≠ 0
            result[q] = None  # immune to this prime
            continue
        killed = None
        pow2 = 1
        for j in range(19):  # check elements 0..18
            val = (pow2 * pp1 - 1) % q
            if val == 0:
                killed = j
                break
            pow2 = (pow2 * 2) % q
        result[q] = killed
    return result


# ---------------------------------------------------------------------------
# Decay curve fitting (least squares on log(count))
# ---------------------------------------------------------------------------

def fit_exponential(xs, ys):
    """Fit log(y) = a*x + b via least squares. Returns (slope, intercept, r_squared)."""
    n = len(xs)
    if n < 2:
        return 0, 0, 0
    log_ys = [math.log(y) for y in ys if y > 0]
    xs_filt = [x for x, y in zip(xs, ys) if y > 0]
    n = len(xs_filt)
    if n < 2:
        return 0, 0, 0

    sx = sum(xs_filt)
    sy = sum(log_ys)
    sxx = sum(x * x for x in xs_filt)
    sxy = sum(x * y for x, y in zip(xs_filt, log_ys))
    syy = sum(y * y for y in log_ys)

    denom = n * sxx - sx * sx
    if abs(denom) < 1e-15:
        return 0, 0, 0

    a = (n * sxy - sx * sy) / denom
    b = (sy - a * sx) / n

    # R²
    ss_res = sum((ly - (a * x + b)) ** 2 for x, ly in zip(xs_filt, log_ys))
    mean_ly = sy / n
    ss_tot = sum((ly - mean_ly) ** 2 for ly in log_ys)
    r_sq = 1 - ss_res / ss_tot if ss_tot > 0 else 0

    return a, b, r_sq


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def build_dashboard_data(data_dir):
    """Build all dashboard data structures from CC*_*.txt files."""

    # Discover files
    pattern = os.path.join(data_dir, "CC*_*.txt")
    files = sorted(glob.glob(pattern))
    if not files:
        print(f"ERROR: No CC*_*.txt files found in {data_dir}", file=sys.stderr)
        sys.exit(1)

    # ---------- Pass 1: read all files, build index ----------
    # all_data[cc][bits] = list of roots
    all_data = defaultdict(lambda: defaultdict(list))
    file_count = 0

    for path in files:
        cc, bits = parse_filename(path)
        if cc is None:
            continue
        roots = read_roots(path)
        if roots:
            all_data[cc][bits].extend(roots)
            file_count += 1

    print(f"Read {file_count} files from {data_dir}")

    # ---------- Summary ----------
    total_roots = 0
    cc_counts = {}
    all_cc = sorted(all_data.keys())
    all_bits_set = set()

    for cc in all_cc:
        count = sum(len(r) for r in all_data[cc].values())
        cc_counts[cc] = count
        total_roots += count
        all_bits_set.update(all_data[cc].keys())

    summary = {
        "total_roots": total_roots,
        "total_files": file_count,
        "cc_range": [min(all_cc), max(all_cc)],
        "bits_range": [min(all_bits_set), max(all_bits_set)],
        "cc_counts": {str(k): v for k, v in sorted(cc_counts.items())},
    }
    print(f"Total roots: {total_roots:,}")

    # ---------- Heatmap ----------
    heatmap = []
    for cc in all_cc:
        for bits in sorted(all_data[cc].keys()):
            count = len(all_data[cc][bits])
            heatmap.append({"cc": cc, "bits": bits, "count": count})

    # ---------- Decay curve (89–91 bit band) ----------
    # Collect all CC lengths with data in the 89-91 band
    decay_band_all = {}
    for cc in all_cc:
        band_count = 0
        for bits in [89, 90, 91]:
            band_count += len(all_data[cc].get(bits, []))
        if band_count > 0:
            decay_band_all[cc] = band_count

    decay_cc_all = sorted(decay_band_all.keys())
    decay_counts_all = [decay_band_all[c] for c in decay_cc_all]

    # Fit only on CC12–CC16 (monotonically decreasing, same search depth)
    # CC7-CC11 have variable search intensity creating non-monotonic bumps
    fit_cc = [c for c in decay_cc_all if 12 <= c <= 16]
    fit_counts = [decay_band_all[c] for c in fit_cc]

    slope, intercept, r_sq = fit_exponential(fit_cc, fit_counts)

    # Extrapolate CC17, CC18
    cc17_est = math.exp(slope * 17 + intercept) if slope != 0 else 0
    cc18_est = math.exp(slope * 18 + intercept) if slope != 0 else 0
    ratio = math.exp(-slope) if slope != 0 else 0  # count(n)/count(n+1)

    decay_89_91 = {
        "cc_lengths": decay_cc_all,
        "counts": decay_counts_all,
        "fit_cc": fit_cc,
        "fit_counts": fit_counts,
        "slope": round(slope, 6),
        "intercept": round(intercept, 4),
        "r_squared": round(r_sq, 6),
        "ratio_per_level": round(ratio, 2),
        "cc17_estimate": round(cc17_est, 3),
        "cc18_estimate": round(cc18_est, 4),
    }
    print(f"Decay fit (CC12-CC16): ratio={ratio:.2f}:1 per CC level, R²={r_sq:.4f}")
    print(f"  CC17 estimate: {cc17_est:.1f}, CC18 estimate: {cc18_est:.3f}")

    # ---------- Prefix analysis (densest CC/bits combo) ----------
    # Find densest bucket
    densest_cc, densest_bits, densest_count = 0, 0, 0
    for cc in all_cc:
        for bits in all_data[cc]:
            c = len(all_data[cc][bits])
            if c > densest_count:
                densest_cc, densest_bits, densest_count = cc, bits, c

    prefix_analysis = _build_prefix_analysis(
        all_data[densest_cc][densest_bits], densest_cc, densest_bits
    )

    # Also build for CC10/90 if different from densest
    if (densest_cc, densest_bits) != (10, 90) and len(all_data.get(10, {}).get(90, [])) > 100:
        prefix_analysis["cc10_90"] = _build_prefix_analysis(
            all_data[10][90], 10, 90
        )

    # ---------- CC14+ crown jewels ----------
    cc14_plus = []
    for cc in sorted(all_data.keys()):
        if cc < 14:
            continue
        for bits in sorted(all_data[cc].keys()):
            for root in all_data[cc][bits]:
                imm = immunization_depth(root)
                hex_str = hex(root)
                actual_bits = root.bit_length()
                cc14_plus.append({
                    "cc": cc,
                    "bits": actual_bits,
                    "file_bits": bits,
                    "hex": hex_str,
                    "immunization": {str(q): v for q, v in imm.items()},
                    "min_kill": min((v for v in imm.values() if v is not None), default=None),
                })

    # Sort: longest chain first, then by bits descending
    cc14_plus.sort(key=lambda x: (-x["cc"], -x["bits"]))
    print(f"CC14+ roots: {len(cc14_plus)}")

    # ---------- High-bit experimental (>128 bits) ----------
    highbit = []
    for cc in all_cc:
        for bits in sorted(all_data[cc].keys()):
            if bits <= 128:
                continue
            count = len(all_data[cc][bits])
            highbit.append({
                "cc": cc,
                "bits": bits,
                "count": count,
                "sample_hex": hex(all_data[cc][bits][0]) if all_data[cc][bits] else None,
            })

    highbit.sort(key=lambda x: (-x["bits"], -x["cc"]))
    print(f"High-bit entries (>128 bit): {len(highbit)}")

    return {
        "summary": summary,
        "heatmap": heatmap,
        "decay_89_91": decay_89_91,
        "prefix_analysis": prefix_analysis,
        "cc14_plus": cc14_plus,
        "highbit": highbit,
    }


def _build_prefix_analysis(roots, cc, bits):
    """Build prefix histogram and Poisson stats for a set of roots."""
    # Use top 20 bits as prefix bucket
    prefix_bits = 20
    prefix_counts = Counter()
    for r in roots:
        bl = r.bit_length()
        if bl > prefix_bits:
            prefix = r >> (bl - prefix_bits)
        else:
            prefix = r
        prefix_counts[prefix] += 1

    # Hit distribution: how many prefixes have k hits
    hit_dist = Counter(prefix_counts.values())
    total_prefixes = len(prefix_counts)
    total_hits = sum(prefix_counts.values())
    lam = total_hits / (1 << prefix_bits) if prefix_bits <= 20 else total_hits / total_prefixes

    # Top hotspots
    top = prefix_counts.most_common(10)
    hotspots = [{"prefix": hex(p), "count": c} for p, c in top]

    return {
        "cc": cc,
        "bits": bits,
        "total_roots": len(roots),
        "prefix_bits": prefix_bits,
        "unique_prefixes": total_prefixes,
        "poisson_lambda": round(lam, 4),
        "hit_distribution": {str(k): v for k, v in sorted(hit_dist.items())},
        "hotspots": hotspots,
    }


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="CC Campaign Dashboard Data Pipeline")
    parser.add_argument("data_dir", help="Directory containing CC*_*.txt files")
    parser.add_argument("--output", "-o", default="cc_dashboard_data.json",
                        help="Output JSON file (default: cc_dashboard_data.json)")
    args = parser.parse_args()

    if not os.path.isdir(args.data_dir):
        print(f"ERROR: {args.data_dir} is not a directory", file=sys.stderr)
        sys.exit(1)

    data = build_dashboard_data(args.data_dir)

    with open(args.output, "w") as f:
        json.dump(data, f, indent=1)

    size_kb = os.path.getsize(args.output) / 1024
    print(f"\nWrote {args.output} ({size_kb:.0f} KB)")


if __name__ == "__main__":
    main()
