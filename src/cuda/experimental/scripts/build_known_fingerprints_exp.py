#!/usr/bin/env python3
# Copyright (c) 2026 Nenad Micic
# SPDX-License-Identifier: Apache-2.0
"""build_known_fingerprints_exp.py — exponent-aware fingerprints.

Capture q^e detail, not just the squarefree prime-support set. Two
chains with the same prime-support but different exponents in (p+1) are
structurally different sieves. Variant D = freq >= 10 (drops the
22K-singleton long tail).

Output format: `q^e * q^e * ...` notation, e.g., `2*3^2*5*7*11*13*19`.
The engine's seed-pool parser accepts `q^e` syntax.

Sibling of build_known_fingerprints_all.py; same palette, default
{2,3,5} force-augment (each at exponent ≥ 1). Pass --pin-7 to also
force 7 (legacy --immune-prime 7 mode).
"""

import argparse
import collections
import csv
import math
import os
import sys
from datetime import date

PRIME_PALETTE = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61]
# 2026-05-14: empirically q=7 is only 24.16% immune in CC10+; q=3, q=5, q=11,
# q=13 are 100%. Default floor is {2,3,5} so 7-excluded chains keep their
# own D_V instead of being aliased into a 7-immune bucket. --pin-7 restores
# the old behavior for legacy --immune-prime 7 campaigns.
PINNED_SIEVE_PRIMES_DEFAULT = (2, 3, 5)
PINNED_SIEVE_PRIMES_WITH_7 = (2, 3, 5, 7)


def factor_p1_over_palette(p_hex):
    p1 = int(p_hex, 16) + 1
    exps = {}
    for q in PRIME_PALETTE:
        e = 0
        while p1 % q == 0:
            p1 //= q
            e += 1
        if e > 0:
            exps[q] = e
    return exps


def augment(exps, pinned):
    """Force-augment with `pinned` primes (each at exponent ≥ 1)."""
    a = dict(exps)
    for q in pinned:
        a[q] = max(a.get(q, 0), 1)
    return tuple(sorted(a.items()))


def fmt_pool_line(t, freq):
    """Emit q^e syntax: 2*3^2*5*7*11*13*19."""
    parts = []
    for q, e in t:
        parts.append(f"{q}^{e}" if e > 1 else f"{q}")
    spec = "*".join(parts)
    bits = sum(e * math.log2(q) for q, e in t)
    return f"{spec}  # freq={freq} bits={bits:.2f}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv",
                    default="./data/"
                            "cc10plus_roots_snapshot_2026-03-19.csv")
    ap.add_argument("--out",
                    default="./pools/"
                            "known_fingerprints_exp_d.txt")
    ap.add_argument("--min-freq", type=int, default=10,
                    help="drop shapes with freq < N (variant D: 10)")
    ap.add_argument("--min-cc", type=int, default=10)
    ap.add_argument("--max-cc", type=int, default=99)
    ap.add_argument("--pin-7", action="store_true",
                    help="force-augment every D_V with 7 (legacy --immune-prime 7 mode); "
                         "default omits 7 so 7-excluded chains keep distinct D_V "
                         "(only 24%% of CC10+ raw fingerprints are 7-immune)")
    args = ap.parse_args()

    pinned_primes = (PINNED_SIEVE_PRIMES_WITH_7
                     if args.pin_7 else PINNED_SIEVE_PRIMES_DEFAULT)

    if not os.path.isfile(args.csv):
        print(f"missing {args.csv}", file=sys.stderr)
        return 2

    counts = collections.Counter()
    total = 0

    with open(args.csv) as f:
        rdr = csv.reader(f)
        next(rdr)
        for row in rdr:
            if not row:
                continue
            cc = int(row[0])
            if cc < args.min_cc or cc > args.max_cc:
                continue
            exps = factor_p1_over_palette(row[1])
            t = augment(exps, pinned_primes)
            counts[t] += 1
            total += 1

    # Filter by freq
    raw_n = len(counts)
    filtered = [(t, c) for t, c in counts.items() if c >= args.min_freq]
    filtered.sort(key=lambda x: (-x[1], x[0]))
    n_v = len(filtered)

    bits = [sum(e * math.log2(q) for q, e in t) for t, _ in filtered]
    bit_min, bit_max = min(bits), max(bits)

    # Cumulative coverage
    covered = sum(c for _, c in filtered)
    coverage_pct = 100.0 * covered / total

    today = date.today().isoformat()
    src_rel = os.path.relpath(args.csv,
                              os.path.dirname(os.path.abspath(args.out)))
    header_lines = [
        f"# generated {today} from {src_rel}",
        f"# via scripts/build_known_fingerprints_exp.py "
        f"(predicate: count exponent of q in (p+1) for q in {PRIME_PALETTE})",
        f"# variant D: freq >= {args.min_freq} (drops singleton/rare-pair tail)",
        f"# augmented with pinned sieve set {set(pinned_primes)} "
        f"(--immune-prime {max(pinned_primes)})",
        f"# format: P1[^e]*P2[^e]*... (^e omitted when e=1)",
        f"# total scanned: {total} CC{args.min_cc}-{args.max_cc} roots",
        f"# raw unique exponent shapes: {raw_n}",
        f"# after freq >= {args.min_freq} filter: {n_v} pool entries",
        f"# coverage of CC10+ data: {covered}/{total} = {coverage_pct:.2f}%",
        f"# D_V bit-range: min={bit_min:.2f} max={bit_max:.2f} "
        f"span={(bit_max - bit_min):.2f} bits",
    ]

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as fh:
        for line in header_lines:
            fh.write(line + "\n")
        for t, freq in filtered:
            fh.write(fmt_pool_line(t, freq) + "\n")

    print(f"wrote {args.out}: {n_v} pool entries", file=sys.stderr)
    print(f"  raw unique shapes:    {raw_n}", file=sys.stderr)
    print(f"  after freq >= {args.min_freq} filter: {n_v}", file=sys.stderr)
    print(f"  coverage of CC10+:    {coverage_pct:.2f}%", file=sys.stderr)
    print(f"  D_V bit-range:        {bit_min:.2f} – {bit_max:.2f} "
          f"(span {bit_max-bit_min:.2f})", file=sys.stderr)
    print(f"  memory at depth=20 n_sp=86: "
          f"device {n_v * 20 * 86 * 4 / 1024 / 1024:.1f} MB", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
