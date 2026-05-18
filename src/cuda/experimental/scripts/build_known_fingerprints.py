#!/usr/bin/env python3
# Copyright (c) 2026 Nenad Micic
# SPDX-License-Identifier: Apache-2.0
"""build_known_fingerprints.py — generate a q-iter seed-pool file from the
historical CC10+ root snapshot.

For each root in --csv:
  - decode root_hex -> int p (the chain root, first-kind)
  - immune fingerprint = { q in PRIME_PALETTE : (p+1) % q == 0 }
Then:
  - count fingerprint frequencies across the whole snapshot (or a CC band)
  - take top-N
  - force-include the pinned-sieve must-keep set {2,3,5} in every variant
    so each emitted D_V is a superset of the default --immune-prime 5 sieve
    floor. Use --pin-7 only for legacy 7-immune experiments.
  - emit '*'-separated lines, frequency-sorted descending, with annotation:
       2*3*5*11*13*19  # freq=NNN bits=BB.B
  - check D_V bit-range coherence. If span exceeds 20 bits, emit bucketed
    files keyed by D_V bit-band (e.g. _70-90bit.txt).

Output: ./pools/known_fingerprints.txt
(plus optional bucket variants alongside it).
"""

import argparse
import collections
import csv
import math
import os
import sys
from datetime import date


# Prime palette used for the immunization fingerprint of a first-kind root.
# Includes 37 to cover known CC18 fingerprints and a few small primes beyond.
PRIME_PALETTE = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61]

# Default sieve floor. q=7 is intentionally not pinned: only 24.16% of the
# CC10+ snapshot is 7-immune, so pinning 7 discards most known structures.
PINNED_SIEVE_PRIMES_DEFAULT = (2, 3, 5)
PINNED_SIEVE_PRIMES_WITH_7 = (2, 3, 5, 7)


def fingerprint_root(p_hex: str) -> frozenset[int]:
    """Immune fingerprint of a first-kind chain root.

    For first-kind chains, q immunizes the root when q | (p+1).
    """
    p = int(p_hex, 16)
    p1 = p + 1
    return frozenset(q for q in PRIME_PALETTE if p1 % q == 0)


def dv_bits(primes: list[int]) -> float:
    """log2(prod primes) — useful for D_V bit-range coherence checks."""
    total = 0.0
    for q in primes:
        total += math.log2(q)
    return total


def fmt_line(primes: list[int], freq: int) -> str:
    """Emit '2*3*5*7*11*13*19  # freq=NNN bits=BB.B'."""
    spec = "*".join(str(q) for q in primes)
    bits = dv_bits(primes)
    return f"{spec}  # freq={freq} bits={bits:.2f}"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv",
                    default="./data/"
                            "cc10plus_roots_snapshot_2026-03-19.csv",
                    help="historical CC10+ snapshot CSV")
    ap.add_argument("--out",
                    default="./pools/"
                            "known_fingerprints.txt",
                    help="output seed-pool file")
    ap.add_argument("--top-n", type=int, default=100,
                    help="emit top-N most-frequent fingerprints (default 100)")
    ap.add_argument("--min-cc", type=int, default=10,
                    help="only count roots with CC >= this (default 10)")
    ap.add_argument("--max-cc", type=int, default=99,
                    help="only count roots with CC <= this (default 99)")
    ap.add_argument("--bucket-bits", type=float, default=20.0,
                    help="if D_V bit span > this, emit bucketed files "
                         "(design doc §3.3, default 20)")
    ap.add_argument("--bucket-width", type=float, default=20.0,
                    help="width of each bucket in bits (default 20)")
    ap.add_argument("--pin-7", action="store_true",
                    help="force-augment every D_V with 7 "
                         "(legacy --immune-prime 7 mode)")
    args = ap.parse_args()

    pinned_primes = (PINNED_SIEVE_PRIMES_WITH_7
                     if args.pin_7 else PINNED_SIEVE_PRIMES_DEFAULT)

    if not os.path.isfile(args.csv):
        print(f"missing {args.csv}", file=sys.stderr)
        return 2

    counts: collections.Counter[frozenset[int]] = collections.Counter()
    total = 0
    cc_totals: collections.Counter[int] = collections.Counter()

    with open(args.csv) as f:
        reader = csv.reader(f)
        header = next(reader)
        if header[:2] != ["cc", "root_hex"]:
            print(f"unexpected header: {header}", file=sys.stderr)
            return 2
        for row in reader:
            if not row:
                continue
            cc = int(row[0])
            if cc < args.min_cc or cc > args.max_cc:
                continue
            fp = fingerprint_root(row[1])
            counts[fp] += 1
            cc_totals[cc] += 1
            total += 1

    if total == 0:
        print("no rows matched CC range", file=sys.stderr)
        return 1

    print(f"# scanned {total} CC{args.min_cc}+ roots; "
          f"{len(counts)} distinct fingerprints", file=sys.stderr)
    for cc in sorted(cc_totals):
        print(f"#   CC{cc:>2}: {cc_totals[cc]}", file=sys.stderr)

    # Build top-N: drop the empty fingerprint (no immunization → uninteresting
    # for q-iter), force-add the pinned sieve set, dedup, take top-N by freq.
    pinned = set(pinned_primes)
    augmented: dict[frozenset[int], int] = {}
    for fp, c in counts.items():
        if not fp:
            continue
        v_set = frozenset(fp | pinned)
        augmented[v_set] = augmented.get(v_set, 0) + c

    top = sorted(augmented.items(), key=lambda kv: (-kv[1], sorted(kv[0])))
    top = top[: args.top_n]

    # Bit-range coherence check.
    primesets = [sorted(s) for s, _ in top]
    bits = [dv_bits(ps) for ps in primesets]
    bit_min, bit_max = min(bits), max(bits)
    span = bit_max - bit_min
    print(f"# top-{len(top)} D_V bits: min={bit_min:.2f} max={bit_max:.2f} "
          f"span={span:.2f}", file=sys.stderr)

    src_csv = os.path.relpath(args.csv,
                              os.path.dirname(os.path.abspath(args.out)))
    today = date.today().isoformat()
    header_lines = [
        f"# generated {today} from {src_csv}",
        f"# via scripts/build_known_fingerprints.py "
        f"(first-kind predicate: q in fp iff q|(p+1))",
        f"# top-{len(top)} most-frequent immune fingerprints, "
        f"CC{args.min_cc}+ band ({total} roots)",
        f"# augmented with pinned sieve set {set(pinned_primes)} "
        f"(--immune-prime {max(pinned_primes)})",
        f"# format: P1*P2*...*Pk  # freq=NNN bits=BB.B  (one D_V per line)",
        f"# D_V bit-range: min={bit_min:.2f} max={bit_max:.2f} "
        f"span={span:.2f} bits",
    ]

    files_written = []

    if span <= args.bucket_bits:
        with open(args.out, "w") as fh:
            for line in header_lines:
                fh.write(line + "\n")
            for primes, freq in top:
                fh.write(fmt_line(sorted(primes), freq) + "\n")
        files_written.append(args.out)
    else:
        # Bucket by floor(bits / width).
        buckets: dict[int, list[tuple[list[int], int, float]]] = {}
        for (primes_fs, freq), b in zip(top, bits):
            key = int(b // args.bucket_width)
            buckets.setdefault(key, []).append((sorted(primes_fs), freq, b))
        base, ext = os.path.splitext(args.out)
        for key in sorted(buckets):
            lo = int(key * args.bucket_width)
            hi = int((key + 1) * args.bucket_width)
            bucket_file = f"{base}_{lo}-{hi}bit{ext}"
            bucket_entries = buckets[key]
            with open(bucket_file, "w") as fh:
                for line in header_lines:
                    fh.write(line + "\n")
                fh.write(f"# bucket: D_V in [{lo}, {hi}) bits "
                        f"({len(bucket_entries)} entries)\n")
                for primes, freq, _b in bucket_entries:
                    fh.write(fmt_line(primes, freq) + "\n")
            files_written.append(bucket_file)
        # Also emit the un-bucketed file (so default --seed-pool-file path
        # is still useful for callers that want the full top-N).
        with open(args.out, "w") as fh:
            for line in header_lines:
                fh.write(line + "\n")
            fh.write(f"# UN-BUCKETED — D_V span {span:.2f} bits > "
                    f"{args.bucket_bits} bit coherence limit\n")
            fh.write(f"# see bucket files alongside this one for per-launch "
                    f"compatible subsets\n")
            for primes_fs, freq in top:
                fh.write(fmt_line(sorted(primes_fs), freq) + "\n")
        files_written.append(args.out)

    for fp in files_written:
        n_lines = sum(1 for L in open(fp) if L.strip() and not L.startswith("#"))
        print(f"wrote {fp}: {n_lines} pool entries", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
