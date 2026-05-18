#!/usr/bin/env python3
"""Verify fingerprint distribution against the CC10+ snapshot CSV.

Reads data/cc10plus_roots_snapshot_2026-03-19.csv (cc, root_hex, digits, bits)
and computes the immune fingerprint of each root over a fixed prime set.
A root p is "immune" to prime q iff (p+1) % q == 0.

When run with no args, also emits
  pools/known_fingerprints.txt — a q-iter --seed-pool-file ready list of
  the top-N most-frequent immune fingerprints, force-augmented with the
  default pinned sieve set {2,3,5}. Use --pin-7 only for legacy
  --immune-prime 7 experiments.

Usage:
  python3 scripts/verify_fingerprints.py
  python3 scripts/verify_fingerprints.py --min-cc 14
"""
import argparse
import collections
import csv
import math
import os
import sys
from datetime import date

PRIMES = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31]

# Extended palette for seed-pool-file emission. Includes 37 to cover known
# CC18 fingerprints, plus a few small primes beyond that.
POOL_PALETTE = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61]

# Default sieve floor. q=7 is intentionally not pinned: only 24.16% of the
# CC10+ snapshot is 7-immune, so pinning 7 discards most known structures.
PINNED_SIEVE_PRIMES_DEFAULT = (2, 3, 5)
PINNED_SIEVE_PRIMES_WITH_7 = (2, 3, 5, 7)

POOL_OUT_PATH = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "pools", "known_fingerprints.txt",
)
POOL_TOP_N = 100
POOL_BUCKET_BITS = 20.0   # design §3.3 coherence cap


def fingerprint(root_hex: str) -> frozenset[int]:
    p = int(root_hex, 16)
    p1 = p + 1
    return frozenset(q for q in PRIMES if p1 % q == 0)


def fingerprint_pool(root_hex: str) -> frozenset[int]:
    """Wider palette used for the seed-pool emission."""
    p = int(root_hex, 16)
    p1 = p + 1
    return frozenset(q for q in POOL_PALETTE if p1 % q == 0)


def _dv_bits(primes):
    return sum(math.log2(q) for q in primes)


def _fmt_pool_line(primes, freq):
    spec = "*".join(str(q) for q in primes)
    return f"{spec}  # freq={freq} bits={_dv_bits(primes):.2f}"


def emit_seed_pool_file(csv_path: str, out_path: str, top_n: int,
                        min_cc: int, max_cc: int,
                        bucket_bits: float,
                        pinned_primes: tuple[int, ...]) -> None:
    """Emit a q-iter --seed-pool-file from the snapshot CSV.

    Mirrors scripts/build_known_fingerprints.py.
    """
    pool_counts = collections.Counter()
    pool_total = 0
    with open(csv_path) as f:
        reader = csv.reader(f)
        next(reader)  # header
        for row in reader:
            if not row:
                continue
            cc = int(row[0])
            if cc < min_cc or cc > max_cc:
                continue
            fp = fingerprint_pool(row[1])
            if not fp:
                continue   # no immunity at all — skip
            pool_counts[fp] += 1
            pool_total += 1

    pinned = frozenset(pinned_primes)
    augmented: dict[frozenset, int] = {}
    for fp, c in pool_counts.items():
        v_set = frozenset(fp | pinned)
        augmented[v_set] = augmented.get(v_set, 0) + c

    top = sorted(augmented.items(), key=lambda kv: (-kv[1], sorted(kv[0])))
    top = top[:top_n]

    primesets = [sorted(s) for s, _ in top]
    bits = [_dv_bits(ps) for ps in primesets]
    bit_min, bit_max = min(bits), max(bits)
    span = bit_max - bit_min

    today = date.today().isoformat()
    src_rel = os.path.relpath(csv_path, os.path.dirname(os.path.abspath(out_path)))
    header_lines = [
        f"# generated {today} from {src_rel}",
        f"# via scripts/verify_fingerprints.py "
        f"(first-kind predicate: q in fp iff q|(p+1))",
        f"# top-{len(top)} most-frequent immune fingerprints, "
        f"CC{min_cc}-CC{max_cc} band ({pool_total} roots)",
        f"# augmented with pinned sieve set {set(pinned_primes)} "
        f"(--immune-prime {max(pinned_primes)})",
        f"# format: P1*P2*...*Pk  # freq=NNN bits=BB.B  (one D_V per line)",
        f"# D_V bit-range: min={bit_min:.2f} max={bit_max:.2f} "
        f"span={span:.2f} bits",
    ]

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    written = []

    # Always emit the full top-N file at out_path.
    with open(out_path, "w") as fh:
        for line in header_lines:
            fh.write(line + "\n")
        if span > bucket_bits:
            fh.write(f"# WARN: D_V span {span:.2f} > {bucket_bits} bit "
                     f"coherence limit (design §3.3); see bucketed files\n")
        for primes_fs, freq in top:
            fh.write(_fmt_pool_line(sorted(primes_fs), freq) + "\n")
    written.append(out_path)

    # Bucket by D_V bit band when the full span exceeds the coherence limit.
    if span > bucket_bits:
        buckets: dict[int, list] = {}
        for (primes_fs, freq), b in zip(top, bits):
            key = int(b // bucket_bits)
            buckets.setdefault(key, []).append((sorted(primes_fs), freq, b))
        base, ext = os.path.splitext(out_path)
        for key in sorted(buckets):
            lo = int(key * bucket_bits)
            hi = int((key + 1) * bucket_bits)
            bp = f"{base}_{lo}-{hi}bit{ext}"
            entries = buckets[key]
            with open(bp, "w") as fh:
                for line in header_lines:
                    fh.write(line + "\n")
                fh.write(f"# bucket: D_V in [{lo}, {hi}) bits "
                        f"({len(entries)} entries)\n")
                for primes, freq, _b in entries:
                    fh.write(_fmt_pool_line(primes, freq) + "\n")
            written.append(bp)

    print("\nseed-pool-file emission:", file=sys.stderr)
    for p in written:
        n = sum(1 for L in open(p) if L.strip() and not L.startswith("#"))
        print(f"  wrote {p}: {n} pool entries", file=sys.stderr)
    print(f"  D_V bits: min={bit_min:.2f} max={bit_max:.2f} "
          f"span={span:.2f}", file=sys.stderr)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default="data/cc10plus_roots_snapshot_2026-03-19.csv")
    ap.add_argument("--min-cc", type=int, default=10)
    ap.add_argument("--max-cc", type=int, default=20)
    ap.add_argument("--top", type=int, default=12)
    ap.add_argument("--pin-7", action="store_true",
                    help="force-augment every D_V with 7 "
                         "(legacy --immune-prime 7 mode)")
    args = ap.parse_args()

    pinned_primes = (PINNED_SIEVE_PRIMES_WITH_7
                     if args.pin_7 else PINNED_SIEVE_PRIMES_DEFAULT)

    if not os.path.isfile(args.csv):
        print(f"missing {args.csv}", file=sys.stderr)
        return 2

    counts = collections.Counter()
    total = 0
    cc_totals = collections.Counter()
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
            fp = fingerprint(row[1])
            counts[fp] += 1
            cc_totals[cc] += 1
            total += 1

    if total == 0:
        print("no rows matched", file=sys.stderr)
        return 1

    print(f"total roots: {total} (CC range [{args.min_cc}..{args.max_cc}])")
    print(f"unique fingerprints: {len(counts)}")
    print(f"per-CC counts:")
    for cc in sorted(cc_totals):
        print(f"  CC{cc:>2}: {cc_totals[cc]}")
    print()
    print(f"top {args.top} fingerprints:")
    for fp, cnt in counts.most_common(args.top):
        members = ",".join(str(q) for q in sorted(fp)) if fp else "(empty)"
        pct = 100.0 * cnt / total
        flags = " <- must-keep base" if fp == frozenset({3, 5, 11, 13, 19}) else ""
        print(f"  {{{members}}}  {cnt:>8d}  {pct:6.2f}%{flags}")

    base = frozenset({3, 5, 11, 13, 19})
    base_count = counts.get(base, 0)
    print()
    print(f"fingerprint exactly {{3,5,11,13,19}}: {base_count} ({100.0 * base_count / total:.2f}%)")
    base_super = sum(c for fp, c in counts.items() if base <= fp)
    print(f"fingerprint ⊇ {{3,5,11,13,19}}: {base_super} ({100.0 * base_super / total:.2f}%)")
    missing_base = total - base_super
    print(f"fingerprint missing some of {{3,5,11,13,19}}: {missing_base} ({100.0 * missing_base / total:.2f}%)")

    # Also emit pools/known_fingerprints.txt as a side effect of the canonical
    # no-arg invocation. Re-scans the CSV with the wider palette; cheap
    # relative to the main verify pass.
    emit_seed_pool_file(
        csv_path=args.csv,
        out_path=POOL_OUT_PATH,
        top_n=POOL_TOP_N,
        min_cc=args.min_cc,
        max_cc=args.max_cc,
        bucket_bits=POOL_BUCKET_BITS,
        pinned_primes=pinned_primes,
    )

    return 0


if __name__ == "__main__":
    sys.exit(main())
