#!/usr/bin/env python3
# Copyright (c) 2026 Nenad Micic
# SPDX-License-Identifier: Apache-2.0
"""build_known_fingerprints_all.py — emit ALL unique immune fingerprints
from the CC10+ snapshot.

Sibling of scripts/build_known_fingerprints.py / verify_fingerprints.py;
shares predicate (q in fp iff q|(p+1)), palette
{3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61}, and force-augment with
the default pinned sieve set {2,3,5}. Use --pin-7 only for legacy
--immune-prime 7 experiments.

UNLIKE the top-N scripts: emits *every* unique D_V (after augment),
frequency-sorted descending, to a separate output file so the existing
top-100 pools/known_fingerprints.txt is left untouched.

Also writes a frequency distribution, a D_V bit-range histogram, and a
device-/host-side memory budget summary to stderr.
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
# q=13 are 100%. Pinning {2,3,5,7} threw away ~76% of chain structures.
# Default floor is {2,3,5} (matches v16 --immune-prime 5). --pin-7 reproduces
# the old 7-immune-only behavior for parity / focused campaigns.
PINNED_SIEVE_PRIMES_DEFAULT = (2, 3, 5)
PINNED_SIEVE_PRIMES_WITH_7 = (2, 3, 5, 7)


def fingerprint_root(p_hex: str) -> frozenset:
    p = int(p_hex, 16)
    p1 = p + 1
    return frozenset(q for q in PRIME_PALETTE if p1 % q == 0)


def dv_bits(primes) -> float:
    return sum(math.log2(q) for q in primes)


def fmt_line(primes, freq: int) -> str:
    spec = "*".join(str(q) for q in primes)
    return f"{spec}  # freq={freq} bits={dv_bits(primes):.2f}"


def freq_bucket_label(freq: int) -> str:
    if freq == 1:
        return "       1x"
    if freq < 10:
        return "    2-9x"
    if freq < 100:
        return "  10-99x"
    if freq < 1_000:
        return "100-999x"
    if freq < 10_000:
        return " 1K-9.9Kx"
    if freq < 100_000:
        return " 10K-99Kx"
    return "    100K+x"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv",
                    default="./data/"
                            "cc10plus_roots_snapshot_2026-03-19.csv")
    ap.add_argument("--out",
                    default="./pools/"
                            "known_fingerprints_all.txt")
    ap.add_argument("--min-cc", type=int, default=10)
    ap.add_argument("--max-cc", type=int, default=99)
    # v16 config knobs (for memory budget report only)
    ap.add_argument("--depth", type=int, default=16,
                    help="kernel depth used in the memory budget report")
    ap.add_argument("--n-sp", type=int, default=43,
                    help="sieve-prime count (--sieve-max 211 -> 43)")
    ap.add_argument("--device-vram-gb", type=float, default=32.0,
                    help="device VRAM cap for budget check (RTX 5090: 32)")
    ap.add_argument("--pin-7", action="store_true",
                    help="force-augment every D_V with 7 (legacy --immune-prime 7 mode); "
                         "default omits 7 so 7-excluded chains stay distinct "
                         "(76%% of CC10+ raw fingerprints are 7-excluded)")
    args = ap.parse_args()

    pinned_primes = (PINNED_SIEVE_PRIMES_WITH_7
                     if args.pin_7 else PINNED_SIEVE_PRIMES_DEFAULT)

    if not os.path.isfile(args.csv):
        print(f"missing {args.csv}", file=sys.stderr)
        return 2

    counts: collections.Counter = collections.Counter()
    total = 0
    cc_totals: collections.Counter = collections.Counter()
    skipped_empty = 0

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
            if not fp:
                skipped_empty += 1
                # Even an empty fingerprint becomes pinned-only after
                # force-augment; keep it so the sieve-floor variant is
                # represented in the all-uniques set.
            counts[fp] += 1
            cc_totals[cc] += 1
            total += 1

    if total == 0:
        print("no rows matched", file=sys.stderr)
        return 1

    pinned = frozenset(pinned_primes)
    augmented: dict = {}
    for fp, c in counts.items():
        v_set = frozenset(fp | pinned)
        augmented[v_set] = augmented.get(v_set, 0) + c

    n_v = len(augmented)
    print(f"# scanned {total} CC{args.min_cc}-{args.max_cc} roots",
          file=sys.stderr)
    print(f"# {len(counts)} raw distinct palette fingerprints "
          f"(of which {skipped_empty} rows had empty fingerprint, mapped "
          f"to pinned-only V after augment)", file=sys.stderr)
    print(f"# {n_v} unique D_V variants after pinned-sieve force-augment "
          f"({set(pinned_primes)})", file=sys.stderr)

    sorted_v = sorted(augmented.items(),
                      key=lambda kv: (-kv[1], sorted(kv[0])))

    # Frequency distribution
    freq_buckets: collections.Counter = collections.Counter()
    for _v, freq in sorted_v:
        freq_buckets[freq_bucket_label(freq)] += 1

    # D_V bit-range histogram (1-bit buckets)
    bits_per_v = [dv_bits(sorted(v)) for v, _ in sorted_v]
    bit_min, bit_max = min(bits_per_v), max(bits_per_v)
    bit_buckets: collections.Counter = collections.Counter()
    for b in bits_per_v:
        bit_buckets[int(math.floor(b))] += 1

    # Memory budget
    per_v_device = args.depth * args.n_sp * 4  # bytes
    total_device_bytes = n_v * per_v_device
    per_v_host = 600  # bytes (v16_pool_seed approx, design doc §10.2)
    total_host_bytes = n_v * per_v_host
    device_mb = total_device_bytes / (1024 * 1024)
    host_mb = total_host_bytes / (1024 * 1024)
    vram_cap_bytes = args.device_vram_gb * (1024 ** 3)

    today = date.today().isoformat()
    src_rel = os.path.relpath(args.csv,
                              os.path.dirname(os.path.abspath(args.out)))
    header_lines = [
        f"# generated {today} from {src_rel}",
        f"# via scripts/build_known_fingerprints_all.py "
        f"(predicate: q in fp iff q|(p+1); palette {PRIME_PALETTE})",
        f"# ALL unique immune fingerprints in CC{args.min_cc}-{args.max_cc} "
        f"band ({total} roots scanned, {skipped_empty} with empty raw fp)",
        f"# augmented with pinned sieve set {set(pinned_primes)} "
        f"(--immune-prime {max(pinned_primes)}) — v16 requires D_V superset of D_sieve",
        f"# format: P1*P2*...*Pk  # freq=NNN bits=BB.B  (one D_V per line)",
        f"# total unique D_V: {n_v}",
        f"# D_V bit-range: min={bit_min:.2f} max={bit_max:.2f} "
        f"span={(bit_max - bit_min):.2f} bits",
    ]

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as fh:
        for line in header_lines:
            fh.write(line + "\n")
        for v_fs, freq in sorted_v:
            fh.write(fmt_line(sorted(v_fs), freq) + "\n")

    n_lines = sum(1 for L in open(args.out) if L.strip() and not L.startswith("#"))
    print(f"\nwrote {args.out}: {n_lines} pool entries\n", file=sys.stderr)

    # --- Report: freq distribution ---
    print("frequency distribution (count of unique D_V per bucket):",
          file=sys.stderr)
    bucket_order = ["    100K+x", " 10K-99Kx", " 1K-9.9Kx",
                    "100-999x", "  10-99x", "    2-9x", "       1x"]
    for label in bucket_order:
        n = freq_buckets.get(label, 0)
        print(f"  {label.strip():>10s}  {n:>7d}", file=sys.stderr)

    # --- Report: D_V bit-range histogram ---
    print("\nD_V bit-range histogram (1-bit buckets):", file=sys.stderr)
    print(f"  min={bit_min:.2f}  max={bit_max:.2f}  "
          f"span={bit_max-bit_min:.2f}", file=sys.stderr)
    max_count = max(bit_buckets.values()) if bit_buckets else 1
    for b in range(int(math.floor(bit_min)), int(math.floor(bit_max)) + 1):
        n = bit_buckets.get(b, 0)
        bar_len = int(40 * n / max_count) if max_count else 0
        print(f"  [{b:>2d},{b+1:>2d})  {n:>6d}  {'#' * bar_len}",
              file=sys.stderr)

    # --- Report: memory budget ---
    print("\nmemory budget (N_V * per-V):", file=sys.stderr)
    print(f"  N_V = {n_v}", file=sys.stderr)
    print(f"  device forbid_Q per V = depth({args.depth}) * "
          f"n_sp({args.n_sp}) * 4 B = {per_v_device} B", file=sys.stderr)
    print(f"  device total = {total_device_bytes} B = {device_mb:.2f} MB "
          f"(VRAM cap {args.device_vram_gb} GB = "
          f"{vram_cap_bytes/1024/1024:.0f} MB, "
          f"{100.0*total_device_bytes/vram_cap_bytes:.4f}%)",
          file=sys.stderr)
    print(f"  host v16_pool_seed per V ~ {per_v_host} B", file=sys.stderr)
    print(f"  host total ~ {total_host_bytes} B = {host_mb:.2f} MB",
          file=sys.stderr)

    current_cap = 8192
    if n_v <= current_cap:
        rec_cap = current_cap
        print(f"  V16_MAX_POOL_SEEDS current={current_cap} -> sufficient "
              f"(N_V={n_v} <= {current_cap})", file=sys.stderr)
    else:
        # Round up to next power of two
        rec_cap = 1
        while rec_cap < n_v:
            rec_cap *= 2
        print(f"  V16_MAX_POOL_SEEDS current={current_cap} -> "
              f"INSUFFICIENT; recommend {rec_cap} (next pow2 >= N_V={n_v})",
              file=sys.stderr)
        print(f"  budget at new cap: device "
              f"{rec_cap*per_v_device/1024/1024:.2f} MB, host "
              f"{rec_cap*per_v_host/1024/1024:.2f} MB", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
