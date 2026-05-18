#!/usr/bin/env python3
# Copyright (c) 2026 Nenad Micic
# SPDX-License-Identifier: Apache-2.0
"""snapshot_replay.py — external prefix test.

For every CC10+ root p_0 in the snapshot CSV, reconstruct the q-iter
decomposition  p_0 + 1 = Q * D_V  where D_V is the palette-fingerprint
of p_0 with the {2,3,5} sieve floor force-augmented (matching the
build_known_fingerprints_exp.py default since 2026-05-14).

Then check whether D_V is present in the loaded pool file. A chain is
"structurally reachable" by the v16 q-iter engine iff:
  (a) D_V is in the pool
  (b) Q bit-length is inside the configured --bits-min / --bits-max band
  (c) the engine's --immune-prime N ≤ min(D_V_primes) so D_V ⊇ D_sieve

This is the python-side ultimatum test. Pair with a cuda single-Q verify
(once v16 grows --verify-q) for end-to-end coverage.
"""

import argparse
import collections
import csv
import math
import sys
from pathlib import Path

PRIME_PALETTE = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61]
PINNED_DEFAULT = (2, 3, 5)


def parse_pool(path):
    """Return dict[tuple[(q,e), ...]] -> (freq, bits) for one pool file."""
    pool = {}
    n_total = 0
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        n_total += 1
        spec = line.split("#", 1)[0].strip()
        # Parse q^e * q^e ... (^e omitted when e=1)
        shape = []
        for tok in spec.split("*"):
            tok = tok.strip()
            if "^" in tok:
                q, e = tok.split("^")
                shape.append((int(q), int(e)))
            else:
                shape.append((int(tok), 1))
        key = tuple(sorted(shape))
        pool[key] = pool.get(key, 0) + 1
    return pool, n_total


def factor_shape(p_plus_1, pinned):
    """Return tuple[(q,e), ...] for the palette factorization, with `pinned`
    augmented to exponent >= 1. p_plus_1 must NOT be divisible by anything
    outside the palette for D_V to be fully captured here; we use whatever
    palette factors we find (matches build_known_fingerprints_exp.py)."""
    n = p_plus_1
    exps = {}
    for q in PRIME_PALETTE:
        e = 0
        while n % q == 0:
            n //= q
            e += 1
        if e > 0:
            exps[q] = e
    for q in pinned:
        exps[q] = max(exps.get(q, 0), 1)
    return tuple(sorted(exps.items()))


def shape_to_int(shape):
    n = 1
    for q, e in shape:
        n *= q ** e
    return n


def fmt_shape(shape):
    return "*".join(f"{q}^{e}" if e > 1 else str(q) for q, e in shape)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--csv",
        default="./data/"
                "cc10plus_roots_snapshot_2026-03-19.csv",
    )
    ap.add_argument(
        "--pool",
        default="./pools/"
                "known_fingerprints_exp_d.txt",
    )
    ap.add_argument("--bits-min", type=int, default=90)
    ap.add_argument("--bits-max", type=int, default=127)
    ap.add_argument("--min-cc", type=int, default=10)
    ap.add_argument("--max-cc", type=int, default=99)
    ap.add_argument("--top-misses", type=int, default=15,
                    help="report top-N D_V shapes that are missing from the pool")
    ap.add_argument("--sample-hits", type=int, default=5,
                    help="emit N sample (cc, p_0, Q, D_V) reachable lines")
    ap.add_argument("--emit-prefixes", type=str, default=None,
                    help="write a TSV of (cc, p_0_hex, D_V, Q_hex) reachable "
                         "prefixes — feed to the engine for end-to-end verify")
    ap.add_argument("--simulate-sieve", action="store_true",
                    help="for each reachable chain, simulate the kernel sieve "
                         "(Q != D_V^-1 mod q for q in palette) — sanity check")
    args = ap.parse_args()

    pool, pool_total = parse_pool(args.pool)
    print(f"loaded pool: {pool_total} entries ({len(pool)} distinct shapes)",
          file=sys.stderr)

    counts = collections.Counter()
    cc_total = collections.Counter()
    cc_in_pool = collections.Counter()
    cc_in_pool_in_band = collections.Counter()
    miss_shapes = collections.Counter()
    sample_hits = []

    bits_min, bits_max = args.bits_min, args.bits_max
    # v16 Q-window formulae (cc16...v16.cu:4255-4257):
    #   Q_min(V) = ceil(2^(N1-1) / D_V)
    #   Q_max(V) = floor((2^N2 - 1) / D_V)
    lo = 1 << (bits_min - 1)
    hi = (1 << bits_max) - 1

    def q_window(DV):
        return (lo + DV - 1) // DV, hi // DV

    total = 0
    n_in_pool = 0
    n_in_pool_in_band = 0
    n_qwindow_ok = 0
    n_sieve_ok = 0

    prefix_fh = open(args.emit_prefixes, "w") if args.emit_prefixes else None
    if prefix_fh:
        prefix_fh.write("# external prefix test: reachable known CC10+ chains\n")
        prefix_fh.write("# columns: cc\tp_0_hex\tQ_hex\tD_V_spec\tp_0_bits\tQ_bits\n")

    with open(args.csv) as f:
        rdr = csv.reader(f)
        hdr = next(rdr)
        if hdr[:2] != ["cc", "root_hex"]:
            print(f"unexpected header: {hdr}", file=sys.stderr)
            return 2
        for row in rdr:
            if not row:
                continue
            cc = int(row[0])
            if cc < args.min_cc or cc > args.max_cc:
                continue
            p0 = int(row[1], 16)
            p1 = p0 + 1
            shape = factor_shape(p1, PINNED_DEFAULT)
            DV = shape_to_int(shape)
            if p1 % DV != 0:
                # Should be impossible by construction, guard anyway.
                continue
            Q = p1 // DV
            Q_bits = Q.bit_length()
            p0_bits = p0.bit_length()
            in_pool = shape in pool
            in_band = bits_min <= p0_bits <= bits_max

            total += 1
            counts[("any", "any")] += 1
            cc_total[cc] += 1
            if in_pool:
                n_in_pool += 1
                cc_in_pool[cc] += 1
                if in_band:
                    n_in_pool_in_band += 1
                    cc_in_pool_in_band[cc] += 1
                    # Engine Q-window check (kernel oracle).
                    q_min, q_max = q_window(DV)
                    if q_min <= Q <= q_max:
                        n_qwindow_ok += 1
                    # Synthetic sieve: a prime p_0 is never divisible by any
                    # palette q (q != p_0), so the engine's sieve survival is
                    # trivially true for any known prime. We assert it.
                    if args.simulate_sieve:
                        ok = all(p0 % q != 0 for q in PRIME_PALETTE)
                        if ok:
                            n_sieve_ok += 1
                    if prefix_fh:
                        prefix_fh.write(
                            f"CC{cc}\t{hex(p0)}\t{hex(Q)}\t{fmt_shape(shape)}\t"
                            f"{p0_bits}\t{Q_bits}\n"
                        )
                    if len(sample_hits) < args.sample_hits:
                        sample_hits.append(
                            (cc, hex(p0), hex(Q), Q_bits, p0_bits, fmt_shape(shape))
                        )
            else:
                miss_shapes[shape] += 1
    if prefix_fh:
        prefix_fh.close()

    # --- report ---
    print()
    print(f"=== Snapshot Replay vs pool ({args.pool.split('/')[-1]}) ===")
    print(f"scanned chains:               {total}")
    print(f"D_V matches pool:             {n_in_pool} "
          f"({100.0*n_in_pool/total:.2f}%)")
    print(f"D_V matches pool AND p_0 in "
          f"[{bits_min},{bits_max}] bits:   {n_in_pool_in_band} "
          f"({100.0*n_in_pool_in_band/total:.2f}%)")
    print(f"  └─ Q in engine [Q_min,Q_max] window:    {n_qwindow_ok} "
          f"({100.0*n_qwindow_ok/max(1,n_in_pool_in_band):.2f}% of reachable)")
    if args.simulate_sieve:
        print(f"  └─ Sieve oracle (p_0 mod q != 0):     {n_sieve_ok} "
              f"({100.0*n_sieve_ok/max(1,n_in_pool_in_band):.2f}% — should be 100%)")
    if args.emit_prefixes:
        print(f"\nWrote external prefix TSV → {args.emit_prefixes}")

    print("\nPer-CC structural coverage (D_V in pool):")
    print(f"  {'CC':>4} {'total':>8} {'in_pool':>10} {'in_band':>10} "
          f"{'pool%':>7} {'band%':>7}")
    for cc in sorted(cc_total):
        n = cc_total[cc]
        np_ = cc_in_pool[cc]
        nb = cc_in_pool_in_band[cc]
        print(f"  CC{cc:<3} {n:>8} {np_:>10} {nb:>10} "
              f"{100.0*np_/n:>6.2f}% {100.0*nb/n:>6.2f}%")

    print(f"\nTop {args.top_misses} missing D_V shapes (chains we still can't find):")
    for shape, n in miss_shapes.most_common(args.top_misses):
        bits = sum(e * math.log2(q) for q, e in shape)
        print(f"  freq={n:>7}  bits={bits:5.2f}  {fmt_shape(shape)}")

    if sample_hits:
        print(f"\nSample {len(sample_hits)} reachable "
              f"(cc, p_0, Q, Q_bits, p_0_bits, D_V):")
        for cc, p0h, Qh, qb, pb, sh in sample_hits:
            print(f"  CC{cc}  p_0={p0h}  Q={Qh}  Q_bits={qb}  "
                  f"p_0_bits={pb}  D_V={sh}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
