#!/usr/bin/env python3
# Copyright (c) 2026 Nenad Mićić <nenad@micic.be>
# SPDX-License-Identifier: Apache-2.0
"""build_known_fingerprints_primorial.py — generate a structural fingerprint
pool from primorial-quotient patterns, using the primorial-quotient seed form
(e.g. `41#/31`, `71#/31`).

Pattern grammar:
    D_V = X#  /  Y1  /  Y2 ... /  Yk
where X ∈ PRIME_CEILINGS, k ∈ DROP_COUNTS, and Y1..Yk ⊂ primes(X#) \\ PINNED.

PINNED is the sieve-floor set kept in every variant ({2, 3, 5} to match
v16's `--immune-prime 5`).  Variants whose D bit-width falls outside
[D_BITS_MIN, D_BITS_MAX] are skipped — those are either too small to
matter (Q-window dominates) or too large to give Q any dynamic range
inside the engine's `--bits-min/--bits-max` band.

Output is in the v16 pool format:
    P1 [^e1] * P2 [^e2] * ...   # primorial=X# drop={Y1,Y2,..} bits=BB.B
(All variants emitted here have e_i = 1 — primorial quotients are squarefree.)

Why this pool exists (separate from the historical-frequency pools):
  - `known_fingerprints_exp_d.txt` filters out singleton-frequency D shapes.
    exotic singleton patterns like `41#/31` appear once in history and
    get dropped from the production pool.
  - This generator is structural / exhaustive, not data-driven: every D in
    the primorial-quotient family inside the bit band gets a variant slot.
    Trade-off: pool size grows, per-V GPU time drops.

Usage:
    python3 scripts/build_known_fingerprints_primorial.py
    python3 scripts/build_known_fingerprints_primorial.py --max-ceiling 53
    python3 scripts/build_known_fingerprints_primorial.py --d-bits 30 80 --drop-max 2
"""

import argparse
from itertools import combinations


# Prime ceilings for X# (the primorial size).  Each cap defines a primorial
# X# = product(primes <= X). This form covers seeds of length 14-17
# historically used in long-chain searches.
DEFAULT_CEILINGS = [41, 43, 47, 53, 59, 61, 67, 71]

# Numbers of primes to drop from X# (the "/Y/Z/W" part).
# Drop=0 emits the unmodified primorial; included for symmetry.
DEFAULT_DROP_COUNTS = [0, 1, 2, 3]

# Sieve-floor primes that must remain in every D (matches --immune-prime 5).
DEFAULT_PINNED = (2, 3, 5)

# Bit-width band for the resulting D.  Too small → Q window is huge but
# variant offers no structural advantage.  Too large → Q range vanishes
# inside `--bits-min/--bits-max`.
DEFAULT_D_BITS_MIN = 30
DEFAULT_D_BITS_MAX = 85


def primes_up_to(n):
    """Sieve of Eratosthenes through n inclusive."""
    if n < 2:
        return []
    sieve = bytearray([1]) * (n + 1)
    sieve[0] = sieve[1] = 0
    for i in range(2, int(n ** 0.5) + 1):
        if sieve[i]:
            for j in range(i * i, n + 1, i):
                sieve[j] = 0
    return [i for i in range(n + 1) if sieve[i]]


def primorial(primes):
    p = 1
    for q in primes:
        p *= q
    return p


def fmt_pool_line(kept_primes, dropped, ceiling):
    """Emit '2*3*5*...' format with annotation."""
    body = "*".join(str(q) for q in kept_primes)
    D = primorial(kept_primes)
    bits = D.bit_length() - 1 + (1 - 1)  # exact integer bit-length
    drop_str = "{" + ",".join(str(q) for q in dropped) + "}" if dropped else "-"
    return f"{body}  # primorial={ceiling}# drop={drop_str} bits={D.bit_length()}"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ceilings", nargs="+", type=int, default=DEFAULT_CEILINGS,
                    help="prime ceilings for X# (default: %(default)s)")
    ap.add_argument("--max-ceiling", type=int, default=None,
                    help="shorthand: keep only ceilings <= this value")
    ap.add_argument("--drop-counts", nargs="+", type=int, default=DEFAULT_DROP_COUNTS,
                    help="number of primes to drop (default: %(default)s)")
    ap.add_argument("--pinned", nargs="+", type=int, default=list(DEFAULT_PINNED),
                    help="primes that must stay in every D (default: %(default)s)")
    ap.add_argument("--d-bits", nargs=2, type=int, metavar=("MIN", "MAX"),
                    default=[DEFAULT_D_BITS_MIN, DEFAULT_D_BITS_MAX],
                    help="D bit-width band (default: %(default)s)")
    ap.add_argument("--out", type=str,
                    default="pools/known_fingerprints_primorial.txt",
                    help="output pool file (default: %(default)s)")
    args = ap.parse_args()

    ceilings = sorted(set(args.ceilings))
    if args.max_ceiling is not None:
        ceilings = [c for c in ceilings if c <= args.max_ceiling]
    drop_counts = sorted(set(args.drop_counts))
    pinned = frozenset(args.pinned)
    d_bits_min, d_bits_max = args.d_bits

    # De-dupe by the kept-prime tuple (same D shape may arise from different
    # (ceiling, drop-set) parameterizations, e.g. 47#/47 = 43#).
    emitted = {}
    raw_count = 0
    skipped_bits = 0
    skipped_dup = 0

    for X in ceilings:
        X_primes = primes_up_to(X)
        droppable = [q for q in X_primes if q not in pinned]
        for k in drop_counts:
            if k > len(droppable):
                continue
            for drop in combinations(droppable, k):
                raw_count += 1
                kept = tuple(q for q in X_primes if q not in drop)
                D = primorial(kept)
                if D.bit_length() < d_bits_min or D.bit_length() > d_bits_max:
                    skipped_bits += 1
                    continue
                if kept in emitted:
                    skipped_dup += 1
                    continue
                emitted[kept] = (X, drop)

    # Sort by D bit-width ascending, then by tuple lex.
    rows = sorted(emitted.items(),
                  key=lambda kv: (primorial(kv[0]).bit_length(), kv[0]))

    # Write
    import os
    out_path = args.out
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    with open(out_path, "w") as fh:
        fh.write(f"# Structural primorial-quotient fingerprint pool.\n")
        fh.write(f"# Generated by scripts/build_known_fingerprints_primorial.py\n")
        fh.write(f"# Ceilings: {ceilings}\n")
        fh.write(f"# Drop counts: {drop_counts}\n")
        fh.write(f"# Pinned (always in D): {sorted(pinned)}\n")
        fh.write(f"# D bit band: [{d_bits_min}, {d_bits_max}]\n")
        fh.write(f"# Total variants emitted: {len(rows)}\n")
        fh.write(f"#\n")
        fh.write(f"# Format: P1*P2*... [# primorial=X# drop={{...}} bits=BB]\n")
        fh.write(f"# (All variants here are squarefree; exponents implicit = 1.)\n")
        fh.write(f"#\n")
        for kept, (ceiling, drop) in rows:
            fh.write(fmt_pool_line(kept, drop, ceiling) + "\n")

    print(f"Generated {len(rows)} unique primorial-quotient variants -> {out_path}")
    print(f"  raw enumerated: {raw_count}")
    print(f"  filtered (D bits outside [{d_bits_min},{d_bits_max}]): {skipped_bits}")
    print(f"  filtered (duplicate D shape): {skipped_dup}")
    if rows:
        first_bits = primorial(rows[0][0]).bit_length()
        last_bits = primorial(rows[-1][0]).bit_length()
        print(f"  D bit range emitted: [{first_bits}, {last_bits}]")


if __name__ == "__main__":
    main()
