#!/usr/bin/env python3
# Copyright (c) 2026 Nenad Mićić <nenad@micic.be>
# SPDX-License-Identifier: Apache-2.0
"""
Ghost Chain Finder — Scan CC10+ roots for chains with many primes beyond official length.
For each root p, compute links: p, 2p+1, 4p+3, ..., 2^k*p + 2^k - 1 up to link 20.
Uses multiprocessing for speed.
"""

import sys
import time
import multiprocessing as mp
from collections import defaultdict

DATA_FILE = "cc10plus_roots_snapshot_2026-03-12.txt"
OUT_FILE  = "ghost_chains_top.txt"

MAX_DEPTH = 20

# For numbers up to 3.317e24 (~82 bits), first 12 primes as witnesses is deterministic.
# Our numbers go up to ~108 bits. Using 15 witnesses is very reliable.
# Fixed witnesses - no random RNG overhead.
WITNESSES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]

# Small primes for trial division (tuple for speed)
def _sieve(limit):
    sieve = bytearray(b'\x01') * (limit + 1)
    sieve[0] = sieve[1] = 0
    for i in range(2, int(limit**0.5) + 1):
        if sieve[i]:
            sieve[i*i::i] = bytearray(len(sieve[i*i::i]))
    return tuple(i for i in range(2, limit + 1) if sieve[i])

SMALL_PRIMES = _sieve(500)  # trial div to 500 is enough to filter most composites fast
TRIAL_PRIMES = _sieve(100000)

def _miller_rabin(n):
    """Miller-Rabin with 15 fixed witnesses. Returns True if probably prime."""
    # Trial division
    for p in SMALL_PRIMES:
        if n == p:
            return True
        if n % p == 0:
            return False

    d = n - 1
    r = 0
    while d & 1 == 0:
        d >>= 1
        r += 1

    for a in WITNESSES:
        if a >= n:
            continue
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        composite = True
        for _ in range(r - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                composite = False
                break
        if composite:
            return False
    return True

def _process_chunk(chunk):
    """Process a chunk of (cc_len, p, hex_val) tuples. Returns list of (prime_count, cc_len, hex_val, pattern_str, p)."""
    results = []
    for cc_len, p, hex_val in chunk:
        pattern = []
        prime_count = 0
        link = p
        for depth in range(MAX_DEPTH):
            if depth > 0:
                link = (link << 1) | 1  # 2*link + 1
            if _miller_rabin(link):
                pattern.append('P')
                prime_count += 1
            else:
                pattern.append('.')
        results.append((prime_count, cc_len, hex_val, ''.join(pattern), p))
    return results

def find_small_factor(n):
    for p in TRIAL_PRIMES:
        if n % p == 0:
            return p
    return 0

def find_factors_for_root(p, pattern):
    small_factors = {}
    link = p
    for depth in range(MAX_DEPTH):
        if depth > 0:
            link = (link << 1) | 1
        if pattern[depth] == '.':
            f = find_small_factor(link)
            if f:
                small_factors[depth] = f
    return small_factors

def main():
    ncpu = mp.cpu_count()
    print(f"Ghost Chain Finder — {ncpu} CPUs", flush=True)
    print(f"Data: {DATA_FILE}", flush=True)
    print(f"Depth: {MAX_DEPTH}, MR witnesses: {len(WITNESSES)}", flush=True)
    print(flush=True)

    # Read all roots
    t0 = time.time()
    roots = []
    with open(DATA_FILE, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            cc_len = int(parts[0][2:])
            hex_val = parts[1]
            p = int(hex_val, 16)
            roots.append((cc_len, p, hex_val))

    t1 = time.time()
    print(f"Loaded {len(roots):,} roots in {t1 - t0:.1f}s", flush=True)

    cc_counts = defaultdict(int)
    for cc_len, _, _ in roots:
        cc_counts[cc_len] += 1
    for k in sorted(cc_counts):
        print(f"  CC{k}: {cc_counts[k]:,}", flush=True)
    print(flush=True)

    # Split into chunks for multiprocessing
    chunk_size = 5000
    chunks = []
    for i in range(0, len(roots), chunk_size):
        chunks.append(roots[i:i+chunk_size])

    print(f"Processing {len(roots):,} roots in {len(chunks)} chunks of ~{chunk_size} using {ncpu} workers...", flush=True)

    results = []
    t_start = time.time()
    done = 0

    with mp.Pool(ncpu) as pool:
        for chunk_result in pool.imap_unordered(_process_chunk, chunks):
            results.extend(chunk_result)
            done += len(chunk_result)
            if done % 50000 < chunk_size:
                elapsed = time.time() - t_start
                rate = done / elapsed if elapsed > 0 else 0
                eta = (len(roots) - done) / rate if rate > 0 else 0
                print(f"  {done:,}/{len(roots):,}  ({100*done/len(roots):.1f}%)  "
                      f"{rate:.0f}/s  elapsed={elapsed:.0f}s  ETA={eta:.0f}s", flush=True)

    t_end = time.time()
    total_time = t_end - t_start
    print(f"\nDone: {len(roots):,} roots in {total_time:.1f}s "
          f"({len(roots)/total_time:.0f} roots/s)\n", flush=True)

    # Sort by prime_count desc, then cc_len desc
    results.sort(key=lambda r: (-r[0], -r[1]))

    # Factor analysis for top results
    print("Finding small factors for top results...", flush=True)
    factors_cache = {}
    for r in results[:500]:
        if r[0] >= 13:
            sf = find_factors_for_root(r[4], r[3])
            if sf:
                factors_cache[r[2]] = sf
    print(f"  Factor analysis done for {len(factors_cache)} entries\n", flush=True)

    # Distribution
    dist = defaultdict(int)
    for r in results:
        dist[r[0]] += 1

    lines = []
    def out(s=""):
        lines.append(s)
        print(s, flush=True)

    def fmt(r, rank=None):
        prime_count, cc_len, hex_val, pattern, p = r
        prefix = f"#{rank:<4}" if rank is not None else "     "
        bits = p.bit_length()
        s = f"{prefix} CC{cc_len:<3} {prime_count}/20  {pattern}  {hex_val}  {bits}bit"
        if hex_val in factors_cache:
            sf = factors_cache[hex_val]
            fstr = ", ".join(f"L{d+1}:{f}" for d, f in sorted(sf.items())[:6])
            s += f"  [{fstr}]"
        return s

    W = 120
    out("=" * W)
    out("GHOST CHAIN ANALYSIS — CC10+ roots, first-kind Cunningham chains extended to depth 20")
    out(f"Total roots: {len(roots):,}   Processed in {total_time:.1f}s")
    out(f"Pattern: P=prime .=composite   L1=p, L2=2p+1, L3=4p+3, ..., L20=2^19*p+2^19-1")
    out("=" * W)

    # (f) Distribution
    out("\n--- DISTRIBUTION: primes out of 20 links ---")
    for k in sorted(dist.keys(), reverse=True):
        if k >= 10:
            pct = 100.0 * dist[k] / len(roots)
            out(f"  {k}/20 primes: {dist[k]:>8,} roots  ({pct:.4f}%)")

    # (a) Top 50 overall
    out("\n--- TOP 50 ROOTS BY TOTAL PRIMES (all CC10+) ---")
    for i in range(min(50, len(results))):
        out(fmt(results[i], i + 1))

    # (b) CC10 top 20
    cc10 = [r for r in results if r[1] == 10]
    out(f"\n--- TOP 20 CC10 ROOTS — hidden gems ({len(cc10):,} total CC10) ---")
    for i in range(min(20, len(cc10))):
        out(fmt(cc10[i], i + 1))

    # (c) CC11 top 10
    cc11 = [r for r in results if r[1] == 11]
    out(f"\n--- TOP 10 CC11 ROOTS ({len(cc11):,} total CC11) ---")
    for i in range(min(10, len(cc11))):
        out(fmt(cc11[i], i + 1))

    # (d) CC12 top 10
    cc12 = [r for r in results if r[1] == 12]
    out(f"\n--- TOP 10 CC12 ROOTS ({len(cc12):,} total CC12) ---")
    for i in range(min(10, len(cc12))):
        out(fmt(cc12[i], i + 1))

    # (e) CC13+ top 10
    cc13p = [r for r in results if r[1] >= 13]
    out(f"\n--- TOP 10 CC13+ ROOTS ({len(cc13p):,} total CC13+) ---")
    for i in range(min(10, len(cc13p))):
        out(fmt(cc13p[i], i + 1))

    # Notable: both L19 and L20 prime
    out("\n--- NOTABLE: BOTH LINK 19 AND 20 PRIME ---")
    far_both = [r for r in results if r[3][18] == 'P' and r[3][19] == 'P']
    out(f"Count: {len(far_both):,}")
    out("Top 10 by total primes:")
    for i in range(min(10, len(far_both))):
        out(fmt(far_both[i], i + 1))

    # Second wind: longest consecutive prime run after first gap
    out("\n--- BEST 'SECOND WIND' CHAINS ---")
    out("(Longest consecutive prime run starting AFTER first composite link)")
    second_wind = []
    for r in results:
        pat = r[3]
        first_dot = pat.find('.')
        if first_dot < 0:
            continue
        best_run = 0
        run = 0
        for ch in pat[first_dot:]:
            if ch == 'P':
                run += 1
                if run > best_run:
                    best_run = run
            else:
                run = 0
        if best_run >= 4:
            second_wind.append((best_run, r))
    second_wind.sort(key=lambda x: (-x[0], -x[1][0]))
    for i in range(min(20, len(second_wind))):
        sw_len, r = second_wind[i]
        out(f"  #{i+1:<4} 2nd-wind={sw_len}  {fmt(r)}")

    out("\n" + "=" * W)

    # Save
    with open(OUT_FILE, 'w') as f:
        f.write('\n'.join(lines) + '\n')
    print(f"\nResults saved to {OUT_FILE}", flush=True)

if __name__ == '__main__':
    main()
