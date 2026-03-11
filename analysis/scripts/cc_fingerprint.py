#!/usr/bin/env python3
# Copyright (c) 2026 Nenad Mićić <nenad@micic.be>
# SPDX-License-Identifier: Apache-2.0
"""
cc_fingerprint.py — Cunningham Chain Factorization Fingerprint Analyzer

Pass 1: For each CC root, produce detailed CSV of factorization structure
        for the root, each chain member, and their p±1 neighborhoods.

Pass 2: Aggregate across many roots to find structural patterns,
        recurring large factors, shadow detection, bit-shell analysis.

Usage:
  # From a list of primes (one per line, decimal or 0x hex):
  python3 cc_fingerprint.py --analyze roots.txt --output fingerprints.csv
  
  # Aggregate existing CSV:
  python3 cc_fingerprint.py --aggregate fingerprints.csv --output summary.txt
  
  # Quick test with known CC13:
  python3 cc_fingerprint.py --test

  # Generate from a range (find CC roots ourselves):
  python3 cc_fingerprint.py --scan --bits 89 --min-chain 8 --count 50

Author: Nenad Micic / Claude
Date: February 2026
"""

import sys
import csv
import os
import argparse
from collections import defaultdict, Counter
from io import StringIO

try:
    import gmpy2
    from gmpy2 import mpz, is_prime as gmp_is_prime, next_prime, iroot
    HAS_GMPY2 = True
except ImportError:
    HAS_GMPY2 = False
    print("WARNING: gmpy2 not found, using slow fallback", file=sys.stderr)


# =============================================================================
# PRIMALITY + FACTORIZATION
# =============================================================================

def is_prime(n):
    """Primality test using gmpy2 or fallback."""
    if n < 2:
        return False
    if HAS_GMPY2:
        return gmp_is_prime(mpz(n)) > 0  # gmpy2 returns 0,1,2
    # Fallback: trial division (slow for big numbers)
    if n < 4:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True


# Small primes for trial division during factorization
SMALL_PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47,
                53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109,
                113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179,
                181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241,
                251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313,
                317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389,
                397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461,
                463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547,
                557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617,
                619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691,
                701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773,
                787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859,
                863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947,
                953, 967, 971, 977, 983, 991, 997]

# Primes used in CRT/wheel/filter (for residue analysis)
CRT_PRIMES = [5, 7, 11, 13, 17, 19]
WHEEL_PRIMES = [23, 29, 31]
FILTER_PRIMES = [37, 41, 43, 47, 53, 59, 61]
SCREEN_PRIMES = [67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113,
                 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179,
                 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239,
                 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307,
                 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373,
                 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439,
                 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503,
                 509, 521, 523, 541, 547]

ALL_ANALYSIS_PRIMES = CRT_PRIMES + WHEEL_PRIMES + FILTER_PRIMES + [67, 71, 73, 79, 83, 89, 97]


def partial_factor(n, trial_limit=10000):
    """
    Factor n using trial division up to trial_limit, then check if remainder is prime.
    Returns list of (prime, exponent) pairs. Last factor may be composite (marked).
    
    For CC analysis we care about:
    - The small prime structure (v2, v3, v5, etc.)
    - Whether the cofactor is prime (smooth neighborhood)
    - The size of the largest prime factor
    """
    if n <= 1:
        return [(n, 1)] if n == 1 else []
    
    n = int(n)
    factors = []
    
    # Trial division by small primes
    for p in SMALL_PRIMES:
        if p * p > n:
            break
        if p > trial_limit:
            break
        e = 0
        while n % p == 0:
            n //= p
            e += 1
        if e > 0:
            factors.append((p, e))
    
    # Continue with odd numbers beyond our prime list
    if SMALL_PRIMES[-1] < trial_limit:
        p = SMALL_PRIMES[-1] + 2
        while p <= trial_limit and p * p <= n:
            e = 0
            while n % p == 0:
                n //= p
                e += 1
            if e > 0:
                factors.append((p, e))
            p += 2
    
    if n > 1:
        if is_prime(n):
            factors.append((n, 1))
        else:
            # Cofactor is composite and large - try Pollard rho
            sub_factors = pollard_rho_factor(n)
            factors.extend(sub_factors)
    
    return sorted(factors)


def pollard_rho_factor(n, max_iter=1000000):
    """Simple Pollard rho for medium-size composites."""
    if n <= 1:
        return []
    if is_prime(n):
        return [(n, 1)]
    if n % 2 == 0:
        e = 0
        while n % 2 == 0:
            n //= 2
            e += 1
        result = [(2, e)]
        if n > 1:
            result.extend(pollard_rho_factor(n))
        return result
    
    if HAS_GMPY2:
        from gmpy2 import gcd as _gcd
    else:
        from math import gcd as _gcd
    
    # Try several c values
    for c in range(1, 20):
        x = 2
        y = 2
        d = 1
        f = lambda v: (v * v + c) % n
        
        iterations = 0
        while d == 1 and iterations < max_iter:
            x = f(x)
            y = f(f(y))
            d = int(_gcd(abs(x - y), n))
            iterations += 1
        
        if d != n and d != 1:
            # Found a factor
            result = []
            for sub_n in [d, n // d]:
                if is_prime(sub_n):
                    result.append((sub_n, 1))
                else:
                    result.extend(pollard_rho_factor(sub_n))
            # Merge duplicates
            merged = defaultdict(int)
            for p, e in result:
                merged[p] += e
            return sorted(merged.items())
    
    # Failed to factor - return as composite
    return [(n, 1)]


# =============================================================================
# CHAIN FOLLOWING
# =============================================================================

def follow_chain_first(p):
    """Follow first-kind chain from root p. Returns list of primes."""
    chain = []
    current = p
    while is_prime(current):
        chain.append(int(current))
        current = 2 * current + 1
    return chain


def follow_chain_second(p):
    """Follow second-kind chain from root p. Returns list of primes."""
    chain = []
    current = p
    while is_prime(current):
        chain.append(int(current))
        current = 2 * current - 1
    return chain


def is_root_first(p):
    """Check if p is a first-kind chain root (predecessor composite)."""
    if not is_prime(p):
        return False
    pred = (p - 1) // 2
    if pred < 2:
        return True
    return not is_prime(pred)


def is_root_second(p):
    """Check if p is a second-kind chain root."""
    if not is_prime(p):
        return False
    pred = (p + 1) // 2
    if pred < 2:
        return True
    return not is_prime(pred)


# =============================================================================
# FINGERPRINT COMPUTATION
# =============================================================================

def v2(n):
    """2-adic valuation."""
    if n == 0:
        return -1
    n = int(n)
    count = 0
    while n % 2 == 0:
        n //= 2
        count += 1
    return count


def vp(n, p):
    """p-adic valuation."""
    if n == 0:
        return -1
    n = int(n)
    count = 0
    while n % p == 0:
        n //= p
        count += 1
    return count


def trailing_ones(n):
    """Count trailing 1-bits in binary."""
    n = int(n)
    count = 0
    while n & 1:
        n >>= 1
        count += 1
    return count


def odd_core(n):
    """Remove all factors of 2."""
    n = int(n)
    while n % 2 == 0:
        n //= 2
    return n


def factors_to_str(factors):
    """Convert factor list to compact string like '2^3 * 3 * 17 * P47'."""
    parts = []
    for p, e in factors:
        bits = p.bit_length() if isinstance(p, int) else int(p).bit_length()
        if bits > 40 and is_prime(p):
            label = f"P{bits}b"  # Large prime, show bit count
        elif bits > 40:
            label = f"C{bits}b"  # Large composite
        else:
            label = str(p)
        
        if e == 1:
            parts.append(label)
        else:
            parts.append(f"{label}^{e}")
    return " * ".join(parts) if parts else "1"


def largest_prime_factor(factors):
    """Return the largest prime factor from a factor list."""
    largest = 1
    for p, e in factors:
        if is_prime(p) and p > largest:
            largest = p
    return largest


def smooth_part(factors, B):
    """Product of prime factors <= B."""
    result = 1
    for p, e in factors:
        if p <= B:
            result *= p ** e
    return result


def compute_residues(n, primes):
    """Compute n mod p for each prime in list."""
    return {p: int(n) % p for p in primes}


def immunizing_residue(p, kind):
    """Return the immunizing residue for prime p and chain kind."""
    if kind == "first":
        return p - 1   # n ≡ p-1 (mod p) → chain never killed by p
    else:
        return 1        # n ≡ 1 (mod p) → chain never killed by p


def kill_position(residue, q, kind, max_pos=32):
    """Compute earliest chain position where q divides a member."""
    r = residue % q
    for i in range(max_pos):
        if r == 0:
            return i
        if kind == "first":
            r = (2 * r + 1) % q
        else:
            r = (2 * r - 1) % q
    return max_pos  # Immunized (or beyond range)


def fingerprint_prime(p, kind="first"):
    """
    Compute full fingerprint of a prime p as a CC member.
    Returns dict with all analysis fields.
    """
    p = int(p)
    bits = p.bit_length()
    
    # Basic structure
    t1 = trailing_ones(p)
    oc = odd_core(p)
    v2_val = v2(p - 1) if p > 1 else 0
    v2_plus = v2(p + 1)
    
    # Factorizations of p-1 and p+1
    pm1_factors = partial_factor(p - 1)
    pp1_factors = partial_factor(p + 1)
    
    # Largest prime factors
    lpf_pm1 = largest_prime_factor(pm1_factors)
    lpf_pp1 = largest_prime_factor(pp1_factors)
    
    # Smooth parts
    smooth_pm1_100 = smooth_part(pm1_factors, 100)
    smooth_pp1_100 = smooth_part(pp1_factors, 100)
    smooth_pm1_1000 = smooth_part(pm1_factors, 1000)
    smooth_pp1_1000 = smooth_part(pp1_factors, 1000)
    
    # Residues mod analysis primes
    residues = compute_residues(p, ALL_ANALYSIS_PRIMES)
    
    # Kill positions for each analysis prime
    kills = {}
    immunized = {}
    for q in ALL_ANALYSIS_PRIMES:
        r = p % q
        kp = kill_position(r, q, kind)
        kills[q] = kp
        immunized[q] = (r == immunizing_residue(q, kind))
    
    # Minimum kill position = chain ceiling from these primes
    min_kill = min(kills[q] for q in ALL_ANALYSIS_PRIMES if kills[q] < 32)
    num_immunized = sum(1 for q in ALL_ANALYSIS_PRIMES if immunized[q])
    
    # p-adic valuations for small primes
    v3_pm1 = vp(p - 1, 3)
    v5_pm1 = vp(p - 1, 5)
    v3_pp1 = vp(p + 1, 3)
    v5_pp1 = vp(p + 1, 5)
    
    return {
        "p": p,
        "bits": bits,
        "trailing_ones": t1,
        "odd_core_bits": oc.bit_length(),
        "v2_pm1": v2_val,
        "v3_pm1": v3_pm1,
        "v5_pm1": v5_pm1,
        "v2_pp1": v2_plus,
        "v3_pp1": v3_pp1,
        "v5_pp1": v5_pp1,
        "pm1_factors": factors_to_str(pm1_factors),
        "pp1_factors": factors_to_str(pp1_factors),
        "lpf_pm1": lpf_pm1,
        "lpf_pm1_bits": lpf_pm1.bit_length(),
        "lpf_pp1": lpf_pp1,
        "lpf_pp1_bits": lpf_pp1.bit_length(),
        "smooth_pm1_100": smooth_pm1_100,
        "smooth_pp1_100": smooth_pp1_100,
        "smooth_pm1_1000": smooth_pm1_1000,
        "smooth_pp1_1000": smooth_pp1_1000,
        "min_kill_pos": min_kill,
        "num_immunized": num_immunized,
        "residues": residues,
        "kills": kills,
        "immunized_primes": [q for q in ALL_ANALYSIS_PRIMES if immunized[q]],
    }


def analyze_chain(root, kind="first"):
    """
    Full analysis of a CC chain starting at root.
    Returns list of fingerprint dicts, one per chain member.
    """
    if kind == "first":
        chain = follow_chain_first(root)
    else:
        chain = follow_chain_second(root)
    
    results = []
    for i, p in enumerate(chain):
        fp = fingerprint_prime(p, kind)
        fp["chain_position"] = i
        fp["chain_length"] = len(chain)
        fp["is_root"] = (i == 0)
        fp["kind"] = kind
        results.append(fp)
    
    return results


# =============================================================================
# CSV OUTPUT (PASS 1)
# =============================================================================

# Core CSV columns (compact)
CSV_COLUMNS = [
    "chain_id", "chain_length", "kind", "chain_position", "is_root",
    "p_hex", "bits", "trailing_ones",
    "v2_pm1", "v3_pm1", "v5_pm1",
    "v2_pp1", "v3_pp1", "v5_pp1",
    "lpf_pm1_bits", "lpf_pp1_bits",
    "min_kill_pos", "num_immunized",
    "pm1_factors", "pp1_factors",
]

# Add residue columns for key primes
RESIDUE_PRIMES_CSV = [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]
for q in RESIDUE_PRIMES_CSV:
    CSV_COLUMNS.append(f"r_mod{q}")

# Add kill-position columns
KILL_PRIMES_CSV = [5, 7, 11, 13, 17, 19, 23, 29, 31]
for q in KILL_PRIMES_CSV:
    CSV_COLUMNS.append(f"kill_{q}")

# Immunized primes (comma-separated list)
CSV_COLUMNS.append("immunized_list")


def chain_to_csv_rows(chain_id, chain_fps):
    """Convert chain fingerprints to CSV rows."""
    rows = []
    for fp in chain_fps:
        row = {
            "chain_id": chain_id,
            "chain_length": fp["chain_length"],
            "kind": fp["kind"],
            "chain_position": fp["chain_position"],
            "is_root": 1 if fp["is_root"] else 0,
            "p_hex": hex(fp["p"]),
            "bits": fp["bits"],
            "trailing_ones": fp["trailing_ones"],
            "v2_pm1": fp["v2_pm1"],
            "v3_pm1": fp["v3_pm1"],
            "v5_pm1": fp["v5_pm1"],
            "v2_pp1": fp["v2_pp1"],
            "v3_pp1": fp["v3_pp1"],
            "v5_pp1": fp["v5_pp1"],
            "lpf_pm1_bits": fp["lpf_pm1_bits"],
            "lpf_pp1_bits": fp["lpf_pp1_bits"],
            "min_kill_pos": fp["min_kill_pos"],
            "num_immunized": fp["num_immunized"],
            "pm1_factors": fp["pm1_factors"],
            "pp1_factors": fp["pp1_factors"],
        }
        
        # Residues
        for q in RESIDUE_PRIMES_CSV:
            row[f"r_mod{q}"] = fp["residues"].get(q, "")
        
        # Kill positions
        for q in KILL_PRIMES_CSV:
            row[f"kill_{q}"] = fp["kills"].get(q, "")
        
        # Immunized list
        row["immunized_list"] = ",".join(str(q) for q in fp["immunized_primes"])
        
        rows.append(row)
    return rows


# =============================================================================
# PASS 2: AGGREGATION
# =============================================================================

def aggregate_csv(csv_path, out=sys.stdout):
    """
    Read fingerprint CSV and produce aggregate analysis.
    Looks for:
    - Recurring large prime factors in p-1, p+1
    - v2/v3/v5 distribution patterns
    - Kill position / immunization patterns
    - Bit-shell differences
    - Residue distributions
    """
    roots = []
    members = []
    
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row.get("is_root") == "1":
                roots.append(row)
            members.append(row)
    
    print("=" * 72, file=out)
    print("  CUNNINGHAM CHAIN FINGERPRINT AGGREGATION", file=out)
    print("=" * 72, file=out)
    print(f"\n  Total chains: {len(roots)}", file=out)
    print(f"  Total members: {len(members)}", file=out)
    
    if not roots:
        print("  No roots found in CSV.", file=out)
        return
    
    # --- Bit-shell distribution ---
    print(f"\n{'─'*72}", file=out)
    print("  BIT-SHELL DISTRIBUTION (roots only)", file=out)
    print(f"{'─'*72}", file=out)
    
    bit_dist = Counter(int(r["bits"]) for r in roots)
    chain_by_bits = defaultdict(list)
    for r in roots:
        chain_by_bits[int(r["bits"])].append(int(r["chain_length"]))
    
    for bits in sorted(bit_dist.keys()):
        lengths = chain_by_bits[bits]
        avg_len = sum(lengths) / len(lengths)
        max_len = max(lengths)
        print(f"  {bits}-bit: {bit_dist[bits]} roots, avg chain={avg_len:.1f}, max={max_len}", file=out)
    
    # --- Chain length distribution ---
    print(f"\n{'─'*72}", file=out)
    print("  CHAIN LENGTH DISTRIBUTION", file=out)
    print(f"{'─'*72}", file=out)
    
    len_dist = Counter(int(r["chain_length"]) for r in roots)
    for length in sorted(len_dist.keys()):
        print(f"  CC{length}: {len_dist[length]} chains", file=out)
    
    # --- v2(p-1) distribution for roots ---
    print(f"\n{'─'*72}", file=out)
    print("  v2(p-1) DISTRIBUTION (roots only)", file=out)
    print(f"{'─'*72}", file=out)
    
    v2_dist = Counter(int(r["v2_pm1"]) for r in roots)
    v2_by_chain = defaultdict(list)
    for r in roots:
        v2_by_chain[int(r["v2_pm1"])].append(int(r["chain_length"]))
    
    for v in sorted(v2_dist.keys()):
        lengths = v2_by_chain[v]
        avg = sum(lengths) / len(lengths)
        print(f"  v2(p-1) = {v}: {v2_dist[v]} roots, avg chain = {avg:.1f}", file=out)
    
    # --- Trailing ones distribution ---
    print(f"\n{'─'*72}", file=out)
    print("  TRAILING ONES (root) — position in 2-adic tree", file=out)
    print(f"{'─'*72}", file=out)
    
    t1_dist = Counter(int(r["trailing_ones"]) for r in roots)
    t1_by_chain = defaultdict(list)
    for r in roots:
        t1_by_chain[int(r["trailing_ones"])].append(int(r["chain_length"]))
    
    for t in sorted(t1_dist.keys()):
        lengths = t1_by_chain[t]
        avg = sum(lengths) / len(lengths)
        print(f"  trailing_ones = {t}: {t1_dist[t]} roots, avg chain = {avg:.1f}", file=out)
    
    # --- Immunization patterns ---
    print(f"\n{'─'*72}", file=out)
    print("  IMMUNIZATION PATTERNS (roots only)", file=out)
    print(f"{'─'*72}", file=out)
    
    imm_dist = Counter(int(r["num_immunized"]) for r in roots)
    imm_by_chain = defaultdict(list)
    for r in roots:
        imm_by_chain[int(r["num_immunized"])].append(int(r["chain_length"]))
    
    for ni in sorted(imm_dist.keys()):
        lengths = imm_by_chain[ni]
        avg = sum(lengths) / len(lengths)
        print(f"  immunized_count = {ni}: {imm_dist[ni]} roots, avg chain = {avg:.1f}", file=out)
    
    # Which primes are most often immunized?
    imm_prime_count = Counter()
    for r in roots:
        imm_list = r.get("immunized_list", "")
        if imm_list:
            for q in imm_list.split(","):
                if q.strip():
                    imm_prime_count[int(q.strip())] += 1
    
    print(f"\n  Most frequently immunized primes (across roots):", file=out)
    for q, count in imm_prime_count.most_common(20):
        pct = 100.0 * count / len(roots)
        print(f"    p={q}: immunized in {count}/{len(roots)} roots ({pct:.1f}%)", file=out)
    
    # --- Kill position analysis ---
    print(f"\n{'─'*72}", file=out)
    print("  KILL POSITION ANALYSIS (what limits chain length?)", file=out)
    print(f"{'─'*72}", file=out)
    
    # For each root, what prime has the minimum kill position?
    limiting_prime = Counter()
    for r in roots:
        min_kp = 999
        min_q = None
        for q in KILL_PRIMES_CSV:
            kp_str = r.get(f"kill_{q}", "32")
            kp = int(kp_str) if kp_str else 32
            if kp < min_kp:
                min_kp = kp
                min_q = q
        if min_q:
            limiting_prime[min_q] += 1
    
    print(f"  Prime most often limiting chain ceiling:", file=out)
    for q, count in limiting_prime.most_common(15):
        pct = 100.0 * count / len(roots)
        print(f"    q={q}: limiting in {count}/{len(roots)} roots ({pct:.1f}%)", file=out)
    
    # --- Residue distribution mod small primes ---
    print(f"\n{'─'*72}", file=out)
    print("  RESIDUE ENRICHMENT (roots vs uniform)", file=out)
    print(f"{'─'*72}", file=out)
    print(f"  (For long chains, certain residues should be enriched)", file=out)
    
    # Split roots into "long" (chain >= 10) and "short" 
    long_roots = [r for r in roots if int(r["chain_length"]) >= 10]
    short_roots = [r for r in roots if int(r["chain_length"]) < 10]
    
    if long_roots and short_roots:
        for q in [5, 7, 11, 13, 17, 19, 23, 29, 31]:
            long_res = Counter(int(r.get(f"r_mod{q}", 0)) for r in long_roots)
            short_res = Counter(int(r.get(f"r_mod{q}", 0)) for r in short_roots)
            
            print(f"\n  mod {q} (long chains ≥10 vs short):", file=out)
            for res in sorted(set(list(long_res.keys()) + list(short_res.keys()))):
                lp = 100.0 * long_res[res] / len(long_roots) if long_roots else 0
                sp = 100.0 * short_res[res] / len(short_roots) if short_roots else 0
                enrichment = lp / sp if sp > 0 else float('inf')
                marker = " <<<" if enrichment > 2.0 else ""
                print(f"    r≡{res}: long={lp:.1f}% short={sp:.1f}% enrichment={enrichment:.2f}x{marker}", file=out)
    else:
        print(f"  (Need both long and short chains for comparison)", file=out)
    
    # --- Largest prime factor analysis ---
    print(f"\n{'─'*72}", file=out)
    print("  LARGEST PRIME FACTOR OF p-1 (roots only)", file=out)
    print(f"{'─'*72}", file=out)
    
    lpf_bits = [(int(r["lpf_pm1_bits"]), int(r["chain_length"]), int(r["bits"])) for r in roots]
    
    # Group by chain length
    lpf_by_chain = defaultdict(list)
    for lpf_b, cl, rb in lpf_bits:
        lpf_by_chain[cl].append(lpf_b)
    
    for cl in sorted(lpf_by_chain.keys()):
        vals = lpf_by_chain[cl]
        avg = sum(vals) / len(vals)
        print(f"  CC{cl}: avg lpf(p-1) = {avg:.1f} bits (from {len(vals)} roots)", file=out)
    
    # --- Recurring large factors ---
    print(f"\n{'─'*72}", file=out)
    print("  RECURRING FACTORS IN p-1 FACTORIZATIONS", file=out)
    print(f"{'─'*72}", file=out)
    
    # Parse factor strings to find recurring primes
    factor_freq = Counter()
    for r in roots:
        fstr = r.get("pm1_factors", "")
        # Extract individual prime values from "2^3 * 3 * 17 * P47b"
        for part in fstr.split(" * "):
            part = part.strip()
            if not part:
                continue
            # Remove exponent
            base = part.split("^")[0]
            if base.startswith("P") or base.startswith("C"):
                continue  # Skip large primes/composites in label form
            try:
                p_val = int(base)
                if p_val > 2:  # Skip 2 (always present)
                    factor_freq[p_val] += 1
            except ValueError:
                pass
    
    print(f"  Small primes appearing in p-1 factorizations:", file=out)
    for p_val, count in factor_freq.most_common(30):
        pct = 100.0 * count / len(roots)
        print(f"    {p_val}: appears in {count}/{len(roots)} roots ({pct:.1f}%)", file=out)
    
    # Same for p+1
    print(f"\n  Small primes appearing in p+1 factorizations:", file=out)
    factor_freq_pp1 = Counter()
    for r in roots:
        fstr = r.get("pp1_factors", "")
        for part in fstr.split(" * "):
            part = part.strip()
            if not part:
                continue
            base = part.split("^")[0]
            if base.startswith("P") or base.startswith("C"):
                continue
            try:
                p_val = int(base)
                if p_val > 2:
                    factor_freq_pp1[p_val] += 1
            except ValueError:
                pass
    
    for p_val, count in factor_freq_pp1.most_common(30):
        pct = 100.0 * count / len(roots)
        print(f"    {p_val}: appears in {count}/{len(roots)} roots ({pct:.1f}%)", file=out)
    
    # --- Per-bit-shell analysis ---
    print(f"\n{'─'*72}", file=out)
    print("  BIT-SHELL COMPARISON (is 89-bit different from 90-bit?)", file=out)
    print(f"{'─'*72}", file=out)
    
    shells = defaultdict(lambda: {"roots": [], "total_chain": 0, "imm_counts": []})
    for r in roots:
        b = int(r["bits"])
        shells[b]["roots"].append(r)
        shells[b]["total_chain"] += int(r["chain_length"])
        shells[b]["imm_counts"].append(int(r["num_immunized"]))
    
    for b in sorted(shells.keys()):
        s = shells[b]
        n_roots = len(s["roots"])
        if n_roots == 0:
            continue
        avg_chain = s["total_chain"] / n_roots
        avg_imm = sum(s["imm_counts"]) / n_roots
        print(f"  {b}-bit shell: {n_roots} roots, avg_chain={avg_chain:.1f}, avg_immunized={avg_imm:.1f}", file=out)
    
    print(f"\n{'='*72}", file=out)
    print("  END OF AGGREGATION", file=out)
    print(f"{'='*72}\n", file=out)


# =============================================================================
# SCANNING MODE — find CC roots in a bit range
# =============================================================================

def scan_for_chains(bits, min_chain, count, kind="first"):
    """Find CC roots by scanning a bit range."""
    if not HAS_GMPY2:
        print("ERROR: scan mode requires gmpy2", file=sys.stderr)
        return []
    
    start = mpz(1) << (bits - 1)
    end = mpz(1) << bits
    
    roots_found = []
    p = int(next_prime(start))
    tested = 0
    
    print(f"Scanning {bits}-bit range for CC{min_chain}+ {kind}-kind roots...", file=sys.stderr)
    
    while p < end and len(roots_found) < count:
        tested += 1
        if tested % 100000 == 0:
            print(f"  tested {tested}, found {len(roots_found)}...", file=sys.stderr)
        
        if kind == "first":
            if p % 6 == 5 and is_root_first(p):
                chain = follow_chain_first(p)
                if len(chain) >= min_chain:
                    roots_found.append((p, len(chain), kind))
                    print(f"  CC{len(chain)} at {hex(p)} ({p.bit_length()} bits)", file=sys.stderr)
        else:
            if p % 6 == 1 and is_root_second(p):
                chain = follow_chain_second(p)
                if len(chain) >= min_chain:
                    roots_found.append((p, len(chain), kind))
                    print(f"  CC{len(chain)} at {hex(p)} ({p.bit_length()} bits)", file=sys.stderr)
        
        p = int(next_prime(mpz(p)))
    
    print(f"Scan complete: tested {tested}, found {len(roots_found)} roots", file=sys.stderr)
    return roots_found


# =============================================================================
# TEST MODE
# =============================================================================

def run_tests():
    """Test with known CC chains."""
    print("=" * 60)
    print("  CC FINGERPRINT ANALYZER — TEST MODE")
    print("=" * 60)
    
    # Known CC chains (first-kind)
    test_chains = [
        (1122659, "CC7"),      # Well-known CC7
        (19099919, "CC8"),     # CC8
        (85864769, "CC9"),     # CC9
        (26089808579, "CC10"), # CC10
    ]
    
    all_rows = []
    chain_id = 0
    
    for root, label in test_chains:
        print(f"\n--- {label} root={root} ---")
        chain = follow_chain_first(root)
        print(f"  Chain length: {len(chain)}")
        print(f"  Chain: {' → '.join(str(p) for p in chain[:6])}{'...' if len(chain) > 6 else ''}")
        
        fps = analyze_chain(root, "first")
        
        # Show root fingerprint
        rfp = fps[0]
        print(f"  Root bits: {rfp['bits']}")
        print(f"  Trailing ones: {rfp['trailing_ones']}")
        print(f"  v2(p-1)={rfp['v2_pm1']}, v3(p-1)={rfp['v3_pm1']}, v5(p-1)={rfp['v5_pm1']}")
        print(f"  p-1 = {rfp['pm1_factors']}")
        print(f"  p+1 = {rfp['pp1_factors']}")
        print(f"  lpf(p-1): {rfp['lpf_pm1_bits']}-bit, lpf(p+1): {rfp['lpf_pp1_bits']}-bit")
        print(f"  Min kill pos: {rfp['min_kill_pos']}")
        print(f"  Immunized ({rfp['num_immunized']}): {rfp['immunized_primes']}")
        
        rows = chain_to_csv_rows(chain_id, fps)
        all_rows.extend(rows)
        chain_id += 1
    
    # Write test CSV
    outpath = "test_fingerprints.csv"
    with open(outpath, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        for row in all_rows:
            writer.writerow({k: row.get(k, "") for k in CSV_COLUMNS})
    
    print(f"\n  CSV written to {outpath}")
    
    # Run aggregation on test data
    print("\n" + "=" * 60)
    print("  AGGREGATION TEST")
    print("=" * 60)
    aggregate_csv(outpath)
    
    print("\n  All tests passed.")


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="CC Chain Factorization Fingerprint Analyzer",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 cc_fingerprint.py --test
  python3 cc_fingerprint.py --analyze roots.txt --output fp.csv
  python3 cc_fingerprint.py --aggregate fp.csv --output summary.txt
  python3 cc_fingerprint.py --scan --bits 89 --min-chain 8 --count 50
        """
    )
    
    parser.add_argument("--test", action="store_true", help="Run tests with known chains")
    parser.add_argument("--analyze", metavar="FILE", help="Analyze roots from file (one per line)")
    parser.add_argument("--aggregate", metavar="CSV", help="Aggregate existing fingerprint CSV")
    parser.add_argument("--scan", action="store_true", help="Scan a bit range for CC roots")
    parser.add_argument("--output", "-o", metavar="FILE", help="Output file")
    parser.add_argument("--bits", type=int, default=89, help="Bit size for scan (default: 89)")
    parser.add_argument("--min-chain", type=int, default=8, help="Minimum chain length (default: 8)")
    parser.add_argument("--count", type=int, default=50, help="Number of roots to find (default: 50)")
    parser.add_argument("--kind", choices=["first", "second"], default="first", help="Chain kind")
    
    args = parser.parse_args()
    
    if args.test:
        run_tests()
        return
    
    if args.aggregate:
        out = open(args.output, 'w') if args.output else sys.stdout
        aggregate_csv(args.aggregate, out)
        if args.output:
            out.close()
            print(f"Aggregation written to {args.output}", file=sys.stderr)
        return
    
    if args.analyze:
        # Read roots from file
        roots = []
        with open(args.analyze, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                if line.startswith("0x") or line.startswith("0X"):
                    roots.append(int(line, 16))
                else:
                    # Try to extract number from various formats
                    # Handle "CC13: 0x1abc..." format
                    parts = line.split()
                    for part in parts:
                        part = part.strip(",:")
                        try:
                            if part.startswith("0x"):
                                roots.append(int(part, 16))
                                break
                            elif part.isdigit():
                                roots.append(int(part))
                                break
                        except ValueError:
                            continue
        
        if not roots:
            print("ERROR: No valid roots found in file", file=sys.stderr)
            return
        
        print(f"Analyzing {len(roots)} roots...", file=sys.stderr)
        
        all_rows = []
        for chain_id, root in enumerate(roots):
            print(f"  [{chain_id+1}/{len(roots)}] root={hex(root)} ({root.bit_length()} bits)", file=sys.stderr)
            fps = analyze_chain(root, args.kind)
            rows = chain_to_csv_rows(chain_id, fps)
            all_rows.extend(rows)
        
        outpath = args.output or "fingerprints.csv"
        with open(outpath, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
            writer.writeheader()
            for row in all_rows:
                writer.writerow({k: row.get(k, "") for k in CSV_COLUMNS})
        
        print(f"CSV written to {outpath} ({len(all_rows)} rows)", file=sys.stderr)
        return
    
    if args.scan:
        found = scan_for_chains(args.bits, args.min_chain, args.count, args.kind)
        
        if not found:
            print("No chains found", file=sys.stderr)
            return
        
        # Analyze all found roots
        all_rows = []
        for chain_id, (root, length, kind) in enumerate(found):
            print(f"  Fingerprinting CC{length} at {hex(root)}...", file=sys.stderr)
            fps = analyze_chain(root, kind)
            rows = chain_to_csv_rows(chain_id, fps)
            all_rows.extend(rows)
        
        outpath = args.output or f"fingerprints_{args.bits}bit.csv"
        with open(outpath, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
            writer.writeheader()
            for row in all_rows:
                writer.writerow({k: row.get(k, "") for k in CSV_COLUMNS})
        
        print(f"CSV written to {outpath} ({len(all_rows)} rows from {len(found)} chains)", file=sys.stderr)
        
        # Auto-aggregate
        print("\n", file=sys.stderr)
        aggregate_csv(outpath)
        return
    
    parser.print_help()


if __name__ == "__main__":
    main()
