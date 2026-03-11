#!/usr/bin/env python3
# Copyright (c) 2026 Nenad Mićić <nenad@micic.be>
# SPDX-License-Identifier: Apache-2.0
"""Validate Cunningham Chain (first kind) files.
For each CC_N file, checks every hex entry starts a chain of length >= N.
Uses deterministic Miller-Rabin (sufficient for numbers < 2^128).
"""

import os, sys, time

def is_prime_miller_rabin(n):
    """Deterministic Miller-Rabin for n < 2^128."""
    if n < 2:
        return False
    if n < 4:
        return True
    if n % 2 == 0:
        return False
    # For n < 3,317,044,064,679,887,385,961,981, these 13 bases suffice
    bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41]
    d = n - 1
    r = 0
    while d % 2 == 0:
        d >>= 1
        r += 1
    for a in bases:
        if a >= n:
            continue
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(r - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    return True

def cc_length(p):
    """Length of Cunningham chain first kind starting at p."""
    q = p
    L = 0
    while is_prime_miller_rabin(q):
        L += 1
        q = 2 * q + 1
    return L

DIR = "."  # directory containing CC_8.txt .. CC_16.txt
FILES = [8, 9, 10, 11, 12, 13, 14, 15, 16]

total_ok = 0
total_fail = 0
total_longer = 0

for expected in FILES:
    fname = os.path.join(DIR, f"all_CC{expected}.txt")
    print(f"=== Validating CC{expected} ===")
    ok = fail = longer = count = 0
    t0 = time.time()

    with open(fname) as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("0x") or s.startswith("0X"):
                s = s[2:]
            if not s:
                continue
            p = int(s, 16)
            if p <= 1:
                continue
            count += 1
            clen = cc_length(p)
            if clen >= expected:
                ok += 1
                if clen > expected:
                    longer += 1
                    print(f"  #{count}: chain={clen} (>{expected}!) p=0x{s}")
            else:
                fail += 1
                print(f"  FAIL #{count}: chain={clen} (expected >={expected}) p=0x{s}")
            if count % 50000 == 0:
                elapsed = time.time() - t0
                print(f"  ... {count} checked ({elapsed:.1f}s)")

    elapsed = time.time() - t0
    print(f"  Results: {count} checked, {ok} OK, {longer} LONGER, {fail} FAIL ({elapsed:.1f}s)\n")
    total_ok += ok
    total_fail += fail
    total_longer += longer

print(f"=== TOTAL: {total_ok} OK, {total_longer} longer-than-expected, {total_fail} FAIL ===")
if total_fail == 0:
    print("ALL CHAINS VALID!")
else:
    print(f"*** {total_fail} FAILURES DETECTED ***")
