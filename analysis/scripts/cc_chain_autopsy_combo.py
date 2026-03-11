#!/usr/bin/env python3
# Copyright (c) 2026 Nenad Mićić <nenad@micic.be>
# SPDX-License-Identifier: Apache-2.0
"""
cc_chain_autopsy_combo.py

Batch chain-autopsy analyzer for Cunningham chains.

What it does per input root:
- Normalizes to true root (optional non-root backward walk)
- Walks chain positions 0..target+lookahead
- Tracks first breaker and counterfactual next breakers
- Checks immunization persistence across positions
- Summarizes p+1 / p-1 construction signals (S-base coverage)
- Computes 2-adic/spine descriptors
- Optional partial factorization of root and breaker neighborhoods

Input supports report/log files and simple root lists.
"""

from __future__ import annotations

import argparse
import csv
import fnmatch
import json
import math
import pathlib
import random
import re
import sys
import time
from collections import Counter
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple


try:
    import gmpy2

    HAS_GMPY2 = True
except Exception:
    HAS_GMPY2 = False


ROOT_PATTERNS = [
    re.compile(r"CC\s*([0-9]{1,2})[ab]?\s*FOUND!?[^\n]*?(0x[0-9a-fA-F]+)", re.IGNORECASE),
    re.compile(r"CC\s*([0-9]{1,2})[ab]?\s*elem\b[^\n]*?root\s*:\s*(0x[0-9a-fA-F]+)", re.IGNORECASE),
    re.compile(r"NON-ROOT\s+CC\s*([0-9]{1,2})[ab]?\s+(0x[0-9a-fA-F]+)", re.IGNORECASE),
    re.compile(r"CC\s*([0-9]{1,2})[ab]?\s*:\s*(0x[0-9a-fA-F]+)", re.IGNORECASE),
]
HEX_RE = re.compile(r"0x[0-9a-fA-F]+")
DEC_RE = re.compile(r"\b[0-9]{6,}\b")


@dataclass
class RootEntry:
    n: int
    cc: Optional[int]
    path: str
    line_no: int


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Batch autopsy analyzer for CC roots")
    p.add_argument("paths", nargs="+", help="Input files/directories")
    p.add_argument("--kind", choices=["first", "second"], default="first")
    p.add_argument("--target", type=int, default=18, help="Target chain length")
    p.add_argument("--lookahead", type=int, default=4, help="Extra positions after target")
    p.add_argument("--next-breakers", type=int, default=3, help="Record up to this many breaker positions")
    p.add_argument("--primes", type=str, default="3,5,7,11,13,17,19,23,29", help="Immunization prime set")
    p.add_argument("--s-base", type=str, default="2,3,5,11,13,19", help="Base construction primes")
    p.add_argument("--min-cc", type=int, default=None, help="Only parse entries with CC >= min")
    p.add_argument("--allow-no-cc", action="store_true", help="Allow parsing plain hex/dec roots without CC tag")
    p.add_argument("--include-name", type=str, default="*", help="Comma-separated filename globs to include")
    p.add_argument("--exclude-name", type=str, default="", help="Comma-separated filename globs to exclude")
    p.add_argument("--max-roots", type=int, default=0, help="Stop after this many unique roots (0=all)")
    p.add_argument("--no-dedupe", action="store_true", help="Do not dedupe roots")
    p.add_argument("--no-normalize", action="store_true", help="Do not walk backward to true root")

    p.add_argument("--factor-depth", type=int, default=0, help="0=off, 1=root +/-1, 2=also breaker")
    p.add_argument("--factor-trial-limit", type=int, default=10000)
    p.add_argument("--factor-time-ms", type=int, default=30, help="Time budget per number")
    p.add_argument("--factor-prime-positions", type=int, default=1, help="Factor p+/-1 for first N prime positions")

    p.add_argument("--neighbors", action="store_true", help="Find nearest lower/upper prime (slower)")
    p.add_argument("--neighbor-max-steps", type=int, default=200000)

    p.add_argument("--tree-n", type=int, default=0, help="2-adic grid height N for column mapping (0=use root bits)")

    p.add_argument("--roots-csv", type=str, default="", help="Write per-root CSV")
    p.add_argument("--positions-csv", type=str, default="", help="Write per-position CSV")
    p.add_argument("--summary-json", type=str, default="", help="Write summary JSON")
    p.add_argument("--quiet", action="store_true")
    return p.parse_args()


def parse_int_tokens(s: str) -> List[int]:
    out: List[int] = []
    for hx in HEX_RE.findall(s):
        try:
            out.append(int(hx, 16))
        except ValueError:
            pass
    for d in DEC_RE.findall(s):
        try:
            out.append(int(d, 10))
        except ValueError:
            pass
    return out


def parse_prime_list(s: str) -> List[int]:
    out: List[int] = []
    for t in s.split(","):
        t = t.strip()
        if not t:
            continue
        q = int(t, 10)
        if q <= 1:
            raise ValueError(f"invalid prime {q}")
        out.append(q)
    if not out:
        raise ValueError("empty prime list")
    return out


def parse_globs(s: str, default_star: bool = False) -> List[str]:
    parts = [x.strip() for x in s.split(",") if x.strip()]
    if parts:
        return parts
    return ["*"] if default_star else []


def name_allowed(name: str, include_globs: List[str], exclude_globs: List[str]) -> bool:
    if not any(fnmatch.fnmatch(name, g) for g in include_globs):
        return False
    if any(fnmatch.fnmatch(name, g) for g in exclude_globs):
        return False
    return True


def iter_input_files(paths: Iterable[str], include_globs: List[str], exclude_globs: List[str]) -> Iterable[pathlib.Path]:
    for raw in paths:
        p = pathlib.Path(raw)
        if not p.exists():
            continue
        if p.is_file():
            if name_allowed(p.name, include_globs, exclude_globs):
                yield p
        elif p.is_dir():
            for f in p.rglob("*"):
                if f.is_file() and name_allowed(f.name, include_globs, exclude_globs):
                    yield f


def extract_entries_from_line(line: str, allow_no_cc: bool) -> List[Tuple[Optional[int], int]]:
    out: List[Tuple[Optional[int], int]] = []
    seen = set()

    for rx in ROOT_PATTERNS:
        for m in rx.finditer(line):
            cc = int(m.group(1))
            n = int(m.group(2), 16)
            key = (cc, n)
            if key in seen:
                continue
            seen.add(key)
            out.append((cc, n))

    if out:
        return out

    if allow_no_cc:
        for n in parse_int_tokens(line):
            out.append((None, n))

    return out


# ---------- primality ----------

def _mr_test(n: int, a: int, d: int, s: int) -> bool:
    x = pow(a, d, n)
    if x == 1 or x == n - 1:
        return True
    for _ in range(s - 1):
        x = (x * x) % n
        if x == n - 1:
            return True
    return False


def is_probable_prime(n: int) -> bool:
    if n < 2:
        return False
    if HAS_GMPY2:
        return bool(gmpy2.is_prime(n) > 0)

    small = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
    for p in small:
        if n == p:
            return True
        if n % p == 0:
            return False

    d = n - 1
    s = 0
    while d % 2 == 0:
        s += 1
        d //= 2

    bases = [2, 325, 9375, 28178, 450775, 9780504, 1795265022]
    if n.bit_length() > 64:
        bases += [3, 5, 7, 11, 13, 17, 19, 23, 29]

    for a in bases:
        a %= n
        if a in (0, 1):
            continue
        if not _mr_test(n, a, d, s):
            return False
    return True


def prev_prime(n: int, max_steps: int) -> Optional[int]:
    if n <= 2:
        return None
    c = n - 1 if (n - 1) & 1 else n - 2
    steps = 0
    while c >= 2 and steps < max_steps:
        if is_probable_prime(c):
            return c
        c -= 2
        steps += 1
    return None


def next_prime(n: int, max_steps: int) -> Optional[int]:
    c = n + 1 if (n + 1) & 1 else n + 2
    steps = 0
    while steps < max_steps:
        if is_probable_prime(c):
            return c
        c += 2
        steps += 1
    return None


# ---------- factorization ----------

_SMALL_PRIME_CACHE: Dict[int, List[int]] = {}


def small_primes_upto(limit: int) -> List[int]:
    limit = max(2, limit)
    if limit in _SMALL_PRIME_CACHE:
        return _SMALL_PRIME_CACHE[limit]
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[0:2] = b"\x00\x00"
    for i in range(2, int(limit ** 0.5) + 1):
        if sieve[i]:
            step = i
            start = i * i
            sieve[start:limit + 1:step] = b"\x00" * (((limit - start) // step) + 1)
    ps = [i for i, b in enumerate(sieve) if b]
    _SMALL_PRIME_CACHE[limit] = ps
    return ps


def pollard_rho_one(n: int, deadline: float, rng: random.Random) -> Optional[int]:
    if n % 2 == 0:
        return 2
    if n % 3 == 0:
        return 3
    if is_probable_prime(n):
        return n

    while time.time() < deadline:
        x = rng.randrange(2, n - 1)
        y = x
        c = rng.randrange(1, n - 1)
        d = 1
        while d == 1 and time.time() < deadline:
            x = (pow(x, 2, n) + c) % n
            y = (pow(y, 2, n) + c) % n
            y = (pow(y, 2, n) + c) % n
            d = math.gcd(abs(x - y), n)
        if d != 1 and d != n:
            return d
    return None


def factor_partial(n: int, trial_limit: int, time_ms: int) -> Dict[str, object]:
    n = abs(int(n))
    if n <= 1:
        return {"factors": [], "cofactor": 1, "cofactor_is_prime": True}

    start = time.time()
    deadline = start + max(0.001, time_ms / 1000.0)
    fac: Counter[int] = Counter()

    for p in small_primes_upto(trial_limit):
        if p * p > n:
            break
        while n % p == 0:
            fac[p] += 1
            n //= p

    if n == 1:
        return {
            "factors": sorted(fac.items()),
            "cofactor": 1,
            "cofactor_is_prime": True,
        }

    if is_probable_prime(n):
        fac[n] += 1
        return {
            "factors": sorted(fac.items()),
            "cofactor": 1,
            "cofactor_is_prime": True,
        }

    rng = random.Random((n ^ (n >> 31)) & ((1 << 64) - 1))
    stack = [n]
    residue = 1

    while stack and time.time() < deadline:
        x = stack.pop()
        if x == 1:
            continue
        if is_probable_prime(x):
            fac[x] += 1
            continue
        d = pollard_rho_one(x, deadline, rng)
        if d is None or d == 1 or d == x:
            residue *= x
            continue
        stack.append(d)
        stack.append(x // d)

    for x in stack:
        residue *= x

    co_is_prime = bool(residue == 1 or is_probable_prime(residue))
    if residue > 1 and co_is_prime:
        fac[residue] += 1
        residue = 1

    return {
        "factors": sorted(fac.items()),
        "cofactor": residue,
        "cofactor_is_prime": co_is_prime,
    }


def format_fact(res: Dict[str, object]) -> str:
    parts: List[str] = []
    for p, e in res["factors"]:
        if e == 1:
            parts.append(str(p))
        else:
            parts.append(f"{p}^{e}")
    c = int(res["cofactor"])
    if c > 1:
        tag = "P" if bool(res["cofactor_is_prime"]) else "C"
        parts.append(f"{tag}{c.bit_length()}b")
    return " * ".join(parts) if parts else "1"


def smallest_factor_hint(n: int, trial_limit: int = 10000) -> Optional[int]:
    n = abs(int(n))
    if n < 2:
        return None
    for p in small_primes_upto(trial_limit):
        if p * p > n:
            break
        if n % p == 0:
            return p
    return None


# ---------- chain helpers ----------

def v2(n: int) -> int:
    n = abs(int(n))
    if n == 0:
        return 0
    c = 0
    while (n & 1) == 0:
        n >>= 1
        c += 1
    return c


def odd_part(n: int) -> int:
    n = abs(int(n))
    while n > 0 and (n & 1) == 0:
        n >>= 1
    return n


def trailing_ones(n: int) -> int:
    return v2(n + 1)


def immune_residue(q: int, kind: str) -> int:
    return q - 1 if kind == "first" else 1


def predecessor(n: int, kind: str) -> Optional[int]:
    if kind == "first":
        if n < 3:
            return None
        if (n - 1) & 1:
            return None
        x = (n - 1) // 2
        return x if x >= 2 else None
    if n < 3:
        return None
    if (n + 1) & 1:
        return None
    x = (n + 1) // 2
    return x if x >= 2 else None


def next_term(n: int, kind: str) -> int:
    return 2 * n + 1 if kind == "first" else 2 * n - 1


def normalize_to_root(n: int, kind: str, max_back: int = 128) -> Tuple[int, int]:
    cur = int(n)
    steps = 0
    while steps < max_back:
        pred = predecessor(cur, kind)
        if pred is None:
            break
        if not is_probable_prime(pred):
            break
        cur = pred
        steps += 1
    return cur, steps


# ---------- analysis ----------

def analyze_root(
    n_input: int,
    cc_input: Optional[int],
    kind: str,
    target: int,
    lookahead: int,
    primes: List[int],
    s_base: List[int],
    normalize: bool,
    next_breakers: int,
    factor_depth: int,
    factor_trial_limit: int,
    factor_time_ms: int,
    factor_prime_positions: int,
    neighbors: bool,
    neighbor_max_steps: int,
    tree_n: int,
) -> Tuple[Dict[str, object], List[Dict[str, object]]]:
    root, back_steps = (normalize_to_root(n_input, kind) if normalize else (int(n_input), 0))

    span = max(target + lookahead, 1)
    terms: List[int] = []
    is_p: List[bool] = []
    t = root
    for i in range(span):
        terms.append(t)
        is_p.append(is_probable_prime(t))
        t = next_term(t, kind)

    chain_len_obs = 0
    for v in is_p:
        if v:
            chain_len_obs += 1
        else:
            break

    first_break_pos: Optional[int] = None
    for i, v in enumerate(is_p):
        if not v:
            first_break_pos = i
            break

    chain_censored = first_break_pos is None
    first_breaker = terms[first_break_pos] if first_break_pos is not None else None

    comp_positions = [i for i, v in enumerate(is_p) if not v]
    next_break_pos = comp_positions[: max(0, next_breakers)]
    next_break_terms = [terms[i] for i in next_break_pos]
    next_break_small_factors = [smallest_factor_hint(terms[i], factor_trial_limit) for i in next_break_pos]

    # Immunization persistence and per-position counts
    immune_counts: List[int] = []
    immune_lists: List[List[int]] = []
    per_q = {}
    for q in primes:
        flags = []
        ir = immune_residue(q, kind)
        for n in terms:
            flags.append((n % q) == ir)
        per_q[q] = flags

    for i in range(span):
        il = [q for q in primes if per_q[q][i]]
        immune_lists.append(il)
        immune_counts.append(len(il))

    root_immune = [q for q in primes if per_q[q][0]]
    persist_q = [q for q in root_immune if all(per_q[q])]
    broken_persist_q = [q for q in root_immune if not all(per_q[q])]

    # Construction coverage per position
    plus_cov = []
    minus_cov = []
    for n in terms:
        plus_cov.append(all((n + 1) % q == 0 for q in s_base))
        minus_cov.append(all((n - 1) % q == 0 for q in s_base))

    root_pm1 = root - 1
    root_pp1 = root + 1

    f_root_pm1 = None
    f_root_pp1 = None
    f_breaker = None

    if factor_depth >= 1:
        f_root_pm1 = factor_partial(root_pm1, factor_trial_limit, factor_time_ms)
        f_root_pp1 = factor_partial(root_pp1, factor_trial_limit, factor_time_ms)
    if factor_depth >= 2 and first_breaker is not None:
        f_breaker = factor_partial(first_breaker, factor_trial_limit, factor_time_ms)

    prime_pos_fact: Dict[int, Dict[str, object]] = {}
    if factor_depth >= 1 and factor_prime_positions > 0:
        cnt = 0
        for i in range(span):
            if is_p[i]:
                prime_pos_fact[i] = {
                    "pm1": factor_partial(terms[i] - 1, factor_trial_limit, factor_time_ms),
                    "pp1": factor_partial(terms[i] + 1, factor_trial_limit, factor_time_ms),
                }
                cnt += 1
                if cnt >= factor_prime_positions:
                    break

    if tree_n <= 0:
        tree_n = root.bit_length()
    row = root.bit_length()
    col = None
    if 1 <= row <= tree_n:
        col = (root & ((1 << (row - 1)) - 1)) << (tree_n - row)

    prev_p = None
    next_p = None
    if neighbors:
        prev_p = prev_prime(root, neighbor_max_steps)
        next_p = next_prime(root, neighbor_max_steps)

    first_break_small = smallest_factor_hint(first_breaker, factor_trial_limit) if first_breaker else None
    first_break_v2_p1 = v2(first_breaker + 1) if first_breaker is not None else None
    first_break_v2_m1 = v2(first_breaker - 1) if first_breaker is not None else None
    first_break_odd_p1 = odd_part(first_breaker + 1) if first_breaker is not None else None
    first_break_odd_m1 = odd_part(first_breaker - 1) if first_breaker is not None else None

    root_row = {
        "input_root_hex": hex(n_input),
        "cc_input": cc_input,
        "kind": kind,
        "normalized_root_hex": hex(root),
        "normalized_bits": root.bit_length(),
        "back_steps": back_steps,
        "target": target,
        "lookahead": lookahead,
        "span": span,
        "chain_len_observed": chain_len_obs,
        "reaches_target": int(chain_len_obs >= target),
        "chain_censored": int(chain_censored),
        "first_break_pos": first_break_pos,
        "first_breaker_hex": hex(first_breaker) if first_breaker is not None else "",
        "first_breaker_bits": first_breaker.bit_length() if first_breaker is not None else "",
        "first_breaker_small_factor": first_break_small if first_break_small is not None else "",
        "first_breaker_v2_p_plus_1": first_break_v2_p1 if first_break_v2_p1 is not None else "",
        "first_breaker_odd_core_p_plus_1": first_break_odd_p1 if first_break_odd_p1 is not None else "",
        "first_breaker_v2_p_minus_1": first_break_v2_m1 if first_break_v2_m1 is not None else "",
        "first_breaker_odd_core_p_minus_1": first_break_odd_m1 if first_break_odd_m1 is not None else "",
        "next_break_positions": ",".join(str(x) for x in next_break_pos),
        "next_break_terms_hex": ",".join(hex(x) for x in next_break_terms),
        "next_break_small_factors": ",".join("" if x is None else str(x) for x in next_break_small_factors),
        "immune_root": ",".join(str(q) for q in root_immune),
        "immune_persist_all_positions": ",".join(str(q) for q in persist_q),
        "immune_persist_broken": ",".join(str(q) for q in broken_persist_q),
        "immune_counts_by_pos": ",".join(str(x) for x in immune_counts),
        "sbase_plus_cover_root": int(plus_cov[0]),
        "sbase_minus_cover_root": int(minus_cov[0]),
        "v2_root_p_plus_1": v2(root + 1),
        "v2_root_p_minus_1": v2(root - 1),
        "odd_core_root_p_plus_1": odd_part(root + 1),
        "odd_core_root_p_minus_1": odd_part(root - 1),
        "spine_k": trailing_ones(root),
        "spine_m": (root + 1) >> trailing_ones(root),
        "tree_n": tree_n,
        "tree_row": row,
        "tree_col": col if col is not None else "",
        "neighbor_prev_prime": prev_p if prev_p is not None else "",
        "neighbor_next_prime": next_p if next_p is not None else "",
        "root_pm1_fact": format_fact(f_root_pm1) if f_root_pm1 else "",
        "root_pp1_fact": format_fact(f_root_pp1) if f_root_pp1 else "",
        "breaker_fact": format_fact(f_breaker) if f_breaker else "",
    }

    positions: List[Dict[str, object]] = []
    for i in range(span):
        n = terms[i]
        rowp = {
            "root_hex": hex(root),
            "cc_input": cc_input if cc_input is not None else "",
            "pos": i,
            "term_hex": hex(n),
            "term_bits": n.bit_length(),
            "is_prime": int(is_p[i]),
            "immune_count": immune_counts[i],
            "immune_list": ",".join(str(q) for q in immune_lists[i]),
            "sbase_plus_cover": int(plus_cov[i]),
            "sbase_minus_cover": int(minus_cov[i]),
            "v2_term_plus_1": v2(n + 1),
            "v2_term_minus_1": v2(n - 1),
        }
        if not is_p[i]:
            sf = smallest_factor_hint(n, factor_trial_limit)
            rowp["small_factor_hint"] = sf if sf is not None else ""
        else:
            rowp["small_factor_hint"] = ""

        if i in prime_pos_fact:
            rowp["term_pm1_fact"] = format_fact(prime_pos_fact[i]["pm1"])
            rowp["term_pp1_fact"] = format_fact(prime_pos_fact[i]["pp1"])
        else:
            rowp["term_pm1_fact"] = ""
            rowp["term_pp1_fact"] = ""

        positions.append(rowp)

    return root_row, positions


def print_summary(rows: List[Dict[str, object]], primes: List[int], s_base: List[int]) -> Dict[str, object]:
    n = len(rows)
    if n == 0:
        print("No roots analyzed.")
        return {"roots": 0}

    cc_vals = [r["cc_input"] for r in rows if isinstance(r["cc_input"], int)]
    chain_dist = Counter(int(r["chain_len_observed"]) for r in rows)
    break_pos_dist = Counter(int(r["first_break_pos"]) for r in rows if r["first_break_pos"] != "" and r["first_break_pos"] is not None)

    persist_hits = Counter()
    broken_hits = Counter()
    for r in rows:
        p_all = str(r["immune_persist_all_positions"]).strip()
        p_broken = str(r["immune_persist_broken"]).strip()
        if p_all:
            for tok in p_all.split(","):
                if tok:
                    persist_hits[int(tok)] += 1
        if p_broken:
            for tok in p_broken.split(","):
                if tok:
                    broken_hits[int(tok)] += 1

    plus_cover = sum(int(r["sbase_plus_cover_root"]) for r in rows)
    minus_cover = sum(int(r["sbase_minus_cover_root"]) for r in rows)

    print("=== Chain Autopsy Combo Summary ===")
    print(f"Roots analyzed: {n}")
    if cc_vals:
        print(f"Input CC range: CC{min(cc_vals)}..CC{max(cc_vals)}")
    print(f"Chain observed distribution: {dict(sorted(chain_dist.items()))}")
    print(f"First-break position distribution: {dict(sorted(break_pos_dist.items()))}")
    print(f"S-base ({s_base}) coverage at root: p+1={plus_cover}/{n} ({100.0*plus_cover/n:.2f}%), p-1={minus_cover}/{n} ({100.0*minus_cover/n:.2f}%)")

    print("Immunization persistence (root-immune and remains immune across scanned positions):")
    for q in primes:
        a = persist_hits[q]
        b = broken_hits[q]
        print(f"  q={q:>3}: persist={a:>6} broken={b:>6}")

    out = {
        "roots": n,
        "cc_min": min(cc_vals) if cc_vals else None,
        "cc_max": max(cc_vals) if cc_vals else None,
        "chain_dist": dict(sorted(chain_dist.items())),
        "first_break_pos_dist": dict(sorted(break_pos_dist.items())),
        "sbase_plus_cover": plus_cover,
        "sbase_minus_cover": minus_cover,
        "persist_hits": dict(sorted(persist_hits.items())),
        "broken_hits": dict(sorted(broken_hits.items())),
    }
    return out


def main() -> int:
    args = parse_args()

    try:
        primes = parse_prime_list(args.primes)
        s_base = parse_prime_list(args.s_base)
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2

    include_globs = parse_globs(args.include_name, default_star=True)
    exclude_globs = parse_globs(args.exclude_name, default_star=False)

    roots: List[RootEntry] = []
    seen = set()
    files_scanned = 0
    lines_scanned = 0

    for f in iter_input_files(args.paths, include_globs, exclude_globs):
        files_scanned += 1
        try:
            with f.open("r", encoding="utf-8", errors="ignore") as fp:
                for ln, line in enumerate(fp, start=1):
                    lines_scanned += 1
                    entries = extract_entries_from_line(line, allow_no_cc=args.allow_no_cc)
                    if args.min_cc is not None:
                        entries = [(cc, n) for (cc, n) in entries if cc is not None and cc >= args.min_cc]
                    if not entries:
                        continue

                    for cc, n in entries:
                        if not args.no_dedupe:
                            if n in seen:
                                continue
                            seen.add(n)
                        roots.append(RootEntry(n=n, cc=cc, path=str(f), line_no=ln))
                        if args.max_roots > 0 and len(roots) >= args.max_roots:
                            break
                    if args.max_roots > 0 and len(roots) >= args.max_roots:
                        break
            if args.max_roots > 0 and len(roots) >= args.max_roots:
                break
        except Exception:
            continue

    if not roots:
        print("No roots found. Check include/exclude filters, min-cc, and paths.")
        return 1

    if not args.quiet:
        print(f"Files scanned: {files_scanned}")
        print(f"Lines scanned: {lines_scanned}")
        print(f"Roots queued:  {len(roots)}")
        print(f"Primality backend: {'gmpy2' if HAS_GMPY2 else 'python-mr'}")

    root_rows: List[Dict[str, object]] = []
    pos_rows: List[Dict[str, object]] = []

    t0 = time.time()
    for i, e in enumerate(roots, start=1):
        rr, pr = analyze_root(
            n_input=e.n,
            cc_input=e.cc,
            kind=args.kind,
            target=args.target,
            lookahead=args.lookahead,
            primes=primes,
            s_base=s_base,
            normalize=not args.no_normalize,
            next_breakers=args.next_breakers,
            factor_depth=args.factor_depth,
            factor_trial_limit=args.factor_trial_limit,
            factor_time_ms=args.factor_time_ms,
            factor_prime_positions=args.factor_prime_positions,
            neighbors=args.neighbors,
            neighbor_max_steps=args.neighbor_max_steps,
            tree_n=args.tree_n,
        )
        rr["source_file"] = e.path
        rr["source_line"] = e.line_no
        root_rows.append(rr)
        pos_rows.extend(pr)

        if (i % 200 == 0 or i == len(roots)) and not args.quiet:
            dt = time.time() - t0
            print(f"Processed {i}/{len(roots)} roots in {dt:.1f}s")

    summary = print_summary(root_rows, primes, s_base)

    if args.roots_csv:
        fieldnames = sorted({k for r in root_rows for k in r.keys()})
        with open(args.roots_csv, "w", newline="", encoding="utf-8") as fp:
            w = csv.DictWriter(fp, fieldnames=fieldnames)
            w.writeheader()
            for r in root_rows:
                w.writerow(r)
        if not args.quiet:
            print(f"Wrote roots CSV: {args.roots_csv}")

    if args.positions_csv:
        fieldnames = sorted({k for r in pos_rows for k in r.keys()})
        with open(args.positions_csv, "w", newline="", encoding="utf-8") as fp:
            w = csv.DictWriter(fp, fieldnames=fieldnames)
            w.writeheader()
            for r in pos_rows:
                w.writerow(r)
        if not args.quiet:
            print(f"Wrote positions CSV: {args.positions_csv}")

    if args.summary_json:
        with open(args.summary_json, "w", encoding="utf-8") as fp:
            json.dump(summary, fp, indent=2, sort_keys=True)
        if not args.quiet:
            print(f"Wrote summary JSON: {args.summary_json}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
