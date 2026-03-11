#!/usr/bin/env python3
# Copyright (c) 2026 Nenad Mićić <nenad@micic.be>
# SPDX-License-Identifier: Apache-2.0
"""
Analyze immunization patterns from CC logs.

Primary dataset (`paths`) is typically positives (e.g. CC10+ roots).
Optional control dataset (`--control-paths`) is typically all survivors.

The script reports:
- overall immunization statistics
- per-CC statistics
- per-bit-length statistics
- optional control-vs-primary lift tables
"""

from __future__ import annotations

import argparse
import fnmatch
import pathlib
import re
import sys
from collections import Counter
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Tuple


HEX_RE = re.compile(r"0x[0-9a-fA-F]+")
ROOT_PATTERNS = [
    # Example: "*** CC7a FOUND! 0x..."
    re.compile(
        r"CC\s*([0-9]{1,2})[ab]?\s*FOUND!?[^\n]*?(0x[0-9a-fA-F]+)",
        re.IGNORECASE,
    ),
    # Example: "[NON-ROOT] CC7a elem, root: 0x..."
    re.compile(
        r"CC\s*([0-9]{1,2})[ab]?\s*elem\b[^\n]*?root\s*:\s*(0x[0-9a-fA-F]+)",
        re.IGNORECASE,
    ),
    # Example (machine-parsable): "NON-ROOT CC17 0x... 72"
    re.compile(
        r"NON-ROOT\s+CC\s*([0-9]{1,2})[ab]?\s+(0x[0-9a-fA-F]+)",
        re.IGNORECASE,
    ),
    # Example: "CC 8a: 0x..."
    re.compile(
        r"CC\s*([0-9]{1,2})[ab]?\s*:\s*(0x[0-9a-fA-F]+)",
        re.IGNORECASE,
    ),
]


@dataclass
class GroupStats:
    per_prime_hits: List[int]
    count: int = 0
    imm_sum: int = 0
    bit_min: Optional[int] = None
    bit_max: Optional[int] = None

    def add(self, bits: int, imm: int, mask: int) -> None:
        self.count += 1
        self.imm_sum += imm
        if self.bit_min is None or bits < self.bit_min:
            self.bit_min = bits
        if self.bit_max is None or bits > self.bit_max:
            self.bit_max = bits
        for i in range(len(self.per_prime_hits)):
            if (mask >> i) & 1:
                self.per_prime_hits[i] += 1


@dataclass
class AnalysisResult:
    files_scanned: int = 0
    lines_scanned: int = 0
    lines_kept: int = 0
    roots_seen: int = 0
    roots_kept: int = 0
    dup_skipped: int = 0
    cc_min_seen: Optional[int] = None
    cc_max_seen: Optional[int] = None
    per_prime_hits: List[int] = field(default_factory=list)
    immunized_count_dist: Counter[int] = field(default_factory=Counter)
    missing_combo_dist: Counter[int] = field(default_factory=Counter)
    cc_dist: Counter[int] = field(default_factory=Counter)
    bit_dist: Counter[int] = field(default_factory=Counter)
    by_cc: Dict[int, GroupStats] = field(default_factory=dict)
    by_bits: Dict[int, GroupStats] = field(default_factory=dict)
    by_cc_bits: Dict[Tuple[int, int], GroupStats] = field(default_factory=dict)
    require_pass: int = 0


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Analyze CC immunization patterns from logs.")
    p.add_argument(
        "paths",
        nargs="+",
        help="Primary input files or directories (recursively scanned)",
    )
    p.add_argument(
        "--control-paths",
        nargs="+",
        default=None,
        help="Optional control dataset (e.g. all survivors)",
    )
    p.add_argument(
        "--min-cc",
        type=int,
        default=10,
        help="Primary dataset: only keep lines with CC >= this value (default: 10)",
    )
    p.add_argument(
        "--control-min-cc",
        type=int,
        default=None,
        help="Control dataset: if set, only keep lines with CC >= this value",
    )
    p.add_argument(
        "--primes",
        type=str,
        default="3,5,7,11,13,17,19,23,29",
        help="Comma-separated prime set for immunization analysis",
    )
    p.add_argument(
        "--require-immunized",
        type=int,
        default=None,
        help="Report pass-rate for entries with at least this many immunized primes",
    )
    p.add_argument(
        "--top-combos",
        type=int,
        default=15,
        help="Show top N missing-prime combos (default: 15)",
    )
    p.add_argument(
        "--allow-no-cc",
        action="store_true",
        help="Primary dataset: also parse lines with hex values but no CC tag",
    )
    p.add_argument(
        "--no-dedupe",
        action="store_true",
        help="Primary dataset: do not dedupe repeated roots",
    )
    p.add_argument(
        "--control-no-dedupe",
        action="store_true",
        help="Control dataset: do not dedupe repeated roots",
    )
    p.add_argument(
        "--group-min-count",
        type=int,
        default=20,
        help="Minimum count for per-group rows (default: 20)",
    )
    p.add_argument(
        "--cc-top",
        type=int,
        default=40,
        help="Max rows in per-CC section (default: 40)",
    )
    p.add_argument(
        "--bits-top",
        type=int,
        default=40,
        help="Max rows in per-bit section (default: 40)",
    )
    p.add_argument(
        "--ccbit-top",
        type=int,
        default=30,
        help="Max rows in per-(CC,bit) section (default: 30)",
    )
    p.add_argument(
        "--include-name",
        type=str,
        default="*",
        help="Comma-separated filename globs to include (default: *)",
    )
    p.add_argument(
        "--exclude-name",
        type=str,
        default="",
        help="Comma-separated filename globs to exclude",
    )
    return p.parse_args()


def parse_prime_list(s: str) -> List[int]:
    out: List[int] = []
    for tok in s.split(","):
        tok = tok.strip()
        if not tok:
            continue
        try:
            q = int(tok, 10)
        except ValueError:
            raise ValueError(f"Invalid prime token: {tok!r}") from None
        if q <= 1:
            raise ValueError(f"Prime must be > 1: {q}")
        out.append(q)
    if not out:
        raise ValueError("Prime list is empty")
    return out


def parse_glob_list(s: str) -> List[str]:
    out = [x.strip() for x in s.split(",") if x.strip()]
    return out if out else ["*"]


def name_allowed(name: str, include_globs: List[str], exclude_globs: List[str]) -> bool:
    if not any(fnmatch.fnmatch(name, pat) for pat in include_globs):
        return False
    if any(fnmatch.fnmatch(name, pat) for pat in exclude_globs):
        return False
    return True


def iter_input_files(
    paths: Iterable[str], include_globs: List[str], exclude_globs: List[str]
) -> Iterable[pathlib.Path]:
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


def format_missing_combo(primes: List[int], mask: int) -> str:
    missing = [str(primes[i]) for i in range(len(primes)) if (mask >> i) & 1 == 0]
    return "{" + ",".join(missing) + "}" if missing else "{}"


def bitcount(x: int) -> int:
    return x.bit_count()


def extract_entries_from_line(line: str, allow_no_cc: bool) -> List[Tuple[Optional[int], str]]:
    out: List[Tuple[Optional[int], str]] = []
    seen = set()

    for rx in ROOT_PATTERNS:
        for m in rx.finditer(line):
            cc = int(m.group(1))
            hx = m.group(2).lower()
            key = (cc, hx)
            if key in seen:
                continue
            seen.add(key)
            out.append((cc, hx))

    if out:
        return out

    if allow_no_cc:
        for hx in HEX_RE.findall(line):
            out.append((None, hx.lower()))

    return out


def get_group(d: Dict, key, nprimes: int) -> GroupStats:
    g = d.get(key)
    if g is None:
        g = GroupStats(per_prime_hits=[0] * nprimes)
        d[key] = g
    return g


def analyze_dataset(
    paths: Iterable[str],
    primes: List[int],
    min_cc: Optional[int],
    allow_no_cc: bool,
    dedupe: bool,
    require_k: Optional[int],
    include_globs: List[str],
    exclude_globs: List[str],
) -> AnalysisResult:
    res = AnalysisResult(per_prime_hits=[0] * len(primes))
    seen_roots = set()

    for f in iter_input_files(paths, include_globs, exclude_globs):
        res.files_scanned += 1
        try:
            with f.open("r", encoding="utf-8", errors="ignore") as fp:
                for line in fp:
                    res.lines_scanned += 1
                    entries = extract_entries_from_line(line, allow_no_cc=allow_no_cc)
                    if not entries:
                        continue

                    if min_cc is not None:
                        entries = [(cc, hx) for (cc, hx) in entries if cc is not None and cc >= min_cc]
                    if not entries:
                        continue

                    res.lines_kept += 1

                    for cc_line, hx in entries:
                        res.roots_seen += 1
                        n = int(hx, 16)

                        if dedupe:
                            if n in seen_roots:
                                res.dup_skipped += 1
                                continue
                            seen_roots.add(n)

                        res.roots_kept += 1

                        n1 = n + 1
                        mask = 0
                        for i, q in enumerate(primes):
                            if n1 % q == 0:
                                mask |= (1 << i)
                                res.per_prime_hits[i] += 1

                        imm = bitcount(mask)
                        bits = n.bit_length()

                        res.immunized_count_dist[imm] += 1
                        res.missing_combo_dist[mask] += 1
                        res.bit_dist[bits] += 1

                        if cc_line is not None:
                            res.cc_min_seen = cc_line if res.cc_min_seen is None else min(res.cc_min_seen, cc_line)
                            res.cc_max_seen = cc_line if res.cc_max_seen is None else max(res.cc_max_seen, cc_line)
                            res.cc_dist[cc_line] += 1
                            get_group(res.by_cc, cc_line, len(primes)).add(bits, imm, mask)
                            get_group(res.by_cc_bits, (cc_line, bits), len(primes)).add(bits, imm, mask)

                        get_group(res.by_bits, bits, len(primes)).add(bits, imm, mask)

                        if require_k is not None and imm >= require_k:
                            res.require_pass += 1
        except (OSError, UnicodeError):
            continue

    return res


def pct(a: int, b: int) -> float:
    return 0.0 if b <= 0 else 100.0 * a / b


def safe_ratio(a: float, b: float) -> str:
    if b == 0.0:
        return "inf" if a > 0.0 else "-"
    return f"{a / b:.2f}x"


def choose_variable_prime_indices(res: AnalysisResult, primes: List[int]) -> List[int]:
    if res.roots_kept == 0:
        return list(range(min(3, len(primes))))
    idx = []
    for i in range(len(primes)):
        r = res.per_prime_hits[i] / res.roots_kept
        if 0.001 < r < 0.999:
            idx.append(i)
    if idx:
        return idx
    return list(range(min(3, len(primes))))


def format_prime_rates(per_prime_hits: List[int], total: int, primes: List[int], idxs: List[int]) -> str:
    parts = []
    for i in idxs:
        parts.append(f"q{primes[i]}={pct(per_prime_hits[i], total):.2f}%")
    return " ".join(parts)


def print_overall(
    title: str,
    res: AnalysisResult,
    primes: List[int],
    min_cc: Optional[int],
    dedupe: bool,
    require_k: Optional[int],
    top_combos: int,
) -> None:
    print(f"=== {title} ===")
    print(f"Files scanned:         {res.files_scanned}")
    print(f"Lines scanned:         {res.lines_scanned}")
    print(f"Lines kept:            {res.lines_kept}")
    print(f"Roots parsed:          {res.roots_seen}")
    print(f"Roots analyzed:        {res.roots_kept}")
    if dedupe:
        print(f"Duplicate roots skip:  {res.dup_skipped}")
    if min_cc is None:
        print("Min CC filter:         disabled")
    else:
        print(f"Min CC filter:         CC{min_cc}+")
    print(f"Prime set:             {primes}")
    if res.cc_min_seen is not None and res.cc_max_seen is not None:
        print(f"CC seen in input:      min=CC{res.cc_min_seen} max=CC{res.cc_max_seen}")

    if res.roots_kept == 0:
        print("No matching roots found.")
        return

    print("\nPer-prime immunization rate (q | n+1):")
    for i, q in enumerate(primes):
        print(
            f"  q={q:>3}: {res.per_prime_hits[i]:>10} / {res.roots_kept:<10} "
            f"({pct(res.per_prime_hits[i], res.roots_kept):6.2f}%)"
        )

    print("\nImmunized-count distribution:")
    for imm in sorted(res.immunized_count_dist):
        cnt = res.immunized_count_dist[imm]
        print(
            f"  {imm:>2}/{len(primes)} immunized: {cnt:>10} "
            f"({pct(cnt, res.roots_kept):6.2f}%)"
        )

    if require_k is not None:
        print(
            f"\nAt least {require_k}/{len(primes)} immunized: "
            f"{res.require_pass}/{res.roots_kept} ({pct(res.require_pass, res.roots_kept):.2f}%)"
        )

    print(f"\nTop {top_combos} missing-prime combos ({{}} means none missing):")
    for mask, cnt in res.missing_combo_dist.most_common(max(1, top_combos)):
        print(
            f"  {format_missing_combo(primes, mask):<36} {cnt:>10} "
            f"({pct(cnt, res.roots_kept):6.2f}%)"
        )


def print_group_stats(
    title: str,
    groups: List[Tuple[str, GroupStats]],
    primes: List[int],
    var_idxs: List[int],
    limit: int,
    min_count: int,
) -> None:
    rows = [(name, g) for name, g in groups if g.count >= min_count]
    if not rows:
        return

    if limit > 0:
        rows = rows[:limit]

    print(f"\n{title} (min count {min_count}):")
    for name, g in rows:
        avg_imm = g.imm_sum / g.count
        rate_str = format_prime_rates(g.per_prime_hits, g.count, primes, var_idxs)
        bmin = g.bit_min if g.bit_min is not None else 0
        bmax = g.bit_max if g.bit_max is not None else 0
        print(
            f"  {name:<14} n={g.count:<8} avgImm={avg_imm:5.2f} "
            f"bits=[{bmin},{bmax}] {rate_str}"
        )


def print_breakdowns(
    res: AnalysisResult,
    primes: List[int],
    cc_top: int,
    bits_top: int,
    ccbit_top: int,
    min_count: int,
) -> None:
    if res.roots_kept == 0:
        return

    var_idxs = choose_variable_prime_indices(res, primes)

    cc_rows = [(f"CC{cc}", g) for cc, g in sorted(res.by_cc.items(), key=lambda x: x[0])]
    print_group_stats("Per-CC summary", cc_rows, primes, var_idxs, cc_top, min_count)

    bit_rows = [(f"{bits}-bit", g) for bits, g in sorted(res.by_bits.items(), key=lambda x: x[0])]
    print_group_stats("Per-bit summary", bit_rows, primes, var_idxs, bits_top, min_count)

    ccbit_rows = []
    for (cc, bits), g in sorted(res.by_cc_bits.items(), key=lambda kv: kv[1].count, reverse=True):
        ccbit_rows.append((f"CC{cc}@{bits}b", g))
    print_group_stats("Top (CC,bit) cells", ccbit_rows, primes, var_idxs, ccbit_top, min_count)


def print_comparison(pos: AnalysisResult, ctrl: AnalysisResult, primes: List[int]) -> None:
    if pos.roots_kept == 0 or ctrl.roots_kept == 0:
        return

    print("\n=== Primary vs Control Comparison ===")
    print(f"Primary roots: {pos.roots_kept}")
    print(f"Control roots: {ctrl.roots_kept}")

    print("\nPer-prime rate and lift (primary/control):")
    for i, q in enumerate(primes):
        rp = pos.per_prime_hits[i] / pos.roots_kept
        rc = ctrl.per_prime_hits[i] / ctrl.roots_kept
        print(
            f"  q={q:>3}: pos={100.0*rp:6.2f}% ctrl={100.0*rc:6.2f}% "
            f"lift={safe_ratio(rp, rc)}"
        )

    print("\nImmunized-count distribution lift:")
    all_k = sorted(set(pos.immunized_count_dist.keys()) | set(ctrl.immunized_count_dist.keys()))
    for k in all_k:
        p_cnt = pos.immunized_count_dist.get(k, 0)
        c_cnt = ctrl.immunized_count_dist.get(k, 0)
        rp = p_cnt / pos.roots_kept
        rc = c_cnt / ctrl.roots_kept
        print(
            f"  {k:>2}/{len(primes)}: pos={100.0*rp:6.2f}% ctrl={100.0*rc:6.2f}% "
            f"lift={safe_ratio(rp, rc)}"
        )

    print("\nBit-distribution lift (share in primary vs control):")
    all_bits = sorted(set(pos.bit_dist.keys()) | set(ctrl.bit_dist.keys()))
    shown = 0
    for b in all_bits:
        p_cnt = pos.bit_dist.get(b, 0)
        c_cnt = ctrl.bit_dist.get(b, 0)
        if p_cnt == 0 and c_cnt == 0:
            continue
        sp = p_cnt / pos.roots_kept
        sc = c_cnt / ctrl.roots_kept
        print(
            f"  {b:>3}-bit: pos_share={100.0*sp:6.2f}% ctrl_share={100.0*sc:6.2f}% "
            f"lift={safe_ratio(sp, sc)}"
        )
        shown += 1
        if shown >= 32:
            break


def main() -> int:
    args = parse_args()

    try:
        primes = parse_prime_list(args.primes)
    except ValueError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2

    require_k = args.require_immunized
    if require_k is not None and (require_k < 0 or require_k > len(primes)):
        print(
            f"ERROR: --require-immunized must be in [0,{len(primes)}], got {require_k}",
            file=sys.stderr,
        )
        return 2

    include_globs = parse_glob_list(args.include_name)
    exclude_globs = [x.strip() for x in args.exclude_name.split(",") if x.strip()]

    pos = analyze_dataset(
        paths=args.paths,
        primes=primes,
        min_cc=args.min_cc,
        allow_no_cc=args.allow_no_cc,
        dedupe=not args.no_dedupe,
        require_k=require_k,
        include_globs=include_globs,
        exclude_globs=exclude_globs,
    )

    print_overall(
        "Primary Dataset",
        pos,
        primes,
        min_cc=args.min_cc,
        dedupe=not args.no_dedupe,
        require_k=require_k,
        top_combos=args.top_combos,
    )
    print_breakdowns(
        pos,
        primes,
        cc_top=args.cc_top,
        bits_top=args.bits_top,
        ccbit_top=args.ccbit_top,
        min_count=max(1, args.group_min_count),
    )

    if args.control_paths:
        ctrl_allow_no_cc = args.control_min_cc is None
        ctrl = analyze_dataset(
            paths=args.control_paths,
            primes=primes,
            min_cc=args.control_min_cc,
            allow_no_cc=ctrl_allow_no_cc,
            dedupe=not args.control_no_dedupe,
            require_k=None,
            include_globs=include_globs,
            exclude_globs=exclude_globs,
        )

        print()
        print_overall(
            "Control Dataset",
            ctrl,
            primes,
            min_cc=args.control_min_cc,
            dedupe=not args.control_no_dedupe,
            require_k=None,
            top_combos=args.top_combos,
        )
        print_breakdowns(
            ctrl,
            primes,
            cc_top=args.cc_top,
            bits_top=args.bits_top,
            ccbit_top=args.ccbit_top,
            min_count=max(1, args.group_min_count),
        )
        print_comparison(pos, ctrl, primes)

    if pos.roots_kept == 0:
        print("Tip: verify --min-cc, --allow-no-cc, and input paths.")
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
