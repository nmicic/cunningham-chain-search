#!/usr/bin/env python3
# Copyright (c) 2026 Nenad Mićić <nenad@micic.be>
# SPDX-License-Identifier: Apache-2.0
"""
generate_2adic_tree_v5.py — Enhanced 2-Adic Tree (row/column layout)
=====================================================================
Generates the classic bottom-up 2-adic tree where each column is a
"ray" (n, n/2, n/4, …) and adds rich colouring modes from the
p+1 / Cunningham-chain analysis family.

Layout
------
  Bottom row contains 2^(rows-1) consecutive integers.
  Every even number n has n/2 directly above it in the same column.
  Odd numbers have empty cells above them.

  For rows=6:  bottom row = 32…63  (32 columns)
               top cell in col 0:  1  (col 0 = ray of 1)
               top cell in col 2:  17 (col 2 = ray of 17)

Colour modes  (--color / -c)
-------------
  classic   Black=prime, gray=composite, red=Mersenne, blue=Fermat,
            green=power-of-2, orange=1  [default, light background]

  dark      Same classes, dark background (like the HTML explorer)

  factor    For primes: solid colour from smallest of {3,5,7,11,13,17,
            19,23,29} that divides p+1.  White = p+1 is a pure power
            of 2 (Mersenne-like).  Composites: grey.

  modp      All cells coloured by  n mod P  rainbow.
            Choose P with --modprime (default 7).

  chain     Cunningham-chain membership: CC1 start=crimson,
            CC1 member=red, CC2 start=blue, CC2 member=steel-blue,
            safe prime=green, plain prime=dark, composite=grey.

  residues  Each cell gets up to 4 coloured border-lines showing which
            of {3,5,7,11,13} divide n (SVG only; PNG uses tinted cell).

  p1grid    Inside each cell draw a small 3×3 grid of squares showing
            which of the 9 small primes divide p+1  (for primes) or n
            (for composites).  Works best with --cell ≥ 40.

Output formats (--format)
-------------------------
  svg       Scalable vector, labelled numbers, best for rows ≤ 10.
  png       Raster, works for any size; --cell sets pixels per cell.
  both      Both SVG and PNG.

Usage examples
--------------
  python generate_2adic_tree_v5.py                        # rows=6, classic SVG
  python generate_2adic_tree_v5.py 8 --color factor       # p+1 factor coloring
  python generate_2adic_tree_v5.py 7 --color chain        # CC chain coloring
  python generate_2adic_tree_v5.py 6 --color modp --modprime 11 --format both
  python generate_2adic_tree_v5.py 8 --color residues --format svg
  python generate_2adic_tree_v5.py 7 --color p1grid --cell 50 --format png
  python generate_2adic_tree_v5.py 10 --color factor --format png --cell 8
  python generate_2adic_tree_v5.py --all-colors 7 --format both  # all modes
  python generate_2adic_tree_v5.py 7 --color chain --no-chain-lines  # without CC lines
"""

import argparse
import math
import os
import time
from typing import Optional

try:
    from PIL import Image, ImageDraw, ImageFont
    HAS_PIL = True
except ImportError:
    HAS_PIL = False

# ─────────────────────────────────────────────────────────────────────────────
# PALETTE
# ─────────────────────────────────────────────────────────────────────────────

NINE_PRIMES = [3, 5, 7, 11, 13, 17, 19, 23, 29]
NINE_RGB = [
    (255, 140,   0),   # 3   orange
    (255, 215,   0),   # 5   gold
    (  0, 210,  90),   # 7   green
    (  0, 200, 255),   # 11  cyan
    (140, 110, 255),   # 13  violet
    (220,  55, 255),   # 17  purple
    (255,  55, 170),   # 19  pink
    (255,  45,  75),   # 23  red
    ( 70, 255, 195),   # 29  mint
]
NINE_HEX = ['#ff8c00','#ffd700','#00d25a','#00c8ff',
            '#8c6eff','#dc37ff','#ff37aa','#ff2d4b','#46ffc3']
NINE_LABELS = ['3','5','7','11','13','17','19','23','29']

# Five residue-border colours for {3,5,7,11,13}
RES_HEX = ['#f472b6','#a78bfa','#67e8f9','#fcd34d','#86efac']
RES_RGB = [(244,114,182),(167,139,250),(103,232,249),(252,211,77),(134,239,172)]

# ─────────────────────────────────────────────────────────────────────────────
# MATH
# ─────────────────────────────────────────────────────────────────────────────

def is_prime(n: int) -> bool:
    if n < 2: return False
    if n == 2: return True
    if n % 2 == 0: return False
    if n < 9: return True
    if n % 3 == 0: return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i+2) == 0: return False
        i += 6
    return True

def sieve(limit: int) -> bytearray:
    """Boolean sieve; index = number, 1 = prime."""
    s = bytearray([1]) * (limit + 1)
    s[0] = s[1] = 0
    for i in range(2, int(limit**0.5) + 1):
        if s[i]:
            s[i*i::i] = bytearray(len(s[i*i::i]))
    return s

def factorize(n: int) -> list:
    f, d = [], 2
    while d*d <= n:
        while n % d == 0: f.append(d); n //= d
        d += 1
    if n > 1: f.append(n)
    return f

def factor_str(n: int) -> str:
    if n <= 1: return str(n)
    c: dict = {}
    for p in factorize(n): c[p] = c.get(p, 0) + 1
    sup = {2:'²',3:'³',4:'⁴',5:'⁵',6:'⁶',7:'⁷',8:'⁸'}
    return ' × '.join(f"{p}{sup.get(e,'')}" for p,e in sorted(c.items()))

def smallest_nine_idx(n: int) -> int:
    for i, p in enumerate(NINE_PRIMES):
        if n % p == 0: return i
    return -1  # pure power of 2

def small_nine_divs(n: int) -> list:
    return [(n % p == 0) for p in NINE_PRIMES]

def hsl_to_rgb(h, s, l):
    h %= 360; s /= 100; l /= 100
    c = (1 - abs(2*l-1)) * s
    x = c * (1 - abs((h/60) % 2 - 1))
    m = l - c/2
    if   h < 60:  r,g,b = c,x,0
    elif h < 120: r,g,b = x,c,0
    elif h < 180: r,g,b = 0,c,x
    elif h < 240: r,g,b = 0,x,c
    elif h < 300: r,g,b = x,0,c
    else:         r,g,b = c,0,x
    return (int((r+m)*255), int((g+m)*255), int((b+m)*255))

def rgb_hex(r,g,b): return '#{:02x}{:02x}{:02x}'.format(r,g,b)

def mod_rgb(n, P):   return hsl_to_rgb(round(360*(n%P)/P), 78, 45)
def mod_hex(n, P):   return rgb_hex(*mod_rgb(n, P))

# ─────────────────────────────────────────────────────────────────────────────
# CUNNINGHAM CHAIN DETECTION
# ─────────────────────────────────────────────────────────────────────────────

def build_cc_data(all_nums: list, prime_set: set):
    """
    Returns:
      mem   dict: n -> set of role strings
                  {'cc1_start','cc1_member','cc2_start','cc2_member','safe'}
      chains list of dicts: [{'type':1|2, 'chain':[p,…], 'len':N}, …]
    """
    mem: dict = {}
    chains: list = []

    def mark(chain, kind):
        for i, p in enumerate(chain):
            if p not in mem: mem[p] = set()
            mem[p].add(f'{kind}_start' if i == 0 else f'{kind}_member')

    primes = sorted(p for p in all_nums if p in prime_set)

    for p in primes:
        # ── CC first kind  p → 2p+1 → 4p+3 → … ────────────────────────────
        pred = (p - 1) / 2
        if pred > 1 and pred == int(pred) and int(pred) in prime_set:
            pass  # not a start — still track membership via mark below
        else:
            ch = [p]; cur = p
            while True:
                nxt = 2*cur + 1
                if nxt not in prime_set: break
                ch.append(nxt); cur = nxt
            if len(ch) >= 2:
                mark(ch, 'cc1')
                chains.append({'type': 1, 'chain': ch, 'len': len(ch)})

        # ── CC second kind  p → 2p−1 → 4p−3 → … ───────────────────────────
        pred2 = (p + 1) / 2
        if pred2 > 1 and pred2 == int(pred2) and int(pred2) in prime_set:
            pass
        else:
            ch2 = [p]; cur = p
            while True:
                nxt = 2*cur - 1
                if nxt < 2 or nxt not in prime_set: break
                ch2.append(nxt); cur = nxt
            if len(ch2) >= 2:
                mark(ch2, 'cc2')
                chains.append({'type': 2, 'chain': ch2, 'len': len(ch2)})

    # safe primes
    for p in primes:
        sg = (p - 1) // 2
        if (p - 1) % 2 == 0 and sg in prime_set:
            if p not in mem: mem[p] = set()
            mem[p].add('safe')

    return mem, chains


# backwards-compat alias used internally
def build_cc_membership(all_nums, prime_set):
    mem, _ = build_cc_data(all_nums, prime_set)
    return mem

# ─────────────────────────────────────────────────────────────────────────────
# GRID CONSTRUCTION
# ─────────────────────────────────────────────────────────────────────────────

def build_grid(num_rows: int):
    num_cols = 2 ** (num_rows - 1)
    bottom_row = num_rows - 1
    grid = [[None]*num_cols for _ in range(num_rows)]

    for col in range(num_cols):
        grid[bottom_row][col] = 2**(num_rows-1) + col

    for row in range(bottom_row, 0, -1):
        for col in range(num_cols):
            n = grid[row][col]
            if n is not None and n % 2 == 0:
                grid[row-1][col] = n // 2

    return grid, num_cols


# ─────────────────────────────────────────────────────────────────────────────
# POSITION MAP  n → (row, col)
# ─────────────────────────────────────────────────────────────────────────────

def build_pos_map(grid) -> dict:
    """Build reverse lookup: number → (row, col)."""
    pm = {}
    for r, row in enumerate(grid):
        for c, n in enumerate(row):
            if n is not None:
                pm[n] = (r, c)
    return pm


# ─────────────────────────────────────────────────────────────────────────────
# CHAIN LINE DRAWING  (corner-to-corner)
# ─────────────────────────────────────────────────────────────────────────────
#
#  CC1 (1st kind, 2p+1): connect via RIGHT-TOP corner of each cell
#  CC2 (2nd kind, 2p−1): connect via LEFT-TOP corner of each cell
#
#  Simple straight lines, no arrowheads.  Small dots at connection corners.
#  "CCk·N" badge at the chain start cell.

CC1_LINE  = '#ef4444'   # red  (consistent with CC mesh)
CC1_DOT   = '#fca5a5'
CC2_LINE  = '#3b82f6'   # blue (consistent with CC mesh)
CC2_DOT   = '#93c5fd'


def _cell_corner(row, col, cw, ch, kind, title_h=0, inset=2):
    """Return (x, y) of the connection corner for a cell.
    
    CC1 (kind=1): right-top corner
    CC2 (kind=2): left-top corner
    """
    if kind == 1:
        # right-top corner, inset slightly
        x = col * cw + cw - inset
        y = title_h + row * ch + inset
    else:
        # left-top corner, inset slightly
        x = col * cw + inset
        y = title_h + row * ch + inset
    return (x, y)


def _cell_center(row, col, cw, ch, title_h=0):
    return (col * cw + cw // 2, title_h + row * ch + ch // 2)


# ── SVG chain lines ──────────────────────────────────────────────────────────

def svg_chain_lines(chains: list, pos_map: dict,
                    cw: int, ch: int, min_len: int,
                    title_h: int = 0) -> list:
    """Return list of SVG element strings for corner-to-corner CC connections."""
    parts = []

    lw     = max(1.0, min(3.0, cw / 12))
    inset  = max(2, cw // 8)
    dot_r  = max(1.5, min(3, cw / 10))

    for info in sorted(chains, key=lambda x: -x['len']):
        if info['len'] < min_len: continue
        ch_type  = info['type']
        members  = info['chain']
        line_col = CC1_LINE if ch_type == 1 else CC2_LINE
        dot_col  = CC1_DOT  if ch_type == 1 else CC2_DOT

        # Collect corner positions for members that are in the grid
        pts = []
        for n in members:
            if n in pos_map:
                r, c = pos_map[n]
                pts.append((_cell_corner(r, c, cw, ch, ch_type, title_h, inset), n))
            else:
                pts.append((None, n))

        # Draw line segments between consecutive visible members
        for i in range(len(pts) - 1):
            (p0, _), (p1, _) = pts[i], pts[i+1]
            if p0 is None or p1 is None: continue
            x0, y0 = p0; x1, y1 = p1
            parts.append(
                f'<line x1="{x0}" y1="{y0}" x2="{x1}" y2="{y1}" '
                f'stroke="{line_col}" stroke-width="{lw:.1f}" '
                f'stroke-opacity="0.8" stroke-linecap="round"/>'
            )

        # Draw small dots at each corner connection point
        for (pt, n) in pts:
            if pt is None: continue
            x, y = pt
            parts.append(
                f'<circle cx="{x}" cy="{y}" r="{dot_r:.1f}" '
                f'fill="{dot_col}" stroke="{line_col}" stroke-width="{max(0.5,lw*0.4):.1f}" '
                f'opacity="0.9"/>'
            )

        # Badge at start: "CC1·N" or "CC2·N"
        if pts[0][0] is not None:
            cx, cy = _cell_center(*pos_map[members[0]], cw, ch, title_h)
            badge_txt = f'CC{ch_type}·{info["len"]}'
            fsize = max(6, min(11, cw // 3))
            bpad = 3
            tw = len(badge_txt) * fsize * 0.55
            bx = cx - tw/2 - bpad
            by = cy - ch//2 - fsize - bpad*2
            parts += [
                f'<rect x="{bx:.0f}" y="{by:.0f}" '
                f'width="{tw + bpad*2:.0f}" height="{fsize + bpad*2:.0f}" '
                f'fill="{line_col}" rx="3" opacity="0.88"/>',
                f'<text x="{cx}" y="{by + fsize/2 + bpad:.0f}" '
                f'font-size="{fsize}" fill="#fff" font-weight="700" '
                f'text-anchor="middle" dominant-baseline="middle" '
                f'font-family="monospace">{badge_txt}</text>',
            ]

    return parts


# ── PNG chain lines ──────────────────────────────────────────────────────────

def png_chain_lines(draw, chains: list, pos_map: dict,
                    cw: int, ch: int, min_len: int,
                    title_h: int = 0):
    """Draw corner-to-corner chain connection lines on a PIL ImageDraw."""

    lw     = max(1, cw // 10)
    inset  = max(2, cw // 8)
    dot_r  = max(1, cw // 10)

    for info in sorted(chains, key=lambda x: -x['len']):
        if info['len'] < min_len: continue
        ch_type  = info['type']
        members  = info['chain']
        line_col = (239,68,68) if ch_type == 1 else (59,130,246)
        dot_col  = (252,165,165) if ch_type == 1 else (147,197,253)

        pts = []
        for n in members:
            if n in pos_map:
                r, c = pos_map[n]
                pts.append(_cell_corner(r, c, cw, ch, ch_type, title_h, inset))
            else:
                pts.append(None)

        # Line segments
        for i in range(len(pts) - 1):
            p0, p1 = pts[i], pts[i+1]
            if p0 is None or p1 is None: continue
            draw.line([p0, p1], fill=line_col, width=lw)

        # Corner dots
        for pt in pts:
            if pt is None: continue
            x, y = pt
            draw.ellipse([x-dot_r, y-dot_r, x+dot_r, y+dot_r],
                         fill=dot_col, outline=line_col, width=max(1, lw//2))

        # Badge at start
        if pts[0] is not None and members[0] in pos_map:
            cx, cy = _cell_center(*pos_map[members[0]], cw, ch, title_h)
            badge = f'CC{ch_type}·{info["len"]}'
            fh = max(7, cw // 4)
            bw2 = len(badge) * fh // 2 + 4
            draw.rectangle([cx-bw2, cy-ch//2-fh-4, cx+bw2, cy-ch//2-2],
                           fill=line_col)
            try:
                fn = try_font(fh)
                draw.text((cx, cy-ch//2-fh//2-3), badge,
                          fill=(255,255,255), font=fn, anchor='mm')
            except Exception:
                pass

# ─────────────────────────────────────────────────────────────────────────────
# CELL STYLE — returns (fill_rgb, fill_hex, border_extras, label_color_hex)
# border_extras: list of (side, hex_color) for SVG residue borders
# ─────────────────────────────────────────────────────────────────────────────

class CellStyle:
    __slots__ = ['fill_rgb','fill_hex','label_hex','border_hex',
                 'res_borders','glow']
    def __init__(self, fill_rgb=(255,255,255), fill_hex='#ffffff',
                 label_hex='#000000', border_hex='#cccccc',
                 res_borders=None, glow=False):
        self.fill_rgb   = fill_rgb
        self.fill_hex   = fill_hex
        self.label_hex  = label_hex
        self.border_hex = border_hex
        self.res_borders = res_borders or []   # [(side, hex)] for residue mode
        self.glow = glow

EMPTY_STYLE = CellStyle((245,245,245),'#f5f5f5','#ccc','#eeeeee')
EMPTY_DARK  = CellStyle((10,10,22),'#0a0a16','#333','#222244')

def classic_style(n: int, dark: bool) -> CellStyle:
    if dark:
        if n == 1:               return CellStyle((255,140,0),'#ff8c00','#000','#ff8c00')
        if (n & (n-1)) == 0:     return CellStyle((0,180,60),'#00b43c','#000','#00b43c')
        if n in {3,5,17,257,65537}: return CellStyle((50,100,255),'#3264ff','#fff','#3264ff')
        m = n+1
        if m > 0 and (m&(m-1))==0 and is_prime(n):
            return CellStyle((200,30,30),'#c81e1e','#fff','#c81e1e')
        if is_prime(n):
            return CellStyle((30,30,30),'#1e1e1e','#eee','#444')
        return CellStyle((40,40,60),'#28283c','#777','#333355')
    else:
        if n == 1:               return CellStyle((255,200,80),'#ffc850','#000','#ccaa00')
        if (n & (n-1)) == 0:     return CellStyle((150,255,150),'#96ff96','#000','#44aa44')
        if n in {3,5,17,257,65537}: return CellStyle((100,160,255),'#64a0ff','#000','#2255cc')
        m = n+1
        if m > 0 and (m&(m-1))==0 and is_prime(n):
            return CellStyle((255,80,80),'#ff5050','#fff','#cc0000')
        if is_prime(n):
            return CellStyle((30,30,30),'#1e1e1e','#fff','#555')
        return CellStyle((200,200,200),'#c8c8c8','#555','#aaa')

def factor_style(n: int, is_p: bool, modprime: int) -> CellStyle:
    if is_p:
        p1 = n + 1
        i = smallest_nine_idx(p1)
        if i >= 0:
            r,g,b = NINE_RGB[i]
            return CellStyle((r,g,b), NINE_HEX[i],
                             '#000' if (r+g+b)>400 else '#fff', NINE_HEX[i], glow=(i<3))
        # Mersenne-like: p+1 is pure power of 2
        return CellStyle((255,255,255),'#ffffff','#000','#ffd700', glow=True)
    # composite
    r,g,b = mod_rgb(n, modprime)
    return CellStyle((r,g,b), rgb_hex(r,g,b), '#000' if (r+g+b)>400 else '#fff',
                     rgb_hex(r//2,g//2,b//2))

def modp_style(n: int, P: int, is_p: bool) -> CellStyle:
    r,g,b = mod_rgb(n, P)
    lum = 0.299*r + 0.587*g + 0.114*b
    lc = '#000' if lum > 140 else '#fff'
    bord = '#ffffff88' if is_p else rgb_hex(r//2,g//2,b//2)
    return CellStyle((r,g,b), rgb_hex(r,g,b), lc, bord)

def chain_style(n: int, is_p: bool, cc_mem: dict) -> CellStyle:
    mem = cc_mem.get(n, set())
    if 'cc1_start'  in mem: return CellStyle((180,20,20),'#b41414','#fff','#ff4444')
    if 'cc1_member' in mem: return CellStyle((220,60,60),'#dc3c3c','#fff','#ff8888')
    if 'cc2_start'  in mem: return CellStyle((20,80,200),'#1450c8','#fff','#4488ff')
    if 'cc2_member' in mem: return CellStyle((60,120,230),'#3c78e6','#fff','#88aaff')
    if 'safe'       in mem: return CellStyle((20,150,60),'#14963c','#fff','#44ff88')
    if is_p:                return CellStyle((30,30,30),'#1e1e1e','#eee','#555')
    return CellStyle((45,45,65),'#2d2d41','#666','#333355')

def residues_style(n: int, is_p: bool, dark: bool) -> CellStyle:
    """Cell fill from classic; residue borders added."""
    base = classic_style(n, dark)
    # borders: left=÷3, top=÷5, right=÷7, bottom=÷11, inner-border=÷13
    sides = ['left','top','right','bottom','inner']
    borders = []
    for i,(p,side) in enumerate(zip([3,5,7,11,13], sides)):
        if n % p == 0:
            borders.append((side, RES_HEX[i]))
    base.res_borders = borders
    return base

def p1grid_style(n: int, is_p: bool, dark: bool) -> CellStyle:
    """Dark base, mini-grid drawn separately."""
    if dark:
        bg = (10,10,24) if is_p else (5,5,15)
        return CellStyle(bg, rgb_hex(*bg), '#aaa', '#222')
    else:
        bg = (230,240,255) if is_p else (245,245,245)
        return CellStyle(bg, rgb_hex(*bg), '#333', '#ccc')

def get_style(n: int, color: str, dark: bool, modprime: int,
              cc_mem: dict, prime_set: set) -> CellStyle:
    is_p = (n in prime_set)
    if color in ('classic','dark'):
        return classic_style(n, dark or color=='dark')
    if color == 'factor':
        return factor_style(n, is_p, modprime)
    if color == 'modp':
        return modp_style(n, modprime, is_p)
    if color == 'chain':
        return chain_style(n, is_p, cc_mem)
    if color == 'residues':
        return residues_style(n, is_p, dark)
    if color == 'p1grid':
        return p1grid_style(n, is_p, dark)
    return classic_style(n, dark)

# ─────────────────────────────────────────────────────────────────────────────
# FONT HELPER
# ─────────────────────────────────────────────────────────────────────────────

def try_font(size: int):
    from PIL import ImageFont
    for p in ["/usr/share/fonts/truetype/dejavu/DejaVuSansMono-Bold.ttf",
              "/usr/share/fonts/truetype/liberation/LiberationMono-Bold.ttf",
              "/usr/share/fonts/truetype/freefont/FreeMono.ttf"]:
        if os.path.exists(p):
            try: return ImageFont.truetype(p, size)
            except: pass
    return ImageFont.load_default()

# ─────────────────────────────────────────────────────────────────────────────
# ════════════════  SVG RENDERER  ════════════════════════════════════════════
# ─────────────────────────────────────────────────────────────────────────────

def render_svg(grid, num_rows: int, num_cols: int,
               color: str, dark: bool, modprime: int,
               cw: int, ch: int, cc_mem: dict, prime_set: set,
               show_numbers: bool,
               chains: list = None, pos_map: dict = None,
               show_chains: bool = True, min_chain: int = 2) -> str:

    bg_page = '#0e0e1c' if (dark or color=='dark') else '#ffffff'
    LEGEND_H = 50
    total_w = num_cols * cw
    total_h = num_rows * ch + LEGEND_H

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'viewBox="0 0 {total_w} {total_h}" width="{total_w}" height="{total_h}">',
        f'<rect width="{total_w}" height="{total_h}" fill="{bg_page}"/>',
    ]

    text_fsize = max(6, min(11, int(cw * 0.28)))

    for row in range(num_rows):
        for col in range(num_cols):
            n = grid[row][col]
            x, y = col * cw, row * ch

            if n is None:
                es = EMPTY_DARK if (dark or color=='dark') else EMPTY_STYLE
                parts.append(f'<rect x="{x}" y="{y}" width="{cw}" height="{ch}" '
                             f'fill="{es.fill_hex}" stroke="{es.border_hex}" stroke-width="0.5"/>')
                continue

            st = get_style(n, color, dark, modprime, cc_mem, prime_set)

            # Base cell
            parts.append(f'<rect x="{x+0.5}" y="{y+0.5}" width="{cw-1}" height="{ch-1}" '
                         f'fill="{st.fill_hex}" stroke="{st.border_hex}" stroke-width="1"/>')

            # Residue borders (thick coloured lines on cell sides)
            if st.res_borders and cw >= 8:
                bw = max(2, cw // 8)
                for side, hx in st.res_borders:
                    if side == 'left':
                        parts.append(f'<line x1="{x+bw//2}" y1="{y}" x2="{x+bw//2}" y2="{y+ch}" '
                                     f'stroke="{hx}" stroke-width="{bw}"/>')
                    elif side == 'top':
                        parts.append(f'<line x1="{x}" y1="{y+bw//2}" x2="{x+cw}" y2="{y+bw//2}" '
                                     f'stroke="{hx}" stroke-width="{bw}"/>')
                    elif side == 'right':
                        parts.append(f'<line x1="{x+cw-bw//2}" y1="{y}" x2="{x+cw-bw//2}" y2="{y+ch}" '
                                     f'stroke="{hx}" stroke-width="{bw}"/>')
                    elif side == 'bottom':
                        parts.append(f'<line x1="{x}" y1="{y+ch-bw//2}" x2="{x+cw}" y2="{y+ch-bw//2}" '
                                     f'stroke="{hx}" stroke-width="{bw}"/>')

            # p1grid mini-squares
            if color == 'p1grid' and cw >= 20:
                is_p = (n in prime_set)
                divs = small_nine_divs(n+1) if is_p else small_nine_divs(n)
                pad = max(1, int(cw*0.08)); avail = cw - 2*pad
                gap = max(1, int(cw*0.03)); sq = max(1,(avail-gap*2)//3)
                for gi in range(9):
                    gr, gc = divmod(gi, 3)
                    sx = x + pad + gc*(sq+gap)
                    sy = y + pad + gr*(sq+gap)
                    if sy + sq > y + ch - 1: continue
                    if divs[gi]:
                        parts.append(f'<rect x="{sx}" y="{sy}" width="{sq}" height="{sq}" '
                                     f'fill="{NINE_HEX[gi]}"/>')
                    else:
                        parts.append(f'<rect x="{sx}" y="{sy}" width="{sq}" height="{sq}" '
                                     f'fill="#0c0c20" stroke="#1c1c38" stroke-width="0.5"/>')
                    if sq >= 8:
                        lc = '#000' if divs[gi] else '#22224a'
                        fs3 = max(4, int(sq*0.48))
                        parts.append(f'<text x="{sx+sq//2}" y="{sy+sq//2}" font-size="{fs3}" '
                                     f'fill="{lc}" text-anchor="middle" dominant-baseline="middle" '
                                     f'font-family="monospace" font-weight="700">{NINE_LABELS[gi]}</text>')

            # Number label
            if show_numbers and cw >= 16:
                fsize = text_fsize
                if n >= 10000: fsize = max(5, fsize-2)
                elif n >= 1000: fsize = max(5, fsize-1)
                cx2, cy2 = x+cw//2, y+ch//2
                # Shadow
                parts.append(f'<text x="{cx2+0.5}" y="{cy2+0.5}" font-size="{fsize}" '
                             f'fill="#00000066" text-anchor="middle" dominant-baseline="middle" '
                             f'font-family="monospace">{n}</text>')
                parts.append(f'<text x="{cx2}" y="{cy2}" font-size="{fsize}" '
                             f'fill="{st.label_hex}" text-anchor="middle" dominant-baseline="middle" '
                             f'font-family="monospace" font-weight="bold">{n}</text>')

    # ── Chain relationship lines ──────────────────────────────────────────────
    if show_chains and chains and pos_map:
        parts += svg_chain_lines(chains, pos_map, cw, ch, min_chain, title_h=0)

    # ── Legend ───────────────────────────────────────────────────────────────
    ly = num_rows * ch
    lbg = '#0a0a16' if (dark or color=='dark') else '#f0f0f0'
    ltxt = '#aaaacc' if (dark or color=='dark') else '#333'
    parts.append(f'<rect x="0" y="{ly}" width="{total_w}" height="{LEGEND_H}" fill="{lbg}"/>')
    parts.append(f'<line x1="0" y1="{ly}" x2="{total_w}" y2="{ly}" stroke="#33335a" stroke-width="1"/>')

    legend_items = _legend_items(color, modprime)
    sq_l = 12; xc = 10; yc = ly + (LEGEND_H - sq_l) // 2
    fs_l = 9
    for label, hx in legend_items:
        parts.append(f'<rect x="{xc}" y="{yc}" width="{sq_l}" height="{sq_l}" fill="{hx}" rx="2"/>')
        parts.append(f'<text x="{xc+sq_l+3}" y="{yc+sq_l//2}" font-size="{fs_l}" '
                     f'fill="{ltxt}" dominant-baseline="middle" font-family="sans-serif">{label}</text>')
        xc += sq_l + 5 + len(label)*5 + 4
        if xc > total_w - 80: break

    # Mode label right
    mode_name = _mode_name(color, modprime)
    parts.append(f'<text x="{total_w-8}" y="{ly+LEGEND_H//2}" font-size="10" '
                 f'fill="{ltxt}" text-anchor="end" dominant-baseline="middle" '
                 f'font-family="monospace">{mode_name} · rows={num_rows}</text>')

    parts.append('</svg>')
    return '\n'.join(parts)

# ─────────────────────────────────────────────────────────────────────────────
# ════════════════  PNG RENDERER  ════════════════════════════════════════════
# ─────────────────────────────────────────────────────────────────────────────

def render_png(grid, num_rows: int, num_cols: int,
               color: str, dark: bool, modprime: int,
               cw: int, ch: int, cc_mem: dict, prime_set: set,
               show_numbers: bool,
               chains: list = None, pos_map: dict = None,
               show_chains: bool = True, min_chain: int = 2) -> 'Image':
    from PIL import Image, ImageDraw

    LEGEND_H = max(30, ch + 10)
    TITLE_H  = max(16, ch // 2)
    iw = num_cols * cw
    ih = TITLE_H + num_rows * ch + LEGEND_H
    bg_col = (10,10,20) if (dark or color=='dark') else (255,255,255)

    img  = Image.new('RGB', (iw, ih), bg_col)
    draw = ImageDraw.Draw(img)
    fn   = try_font(max(7, int(cw * 0.30)))
    fs   = try_font(max(7, min(10, int(cw * 0.18))))

    for row in range(num_rows):
        for col in range(num_cols):
            n = grid[row][col]
            px, py = col * cw, TITLE_H + row * ch

            if n is None:
                ec = (10,10,20) if (dark or color=='dark') else (245,245,245)
                draw.rectangle([px, py, px+cw-1, py+ch-1], fill=ec)
                continue

            st = get_style(n, color, dark, modprime, cc_mem, prime_set)
            r1,c1 = px, py; r2,c2 = px+cw-2, py+ch-2

            if st.glow and cw >= 8:
                gv = tuple(v//5 for v in st.fill_rgb)
                draw.rectangle([r1-1, c1-1, r2+1, c2+1], fill=gv)
            draw.rectangle([r1, c1, r2, c2], fill=st.fill_rgb)

            # Residue mode: coloured side bars
            if color == 'residues' and st.res_borders and cw >= 6:
                bw = max(2, cw // 8)
                for side, hx in st.res_borders:
                    rc = tuple(int(hx[i:i+2],16) for i in (1,3,5))
                    if side == 'left':
                        draw.rectangle([r1,c1,r1+bw,c2], fill=rc)
                    elif side == 'top':
                        draw.rectangle([r1,c1,r2,c1+bw], fill=rc)
                    elif side == 'right':
                        draw.rectangle([r2-bw,c1,r2,c2], fill=rc)
                    elif side == 'bottom':
                        draw.rectangle([r1,c2-bw,r2,c2], fill=rc)

            # p1grid mini squares
            if color == 'p1grid' and cw >= 18:
                is_p = (n in prime_set)
                divs = small_nine_divs(n+1) if is_p else small_nine_divs(n)
                pad = max(1, int(cw*0.08)); avail = cw - 2*pad
                gap = max(1, int(cw*0.03)); sq = max(1,(avail-gap*2)//3)
                for gi in range(9):
                    gr2, gc2 = divmod(gi, 3)
                    sx = px + pad + gc2*(sq+gap)
                    sy = py + pad + gr2*(sq+gap)
                    if sy + sq > py + ch - 1: continue
                    if divs[gi]:
                        draw.rectangle([sx,sy,sx+sq-1,sy+sq-1], fill=NINE_RGB[gi])
                    else:
                        draw.rectangle([sx,sy,sx+sq-1,sy+sq-1], fill=(12,12,32))

            # Number label
            if show_numbers and cw >= 18:
                label = str(n)
                cx2, cy2 = px + cw//2, py + ch//2
                lum = 0.299*st.fill_rgb[0]+0.587*st.fill_rgb[1]+0.114*st.fill_rgb[2]
                tc = (20,20,30) if lum > 140 else (210,210,225)
                if color == 'p1grid': tc = (200,200,220) if (dark or color=='dark') else (30,30,50)
                draw.text((cx2+1,cy2+1), label, fill=(0,0,0), font=fn, anchor='mm')
                draw.text((cx2, cy2),    label, fill=tc,      font=fn, anchor='mm')

    # ── Chain relationship lines ──────────────────────────────────────────────
    if show_chains and chains and pos_map:
        png_chain_lines(draw, chains, pos_map, cw, ch, min_chain, title_h=TITLE_H)

    # Title bar
    mode_name = _mode_name(color, modprime)
    bot_start = 2**(num_rows-1); bot_end = 2**num_rows - 1
    title = f'2-Adic Tree  rows={num_rows}  ·  {mode_name}  ·  bottom: {bot_start}–{bot_end}'
    tc_bar = (12,12,26) if (dark or color=='dark') else (240,242,250)
    draw.rectangle([0,0,iw,TITLE_H], fill=tc_bar)
    draw.text((8, TITLE_H//2), title, fill=(160,165,210), font=fs, anchor='lm')

    # Legend bar
    leg_y = TITLE_H + num_rows * ch
    lbg = (10,10,20) if (dark or color=='dark') else (240,240,248)
    draw.rectangle([0, leg_y, iw, ih], fill=lbg)
    draw.line([0, leg_y, iw, leg_y], fill=(35,35,65), width=1)

    sq_l = min(12, LEGEND_H-8); xc = 10; yc = leg_y + (LEGEND_H - sq_l)//2
    for label, hx in _legend_items(color, modprime):
        rc = tuple(int(hx[i:i+2],16) for i in (1,3,5))
        draw.rectangle([xc,yc,xc+sq_l,yc+sq_l], fill=rc, outline=(60,60,80))
        draw.text((xc+sq_l+3, yc+sq_l//2), label, fill=(140,140,170), font=fs, anchor='lm')
        xc += sq_l + 4 + len(label)*5 + 4
        if xc > iw - 90: break
    draw.text((iw-8, leg_y+LEGEND_H//2), mode_name.upper(), fill=(160,160,200), font=fs, anchor='rm')

    return img

# ─────────────────────────────────────────────────────────────────────────────
# LEGEND + MODE NAME
# ─────────────────────────────────────────────────────────────────────────────

def _mode_name(color: str, modprime: int) -> str:
    return {
        'classic':  'Classic (light)',
        'dark':     'Classic (dark)',
        'factor':   'p+1 Factor Colour',
        'modp':     f'n mod {modprime}',
        'chain':    'Cunningham Chains',
        'residues': 'Residues mod {3,5,7,11,13}',
        'p1grid':   '3×3 p+1 Mini-Grid',
    }.get(color, color)

def _legend_items(color: str, modprime: int) -> list:
    """Returns [(label, '#rrggbb')]."""
    if color in ('classic','dark'):
        return [('Prime','#1e1e1e'),('Composite','#888888'),
                ('Mersenne','#c81e1e'),('Fermat','#3264ff'),
                ('2^n','#00b43c'),('1','#ff8c00')]
    if color == 'factor':
        return list(zip(NINE_LABELS, NINE_HEX)) + [('pure-2^k p+1','#ffffff')]
    if color == 'modp':
        P = modprime
        return [(str(r), rgb_hex(*hsl_to_rgb(round(360*r/P),78,45))) for r in range(P)]
    if color == 'chain':
        return [('CC1 start','#b41414'),('CC1 member','#dc3c3c'),
                ('CC2 start','#1450c8'),('CC2 member','#3c78e6'),
                ('Safe prime','#14963c'),('Prime','#1e1e1e'),
                ('CC1 line (→corner)','#ef4444'),('CC2 line (←corner)','#3b82f6')]
    if color == 'residues':
        return [('÷3 left','#f472b6'),('÷5 top','#a78bfa'),
                ('÷7 right','#67e8f9'),('÷11 bottom','#fcd34d'),('÷13 inner','#86efac')]
    if color == 'p1grid':
        return list(zip([f'p+1÷{l}' for l in NINE_LABELS], NINE_HEX))
    return []

# ─────────────────────────────────────────────────────────────────────────────
# SAVE HELPERS
# ─────────────────────────────────────────────────────────────────────────────

def save(grid, num_rows, num_cols, color, dark, modprime, cw_png, ch_png,
         cw_svg, ch_svg, cc_mem, prime_set, fmt, show_numbers, output_dir,
         chains=None, pos_map=None, show_chains=True, min_chain=2):

    stem = f'2adic_rows{num_rows}_{color}'
    if color == 'modp': stem += f'_mod{modprime}'
    if dark and color not in ('dark',): stem += '_dark'

    if fmt in ('png','both'):
        if not HAS_PIL:
            print('  ⚠  Pillow not installed (pip install Pillow)')
        else:
            print(f'  PNG  {stem}.png …', end=' ', flush=True)
            img = render_png(grid, num_rows, num_cols, color, dark, modprime,
                             cw_png, ch_png, cc_mem, prime_set, show_numbers,
                             chains, pos_map, show_chains, min_chain)
            path = os.path.join(output_dir, stem+'.png')
            img.save(path, 'PNG')
            kb = os.path.getsize(path)//1024
            print(f'✓  ({img.width}×{img.height}, {kb}KB)')

    if fmt in ('svg','both'):
        print(f'  SVG  {stem}.svg …', end=' ', flush=True)
        svg = render_svg(grid, num_rows, num_cols, color, dark, modprime,
                         cw_svg, ch_svg, cc_mem, prime_set, show_numbers,
                         chains, pos_map, show_chains, min_chain)
        path = os.path.join(output_dir, stem+'.svg')
        with open(path,'w',encoding='utf-8') as f: f.write(svg)
        kb = os.path.getsize(path)//1024
        print(f'✓  ({kb}KB)')

# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

ALL_COLORS = ['classic','dark','factor','modp','chain','residues','p1grid']

def main():
    parser = argparse.ArgumentParser(
        description='2-Adic Tree (row layout) — enhanced colouring',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)

    parser.add_argument('rows', type=int, nargs='?', default=6,
                        help='Number of rows (default 6 → 32 columns)')
    parser.add_argument('--color', '-c', default='classic', choices=ALL_COLORS,
                        help='Colour mode (default: classic)')
    parser.add_argument('--format', '-f', default='svg', dest='fmt',
                        choices=['svg','png','both'],
                        help='Output format: svg | png | both (default: svg)')
    parser.add_argument('--dark', action='store_true',
                        help='Dark background for classic/residues modes')
    parser.add_argument('--modprime', type=int, default=7, choices=NINE_PRIMES,
                        help='Modulus prime for modp mode (default 7)')
    parser.add_argument('--cell', type=int, default=0, metavar='PX',
                        help='Cell size in px for both outputs (0 = auto)')
    parser.add_argument('--cell-w', type=int, default=0, metavar='PX',
                        help='Override cell width (SVG)')
    parser.add_argument('--cell-h', type=int, default=0, metavar='PX',
                        help='Override cell height (SVG)')
    parser.add_argument('--cell-png', type=int, default=0, metavar='PX',
                        help='Override square cell size for PNG')
    parser.add_argument('--no-numbers', action='store_true',
                        help='Omit number labels')
    parser.add_argument('--no-chain-lines', action='store_true',
                        help='Disable Cunningham chain relationship lines')
    parser.add_argument('--min-chain', type=int, default=2, metavar='N',
                        help='Minimum chain length to draw lines for (default 2)')
    parser.add_argument('--all-colors', action='store_true',
                        help='Render all colour modes')
    parser.add_argument('--output', '-o', default='.',
                        help='Output directory (default: .)')

    args = parser.parse_args()
    os.makedirs(args.output, exist_ok=True)

    num_rows = args.rows
    num_cols = 2**(num_rows-1)
    bot_start = 2**(num_rows-1); bot_end = 2**num_rows - 1

    print(f'\n2-Adic Tree  rows={num_rows}  cols={num_cols}  bottom: {bot_start}–{bot_end}')

    # ── auto cell sizes ──────────────────────────────────────────────────────
    # SVG: target ~1200px wide
    auto_svg_w = max(8, 1200 // num_cols)
    auto_svg_h = max(6, int(auto_svg_w * 0.70))
    cw_svg = args.cell_w or args.cell or auto_svg_w
    ch_svg = args.cell_h or args.cell or auto_svg_h

    # PNG: target ~1024px wide
    auto_png = max(1, 1024 // num_cols)
    cs_png = args.cell_png or args.cell or auto_png
    cw_png = cs_png; ch_png = max(1, int(cs_png * 0.70))

    print(f'  cell(SVG) = {cw_svg}×{ch_svg}px  cell(PNG) = {cw_png}×{ch_png}px')

    # ── build grid ───────────────────────────────────────────────────────────
    t0 = time.time()
    print('  Building grid …', end=' ', flush=True)
    grid, _ = build_grid(num_rows)
    print(f'done ({time.time()-t0:.2f}s)')

    # ── prime set + CC membership ────────────────────────────────────────────
    print('  Sieve …', end=' ', flush=True)
    sv = sieve(bot_end)
    prime_set = set(n for row in grid for n in row if n is not None and sv[n])
    print(f'{len(prime_set)} primes found')

    cc_mem: dict = {}
    cc_chains: list = []
    # Always compute CC data — needed for chain lines on all modes
    print('  Computing Cunningham chains …', end=' ', flush=True)
    all_nums = [n for row in grid for n in row if n is not None]
    cc_mem, cc_chains = build_cc_data(all_nums, prime_set)
    print(f'{len(cc_chains)} chains, {len(cc_mem)} members')

    pos_map = build_pos_map(grid)
    show_numbers = not args.no_numbers
    show_chains  = not args.no_chain_lines
    min_chain    = args.min_chain

    # ── render ───────────────────────────────────────────────────────────────
    if args.all_colors:
        print(f'\nRendering all {len(ALL_COLORS)} colour modes …')
        for c in ALL_COLORS:
            print(f'\n[{c}]')
            save(grid, num_rows, num_cols, c, args.dark, args.modprime,
                 cw_png, ch_png, cw_svg, ch_svg, cc_mem, prime_set,
                 args.fmt, show_numbers, args.output,
                 cc_chains, pos_map, show_chains, min_chain)
    else:
        save(grid, num_rows, num_cols, args.color, args.dark, args.modprime,
             cw_png, ch_png, cw_svg, ch_svg, cc_mem, prime_set,
             args.fmt, show_numbers, args.output,
             cc_chains, pos_map, show_chains, min_chain)

    print('\nDone.')


if __name__ == '__main__':
    main()
