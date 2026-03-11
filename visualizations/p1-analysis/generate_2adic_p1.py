#!/usr/bin/env python3
# Copyright (c) 2026 Nenad Mićić <nenad@micic.be>
# SPDX-License-Identifier: Apache-2.0
"""
generate_2adic_p1.py — 2-Adic Tree · p+1 Factorisation Analysis
================================================================
Outputs SVG and/or PNG of the 2-adic top-down landscape with
optional p+1 colouring overlays.

Usage
-----
  python generate_2adic_p1.py                           # D=16, depth, PNG
  python generate_2adic_p1.py --D 32 --format svg       # SVG only
  python generate_2adic_p1.py --D 64 --format both      # PNG + SVG
  python generate_2adic_p1.py --D 1024 --mode factor    # large PNG
  python generate_2adic_p1.py --D 256 --mode grid3      # 3x3 mini-grid
  python generate_2adic_p1.py --D 128 --mode modp --modprime 11 --format svg
  python generate_2adic_p1.py --D 32 --all --format both

Modes
-----
  depth    3D depth (Y-level brightness, red = prime)
  prime    Primes vs composites
  chain    Cunningham chain length heat-map
  factor   Solid colour = smallest of {3..29} dividing p+1
  grid2    2x2 mosaic: does p+1 share {3,5,7,11}?
  grid3    3x3 mosaic: all 9 small primes in p+1
  modp     n mod P rainbow (pick --modprime)

Format
------
  --format png     PNG only (default)
  --format svg     SVG only
  --format both    PNG + SVG

D limits
--------
  SVG is skipped automatically for D > 256 (too large).
  PNG works up to D=1024; uses 1-px cells at that scale.
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

NINE_PRIMES    = [3, 5, 7, 11, 13, 17, 19, 23, 29]
NINE_COLORS_RGB = [
    (255, 140,   0),  # 3  orange
    (255, 215,   0),  # 5  gold
    (  0, 220, 100),  # 7  green
    (  0, 200, 255),  # 11 cyan
    (140, 110, 255),  # 13 violet
    (220,  55, 255),  # 17 purple
    (255,  55, 170),  # 19 pink
    (255,  45,  75),  # 23 red
    ( 70, 255, 195),  # 29 mint
]
NINE_HEX    = ['#ff8c00','#ffd700','#00dc64','#00c8ff',
               '#8c6eff','#dc37ff','#ff37aa','#ff2d4b','#46ffc3']
NINE_LABELS = ['3','5','7','11','13','17','19','23','29']

BG_RGB          = (  8,   8,  16)
CELL_EMPTY_RGB  = ( 10,  10,  22)
CHAIN_GOLD_RGB  = (255, 215,   0)
CHAIN_ORG_RGB   = (255, 136,   0)
PRIME_RED_RGB   = (255,  68, 100)
WHITE_RGB       = (255, 255, 255)
TEXT_DARK_RGB   = ( 20,  20,  30)
TEXT_LIGHT_RGB  = (200, 200, 220)
TEXT_DIM_RGB    = ( 90,  90, 130)
BG_HEX          = '#08080f'
CHAIN_GOLD_HEX  = '#ffd700'
CHAIN_ORG_HEX   = '#ff8800'
PRIME_RED_HEX   = '#ff4464'

ALL_MODES  = ['depth','prime','chain','factor','grid2','grid3','modp']
VALID_D    = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]
SVG_MAX_D  = 256   # larger D -> SVG skipped

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

def factorize(n: int) -> list:
    factors, d = [], 2
    while d*d <= n:
        while n % d == 0:
            factors.append(d); n //= d
        d += 1
    if n > 1: factors.append(n)
    return factors

def factor_str(n: int) -> str:
    if n <= 1: return str(n)
    counts: dict = {}
    for f in factorize(n):
        counts[f] = counts.get(f, 0) + 1
    sup = {2:'²',3:'³',4:'⁴',5:'⁵',6:'⁶',7:'⁷',8:'⁸'}
    return ' × '.join(f"{p}{sup.get(e,'')}" for p,e in sorted(counts.items()))

def small_prime_divs(n: int) -> list:
    return [(n % p == 0) for p in NINE_PRIMES]

def smallest_nine_idx(n: int) -> int:
    for i, p in enumerate(NINE_PRIMES):
        if n % p == 0: return i
    return -1

def get_3d_pos(n: int, D: int) -> Optional[tuple]:
    if n <= 0: return None
    D2 = D * D
    bottom, level = n, 0
    while bottom < D2:  bottom *= 2;  level += 1
    while bottom >= 2*D2: bottom //= 2; level -= 1
    if D2 <= bottom < 2*D2:
        idx = bottom - D2
        return (idx % D, level, idx // D, bottom)
    return None

def follow_chain(start: int, max_len: int = 60) -> list:
    chain, cur = [], start
    while is_prime(cur) and len(chain) < max_len:
        chain.append(cur); cur = 2*cur + 1
    return chain

def is_chain_start(n: int) -> bool:
    if not is_prime(n): return False
    if n <= 2: return n == 2
    pred = (n - 1) / 2
    return not pred.is_integer() or pred < 2 or not is_prime(int(pred))

# ─────────────────────────────────────────────────────────────────────────────
# GRID
# ─────────────────────────────────────────────────────────────────────────────

class Cell:
    __slots__ = ['X','Z','top_num','Y','bottom_num','is_prime',
                 'chain_len','is_chain_start','is_chain_member',
                 'p1_val','p1_str','p1_divs','p1_idx','n_divs']
    def __init__(self, x, z, bn):
        self.X=x; self.Z=z; self.top_num=None; self.Y=0; self.bottom_num=bn
        self.is_prime=False; self.chain_len=0
        self.is_chain_start=False; self.is_chain_member=False
        self.p1_val=None; self.p1_str=''; self.p1_divs=[False]*9
        self.p1_idx=-1; self.n_divs=[False]*9

def sieve(limit: int) -> list:
    """Sieve of Eratosthenes; returns bool array indexed 0..limit."""
    s = bytearray([1]) * (limit + 1)
    s[0] = s[1] = 0
    for i in range(2, int(limit**0.5) + 1):
        if s[i]:
            s[i*i::i] = bytearray(len(s[i*i::i]))
    return s

def build_grid(D: int, verbose: bool = True):
    D2 = D*D; bot_start, bot_end = D2, 2*D2-1
    t0 = time.time()
    if verbose:
        print(f"  Building D={D} grid (range {bot_start}–{bot_end}) …", flush=True)

    # For large D use a sieve — much faster than trial division per number
    use_sieve = (D >= 64)
    prime_sieve = sieve(bot_end) if use_sieve else None
    def prime_check(n):
        return bool(prime_sieve[n]) if use_sieve else is_prime(n)

    grid = [[Cell(x, z, bot_start+z*D+x) for x in range(D)] for z in range(D)]
    max_y = 0

    for n in range(2, bot_end+1):
        if n.bit_length() < 2: continue
        pos = get_3d_pos(n, D)
        if pos is None: continue
        x, y, z, _ = pos
        if 0 <= x < D and 0 <= z < D:
            cell = grid[z][x]
            if cell.top_num is None or y > cell.Y:
                cell.top_num = n; cell.Y = y; cell.is_prime = prime_check(n)
            if y > max_y: max_y = y

    all_chains, chain_set = [], set()
    for z in range(D):
        for x in range(D):
            cell = grid[z][x]
            if cell.top_num and is_chain_start(cell.top_num):
                cell.is_chain_start = True
                ch = follow_chain(cell.top_num)
                cell.chain_len = len(ch)
                if len(ch) >= 2:
                    positions = []
                    for num in ch:
                        p = get_3d_pos(num, D)
                        positions.append((p[0], p[2]) if p else (-1,-1))
                    all_chains.append({'chain':ch,'positions':positions})
                    chain_set.update(ch)

    for z in range(D):
        for x in range(D):
            cell = grid[z][x]
            if not cell.top_num: continue
            if cell.top_num in chain_set: cell.is_chain_member = True
            p1 = cell.top_num + 1
            cell.p1_val = p1; cell.p1_str = factor_str(p1)
            cell.p1_divs = small_prime_divs(p1)
            cell.p1_idx  = smallest_nine_idx(p1)
            cell.n_divs  = small_prime_divs(cell.top_num)

    if verbose:
        primes = sum(1 for row in grid for c in row if c.is_prime)
        print(f"  Done in {time.time()-t0:.1f}s · {primes} primes · {len(all_chains)} CC starts")
    return grid, all_chains, max_y

# ─────────────────────────────────────────────────────────────────────────────
# COLOUR HELPERS
# ─────────────────────────────────────────────────────────────────────────────

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

def rgb_hex(rgb):
    return '#{:02x}{:02x}{:02x}'.format(*rgb)

def mod_color_rgb(n, P):
    return hsl_to_rgb(round(360*(n%P)/P), 78, 45)

def depth_color_rgb(Y, max_y, prime):
    ny = Y/max_y if max_y else 0
    return hsl_to_rgb(0 if prime else 220, 50 if prime else 25, 12+ny*50)

def get_cell_color(cell: Cell, mode: str, max_y: int,
                   highlight_chains: bool, min_chain: int, modprime: int):
    """Returns (bg_rgb, border_rgb_or_None, glow_bool)."""
    if not cell.top_num: return CELL_EMPTY_RGB, None, False
    if highlight_chains and cell.is_chain_member:
        if cell.is_chain_start and cell.chain_len >= min_chain:
            return CHAIN_GOLD_RGB, WHITE_RGB, True
        return CHAIN_ORG_RGB, CHAIN_GOLD_RGB, False

    if mode == 'depth':
        return depth_color_rgb(cell.Y, max_y, cell.is_prime), \
               (PRIME_RED_RGB if cell.is_prime else None), False
    if mode == 'chain':
        if cell.is_chain_start and cell.chain_len >= min_chain:
            return hsl_to_rgb(max(0,60-cell.chain_len*8),90,50), WHITE_RGB, True
        if cell.is_chain_start: return (136,85,0), None, False
        if cell.is_prime:       return (68,34,34), None, False
        return (20,20,40), None, False
    if mode == 'prime':
        if cell.is_prime: return PRIME_RED_RGB, None, False
        ny = cell.Y/max_y if max_y else 0
        return hsl_to_rgb(220, 20, int(10+ny*25)), None, False
    if mode == 'factor':
        if cell.is_prime:
            i = cell.p1_idx
            return (NINE_COLORS_RGB[i] if i>=0 else WHITE_RGB), None, (i<0)
        return mod_color_rgb(cell.top_num, modprime), None, False
    if mode == 'modp':
        h = round(360*(cell.top_num%modprime)/modprime)
        l = 55 if cell.is_prime else 35
        return hsl_to_rgb(h,75,l), ((200,200,200) if cell.is_prime else None), False
    if mode in ('grid2','grid3'):
        return ((12,12,30) if cell.is_prime else (6,6,18)), None, False
    return CELL_EMPTY_RGB, None, False

# ─────────────────────────────────────────────────────────────────────────────
# MINI-GRID (shared logic; callers draw via callback)
# ─────────────────────────────────────────────────────────────────────────────

def mini_grid_rects(px, py, cs, cell: Cell, mode: str):
    """Yield (sx,sy,sq,idx,lit) for each mini-grid square."""
    gn = 3 if mode == 'grid3' else 2
    divs = cell.p1_divs if cell.is_prime else cell.n_divs
    pad  = max(1, int(cs*0.08)); avail = cs - 2*pad
    gap  = max(1, int(cs*0.03)); sq = max(1, (avail - gap*(gn-1)) // gn)
    for row in range(gn):
        for col in range(gn):
            idx = row*gn + col
            if idx >= 9: continue
            yield (px+pad+col*(sq+gap), py+pad+row*(sq+gap), sq, idx, divs[idx])

# ─────────────────────────────────────────────────────────────────────────────
# PNG RENDERER
# ─────────────────────────────────────────────────────────────────────────────

def try_font(size):
    for p in ["/usr/share/fonts/truetype/dejavu/DejaVuSansMono-Bold.ttf",
              "/usr/share/fonts/truetype/liberation/LiberationMono-Bold.ttf",
              "/usr/share/fonts/truetype/freefont/FreeMono.ttf"]:
        if os.path.exists(p):
            try:
                return ImageFont.truetype(p, size)
            except: pass
    return ImageFont.load_default()

def auto_cell_png(D: int) -> int:
    """Target ~1024px canvas."""
    return max(1, 1024 // D)

def render_png(D, mode, cs, modprime, min_chain, show_chains, show_numbers,
               grid, all_chains, max_y):
    LEGEND_H = max(40, cs*2 + 18)
    TITLE_H  = max(16, cs)
    iw, ih   = D*cs, TITLE_H + D*cs + LEGEND_H

    img  = Image.new('RGB', (iw, ih), BG_RGB)
    draw = ImageDraw.Draw(img)
    fn   = try_font(max(7, int(cs*0.32)))
    fs   = try_font(max(7, min(11, int(cs*0.18))))

    for z in range(D):
        for x in range(D):
            cell = grid[z][x]
            bg, border, glow = get_cell_color(cell, mode, max_y, show_chains, min_chain, modprime)
            px, py = x*cs, TITLE_H + z*cs

            # Clamp to at least 1px width/height (cs=1 would make x1<x0)
            rx0, ry0 = px, py
            rx1, ry1 = max(px, px+cs-2), max(py, py+cs-2)

            if glow and cs >= 6:
                draw.rectangle([rx0-1,ry0-1,rx1+1,ry1+1], fill=tuple(v//4 for v in bg))
            draw.rectangle([rx0, ry0, rx1, ry1], fill=bg)
            if border and cs >= 4:
                lw = max(1, cs//12)
                draw.rectangle([rx0,ry0,rx1,ry1], outline=border, width=lw)

            if mode in ('grid2','grid3') and cell.top_num and cs >= 5:
                for sx,sy,sq,idx,lit in mini_grid_rects(px, py, cs, cell, mode):
                    if lit:
                        draw.rectangle([sx,sy,sx+sq-1,sy+sq-1], fill=NINE_COLORS_RGB[idx])
                    else:
                        draw.rectangle([sx,sy,sx+sq-1,sy+sq-1], fill=(12,12,32),
                                       outline=(26,26,52), width=1)

            if show_numbers and cell.top_num and cs >= 18:
                label = str(cell.top_num)
                cx, cy = px+cs//2, py+cs//2
                draw.text((cx+1,cy+1), label, fill=(0,0,0), font=fn, anchor='mm')
                tc = TEXT_DARK_RGB if sum(bg) > 350 else TEXT_LIGHT_RGB
                if mode in ('grid2','grid3'): tc = TEXT_LIGHT_RGB
                draw.text((cx,cy), label, fill=tc, font=fn, anchor='mm')

    # chain lines
    if show_chains and mode not in ('grid2','grid3') and cs >= 3:
        for info in all_chains:
            if len(info['chain']) < min_chain: continue
            pts = [(x*cs+cs//2, TITLE_H+z*cs+cs//2)
                   for x,z in info['positions'] if 0<=x<D and 0<=z<D]
            if len(pts) < 2: continue
            lw = max(1, cs//10)
            for i in range(len(pts)-1):
                draw.line([pts[i],pts[i+1]], fill=CHAIN_GOLD_RGB, width=lw)
            dr = max(2, cs//7)
            for i,pt in enumerate(pts):
                c = WHITE_RGB if i==0 else CHAIN_GOLD_RGB
                draw.ellipse([pt[0]-dr,pt[1]-dr,pt[0]+dr,pt[1]+dr], fill=c)

    # title
    mode_label = _mode_name(mode, modprime)
    title = f'2-Adic Tree  D={D}  ·  {mode_label}  ·  {D*D}–{2*D*D-1}'
    draw.rectangle([0,0,iw,TITLE_H], fill=(12,12,28))
    draw.text((8, TITLE_H//2), title, fill=(180,180,220), font=fs, anchor='lm')

    # legend
    ly0 = TITLE_H + D*cs
    draw.rectangle([0,ly0,iw,ih], fill=(10,10,22))
    draw.line([0,ly0,iw,ly0], fill=(32,32,62), width=1)
    sq = min(13, LEGEND_H-8); xc = 10; yc = ly0 + (LEGEND_H-sq)//2

    items = _legend_items(mode, modprime, min_chain)
    for label, color in items:
        draw.rectangle([xc,yc,xc+sq,yc+sq], fill=color, outline=(50,50,70))
        draw.text((xc+sq+3, yc+sq//2), label, fill=TEXT_DIM_RGB, font=fs, anchor='lm')
        xc += sq + 4 + len(label)*5 + 4
        if xc > iw - 90: break
    draw.text((iw-8, ly0+LEGEND_H//2), mode_label.upper(),
              fill=(170,170,210), font=fs, anchor='rm')
    return img

# ─────────────────────────────────────────────────────────────────────────────
# SVG RENDERER
# ─────────────────────────────────────────────────────────────────────────────

def auto_cell_svg(D: int) -> int:
    """Target ~1200px wide for SVG."""
    return max(4, 1200 // D)

def render_svg(D, mode, cs, modprime, min_chain, show_chains, show_numbers,
               grid, all_chains, max_y):
    LEGEND_H = max(36, cs*2 + 14)
    TITLE_H  = max(16, cs)
    iw = D*cs; ih = TITLE_H + D*cs + LEGEND_H
    mode_label = _mode_name(mode, modprime)
    title = f'2-Adic Tree  D={D}  ·  {mode_label}  ·  {D*D}–{2*D*D-1}'

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {iw} {ih}" width="{iw}" height="{ih}">',
        f'<rect width="{iw}" height="{ih}" fill="{BG_HEX}"/>',
        f'<rect x="0" y="0" width="{iw}" height="{TITLE_H}" fill="#0c0c1c"/>',
        f'<text x="8" y="{TITLE_H//2}" font-size="{max(8,TITLE_H-5)}" '
        f'fill="#b4b4dc" font-family="monospace" dominant-baseline="middle">{title}</text>',
    ]

    # cells
    for z in range(D):
        for x in range(D):
            cell = grid[z][x]
            bg, border, glow = get_cell_color(cell, mode, max_y, show_chains, min_chain, modprime)
            px = x*cs; py = TITLE_H + z*cs
            bh = rgb_hex(bg)

            if glow:
                parts.append(f'<rect x="{px-1}" y="{py-1}" width="{cs+2}" height="{cs+2}" '
                              f'fill="{bh}" opacity="0.25"/>')
            bstr = (f' stroke="{rgb_hex(border)}" stroke-width="1.5"' if border else '')
            parts.append(f'<rect x="{px+0.5}" y="{py+0.5}" width="{cs-1}" height="{cs-1}" '
                         f'fill="{bh}"{bstr}/>')

            if mode in ('grid2','grid3') and cell.top_num:
                for sx,sy,sq,idx,lit in mini_grid_rects(px, py, cs, cell, mode):
                    if lit:
                        parts.append(f'<rect x="{sx}" y="{sy}" width="{sq}" height="{sq}" '
                                     f'fill="{NINE_HEX[idx]}"/>')
                    else:
                        parts.append(f'<rect x="{sx}" y="{sy}" width="{sq}" height="{sq}" '
                                     f'fill="#0c0c20" stroke="#1a1a38" stroke-width="0.5"/>')
                    if sq >= 9:
                        lc = '#000' if lit else '#20203a'
                        fs2 = max(5, int(sq*0.45))
                        parts.append(f'<text x="{sx+sq//2}" y="{sy+sq//2}" font-size="{fs2}" '
                                     f'fill="{lc}" font-family="monospace" font-weight="700" '
                                     f'text-anchor="middle" dominant-baseline="middle">'
                                     f'{NINE_LABELS[idx]}</text>')

            if show_numbers and cell.top_num and cs >= 16:
                tc = '#141418' if sum(bg) > 350 else '#c8c8dc'
                if mode in ('grid2','grid3'): tc = '#c0c0d8'
                fsize = max(7, int(cs*0.28)) if cell.top_num < 1000 else max(5, int(cs*0.22))
                cx2, cy2 = px+cs//2, py+cs//2
                parts.append(f'<text x="{cx2+0.5}" y="{cy2+0.5}" font-size="{fsize}" '
                              f'fill="#000" font-family="monospace" font-weight="600" '
                              f'text-anchor="middle" dominant-baseline="middle">{cell.top_num}</text>')
                parts.append(f'<text x="{cx2}" y="{cy2}" font-size="{fsize}" '
                              f'fill="{tc}" font-family="monospace" font-weight="600" '
                              f'text-anchor="middle" dominant-baseline="middle">{cell.top_num}</text>')

    # chain lines
    if show_chains and mode not in ('grid2','grid3'):
        for info in all_chains:
            if len(info['chain']) < min_chain: continue
            pts = [(x*cs+cs//2, TITLE_H+z*cs+cs//2)
                   for x,z in info['positions'] if 0<=x<D and 0<=z<D]
            if len(pts) < 2: continue
            d_attr = 'M ' + ' L '.join(f'{px},{py}' for px,py in pts)
            lw = max(1, cs//8)
            dash = f'{lw*3},{lw*2}'
            parts.append(f'<path d="{d_attr}" fill="none" stroke="#ffd700" '
                         f'stroke-width="{lw}" stroke-opacity="0.7" stroke-dasharray="{dash}"/>')
            dr = max(2, cs//7)
            for i,(px2,py2) in enumerate(pts):
                fc = '#ffffff' if i==0 else '#ffd700'
                parts.append(f'<circle cx="{px2}" cy="{py2}" r="{dr}" fill="{fc}"/>')

    # legend
    ly0 = TITLE_H + D*cs
    parts.append(f'<rect x="0" y="{ly0}" width="{iw}" height="{LEGEND_H}" fill="#0a0a16"/>')
    parts.append(f'<line x1="0" y1="{ly0}" x2="{iw}" y2="{ly0}" stroke="#202040" stroke-width="1"/>')
    sq = min(13, LEGEND_H-8); xc = 10; yc = ly0 + (LEGEND_H-sq)//2
    fsl = max(7, sq-2)

    for label, color in _legend_items_hex(mode, modprime, min_chain):
        parts.append(f'<rect x="{xc}" y="{yc}" width="{sq}" height="{sq}" fill="{color}" rx="2"/>')
        parts.append(f'<text x="{xc+sq+3}" y="{yc+sq//2}" font-size="{fsl}" '
                     f'fill="#7878aa" font-family="monospace" dominant-baseline="middle">'
                     f'{label}</text>')
        xc += sq + 4 + len(label)*fsl//2 + 6
        if xc > iw - 90: break

    parts.append(f'<text x="{iw-8}" y="{ly0+LEGEND_H//2}" font-size="{fsl+1}" '
                 f'fill="#aaaacc" font-family="monospace" dominant-baseline="middle" '
                 f'text-anchor="end">{mode_label.upper()}</text>')
    parts.append('</svg>')
    return '\n'.join(parts)

# ─────────────────────────────────────────────────────────────────────────────
# COMPARISON SHEET
# ─────────────────────────────────────────────────────────────────────────────

def render_comparison_png(D, cs, modprime, min_chain, output_dir, grid, all_chains, max_y):
    modes = ['depth','prime','factor','grid3','modp','chain','grid2']
    images = []
    for m in modes:
        print(f"    {m} …", end=' ', flush=True)
        images.append(render_png(D, m, cs, modprime, min_chain, True,
                                 cs >= 18, grid, all_chains, max_y))
        print('ok')

    cols = 2; rows = math.ceil(len(modes)/2)
    W,H = images[0].width, images[0].height
    PAD = 5
    sheet_w = cols*W + (cols+1)*PAD
    sheet_h = 40 + rows*H + (rows+1)*PAD

    sheet = Image.new('RGB', (sheet_w, sheet_h), (6,6,14))
    draw  = ImageDraw.Draw(sheet)
    fnt   = try_font(15)
    draw.text((sheet_w//2, 22), f'2-Adic Tree — All Modes  (D={D})',
              fill=(255,215,0), font=fnt, anchor='mm')
    for i, im in enumerate(images):
        r,c = divmod(i, cols)
        sheet.paste(im, (PAD + c*(W+PAD), 40 + PAD + r*(H+PAD)))

    path = os.path.join(output_dir, f'2adic_comparison_D{D}.png')
    sheet.save(path, 'PNG')
    print(f"  Saved comparison sheet → {path}")
    return path

# ─────────────────────────────────────────────────────────────────────────────
# HELPERS
# ─────────────────────────────────────────────────────────────────────────────

def _mode_name(mode, modprime):
    return {'depth':'3D Depth','prime':'Primes','chain':'Cunningham Chains',
            'factor':'p+1 Factor Color','grid2':'2×2 Mini-Grid',
            'grid3':'3×3 Mini-Grid','modp':f'n mod {modprime}'}.get(mode, mode)

def _legend_items(mode, modprime, min_chain):
    """Return [(label, rgb_tuple)]."""
    if mode == 'factor':
        items = list(zip(NINE_LABELS, NINE_COLORS_RGB)) + [('?p+1', WHITE_RGB)]
    elif mode == 'modp':
        items = [(str(r), hsl_to_rgb(round(360*r/modprime),75,45)) for r in range(modprime)]
    elif mode in ('grid2','grid3'):
        n = 4 if mode=='grid2' else 9
        items = [(f'÷{NINE_LABELS[i]}', NINE_COLORS_RGB[i]) for i in range(n)]
    else:
        items = [('prime', PRIME_RED_RGB), ('composite',(30,30,90)),
                 (f'CC≥{min_chain}', CHAIN_GOLD_RGB), ('CC', CHAIN_ORG_RGB)]
    return items

def _legend_items_hex(mode, modprime, min_chain):
    return [(lbl, rgb_hex(col)) for lbl, col in _legend_items(mode, modprime, min_chain)]

# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='2-Adic Tree · p+1 Analysis — PNG and SVG output',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)

    parser.add_argument('--D', type=int, default=16, metavar='D',
                        choices=VALID_D,
                        help=f'Grid dimension. Choices: {VALID_D}. Default 16.')
    parser.add_argument('--mode', default='depth', choices=ALL_MODES,
                        help='Colour mode. Default: depth')
    parser.add_argument('--format', dest='fmt', default='png',
                        choices=['png','svg','both'],
                        help='Output format: png | svg | both. Default: png')
    parser.add_argument('--cell', type=int, default=0, metavar='PX',
                        help='Cell size in px for both formats (0 = auto).')
    parser.add_argument('--cell-png', type=int, default=0, metavar='PX',
                        help='Override cell size for PNG only.')
    parser.add_argument('--cell-svg', type=int, default=0, metavar='PX',
                        help='Override cell size for SVG only.')
    parser.add_argument('--modprime', type=int, default=7, choices=NINE_PRIMES,
                        help='Modulus for modp mode. Default 7.')
    parser.add_argument('--chains', type=int, default=3, metavar='N',
                        help='Min Cunningham chain length. Default 3.')
    parser.add_argument('--no-chains', action='store_true',
                        help='Disable chain overlay.')
    parser.add_argument('--no-numbers', action='store_true',
                        help='Skip number labels on cells.')
    parser.add_argument('--all', action='store_true',
                        help='Render all 7 modes + comparison sheet.')
    parser.add_argument('--output', '-o', default='.',
                        help='Output directory. Default: current dir.')

    args = parser.parse_args()
    D = args.D
    os.makedirs(args.output, exist_ok=True)

    # auto cell sizes
    cs_png = args.cell_png or args.cell or auto_cell_png(D)
    cs_svg = args.cell_svg or args.cell or auto_cell_svg(D)

    show_chains  = not args.no_chains
    show_numbers = not args.no_numbers

    print(f"\n2-Adic Tree  D={D}  format={args.fmt}  "
          f"cell_png={cs_png}px  cell_svg={cs_svg}px")
    if args.fmt in ('svg','both') and D > SVG_MAX_D:
        print(f"  ⚠  SVG output will be skipped for D > {SVG_MAX_D} (file too large).")
    if not HAS_PIL and args.fmt in ('png','both'):
        print("  ⚠  Pillow not installed.  Install with: pip install Pillow")

    # build grid once
    grid, all_chains, max_y = build_grid(D)

    def save_mode(mode):
        stem = f'2adic_D{D}_{mode}' + (f'_mod{args.modprime}' if mode=='modp' else '')

        if args.fmt in ('png','both') and HAS_PIL:
            print(f"  PNG  {stem}.png …", end=' ', flush=True)
            img = render_png(D, mode, cs_png, args.modprime, args.chains,
                             show_chains, show_numbers, grid, all_chains, max_y)
            path = os.path.join(args.output, stem+'.png')
            img.save(path, 'PNG')
            kb = os.path.getsize(path)//1024
            print(f'✓  ({img.width}×{img.height}, {kb}KB)')

        if args.fmt in ('svg','both'):
            if D > SVG_MAX_D:
                print(f"  SVG  skipped (D={D} > {SVG_MAX_D})")
            else:
                print(f"  SVG  {stem}.svg …", end=' ', flush=True)
                svg = render_svg(D, mode, cs_svg, args.modprime, args.chains,
                                 show_chains, show_numbers, grid, all_chains, max_y)
                path = os.path.join(args.output, stem+'.svg')
                with open(path,'w',encoding='utf-8') as f: f.write(svg)
                kb = os.path.getsize(path)//1024
                print(f'✓  ({kb}KB)')

    if args.all:
        print(f"\nRendering all {len(ALL_MODES)} modes …")
        for m in ALL_MODES:
            print(f"\n[{m}]")
            save_mode(m)
        if HAS_PIL and args.fmt in ('png','both'):
            print("\nBuilding comparison sheet …")
            render_comparison_png(D, cs_png, args.modprime, args.chains,
                                  args.output, grid, all_chains, max_y)
    else:
        save_mode(args.mode)

    print("\nDone.")

if __name__ == '__main__':
    main()
