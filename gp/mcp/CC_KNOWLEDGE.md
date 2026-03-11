# Cunningham Chain Mathematics — Agent Knowledge Base

Use this knowledge when working with the cc-math MCP server tools.

## What is a Cunningham Chain?

A sequence of primes where each is derived from the previous:
- **First kind (CC1a):** p₁, p₂=2p₁+1, p₃=2p₂+1, ... (each is 2×previous+1)
- **Second kind (CC1b):** p₁, p₂=2p₁-1, p₃=2p₂-1, ...

Example: 2 → 5 → 11 → 23 → 47 is a CC5a (first-kind, length 5).
The chain breaks when the next number (2×47+1=95=5×19) is composite.

**Root:** the first prime in a chain (no predecessor). 2 is the root of the chain above.
**Breaker:** the first composite that ends the chain (95 above).

## Key Concepts

### Chain Length
CC7a means a first-kind chain of length 7. The famous root 1122659 gives:
1122659 → 2245319 → 4490639 → 8981279 → 17962559 → 35925119 → 71850239

### Safe Primes (CC2a special case)
A safe prime is p=2q+1 where both p and q are prime. This is just a CC2a chain.
Example: q=11, p=23 → the chain 11→23 is CC2a.

### Shadow Spiral / Immunization
For each small prime q, the root p has a residue r = p mod q.
- If r = q-1, the prime q can NEVER kill any member of the chain. We call this **immune**.
- Otherwise, q kills at a specific position (the **kill position**).

Long chains need many immune primes — it's like threading a needle through
multiple sieves simultaneously.

### The Search Pipeline
Finding long CC chains uses a multi-stage filter:
1. **CRT (Chinese Remainder Theorem):** construct candidates that avoid small-prime kills
2. **Sieve:** check against medium primes cheaply
3. **PRP (Probable Prime):** Miller-Rabin/BPSW on survivors
4. **Chain walk:** follow the chain to measure actual length

Each stage is cheaper but less selective than the next.

## Available Tools

### `cc_check(number)` — Quick triage
"Is 1122659 a Cunningham Chain root?" → Yes, CC7a root.
"What about 47?" → Part of a CC5a chain, but not a root (root is 2).

### `cc_analyze(number)` — Deep analysis
Shows chain walk, shadow spiral (immunization), and classification.
Use for understanding WHY a chain is as long as it is.

### `cc_immunization(number)` — Residue analysis
Shows which small primes are immune vs vulnerable, and where vulnerable
primes kill. This explains the "anatomy" of a chain.

### `cc_search(bits, target)` — Find chains
Search for CC roots at a given bit size. Examples:
- 30-bit CC4+ roots: lots of hits
- 40-bit CC5+ roots: moderate hits
- 60-bit CC7+ roots: rare
- 89-bit CC10+ roots: very rare (CC18 is the world-record frontier)

### `cc_safe_prime(bits)` — Generate safe primes
Uses delta-sieve algorithm (beats OpenSSL). Works up to ~1024 bits in GP
(larger sizes better in the C implementation).

### `gp_eval(code)` — Raw GP access
For anything not covered by the above tools. Full cc_lib_v10 is loaded with
37 functions across 10 layers. Key functions:

```
cc_walk(p)              — display chain members
cc_shadow(p)            — small-prime residue analysis
cc_autopsy(p)           — deep analysis of the breaker
cc_full(p)              — everything at once
cc_search(bits, target) — CRT+sieve+PRP search
cc_construct(bits, target, S_limit) — constructor search (p = S*R - 1)
cc_sp_search(bits)      — safe prime search
cc_sp_delta_sieve(bits) — OpenSSL-style delta-sieve
cc_x_forbidden_mask(q, depth)  — bitmask sieve screening
cc_x_line_filter(p, depth)     — sieve-only pre-screen
cc_x_bitwin_walk(n)            — BiTwin chain analysis
cc_x_primorial_scan(bits, prim_n, target) — algebraic form search
cc_x_root_depth(p)             — steps back to root
```

## Interesting Numbers to Explore

| Number | What it is |
|--------|-----------|
| 1122659 | CC7a root (good teaching example) |
| 2 | Root of CC5a: 2→5→11→23→47 |
| 0x3007d09c969e63d7e80a777 | CC16a root (90 bits) |
| 89 | Sophie Germain prime (89→179, CC2a) |
| 11 | Safe prime pair: q=5, p=11 OR q=11, p=23 |
| 6 | BiTwin center: (5,7) both prime, chain continues |

## Teaching Patterns

When explaining CC to someone:
1. Start with the chain 2→5→11→23→47 (CC5a) — small, easy to verify by hand
2. Show that 95=5×19 breaks the chain (the breaker)
3. Introduce immunization: why does prime 3 never kill this chain?
   Because 2 mod 3 = 2 = 3-1 (immune!)
4. Show a longer chain (1122659, CC7a) and its shadow spiral
5. Explain that CC18 is the frontier — nobody has found one yet at the target bit sizes
