# cc_lib_v10.gp — Canned Sessions

Companion to [HOWTO_cc_lib_v10.md](HOWTO_cc_lib_v10.md).

Copy-paste workbook: short GP sessions for the most common tasks.

---

## 1. Load and Smoke Check

```gp
\r cc_lib_v10.gp

cc_chain_test(1122659)      \\ expected: 7
cc_is_root(1122659)         \\ expected: 1
cc_sp_verify(11)            \\ expected: [11, 23, 5]
```

---

## 2. Analyze a Known Root

```gp
p = 1122659;
cc_walk(p)
cc_shadow(p)
cc_autopsy(p)
cc_full(p)
```

---

## 3. Analyze a Non-Root Member

```gp
root = 1122659;
m1 = 2*(root + 1) - 1;
cc_walk(m1)
cc_candidate_profile(m1, 1, 0)
cc_candidate_profile(m1, 1, 1)
```

---

## 4. Analyze the Breaker Composite

```gp
v = cc_walk_quiet(1122659);
breaker = v[4];
cc_candidate_profile(breaker, 1, 0)
cc_candidate_positions(breaker, 8, 1, 0)
```

---

## 5. CRT Search and Inspect First Hit

```gp
hits = cc_search(40, 5, 1, 0, 0, 20000, 0);
#hits
if(#hits > 0,
  p = hits[1][1];
  cc_walk(p);
  cc_chain_verify(p);
);
```

---

## 6. Constructor Search and Inspect First Hit

```gp
hits = cc_construct(40, 4, 13, "B", 0, 0, 20000, 0);
#hits
if(#hits > 0,
  p = hits[1][1];
  cc_walk(p);
  cc_chain_verify(p);
);
```

---

## 7. Safe Prime — Quick Generation

```gp
\\ 512-bit safe prime (~instant)
cc_sp_search(512)

\\ 1024-bit safe prime (~1 second)
cc_sp_search(1024)
```

---

## 8. Safe Prime — Search Space Analysis (POC)

This session shows why combined trial division matters. Use this to
explain the v35-03 optimization to stakeholders or to verify the math.

```gp
\\ Analyze 4096-bit safe prime search space
cc_sp_info(4096)

\\ Key output:
\\   Trial primes: 9805 (2..102397)
\\   Est trial survive: 0.625% (1 in 160)
\\   Speedup: ~160x from trial division alone
```

---

## 9. Safe Prime — Verify the Core Algorithm

This is the POC: demonstrate that cc_sp_trial_check correctly identifies
which candidates need BPSW and which can be skipped.

```gp
\\ Known safe prime: q=11, p=23
cc_sp_trial_check(11, [3, 5, 7])    \\ 1 — both survive trial
cc_sp_verify(11)                     \\ [11, 23, 5] — proven

\\ Not a safe prime: q=13, 2q+1=27=3^3
cc_sp_trial_check(13, [3, 5, 7])    \\ 0 — caught by trial (27 % 3 == 0)
cc_sp_verify(13)                     \\ 0

\\ q=29, p=59 — safe prime
cc_sp_trial_check(29, [3, 5, 7, 11, 13])  \\ 1
cc_sp_verify(29)                           \\ [29, 59, 6]

\\ q=41, 2q+1=83 — safe prime
cc_sp_trial_check(41, cc_sp_trial_primes(3, 100))  \\ 1
cc_sp_verify(41)                                     \\ [41, 83, 7]
```

---

## 10. Safe Prime — 4096-bit Full Workflow

```gp
\\ Search
result = cc_sp_search(4096)

\\ Extract q and p
q = result[1];
p = result[2];

\\ Verify correctness
ispseudoprime(q)                \\ should be 1
ispseudoprime(p)                \\ should be 1
p == 2*q + 1                    \\ should be 1
#binary(p)                      \\ should be 4096

\\ Analyze as CC2 chain
cc_walk(q)

\\ Export hex for use with OpenSSL
printf("q = %s\n", Str(q));
printf("p = %s\n", Str(p));
```

---

## 11. Safe Prime — Benchmark Trial Division Effect

```gp
\\ Small bits: both methods instant, no measurable speedup
cc_sp_bench(64)

\\ Medium bits: trial speedup becomes visible
cc_sp_bench(256)

\\ At 512+ bits, trial division provides significant speedup
\\ because BPSW cost grows cubically while trial stays linear
cc_sp_bench(512, 50000)
```

---

## 12. Safe Prime — Silent Batch Generation

```gp
\\ Generate 5 safe primes at 2048 bits, silently
for(i = 1, 5, \
  my(r = cc_sp_search(2048, 0, 10000000, 0)); \
  if(r, printf("[%d] %d-bit safe prime: q=%s\n", i, r[3], cc_fmt_hex(r[1], 40))); \
)
```

---

## 13. Safe Prime — Cross-Validate with C Code

After generating a safe prime with the C code (`cc_v35_03 --safeprime`),
verify it in GP:

```gp
\\ Paste the hex q value from C output
q = 0xYOUR_HEX_Q_HERE;

\\ Verify
cc_sp_verify(q)          \\ should return [q, 2q+1, bits]
cc_walk(q)               \\ walk the CC2 chain
cc_classify(q)           \\ should include "Sophie Germain"
cc_classify(2*q+1)       \\ should include "Safe prime"
```

---

## 14. Compare Two Roots Side by Side

```gp
p1 = 1122659;
p2 = 0x3007d09c969e63d7e80a777;
cc_compare(p1, p2)
```

---

## 15. Second-Kind Quick Check

```gp
p = 1122659;
cc_walk(p, 2)
cc_shadow(p, 2)
```

---

## 16. Analyze a Real GPU Hit from Log Paste

```gp
p = 0x4300E6A5560F71F54B7E593;
R = cc_root(p);
root = R[1];

cc_profile(root, 1, 1)
cc_shadow(root)
cc_chain_verify(root)
```

---

## 17. Minimal Manual Triage Workflow

```gp
p = 0x4300E6A5560F71F54B7E593;
cc_candidate_profile(p, 1, 0)
cc_motif(p, 1, 0)
cc_candidate_positions(p, 10, 1, 0)
```

---

## 18. OpenSSL Trial Prime Analysis

See how many BPSW tests we save by using more trial primes than OpenSSL.

```gp
\\ At 4096 bits: OpenSSL uses 1024 primes, we use ~9805
cc_sp_trial_analysis(4096)

\\ Key output:
\\   OpenSSL: 1 in 98 survive trial → most candidates hit BPSW
\\   Ours:    1 in 160 survive trial → 39% fewer BPSW tests
```

---

## 19. OpenSSL Delta-Sieve — Head-to-Head Comparison

```gp
\\ Quick comparison at 512 bits (5 runs each)
cc_sp_compare(512, 5)

\\ Longer comparison at 1024 bits
cc_sp_compare(1024, 3)
```

---

## 20. OpenSSL-Style Delta-Sieve — Direct Use

```gp
\\ Delta-sieve with OpenSSL's default trial count
cc_sp_delta_sieve(512)

\\ Delta-sieve with 1000 trial primes (our recommendation)
cc_sp_delta_sieve(512, 1000)

\\ Silent mode for scripting
result = cc_sp_delta_sieve(512, 0, 10000000, 0, 0);
if(result, printf("q = %s\n", Str(result[1])));
```

---

## 21. Bitmask Sieve — Fast Candidate Screening

Build bitmask tables once, test candidates with O(1) `bittest()` per prime.

```gp
\\ Build masks for CC18 screening
sp = [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47];
masks = vector(#sp, i, cc_x_forbidden_mask(sp[i], 18));

\\ Screen a known CC7a root — should pass for depth 7
\\ (may fail for depth 18 since it's only CC7)
p = 1122659;
survived = 1;
for(i = 1, #sp, if(bittest(masks[i], p % sp[i]), survived = 0; break));
printf("1122659 passes CC18 bitmask screen: %s\n", if(survived, "YES", "NO"))

\\ Compare vector vs bitmask approach
mask5 = cc_x_forbidden_mask(5, 3);
vec5 = cc_forbidden_res(5, 3, 1);
printf("Mask: %d, Vector: %s\n", mask5, Str(vec5))
\\ Both encode the same information: residues {0, 2, 3} are forbidden
```

---

## 22. Periodic Forbidden Tables — Short Periods

Some primes have surprisingly short multiplicative order of 2,
meaning the forbidden residue pattern repeats quickly.

```gp
\\ Survey all primes 3..100: find the ones with short periods
printf("Prime  ord_2  (shorter period = easier to handle)\n");
forprime(q = 3, 100,
  my(pt = cc_x_periodic_table(q));
  printf("  %3d  %3d%s\n", q, pt[2], if(pt[2] <= 10, "  ← short!", ""))
)

\\ Prime 31 has period 5 — only 5 unique forbidden residues
pt31 = cc_x_periodic_table(31);
printf("Prime 31: period=%d, table=%s\n", pt31[2], Str(pt31[1]))
\\ For CC18, positions 0..17 mod 5 cycle through the same 5 residues

\\ Prime 73 has period 9
pt73 = cc_x_periodic_table(73);
printf("Prime 73: period=%d\n", pt73[2])
```

---

## 23. Sieve-Only Line Filter — Pre-Screen Before BPSW

`cc_x_line_filter` follows the chain using only trial division.
Much cheaper than `isprime()` at each step.

```gp
\\ Known CC7a root: survives 7 sieve steps
cc_x_line_filter(1122659, 7)     \\ 7

\\ Try deeper: does it survive 10 sieve steps?
cc_x_line_filter(1122659, 10)    \\ probably < 10 (chain is only 7 long)

\\ Use as pre-filter in a search loop
count_tested = 0; count_survived = 0;
for(k = 1, 10000,
  my(p = 6*k + 5);          \\ candidates ≡ 5 mod 6
  if(cc_x_line_filter(p, 5) >= 5,
    count_survived++;
    \\ only NOW spend time on isprime()
  );
  count_tested++;
);
printf("Tested %d, survived sieve: %d (%.2f%%)\n",
       count_tested, count_survived, 100.0*count_survived/count_tested)
```

---

## 24. BiTwin Chain Discovery

BiTwin chains require BOTH a first-kind chain from n-1 AND a second-kind
chain from n+1 to be simultaneously prime.

```gp
\\ Quick check: center 6
cc_x_bitwin_walk(6)    \\ [4, 2, 2] — CC4a from 5, CC2b from 7

\\ Scan for best BiTwin in range
best = [0, 0];
forstep(n = 2, 50000, 2,
  my(bt = cc_x_bitwin_walk(n));
  if(bt[3] > best[2], best = [n, bt[3]];
    printf("  New best: center=%d BiTwin=%d (CC%da + CC%db)\n",
           n, bt[3], bt[1], bt[2])
  )
);
printf("\nBest BiTwin center in [2, 50000]: %d (BT%d)\n", best[1], best[2])

\\ Deep analysis of the best find
bt = cc_x_bitwin_walk(best[1]);
printf("\nFirst-kind (CC%da) from %d:\n", bt[1], best[1]-1);
cc_walk(best[1] - 1, 1);
printf("\nSecond-kind (CC%db) from %d:\n", bt[2], best[1]+1);
cc_walk(best[1] + 1, 2);
```

---

## 25. BiTwin Mask — Constraint Analysis

How constrained is the BiTwin search space?

```gp
\\ For each small prime, how many residues survive for BiTwin depth 3?
forprime(q = 3, 47, my(mask = cc_x_bitwin_mask(q, 3)); my(surv = q - hammingweight(mask)); printf("  prime %2d: %d/%d survive (%.0f%%)\n", q, surv, q, 100.0*surv/q))

\\ Compare with single-kind CC (less constrained)
printf("\nBiTwin vs CC1 at depth 5:\n");
forprime(q = 3, 23, my(bt = hammingweight(cc_x_bitwin_mask(q, 5))); my(cc = hammingweight(cc_x_forbidden_mask(q, 5))); printf("  prime %2d: BiTwin forbids %d/%d, CC1 forbids %d/%d\n", q, bt, q, cc, q))
```

---

## 26. Primorial Algebraic Scan — Record-Style Search

Historical CC records used the form `p = k * P# * 2^m - 1`.

```gp
\\ Quick scan at 20 bits with 7# = 210
results = cc_x_primorial_scan(20, 4, 3);
printf("Found %d chains of CC3+ at 20 bits\n", #results);
for(i = 1, min(5, #results),
  printf("  p=%d CC%d\n", results[i][1], results[i][2])
)

\\ Longer scan at 30 bits with 11# = 2310
results = cc_x_primorial_scan(30, 5, 4, 50000);
printf("Found %d chains of CC4+ at 30 bits\n", #results);

\\ Analyze the best hit
if(#results > 0,
  my(best = results[1]);
  for(i = 2, #results, if(results[i][2] > best[2], best = results[i]));
  printf("Best: CC%d root %d\n", best[2], best[1]);
  cc_walk(best[1])
)
```

---

## 27. Root Depth — Fragment vs True Root Census

Every prime in a chain has a "depth" — how many steps back to the root.
Only depth-0 primes are true roots.

```gp
\\ Trace the chain 2 → 5 → 11 → 23 → 47
foreach([2, 5, 11, 23, 47], p,
  printf("  p=%2d  depth=%d  is_root=%d\n",
         p, cc_x_root_depth(p), cc_is_root(p))
)

\\ Census: search for CC5+ hits and classify
hits = cc_search(35, 5, 1, 0, 0, 100000, 0);
my(depths = vector(10));
for(i = 1, #hits,
  my(d = cc_x_root_depth(hits[i][1]));
  if(d >= 0 && d < 10, depths[d+1]++)
);
printf("CC5+ depth distribution (35-bit):\n");
for(d = 0, 9, if(depths[d+1] > 0,
  printf("  depth %d: %d hits%s\n", d, depths[d+1], if(d==0, " (true roots)", ""))
))
```

---

## 28. Combined Pipeline — Line Filter + Bitmask + Verify

The full Layer 10 pipeline: bitmask screen → line filter → proven verify.

```gp
\\ Setup: build bitmask table for CC7 screening
sp = [5, 7, 11, 13, 17, 19, 23];
masks = vector(#sp, i, cc_x_forbidden_mask(sp[i], 7));

\\ Three-stage pipeline
stage1 = 0; stage2 = 0; stage3 = 0;
for(k = 1, 100000, my(p = 6*k + 5); my(pass = 1); for(i = 1, #sp, if(bittest(masks[i], p % sp[i]), pass = 0; break)); if(!pass, next); stage1++; if(cc_x_line_filter(p, 7) < 7, next); stage2++; my(len = cc_chain_test(p)); if(len >= 7, stage3++; printf("  CC%d root: %d (depth=%d)\n", len, p, cc_x_root_depth(p))))
printf("\nPipeline: 100K candidates -> %d bitmask -> %d line-filter -> %d CC7+\n", stage1, stage2, stage3)
```

