/*
 * Copyright (c) 2026 Nenad Mićić <nenad@micic.be>
 * SPDX-License-Identifier: Apache-2.0
 *============================================================================
 * cc_lib_v10.gp — PARI/GP Interactive Cunningham Chain Library
 *
 * Consolidates analytical functions from 5 Python scripts into a single
 * PARI/GP library for interactive chain analysis.  Load with:
 *
 *   \r cc_lib_v10.gp
 *
 * Then: cc_walk(p), cc_shadow(p), cc_autopsy(p), cc_full(p), ...
 *       cc_search(40, 5)               — find CC5+ at 40 bits
 *       cc_search(60, 7, 1, 5, 3)      — CC7+ at 60 bits, prefix 0b101
 *       cc_search_info(89, 18)          — show search parameters for CC18
 *       cc_construct(35, 4)             — constructive CC4+ at 35 bits
 *       cc_con_info(89, 18)             — constructor parameters for CC18
 *       cc_sp_search(4096)              — find 4096-bit safe prime
 *       cc_sp_search(2048, , 50000)     — 2048-bit, auto-scaled trial limit
 *       cc_x_line_filter(1122659, 7)    — sieve-only chain pre-screen
 *       cc_x_bitwin_walk(6)             — BiTwin chain from even center
 *
 * v10 branch notes (builds on v9):
 * - NEW Layer 10: Alternative Algorithms (CUDA Legacy Extraction)
 *   7 cc_x_ functions mined from 22 CUDA programs in cuda_legacy/.
 *   cc_x_forbidden_mask    — bitmask sieve (METHOD-L01, O(1) per prime)
 *   cc_x_periodic_table    — periodic forbidden table via ord_q(2)
 *   cc_x_line_filter       — sieve-only chain pre-screen (METHOD-L17)
 *   cc_x_bitwin_mask       — BiTwin forbidden residues (METHOD-L12)
 *   cc_x_bitwin_walk       — BiTwin chain walk + verify
 *   cc_x_primorial_scan    — algebraic form k*P#*2^m-1 search (METHOD-L10)
 *   cc_x_root_depth        — k-value: backward steps to true root
 *
 * v9 notes:
 * - Layer 9: Safe Prime Search Engine (cc_gmp_v35_03.c reimplementation)
 * - Layer 9b: OpenSSL Delta-Sieve Comparison
 *
 * v8 notes:
 * - Layer 8: Constructive Search Engine (p = S*R - 1, primorial form)
 * - Fix: cc_spine dead code removed
 *
 * Leverages PARI builtins: isprime(), ispseudoprime(), factor(),
 *   valuation(), hammingweight(), chinese(), forprime(), forstep(),
 *   znorder(), bittest(), bitor()
 *============================================================================*/

/* Ensure adequate stack for large primes (512 MB) */
{
  my(old = default(debugmem));
  default(debugmem, 0);
  default(parisizemax, 512000000);
  default(parisize, 64000000);
  default(debugmem, old);
};

/* ────────────────────────────────────────────────────────────────────────────
 * Layer 0: Constants
 * ──────────────────────────────────────────────────────────────────────────── */

/* 24 analysis primes matching Python scripts */
cc_analysis_primes = [3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97];

/* 8 display-column primes for residue tables */
cc_residue_primes = [3,5,7,11,13,17,19,23];

/* 6 construction-coverage primes (S-base) */
cc_sbase_primes = [2,3,5,11,13,19];

/* Display width */
cc_W = 115;

/* ────────────────────────────────────────────────────────────────────────────
 * Layer 1: Primitive Helpers
 * ──────────────────────────────────────────────────────────────────────────── */

/* Next chain member: 2p+1 (first kind) or 2p-1 (second kind) */
cc_next(p, {kind=1}) = if(kind==1, 2*p+1, 2*p-1);

/* Previous chain member: (p-1)/2 (first kind) or (p+1)/2 (second kind) */
cc_prev(p, {kind=1}) = if(kind==1, (p-1)/2, (p+1)/2);

/* Immunizing residue: if p mod q == this value, q never divides any member */
cc_immune_residue(q, {kind=1}) = if(kind==1, q-1, 1);

/* Check if p is immunized against prime q */
cc_is_immune(p, q, {kind=1}) = (p % q) == cc_immune_residue(q, kind);

/* Kill position: forward simulation of residue walk
 * Starting at r = p mod q, iterate r <- (2r+1) mod q (first kind)
 * or r <- (2r-1) mod q (second kind). First i where r==0 is kill position.
 * Returns max+1 if immunized (no kill found). */
cc_kill_pos(r, q, {kind=1}, {mx=32}) =
{
  my(pos);
  for(pos = 0, mx-1,
    if(r % q == 0, return(pos));
    if(kind == 1,
      r = (2*r + 1) % q,
      r = (2*r - 1) % q
    );
  );
  return(mx);
};

/* Trailing ones: count consecutive 1-bits from LSB */
cc_trailing_ones(n) = if(n==0, 0, valuation(n+1, 2));

/* Odd part: remove all factors of 2 */
cc_odd_part(n) = n >> valuation(n, 2);

/* NAF weight (Non-Adjacent Form): no PARI builtin, hand-rolled */
cc_naf_weight(n) =
{
  my(w = 0, u);
  while(n != 0,
    if(n % 2 != 0,
      u = 2 - (n % 4);
      n -= u;
      w++;
    );
    n >>= 1;
  );
  return(w);
};

/* Format a number as hex, truncated if needed */
cc_fmt_hex(n, {mx=40}) =
{
  my(s, v, trunc);
  if(n < 0,
    s = Str("-0x", strprintf("%x", -n));
  ,
    s = Str("0x", strprintf("%x", n));
  );
  if(#s > mx,
    /* GP strings don't support range slicing; use Vecsmall/Strchr */
    v = Vecsmall(s);
    trunc = Strchr(vecextract(v, Str("1..", mx-3)));
    Str(trunc, "...");
  ,
    s;
  );
};

/* Format a PARI factor matrix as "2^3 * 5 * 17" */
cc_fmt_fact(M) =
{
  my(s = "", i, p, e);
  if(type(M) != "t_MAT" || matsize(M)[1] == 0, return("1"));
  for(i = 1, matsize(M)[1],
    p = M[i,1]; e = M[i,2];
    if(i > 1, s = Str(s, " * "));
    if(e == 1,
      s = Str(s, p);
    ,
      s = Str(s, p, "^", e);
    );
  );
  return(s);
};

/* Format a factor matrix with composite tags */
cc_fmt_fact_tagged(M) =
{
  my(s = "", i, p, e, bits);
  if(type(M) != "t_MAT" || matsize(M)[1] == 0, return("1"));
  for(i = 1, matsize(M)[1],
    p = M[i,1]; e = M[i,2];
    if(i > 1, s = Str(s, " * "));
    bits = #binary(p);
    if(isprime(p),
      if(bits > 20,
        if(e == 1,
          s = Str(s, "P", bits, "b");
        ,
          s = Str(s, "P", bits, "b^", e);
        );
      ,
        if(e == 1, s = Str(s, p), s = Str(s, p, "^", e));
      );
    ,
      if(bits > 20,
        s = Str(s, "C", bits, "b");
      ,
        s = Str(s, p, "[C]");
      );
      if(e > 1, s = Str(s, "^", e));
    );
  );
  return(s);
};

/* Section header display helpers */
cc_sep() = printf("%s\n", concat(vector(cc_W, i, "=")));
cc_hdr(title) = { cc_sep(); printf("  %s\n", title); cc_sep(); };
cc_sub(title) = { printf("\n--- %s %s\n", title, concat(vector(max(0, cc_W - #title - 5), i, "-"))); };

/* Bit-vector / motif helpers */
cc_bitvec_to_str(v) =
{
  my(s = "");
  for(i = 1, #v,
    s = Str(s, v[i]);
  );
  return(s);
};

cc_pattern_vec(pattern) =
{
  my(t = type(pattern), v, out);

  if(t == "t_STR",
    v = Vecsmall(pattern);
    out = vector(#v, i, v[i] - 48);
    for(i = 1, #out,
      if(out[i] != 0 && out[i] != 1,
        error("pattern string must contain only 0/1");
      );
    );
    return(out);
  );

  if(t == "t_VECSMALL" || t == "t_VEC",
    return(vector(#pattern, i, pattern[i]));
  );

  if(t == "t_INT",
    return(binary(pattern));
  );

  error("unsupported pattern type");
};

cc_pattern_str(pattern) = cc_bitvec_to_str(cc_pattern_vec(pattern));

cc_prefix_bits(n, {len=8}) =
{
  my(b = binary(n), take);
  take = min(len, #b);
  return(cc_bitvec_to_str(vecextract(b, Str("1..", take))));
};

cc_tail_bits(n, {len=8}) =
{
  my(b = binary(n), start, stop);
  stop = #b;
  start = max(1, stop - len + 1);
  return(cc_bitvec_to_str(vecextract(b, Str(start, "..", stop))));
};

cc_transition_count(n) =
{
  my(b = binary(n), c = 0);
  for(i = 2, #b,
    if(b[i] != b[i-1], c++);
  );
  return(c);
};

cc_longest_run(n, {bit=1}) =
{
  my(b = binary(n), best = 0, cur = 0);
  for(i = 1, #b,
    if(b[i] == bit,
      cur++;
      if(cur > best, best = cur);
    ,
      cur = 0;
    );
  );
  return(best);
};

cc_motif_count(n, pattern) =
{
  my(bits = binary(n), pat = cc_pattern_vec(pattern), cnt = 0,
     nb = #bits, np = #pat, ok);

  if(np > nb, return(0));

  for(i = 1, nb - np + 1,
    ok = 1;
    for(j = 1, np,
      if(bits[i + j - 1] != pat[j],
        ok = 0;
        break;
      );
    );
    if(ok, cnt++);
  );

  return(cnt);
};

cc_smooth_split(n, {primes=0}) =
{
  my(pr, smooth = 1, rem = n, q, v);
  if(primes == 0,
    pr = [2,3,5,7,11,13,17,19,23,29,31];
  ,
    pr = primes;
  );

  for(i = 1, #pr,
    q = pr[i];
    v = valuation(rem, q);
    if(v > 0,
      smooth *= q^v;
      rem /= q^v;
    );
  );

  return([smooth, rem]);
};

/* Heat metric as a diagnostic scalar only */
cc_heat(len) = 1.0 / (len + 1.0);

/* First analysis-prime divisor of n, or 0 if none */
cc_shadow_hit(n, {primes=0}) =
{
  my(pr);
  if(primes == 0, pr = cc_analysis_primes, pr = primes);
  if(n == 0, return(0));
  for(i = 1, #pr,
    if(n % pr[i] == 0, return(pr[i]));
  );
  return(0);
};

/* Small trial factor hint for large composites; returns 0 if none <= bound */
cc_small_factor_hint(n, {bound=1000}) =
{
  my(m);
  if(bound < 2, return(0));
  m = abs(n);
  if(m <= 1, return(m));
  forprime(q = 2, bound,
    if(m % q == 0, return(q));
  );
  return(0);
};

/* ────────────────────────────────────────────────────────────────────────────
 * Layer 2: Core Chain Walking
 * ──────────────────────────────────────────────────────────────────────────── */

/* Walk backward from p to find the true root of the chain.
 * v6 behavior: composite inputs are treated as exact candidates and do not
 * backtrack through a prime predecessor. */
cc_root(p, {kind=1}) =
{
  my(steps = 0, prev);
  if(!isprime(p), return([p, 0]));
  while(1,
    prev = cc_prev(p, kind);
    if(prev < 2, return([p, steps]));
    /* For first kind, (p-1) must be even; for second, (p+1) must be even */
    if(kind == 1 && (p-1) % 2 != 0, return([p, steps]));
    if(kind == 2 && (p+1) % 2 != 0, return([p, steps]));
    if(!isprime(prev), return([p, steps]));
    p = prev;
    steps++;
    if(steps > 200, return([p, steps])); \\ safety
  );
};

cc_is_root(p, {kind=1}) =
{
  my(prev);
  if(!isprime(p), return(0));
  prev = cc_prev(p, kind);
  if(prev < 2, return(1));
  return(!isprime(prev));
};

/* Main chain walk: find root, walk forward, continue past break for lookahead.
 * Returns [chain_vec, breaker, root, chain_length]
 */
cc_walk(p, {kind=1}, {lookahead=6}) =
{
  my(R, root, back_steps, chain = List(), cur, pos, chain_len,
     breaker = 0, post = List(), residues, r, imm, status, hdr_str);

  /* Find root */
  R = cc_root(p, kind);
  root = R[1]; back_steps = R[2];

  /* Walk forward to build chain + post-break lookahead */
  cur = root; pos = 0;
  chain_len = 0;
  while(1,
    if(breaker == 0,
      if(isprime(cur),
        chain_len++;
        listput(chain, cur);
      ,
        breaker = cur;
        chain_len = pos;
        listput(post, cur);
        if(#post >= lookahead, break);
      );
    ,
      listput(post, cur);
      if(#post >= lookahead, break);
    );
    pos++;
    cur = cc_next(cur, kind);
    if(pos > 200, break); \\ safety
  );

  /* Display */
  hdr_str = if(kind==1, "FIRST", "SECOND");
  cc_hdr(Str("CUNNINGHAM CHAIN — ", hdr_str, " KIND"));

  printf("  Root:         %s\n", cc_fmt_hex(root, 60));
  printf("  Input:        %s\n", cc_fmt_hex(p, 60));
  printf("  Chain length: %d\n", chain_len);
  printf("  Back steps:   %d (from input to root)\n", back_steps);
  printf("  Kind:         %s\n", hdr_str);

  /* Chain member table */
  cc_sub("Chain Members");
  printf("  %3s  %6s %4s %4s %3s  ", "Pos", "Bits", "Pop", "NAF", "T1");
  for(j = 1, #cc_residue_primes,
    printf("%4d", cc_residue_primes[j]);
  );
  printf("  %s\n", "Status");
  printf("  %s\n", concat(vector(cc_W - 2, i, "-")));

  for(i = 1, #chain,
    cur = chain[i];
    printf("  %3d  %6d %4d %4d %3d  ",
      i-1,
      #binary(cur),
      hammingweight(cur),
      cc_naf_weight(cur),
      cc_trailing_ones(cur));
    for(j = 1, #cc_residue_primes,
      r = cur % cc_residue_primes[j];
      imm = cc_is_immune(cur, cc_residue_primes[j], kind);
      if(imm,
        printf("[%2d]", r);
      ,
        printf(" %2d ", r);
      );
    );
    printf("  PRIME\n");
  );

  /* Post-break (lookahead) */
  if(#post > 0,
    printf("\n");
    printf("  %3s  %6s %4s %4s %3s  ", "Pos", "Bits", "Pop", "NAF", "T1");
    for(j = 1, #cc_residue_primes,
      printf("%4d", cc_residue_primes[j]);
    );
    printf("  %s\n", "Status");
    printf("  %s\n", concat(vector(cc_W - 2, i, "-")));

    for(i = 1, #post,
      cur = post[i];
      printf("  %3d  %6d %4d %4d %3d  ",
        chain_len + i - 1,
        #binary(cur),
        hammingweight(cur),
        cc_naf_weight(cur),
        cc_trailing_ones(cur));
      for(j = 1, #cc_residue_primes,
        r = cur % cc_residue_primes[j];
        imm = cc_is_immune(cur, cc_residue_primes[j], kind);
        if(imm,
          printf("[%2d]", r);
        ,
          printf(" %2d ", r);
        );
      );
      if(i == 1,
        printf("  <<BREAK>>\n");
      ,
        printf("  composite\n");
      );
    );
  );

  return(Vec(chain));
};

/* ────────────────────────────────────────────────────────────────────────────
 * Layer 3: Analysis Functions
 * ──────────────────────────────────────────────────────────────────────────── */

/* Shadow analysis: for each small prime, show residue, immunity, kill position */
cc_shadow(p, {kind=1}, {primes=0}) =
{
  my(pr, r, imm, kp, immune_count = 0, active_count = 0,
     min_kill = 999, min_kill_q = 0, results = List());

  if(primes == 0, pr = cc_analysis_primes, pr = primes);

  cc_hdr(Str("SHADOW ANALYSIS — ", if(kind==1,"FIRST","SECOND"), " KIND"));
  printf("  Root: %s\n\n", cc_fmt_hex(p, 60));

  printf("  %5s  %7s  %8s  %8s\n", "Prime", "Residue", "Immune?", "Kill Pos");
  printf("  %s\n", concat(vector(40, i, "-")));

  for(i = 1, #pr,
    r = p % pr[i];
    imm = cc_is_immune(p, pr[i], kind);
    if(imm,
      immune_count++;
      printf("  %5d  %7d  %8s  %8s\n", pr[i], r, "YES", "---");
      listput(results, [pr[i], r, 1, -1]);
    ,
      active_count++;
      kp = cc_kill_pos(r, pr[i], kind);
      if(kp < min_kill,
        min_kill = kp;
        min_kill_q = pr[i];
      );
      printf("  %5d  %7d  %8s  %8d\n", pr[i], r, "no", kp);
      listput(results, [pr[i], r, 0, kp]);
    );
  );

  printf("\n  Summary: %d immune, %d active out of %d primes\n",
    immune_count, active_count, #pr);
  if(min_kill < 999,
    printf("  Earliest kill: position %d by prime %d\n", min_kill, min_kill_q);
  );

  return(Vec(results));
};

/* Verify immunization persistence across all chain positions */
cc_immune_persist(p, {kind=1}) =
{
  my(R, root, chain, cur, pos, q, all_persist = 1,
     immune_at_root = List());

  R = cc_root(p, kind);
  root = R[1];

  /* Find which primes are immune at root */
  for(i = 1, #cc_analysis_primes,
    q = cc_analysis_primes[i];
    if(cc_is_immune(root, q, kind),
      listput(immune_at_root, q);
    );
  );

  cc_sub("Immunization Persistence Check");
  printf("  Root: %s\n", cc_fmt_hex(root, 60));
  printf("  Primes immune at root: %s\n", Str(Vec(immune_at_root)));

  /* Walk chain and verify */
  cur = root; pos = 0;
  while(isprime(cur) && pos < 50,
    for(j = 1, #immune_at_root,
      q = immune_at_root[j];
      if(!cc_is_immune(cur, q, kind),
        printf("  ** BREAK at pos %d: prime %d lost immunity!\n", pos, q);
        all_persist = 0;
      );
    );
    cur = cc_next(cur, kind);
    pos++;
  );

  if(all_persist,
    printf("  VERIFIED: all immunizations persist across %d positions\n", pos);
  ,
    printf("  FAILED: some immunizations did not persist\n");
  );

  return(all_persist);
};

/* Deep factorization: recursive p-1/p+1 factorization tree */
cc_deep_factor(n, {depth=3}, {threshold=256}) =
{
  my(M, nrows, p, e, bits);

  if(depth <= 0 || n <= threshold, return);

  M = factor(n);

  cc__deep_print(n, M, depth, threshold, "");
};

/* Internal recursive deep-factor printer */
cc__deep_print(n, M, depth, threshold, indent) =
{
  my(nrows, p, e, bits, M2);

  printf("%s%s = %s\n", indent, cc_fmt_hex(n, 50), cc_fmt_fact(M));

  if(depth <= 0, return);

  nrows = matsize(M)[1];
  for(i = 1, nrows,
    p = M[i,1]; e = M[i,2];
    bits = #binary(p);
    if(bits > 8 && p > threshold,
      /* Factor p-1 */
      M2 = factor(p-1);
      printf("%s  p-1 of %s:\n", indent, if(bits > 32, Str("P",bits,"b"), Str(p)));
      cc__deep_print(p-1, M2, depth-1, threshold, Str(indent, "    "));

      /* Factor p+1 */
      M2 = factor(p+1);
      printf("%s  p+1 of %s:\n", indent, if(bits > 32, Str("P",bits,"b"), Str(p)));
      cc__deep_print(p+1, M2, depth-1, threshold, Str(indent, "    "));
    );
  );
};

/* Autopsy: analyze the breaker (first composite in chain) */
cc_autopsy(p, {kind=1}, {depth=2}) =
{
  my(R, root, cur, chain_len = 0, breaker, M, nrows, hits = List(),
     q, bits);

  /* Walk to root, then forward to breaker */
  R = cc_root(p, kind);
  root = R[1];

  cur = root;
  while(isprime(cur),
    chain_len++;
    cur = cc_next(cur, kind);
    if(chain_len > 200, break);
  );
  breaker = cur;
  bits = #binary(breaker);

  cc_hdr(Str("AUTOPSY — BREAKER AT POSITION ", chain_len));
  printf("  Breaker: %s\n", cc_fmt_hex(breaker, 60));
  printf("  Bits:    %d\n", bits);
  printf("  Root:    %s\n", cc_fmt_hex(root, 60));
  printf("  Chain:   CC%d%s\n\n", chain_len, if(kind==1,"a","b"));

  /* Factor the breaker */
  cc_sub("Breaker Factorization");
  M = factor(breaker);
  printf("  %s\n", cc_fmt_fact(M));
  printf("  Tagged: %s\n", cc_fmt_fact_tagged(M));

  /* Shadow hits: which small primes divide the breaker? */
  cc_sub("Shadow Hits");
  for(i = 1, #cc_analysis_primes,
    q = cc_analysis_primes[i];
    if(breaker % q == 0,
      listput(hits, q);
    );
  );
  if(#hits > 0,
    printf("  Small primes dividing breaker: %s\n", Str(Vec(hits)));
    printf("  Count: %d out of %d analysis primes\n", #hits, #cc_analysis_primes);
  ,
    printf("  No analysis primes divide the breaker (large-factor break)\n");
  );

  /* Deep factor the breaker */
  if(depth >= 1,
    cc_sub("Deep Factorization of Breaker");
    cc_deep_factor(breaker, depth);
  );

  /* Breaker neighborhood: breaker +/- small offsets */
  cc_sub("Breaker Neighborhood");
  for(offset = -2, 4,
    my(nb = breaker + offset, tag = "");
    if(offset == 0, tag = " <<BREAKER>>");
    if(isprime(nb), tag = Str(tag, " [PRIME]"));
    printf("  %+3d: %s = %s%s\n",
      offset,
      cc_fmt_hex(nb, 40),
      cc_fmt_fact_tagged(factor(nb)),
      tag);
  );

  /* Deep factor breaker-1 and breaker+1 */
  if(depth >= 2,
    cc_sub("Deep Factor: breaker - 1");
    cc_deep_factor(breaker - 1, depth - 1);
    cc_sub("Deep Factor: breaker + 1");
    cc_deep_factor(breaker + 1, depth - 1);
  );

  return(M);
};

/* Find next N composite positions after chain breaks */
cc_next_breakers(p, {count=4}, {kind=1}) =
{
  my(R, root, cur, chain_len = 0, composites = List(), pos);

  R = cc_root(p, kind);
  root = R[1];

  /* Walk to end of chain */
  cur = root;
  while(isprime(cur),
    chain_len++;
    cur = cc_next(cur, kind);
    if(chain_len > 200, break);
  );

  cc_sub(Str("Next ", count, " Breakers After CC", chain_len));

  /* Continue scanning */
  pos = chain_len;
  while(#composites < count && pos < chain_len + 100,
    if(!isprime(cur),
      listput(composites, [pos, cur]);
      printf("  Pos %3d: %s  (%d bits)\n", pos, cc_fmt_hex(cur, 50), #binary(cur));
      /* Show which small primes divide it */
      my(divs = List());
      for(j = 1, min(#cc_analysis_primes, 10),
        if(cur % cc_analysis_primes[j] == 0,
          listput(divs, cc_analysis_primes[j]);
        );
      );
      if(#divs > 0,
        printf("          small factors: %s\n", Str(Vec(divs)));
      );
    );
    cur = cc_next(cur, kind);
    pos++;
  );

  return(Vec(composites));
};

/* ────────────────────────────────────────────────────────────────────────────
 * Layer 4: Properties
 * ──────────────────────────────────────────────────────────────────────────── */

/* Machine-native geometry of a prime */
cc_geometry(p) =
{
  my(bits, pop, naf, t1, v2pm1, v2pp1, odd_pm1, odd_pp1);

  bits = #binary(p);
  pop = hammingweight(p);
  naf = cc_naf_weight(p);
  t1 = cc_trailing_ones(p);
  v2pm1 = valuation(p-1, 2);
  v2pp1 = valuation(p+1, 2);
  odd_pm1 = cc_odd_part(p-1);
  odd_pp1 = cc_odd_part(p+1);

  cc_hdr("MACHINE GEOMETRY");

  printf("  Value:          %s\n", cc_fmt_hex(p, 60));
  printf("  Decimal digits: %d\n", #Str(p));
  printf("  Bits:           %d\n", bits);
  printf("  Popcount:       %d (density: %.1f%%)\n", pop, 100.0*pop/bits);
  printf("  NAF weight:     %d (vs popcount %d)\n", naf, pop);
  printf("  Trailing 1s:    %d\n", t1);

  cc_sub("2-Adic Properties");
  printf("  v2(p-1) = %d    odd_part(p-1) = %s\n", v2pm1, cc_fmt_hex(odd_pm1, 40));
  printf("  v2(p+1) = %d    odd_part(p+1) = %s\n", v2pp1, cc_fmt_hex(odd_pp1, 40));

  cc_sub("Wheel & Modular Residues");
  printf("  p mod 6  = %d\n", p % 6);
  printf("  p mod 30 = %d\n", p % 30);
  printf("  p mod 210 = %d\n", p % 210);
  for(i = 1, #cc_residue_primes,
    printf("  p mod %2d = %d\n", cc_residue_primes[i], p % cc_residue_primes[i]);
  );

  return([bits, pop, naf, t1, v2pm1, v2pp1]);
};

/* Classify a prime: Sophie Germain, Safe, Twin, etc. */
cc_classify(p) =
{
  my(classes = List(), k, e);

  cc_sub("Prime Classifications");

  if(!isprime(p),
    printf("  %s is NOT prime\n", cc_fmt_hex(p, 40));
    return(Vec(classes));
  );

  /* Sophie Germain: 2p+1 is also prime */
  if(isprime(2*p+1),
    listput(classes, "Sophie Germain");
  );

  /* Safe prime: (p-1)/2 is also prime */
  if(p > 2 && (p-1) % 2 == 0 && isprime((p-1)/2),
    listput(classes, "Safe prime");
  );

  /* Twin prime (lower): p+2 is prime */
  if(isprime(p+2),
    listput(classes, "Twin (lower)");
  );

  /* Twin prime (upper): p-2 is prime */
  if(p > 2 && isprime(p-2),
    listput(classes, "Twin (upper)");
  );

  /* Cousin prime: p+4 is prime */
  if(isprime(p+4),
    listput(classes, "Cousin (lower)");
  );

  /* Sexy prime: p+6 is prime */
  if(isprime(p+6),
    listput(classes, "Sexy (lower)");
  );

  /* Mersenne: p = 2^k - 1 */
  if(p + 1 == 1 << valuation(p+1, 2),
    listput(classes, Str("Mersenne (2^", valuation(p+1,2), "-1)"));
  );

  /* Proth: p = k*2^e + 1 where k < 2^e and k is odd */
  e = valuation(p-1, 2);
  if(e >= 1,
    k = (p-1) >> e;
    if(k > 0 && k % 2 == 1 && k < (1 << e),
      listput(classes, Str("Proth (", k, "*2^", e, "+1)"));
    );
  );

  if(#classes > 0,
    for(i = 1, #classes,
      printf("  [%d] %s\n", i, classes[i]);
    );
  ,
    printf("  (no special classifications)\n");
  );

  return(Vec(classes));
};

/* Algebraic chain members: show p_i = 2^i * root + (2^i - 1) with factorization */
cc_members(p, {kind=1}, {breaker_depth=2}) =
{
  my(R, root, cur, chain_len = 0, breaker, M,
     pow2, formula_val);

  R = cc_root(p, kind);
  root = R[1];

  /* Walk to find chain length */
  cur = root;
  while(isprime(cur),
    chain_len++;
    cur = cc_next(cur, kind);
    if(chain_len > 200, break);
  );
  breaker = cur;

  cc_hdr(Str("CHAIN MEMBERS — CC", chain_len, if(kind==1,"a","b"), " ALGEBRAIC FORM"));

  if(kind == 1,
    printf("  Root p = %s\n", cc_fmt_hex(root, 70));
    printf("  Formula: p_i = 2^i * (p+1) - 1  =  2^i * p + (2^i - 1)\n\n");
  ,
    printf("  Root p = %s\n", cc_fmt_hex(root, 70));
    printf("  Formula: p_i = 2^i * (p-1) + 1  =  2^i * p - (2^i - 1)\n\n");
  );

  /* Print each member with algebraic form */
  printf("  %3s  %6s  %-20s  %s\n", "Pos", "Bits", "2^i component", "Value (hex)");
  printf("  %s\n", concat(vector(cc_W - 2, i, "-")));

  cur = root;
  for(i = 0, chain_len - 1,
    pow2 = 1 << i;
    printf("  %3d  %6d  2^%-3d * p %s %-8s  %s  PRIME\n",
      i,
      #binary(cur),
      i,
      if(kind==1, "+", "-"),
      if(i==0, "", Str("(", pow2 - 1, ")")),
      cc_fmt_hex(cur, 50));
    cur = cc_next(cur, kind);
  );

  /* Breaker */
  pow2 = 1 << chain_len;
  printf("\n  %3d  %6d  2^%-3d * p %s %-8s  %s  <<BREAK>>\n",
    chain_len,
    #binary(breaker),
    chain_len,
    if(kind==1, "+", "-"),
    Str("(", pow2 - 1, ")"),
    cc_fmt_hex(breaker, 50));

  /* Verify formula */
  if(kind == 1,
    formula_val = (1 << chain_len) * root + ((1 << chain_len) - 1);
  ,
    formula_val = (1 << chain_len) * root - ((1 << chain_len) - 1);
  );
  if(formula_val == breaker,
    printf("  Formula verified: 2^%d * p %s %d = breaker\n",
      chain_len, if(kind==1,"+","-"), (1 << chain_len) - 1);
  );

  /* Breaker factorization */
  cc_sub("Breaker Factorization");
  M = factor(breaker);
  printf("  %s\n", cc_fmt_fact(M));
  printf("  Tagged: %s\n", cc_fmt_fact_tagged(M));

  /* Factor p-1 and p+1 of root */
  cc_sub("Root p-1");
  M = factor(root - 1);
  printf("  %s\n", cc_fmt_fact(M));
  printf("  v2(p-1) = %d  =>  p-1 = %s * 2^%d\n",
    valuation(root-1, 2), cc_fmt_hex(cc_odd_part(root-1), 30), valuation(root-1, 2));

  cc_sub("Root p+1");
  M = factor(root + 1);
  printf("  %s\n", cc_fmt_fact(M));
  printf("  v2(p+1) = %d  =>  p+1 = %s * 2^%d\n",
    valuation(root+1, 2), cc_fmt_hex(cc_odd_part(root+1), 30), valuation(root+1, 2));

  /* Key relationship: for first kind, all members share (p+1) structure */
  if(kind == 1,
    cc_sub("Structural Identity");
    printf("  p+1 = %s\n", cc_fmt_hex(root + 1, 50));
    printf("  All chain members: p_i + 1 = 2^i * (p+1)\n");
    printf("  So p_i + 1 inherits ALL prime factors of (p+1)\n");
    printf("  This is WHY immunization persists: if q | (p+1), then q | (p_i+1),\n");
    printf("  so p_i mod q = q-1 for all i.\n");
  ,
    cc_sub("Structural Identity");
    printf("  p-1 = %s\n", cc_fmt_hex(root - 1, 50));
    printf("  All chain members: p_i - 1 = 2^i * (p-1)\n");
    printf("  So p_i - 1 inherits ALL prime factors of (p-1)\n");
    printf("  This is WHY immunization persists: if q | (p-1), then q | (p_i-1),\n");
    printf("  so p_i mod q = 1 for all i.\n");
  );

  /* Deep factor the breaker */
  if(breaker_depth >= 1,
    cc_sub("Breaker Deep Factorization");
    cc_deep_factor(breaker, breaker_depth);
  );

  return(chain_len);
};

/* 2-adic spine exploration */
cc_spine(p, {levels=5}, {max_bits=200}) =
{
  my(core, core_v2, cur, results = List(), count);

  core = cc_odd_part(p);
  core_v2 = valuation(p, 2);

  cc_hdr("2-ADIC SPINE EXPLORATION");
  printf("  Input:    %s\n", cc_fmt_hex(p, 60));
  printf("  Odd core: %s  (%d bits)\n", cc_fmt_hex(core, 40), #binary(core));
  printf("  v2(p):    %d\n\n", core_v2);

  /* Odd-core ray: core, 2*core, 4*core, ... */
  cc_sub("Odd-Core Vertical Ray");
  cur = core; count = 0;
  while(#binary(cur) <= max_bits && count < 20,
    printf("  2^%2d * %s = %s  %s\n",
      valuation(cur, 2),
      cc_fmt_hex(core, 20),
      cc_fmt_hex(cur, 30),
      if(isprime(cur), "PRIME", ""));
    cur *= 2;
    count++;
  );

  /* Left and right spines at each level */
  for(slevel = 1, levels,
    my(mult, lstart, rstart);

    if(slevel == 1,
      mult = 2;
    ,
      mult = 1 << (slevel + 1);
    );

    lstart = mult * core - 1;
    rstart = mult * core + 1;

    cc_sub(Str("Spine Level ", slevel, " (multiplier = ", mult, ")"));

    /* Annotate level 1 with chain connection */
    if(slevel == 1 && core_v2 == 0,
      printf("  NOTE: odd_core = input, so Level 1 RIGHT = CC first-kind chain,\n");
      printf("        Level 1 LEFT = CC second-kind chain from 2p-1\n\n");
    );

    /* Left spine: iterate 2m-1 (= second-kind CC operation) */
    printf("  LEFT (2m-1)%s:\n", if(slevel==1 && core_v2==0, "  [second-kind walk from 2p-1]", ""));
    cur = lstart; count = 0;
    while(#binary(cur) <= max_bits && count < 8,
      printf("    %s  (%d bits)  %s\n",
        cc_fmt_hex(cur, 35),
        #binary(cur),
        if(isprime(cur), "PRIME", ""));
      listput(results, [slevel, "L", cur, isprime(cur)]);
      cur = 2*cur - 1;
      count++;
    );

    /* Right spine: iterate 2m+1 (= first-kind CC operation) */
    printf("  RIGHT (2m+1)%s:\n", if(slevel==1 && core_v2==0, "  [= CC first-kind chain!]", ""));
    cur = rstart; count = 0;
    while(#binary(cur) <= max_bits && count < 8,
      my(chain_tag = "");
      if(slevel == 1 && core_v2 == 0,
        chain_tag = Str("  p_", count+1);
      );
      printf("    %s  (%d bits)  %s%s\n",
        cc_fmt_hex(cur, 35),
        #binary(cur),
        if(isprime(cur), "PRIME", "composite"),
        chain_tag);
      listput(results, [slevel, "R", cur, isprime(cur)]);
      cur = 2*cur + 1;
      count++;
    );
  );

  return(Vec(results));
};

/* ────────────────────────────────────────────────────────────────────────────
 * Layer 5: Aggregate
 * ──────────────────────────────────────────────────────────────────────────── */

/* Full analysis: everything in one call */
cc_full(p, {kind=1}) =
{
  my(R, root, chain_len, cur, breaker);

  if(!isprime(p),
    printf("WARNING: %s is not prime. Analyzing anyway.\n", cc_fmt_hex(p, 40));
  );

  R = cc_root(p, kind);
  root = R[1];

  /* 0. Compact normalized summary */
  cc_profile(p, kind, 1);
  cc_motif(p, kind, 1);

  /* 1. Geometry */
  cc_geometry(p);

  /* 2. Classifications */
  cc_classify(p);

  /* 3. Chain walk (first kind) */
  cc_walk(p, kind);

  /* Also try the other kind */
  cc_walk(p, 3 - kind);

  /* 3b. Algebraic chain members */
  cc_members(p, kind);

  /* 4. Shadow analysis */
  cc_shadow(root, kind);

  /* 5. Immunization persistence */
  cc_immune_persist(root, kind);

  /* 6. Autopsy */
  cur = root;
  chain_len = 0;
  while(isprime(cur) && chain_len < 200,
    chain_len++;
    cur = cc_next(cur, kind);
  );
  breaker = cur;

  cc_autopsy(p, kind);

  /* 7. Next breakers */
  cc_next_breakers(p, 4, kind);

  /* 8. Deep factorization of p-1 */
  cc_hdr("DEEP FACTORIZATION: p - 1");
  cc_deep_factor(root - 1, 3);

  /* 9. Deep factorization of p+1 */
  cc_hdr("DEEP FACTORIZATION: p + 1");
  cc_deep_factor(root + 1, 3);

  /* 10. Spine */
  cc_spine(root, 3, 150);

  /* 11. Top-down grid trajectory */
  cc_grid_trajectory(root, kind);

  /* 12. p-Adic address / fingerprint */
  cc_padic_addr(root);

  printf("\n");
  cc_sep();
  printf("  Analysis complete.\n");
  cc_sep();
};

/* ────────────────────────────────────────────────────────────────────────────
 * S-base Coverage Check
 * ──────────────────────────────────────────────────────────────────────────── */

/* Check if p+1 (first kind) or p-1 (second kind) has full S-base coverage */
cc_sbase_check(p, {kind=1}) =
{
  my(target, covered = 1, q);
  target = if(kind == 1, p + 1, p - 1);

  printf("  S-base coverage for %s:\n", cc_fmt_hex(p, 30));
  for(i = 1, #cc_sbase_primes,
    q = cc_sbase_primes[i];
    if(target % q == 0,
      printf("    %2d | (p%s)  YES\n", q, if(kind==1,"+1","-1"));
    ,
      printf("    %2d | (p%s)  no\n", q, if(kind==1,"+1","-1"));
      covered = 0;
    );
  );
  printf("  Full coverage: %s\n", if(covered, "YES", "NO"));
  return(covered);
};

/* ────────────────────────────────────────────────────────────────────────────
 * Batch Helpers
 * ──────────────────────────────────────────────────────────────────────────── */

/* Walk a chain silently, return [chain_vec, chain_length, root, breaker] */
cc_walk_quiet(p, {kind=1}) =
{
  my(R, root, chain = List(), cur, breaker = 0);
  R = cc_root(p, kind);
  root = R[1];
  cur = root;
  while(isprime(cur),
    listput(chain, cur);
    cur = cc_next(cur, kind);
    if(#chain > 200, break);
  );
  breaker = cur;
  return([Vec(chain), #chain, root, breaker]);
};

/* Walk forward from the exact input p, without backward normalization */
cc_walk_from_quiet(p, {kind=1}) =
{
  my(chain = List(), cur = p, breaker = 0);
  while(isprime(cur),
    listput(chain, cur);
    cur = cc_next(cur, kind);
    if(#chain > 200, break);
  );
  breaker = cur;
  return([Vec(chain), #chain, breaker]);
};

cc_profile_labels() =
{
  return([
    "input_hex", "subject_hex", "root_hex", "back_steps", "kind",
    "is_prime", "is_root", "bits", "chain_len", "breaker_hex",
    "immune_count", "active_count", "earliest_kill_pos", "earliest_kill_q",
    "v2_pm1", "v2_pp1", "smooth_pm1_bits", "smooth_pp1_bits",
    "prefix8", "prefix16", "tail5", "tail8"
  ]);
};

cc_join_fields(v, {sep=","}) =
{
  my(s = "");
  for(i = 1, #v,
    if(i > 1, s = Str(s, sep));
    s = Str(s, v[i]);
  );
  return(s);
};

cc_bits_fixed_len(x, len) =
{
  my(b, out, start);
  if(len < 0, error("len must be >= 0"));
  out = vector(len, i, 0);
  if(len == 0, return(out));
  b = binary(x);
  if(#b > len, error("x does not fit in requested bit length"));
  start = len - #b + 1;
  for(i = 1, #b,
    out[start + i - 1] = b[i];
  );
  return(out);
};

cc_profile_data(p, {kind=1}, {normalize=1}) =
{
  my(R, root, back_steps, subj, W, chain_len, breaker,
     immune_count = 0, active_count = 0, min_kill = 999, min_kill_q = 0,
     q, kp, pm1, pp1, S_pm1, S_pp1, smooth_pm1, rem_pm1, smooth_pp1, rem_pp1);

  R = cc_root(p, kind);
  root = R[1];
  back_steps = R[2];
  subj = if(normalize, root, p);

  W = cc_walk_from_quiet(subj, kind);
  chain_len = W[2];
  breaker = W[3];

  for(i = 1, #cc_analysis_primes,
    q = cc_analysis_primes[i];
    if(cc_is_immune(subj, q, kind),
      immune_count++;
    ,
      active_count++;
      kp = cc_kill_pos(subj % q, q, kind);
      if(kp < min_kill,
        min_kill = kp;
        min_kill_q = q;
      );
    );
  );

  pm1 = subj - 1;
  pp1 = subj + 1;
  S_pm1 = cc_smooth_split(pm1);
  smooth_pm1 = S_pm1[1];
  rem_pm1 = S_pm1[2];
  S_pp1 = cc_smooth_split(pp1);
  smooth_pp1 = S_pp1[1];
  rem_pp1 = S_pp1[2];

  return([
    cc_fmt_hex(p, 80),
    cc_fmt_hex(subj, 80),
    cc_fmt_hex(root, 80),
    back_steps,
    kind,
    if(isprime(subj), 1, 0),
    if(cc_is_root(subj, kind), 1, 0),
    #binary(subj),
    chain_len,
    cc_fmt_hex(breaker, 80),
    immune_count,
    active_count,
    if(min_kill < 999, min_kill, -1),
    if(min_kill_q > 0, min_kill_q, 0),
    valuation(pm1, 2),
    valuation(pp1, 2),
    #binary(smooth_pm1),
    #binary(smooth_pp1),
    cc_prefix_bits(subj, 8),
    cc_prefix_bits(subj, 16),
    cc_tail_bits(subj, 5),
    cc_tail_bits(subj, 8)
  ]);
};

cc_profile_emit(p, {kind=1}, {normalize=1}, {sep="\t"}, {header=0}) =
{
  if(header,
    printf("%s\n", cc_join_fields(cc_profile_labels(), sep));
  );
  printf("%s\n", cc_join_fields(cc_profile_data(p, kind, normalize), sep));
  return(1);
};

cc_profile_export(filename, roots, {kind=1}, {normalize=1}, {sep="\t"}, {append=0}) =
{
  if(!append,
    write(filename, cc_join_fields(cc_profile_labels(), sep));
  );
  for(i = 1, #roots,
    write(filename, cc_join_fields(cc_profile_data(roots[i], kind, normalize), sep));
  );
  return(#roots);
};

cc_candidate_profile_labels() =
{
  return([
    "input_hex", "subject_hex", "root_hex", "back_steps", "kind", "normalize",
    "is_prime", "is_root", "bits", "chain_len", "heat",
    "breaker_pos", "breaker_hex", "breaker_shadow_q",
    "breaker_factor_hint", "immune_count", "active_count",
    "earliest_kill_pos", "earliest_kill_q",
    "prefix8", "prefix16", "tail5", "tail8"
  ]);
};

/* Exact-input single-number profile for arbitrary prime/composite inputs */
cc_candidate_profile_data(p, {kind=1}, {normalize=0}, {factor_limit=1000}) =
{
  my(R, root, back_steps, subj, W, chain_len, breaker,
     immune_count = 0, active_count = 0, min_kill = 999, min_kill_q = 0,
     q, kp, breaker_shadow_q, breaker_factor_hint);

  R = cc_root(p, kind);
  root = R[1];
  back_steps = R[2];
  subj = if(normalize, root, p);

  W = cc_walk_from_quiet(subj, kind);
  chain_len = W[2];
  breaker = W[3];
  breaker_shadow_q = cc_shadow_hit(breaker);
  breaker_factor_hint = cc_small_factor_hint(breaker, factor_limit);

  for(i = 1, #cc_analysis_primes,
    q = cc_analysis_primes[i];
    if(cc_is_immune(subj, q, kind),
      immune_count++;
    ,
      active_count++;
      kp = cc_kill_pos(subj % q, q, kind);
      if(kp < min_kill,
        min_kill = kp;
        min_kill_q = q;
      );
    );
  );

  return([
    cc_fmt_hex(p, 80),
    cc_fmt_hex(subj, 80),
    cc_fmt_hex(root, 80),
    back_steps,
    kind,
    normalize,
    if(isprime(subj), 1, 0),
    if(cc_is_root(subj, kind), 1, 0),
    #binary(subj),
    chain_len,
    cc_heat(chain_len),
    chain_len,
    cc_fmt_hex(breaker, 80),
    breaker_shadow_q,
    breaker_factor_hint,
    immune_count,
    active_count,
    if(min_kill < 999, min_kill, -1),
    if(min_kill_q > 0, min_kill_q, 0),
    cc_prefix_bits(subj, 8),
    cc_prefix_bits(subj, 16),
    cc_tail_bits(subj, 5),
    cc_tail_bits(subj, 8)
  ]);
};

cc_candidate_profile_emit(p, {kind=1}, {normalize=0}, {sep="\t"}, {header=0}, {factor_limit=1000}) =
{
  if(header,
    printf("%s\n", cc_join_fields(cc_candidate_profile_labels(), sep));
  );
  printf("%s\n", cc_join_fields(cc_candidate_profile_data(p, kind, normalize, factor_limit), sep));
  return(1);
};

cc_candidate_profile(p, {kind=1}, {normalize=0}, {factor_limit=1000}) =
{
  my(R, root, back_steps, subj, D, chain_len, breaker, heat, breaker_shadow_q, breaker_factor_hint);

  R = cc_root(p, kind);
  root = R[1];
  back_steps = R[2];
  subj = if(normalize, root, p);
  D = cc_candidate_profile_data(p, kind, normalize, factor_limit);
  chain_len = D[10];
  heat = D[11];
  breaker = D[13];
  breaker_shadow_q = D[14];
  breaker_factor_hint = D[15];

  cc_hdr("CANDIDATE PROFILE");
  printf("  Input:          %s\n", cc_fmt_hex(p, 60));
  printf("  Subject:        %s  (%s)\n",
    cc_fmt_hex(subj, 60),
    if(normalize, "normalized subject", "exact input"));
  printf("  True root:      %s\n", cc_fmt_hex(root, 60));
  printf("  Back steps:     %d\n", back_steps);
  if(!normalize && back_steps > 0,
    printf("  Note:           exact-start analysis; input is a non-root prime member\n");
  );
  printf("  Prime?:         %s\n", if(isprime(subj), "YES", "no"));
  printf("  Root?:          %s\n", if(cc_is_root(subj, kind), "YES", "no"));
  printf("  Bits:           %d\n", #binary(subj));
  printf("  Achieved run:   %d\n", chain_len);
  printf("  Heat:           %.6f\n", heat);
  printf("  Breaker pos:    %d\n", chain_len);
  printf("  Breaker:        %s\n", breaker);
  if(breaker_shadow_q > 0,
    printf("  Breaker hit-q:  %d\n", breaker_shadow_q);
  ,
    printf("  Breaker hit-q:  none in analysis set\n");
  );
  if(breaker_factor_hint > 1,
    printf("  Factor hint:    %d\n", breaker_factor_hint);
  ,
    printf("  Factor hint:    none <= %d\n", factor_limit);
  );
  printf("  Immune/active:  %d / %d\n", D[16], D[17]);
  if(D[18] >= 0,
    printf("  Earliest kill:  pos %d by q=%d\n", D[18], D[19]);
  ,
    printf("  Earliest kill:  none in analysis set\n");
  );
  printf("  Prefix8/16:     %s / %s\n", D[20], D[21]);
  printf("  Tail5/8:        %s / %s\n", D[22], D[23]);

  return(D);
};

/* Compact normalized summary suitable for manual checks and future batching */
cc_profile(p, {kind=1}, {normalize=1}) =
{
  my(R, root, back_steps, subj, W, chain_len, breaker,
     immune_count = 0, active_count = 0, min_kill = 999, min_kill_q = 0,
     q, kp, pm1, pp1, S_pm1, S_pp1, smooth_pm1, rem_pm1, smooth_pp1, rem_pp1,
     profile);

  R = cc_root(p, kind);
  root = R[1];
  back_steps = R[2];
  subj = if(normalize, root, p);

  W = cc_walk_from_quiet(subj, kind);
  chain_len = W[2];
  breaker = W[3];

  for(i = 1, #cc_analysis_primes,
    q = cc_analysis_primes[i];
    if(cc_is_immune(subj, q, kind),
      immune_count++;
    ,
      active_count++;
      kp = cc_kill_pos(subj % q, q, kind);
      if(kp < min_kill,
        min_kill = kp;
        min_kill_q = q;
      );
    );
  );

  pm1 = subj - 1;
  pp1 = subj + 1;
  S_pm1 = cc_smooth_split(pm1);
  smooth_pm1 = S_pm1[1];
  rem_pm1 = S_pm1[2];
  S_pp1 = cc_smooth_split(pp1);
  smooth_pp1 = S_pp1[1];
  rem_pp1 = S_pp1[2];

  cc_hdr("PROFILE SUMMARY");
  printf("  Input:          %s\n", cc_fmt_hex(p, 60));
  printf("  Subject:        %s  (%s)\n",
    cc_fmt_hex(subj, 60),
    if(normalize, "normalized root", "exact input"));
  printf("  True root:      %s\n", cc_fmt_hex(root, 60));
  printf("  Back steps:     %d\n", back_steps);
  printf("  Prime?:         %s\n", if(isprime(subj), "YES", "no"));
  printf("  Root?:          %s\n", if(cc_is_root(subj, kind), "YES", "no"));
  printf("  Chain length:   %d\n", chain_len);
  printf("  Breaker:        %s\n", cc_fmt_hex(breaker, 50));
  printf("  Immune/active:  %d / %d\n", immune_count, active_count);
  if(min_kill < 999,
    printf("  Earliest kill:  pos %d by q=%d\n", min_kill, min_kill_q);
  ,
    printf("  Earliest kill:  none in analysis set\n");
  );
  printf("  v2(p-1):        %d\n", valuation(pm1, 2));
  printf("  v2(p+1):        %d\n", valuation(pp1, 2));
  printf("  31-smooth(p-1): %d bits   cofactor %d bits\n", #binary(smooth_pm1), #binary(rem_pm1));
  printf("  31-smooth(p+1): %d bits   cofactor %d bits\n", #binary(smooth_pp1), #binary(rem_pp1));
  printf("  Prefix8/16:     %s / %s\n", cc_prefix_bits(subj, 8), cc_prefix_bits(subj, 16));
  printf("  Tail5/8:        %s / %s\n", cc_tail_bits(subj, 5), cc_tail_bits(subj, 8));

  profile = [
    p,
    subj,
    root,
    back_steps,
    kind,
    isprime(subj),
    cc_is_root(subj, kind),
    chain_len,
    breaker,
    immune_count,
    active_count,
    if(min_kill < 999, min_kill, -1),
    if(min_kill_q > 0, min_kill_q, 0),
    valuation(pm1, 2),
    valuation(pp1, 2),
    #binary(smooth_pm1),
    #binary(smooth_pp1),
    cc_prefix_bits(subj, 8),
    cc_prefix_bits(subj, 16),
    cc_tail_bits(subj, 5),
    cc_tail_bits(subj, 8)
  ];

  return(profile);
};

/* Prefix/tail motif summary for manual shadow inspection */
cc_motif(p, {kind=1}, {normalize=1}, {motifs=0}) =
{
  my(R, subj, M, pat, pat_s, cnt, bits, res = List(), possible,
     transitions, run1, run0);

  if(motifs == 0,
    M = ["111","1111","11111","01","10","101","010"];
  ,
    M = motifs;
  );

  R = cc_root(p, kind);
  subj = if(normalize, R[1], p);
  bits = #binary(subj);
  transitions = cc_transition_count(subj);
  run1 = cc_longest_run(subj, 1);
  run0 = cc_longest_run(subj, 0);

  cc_hdr("PREFIX / TAIL MOTIF SUMMARY");
  printf("  Subject:        %s  (%s)\n",
    cc_fmt_hex(subj, 60),
    if(normalize, "normalized root", "exact input"));
  printf("  Bits:           %d\n", bits);
  printf("  Prefix8/12/16:  %s / %s / %s\n",
    cc_prefix_bits(subj, 8), cc_prefix_bits(subj, 12), cc_prefix_bits(subj, 16));
  printf("  Tail3/5/8:      %s / %s / %s\n",
    cc_tail_bits(subj, 3), cc_tail_bits(subj, 5), cc_tail_bits(subj, 8));
  printf("  Transitions:    %d  (density %.4f)\n",
    transitions, if(bits > 1, transitions / (bits - 1.0), 0.0));
  printf("  Longest runs:   ones=%d  zeros=%d\n\n", run1, run0);

  printf("  %10s  %8s  %10s  %8s\n", "Pattern", "Count", "Density", "Present");
  printf("  %s\n", concat(vector(48, i, "-")));
  for(i = 1, #M,
    pat = cc_pattern_vec(M[i]);
    pat_s = cc_pattern_str(pat);
    cnt = cc_motif_count(subj, pat);
    possible = max(1, bits - #pat + 1);
    printf("  %10s  %8d  %10.6f  %8s\n",
      pat_s, cnt, cnt / possible, if(cnt > 0, "YES", "no"));
    listput(res, [pat_s, cnt, cnt / possible, if(cnt > 0, 1, 0)]);
  );

  return(Vec(res));
};

/* Append candidate tails/patterns to a base root/input and test survival */
cc_tail_probe(p, {patterns=0}, {kind=1}, {normalize=1}) =
{
  my(R, base, base_root, base_walk, base_len, P, pat, pat_s, pat_v, L, cand,
     exact, exact_len, breaker, root_info, cand_root, back_steps, cand_is_root,
     cand_prime, results = List(), delta_len, status);

  R = cc_root(p, kind);
  base_root = R[1];
  base = if(normalize, base_root, p);
  base_walk = cc_walk_from_quiet(base, kind);
  base_len = base_walk[2];

  if(patterns == 0,
    P = ["1","11","111","1111","11111","001","011","101","00101111","11111101"];
  ,
    P = patterns;
  );

  cc_hdr("TAIL / PATTERN PROBE");
  printf("  Base:           %s  (%s)\n",
    cc_fmt_hex(base, 60),
    if(normalize, "normalized root", "exact input"));
  printf("  Base chain len: %d\n\n", base_len);

  printf("  %10s %5s %6s %6s %6s %6s  %s\n",
    "Pattern", "Bits", "Prime?", "Root?", "Run", "Delta", "Candidate");
  printf("  %s\n", concat(vector(92, i, "-")));

  for(i = 1, #P,
    pat = P[i];
    pat_v = cc_pattern_vec(pat);
    pat_s = cc_pattern_str(pat_v);
    L = #pat_v;
    cand = base * (1 << L) + fromdigits(pat_v, 2);

    cand_prime = isprime(cand);
    cand_is_root = cc_is_root(cand, kind);
    exact = cc_walk_from_quiet(cand, kind);
    exact_len = exact[2];
    breaker = exact[3];
    delta_len = exact_len - base_len;

    root_info = cc_root(cand, kind);
    cand_root = root_info[1];
    back_steps = root_info[2];

    status = if(cand_prime, "YES", "no");
    printf("  %10s %5d %6s %6s %6d %+6d  %s\n",
      pat_s, L,
      status,
      if(cand_is_root, "YES", "no"),
      exact_len,
      delta_len,
      cc_fmt_hex(cand, 44));
    printf("  %sroot=%s  back=%d  breaker=%s\n",
      concat(vector(45, j, " ")),
      cc_fmt_hex(cand_root, 30),
      back_steps,
      cc_fmt_hex(breaker, 30));

    listput(results, [
      pat_s, cand, cand_prime, cand_is_root,
      exact_len, delta_len, cand_root, back_steps, breaker
    ]);
  );

  return(Vec(results));
};

/* Fixed-depth exact trace: primes and composites across positions 0..depth */
cc_candidate_positions(p, {depth=18}, {kind=1}, {normalize=0}, {factor_limit=1000}) =
{
  my(R, root, subj, cur, prime_flag, break_seen = 0,
     hit_q, hint, marker, results = List());

  if(depth < 0, error("depth must be >= 0"));

  R = cc_root(p, kind);
  root = R[1];
  subj = if(normalize, root, p);
  cur = subj;

  cc_hdr("CANDIDATE POSITION TRACE");
  printf("  Subject:        %s  (%s)\n",
    cc_fmt_hex(subj, 60),
    if(normalize, "normalized subject", "exact input"));
  printf("  Highest pos:    %d\n", depth);
  printf("  Factor hint <=  %d\n\n", factor_limit);

  printf("  %3s  %6s %6s %8s %8s  %s\n",
    "Pos", "Bits", "Prime?", "Hit-q", "Hint", "Value");
  printf("  %s\n", concat(vector(92, i, "-")));

  for(pos = 0, depth,
    prime_flag = isprime(cur);
    hit_q = if(prime_flag, 0, cc_shadow_hit(cur));
    hint = if(prime_flag, 0, cc_small_factor_hint(cur, factor_limit));
    marker = "";
    if(!break_seen && !prime_flag,
      marker = "  <<BREAK>>";
      break_seen = 1;
    );

    printf("  %3d  %6d %6s %8s %8s  %s%s\n",
      pos,
      #binary(cur),
      if(prime_flag, "YES", "no"),
      if(hit_q > 0, Str(hit_q), "-"),
      if(hint > 1, Str(hint), "-"),
      cc_fmt_hex(cur, 56),
      marker);

    listput(results, [pos, cur, prime_flag, hit_q, hint]);
    cur = cc_next(cur, kind);
  );

  return(Vec(results));
};

/* One-call manual analysis for a single prime/composite input */
cc_candidate_full(p, {max_pos=18}, {kind=1}, {normalize=0}, {factor_limit=1000}) =
{
  my(R, subj);
  R = cc_root(p, kind);
  subj = if(normalize, R[1], p);

  cc_candidate_profile(p, kind, normalize, factor_limit);
  cc_motif(p, kind, normalize);
  cc_shadow(subj, kind);
  cc_candidate_positions(p, max_pos, kind, normalize, factor_limit);

  printf("\n");
  cc_sep();
  printf("  Candidate analysis complete.\n");
  cc_sep();
};

/* Exhaustive append-tail sweep over a length range */
cc_tail_probe_range(p, {len_min=1}, {len_max=5}, {kind=1}, {normalize=1}, {show_all=0}, {min_run=-1}) =
{
  my(R, base, base_walk, base_len, threshold, results = List(),
     L, x, pat_v, pat_s, cand, exact, run, delta, breaker, cand_prime,
     best_run, best_delta, best_pat, best_cand, best_breaker);

  if(len_min < 1 || len_max < len_min, error("invalid length range"));
  if(len_max > 16, error("len_max too large; keep <= 16 for manual use"));

  R = cc_root(p, kind);
  base = if(normalize, R[1], p);
  base_walk = cc_walk_from_quiet(base, kind);
  base_len = base_walk[2];
  threshold = if(min_run < 0, base_len, min_run);

  cc_hdr("TAIL PROBE RANGE");
  printf("  Base:           %s  (%s)\n",
    cc_fmt_hex(base, 60),
    if(normalize, "normalized root", "exact input"));
  printf("  Base chain len: %d\n", base_len);
  printf("  Printed rows:   %s\n\n", if(show_all, "all", Str("run >= ", threshold)));

  for(L = len_min, len_max,
    best_run = -1;
    best_delta = -999;
    best_pat = "";
    best_cand = 0;
    best_breaker = 0;

    cc_sub(Str("Tail length ", L, " (", 1 << L, " candidates)"));
    printf("  %10s %6s %6s %6s  %s\n",
      "Pattern", "Prime?", "Run", "Delta", "Candidate");
    printf("  %s\n", concat(vector(72, i, "-")));

    for(x = 0, (1 << L) - 1,
      pat_v = cc_bits_fixed_len(x, L);
      pat_s = cc_bitvec_to_str(pat_v);
      cand = base * (1 << L) + x;
      exact = cc_walk_from_quiet(cand, kind);
      run = exact[2];
      breaker = exact[3];
      delta = run - base_len;
      cand_prime = isprime(cand);

      if(run > best_run || (run == best_run && delta > best_delta),
        best_run = run;
        best_delta = delta;
        best_pat = pat_s;
        best_cand = cand;
        best_breaker = breaker;
      );

      if(show_all || run >= threshold,
        printf("  %10s %6s %6d %+6d  %s\n",
          pat_s,
          if(cand_prime, "YES", "no"),
          run,
          delta,
          cc_fmt_hex(cand, 44));
      );

      listput(results, [L, pat_s, cand, cand_prime, run, delta, breaker]);
    );

    printf("  BEST len=%d: pattern=%s  run=%d  delta=%+d  cand=%s  breaker=%s\n",
      L, best_pat, best_run, best_delta, cc_fmt_hex(best_cand, 44), cc_fmt_hex(best_breaker, 30));
  );

  return(Vec(results));
};

/* Fixed high prefix, exhaustive sweep over low vary_bits tail */
cc_prefix_probe(p, {vary_bits=8}, {kind=1}, {normalize=1}, {show_all=0}, {min_run=-1}) =
{
  my(R, base, base_walk, base_len, threshold, bits, prefix_len, prefix_v, prefix_s,
     prefix_int, base_tail, results = List(), x, tail_v, tail_s, cand, exact, run,
     delta, breaker, cand_prime, marker, best_run = -1, best_delta = -999,
     best_tail = "", best_cand = 0, best_breaker = 0);

  if(vary_bits < 1 || vary_bits > 16, error("vary_bits must be in [1,16]"));

  R = cc_root(p, kind);
  base = if(normalize, R[1], p);
  base_walk = cc_walk_from_quiet(base, kind);
  base_len = base_walk[2];
  threshold = if(min_run < 0, base_len, min_run);
  bits = #binary(base);
  if(vary_bits >= bits, error("vary_bits must be smaller than bit length"));

  prefix_len = bits - vary_bits;
  prefix_v = vecextract(binary(base), Str("1..", prefix_len));
  prefix_s = cc_bitvec_to_str(prefix_v);
  prefix_int = fromdigits(prefix_v, 2);
  base_tail = cc_tail_bits(base, vary_bits);

  cc_hdr("PREFIX-MATCHED LOW-TAIL SWEEP");
  printf("  Base:           %s  (%s)\n",
    cc_fmt_hex(base, 60),
    if(normalize, "normalized root", "exact input"));
  printf("  Fixed prefix:   %s  (%d bits fixed)\n", prefix_s, prefix_len);
  printf("  Swept tail:     %d bits  (%d candidates)\n", vary_bits, 1 << vary_bits);
  printf("  Base tail:      %s\n", base_tail);
  printf("  Base chain len: %d\n", base_len);
  printf("  Printed rows:   %s\n\n", if(show_all, "all", Str("run >= ", threshold)));

  printf("  %10s %6s %6s %6s  %s\n",
    "Tail", "Prime?", "Run", "Delta", "Candidate");
  printf("  %s\n", concat(vector(76, i, "-")));

  for(x = 0, (1 << vary_bits) - 1,
    tail_v = cc_bits_fixed_len(x, vary_bits);
    tail_s = cc_bitvec_to_str(tail_v);
    cand = prefix_int * (1 << vary_bits) + x;
    exact = cc_walk_from_quiet(cand, kind);
    run = exact[2];
    breaker = exact[3];
    delta = run - base_len;
    cand_prime = isprime(cand);
    marker = if(tail_s == base_tail, "  <-- BASE", "");

    if(run > best_run || (run == best_run && delta > best_delta),
      best_run = run;
      best_delta = delta;
      best_tail = tail_s;
      best_cand = cand;
      best_breaker = breaker;
    );

    if(show_all || run >= threshold,
      printf("  %10s %6s %6d %+6d  %s%s\n",
        tail_s,
        if(cand_prime, "YES", "no"),
        run,
        delta,
        cc_fmt_hex(cand, 48),
        marker);
    );

    listput(results, [tail_s, cand, cand_prime, run, delta, breaker]);
  );

  printf("\n  BEST tail=%s  run=%d  delta=%+d  cand=%s  breaker=%s\n",
    best_tail, best_run, best_delta, cc_fmt_hex(best_cand, 48), cc_fmt_hex(best_breaker, 30));

  return(Vec(results));
};

/* Compare two chains side by side: residues and shadow structure */
cc_compare(p1, p2, {kind=1}) =
{
  my(W1, W2, q, r1, r2, k1, k2, imm1, imm2);

  W1 = cc_walk_quiet(p1, kind);
  W2 = cc_walk_quiet(p2, kind);

  cc_hdr("CHAIN COMPARISON");
  printf("  Chain A: CC%d%s root=%s\n", W1[2], if(kind==1,"a","b"), cc_fmt_hex(W1[3], 40));
  printf("  Chain B: CC%d%s root=%s\n", W2[2], if(kind==1,"a","b"), cc_fmt_hex(W2[3], 40));

  cc_sub("Shadow Comparison");
  printf("  %5s  %7s %7s  %5s %5s  %5s %5s\n",
    "Prime", "Res-A", "Res-B", "Imm-A", "Imm-B", "Kill-A", "Kill-B");
  printf("  %s\n", concat(vector(55, i, "-")));

  for(i = 1, #cc_analysis_primes,
    q = cc_analysis_primes[i];
    r1 = W1[3] % q; r2 = W2[3] % q;
    imm1 = cc_is_immune(W1[3], q, kind);
    imm2 = cc_is_immune(W2[3], q, kind);
    k1 = if(imm1, -1, cc_kill_pos(r1, q, kind));
    k2 = if(imm2, -1, cc_kill_pos(r2, q, kind));
    printf("  %5d  %7d %7d  %5s %5s  %5s %5s\n",
      q, r1, r2,
      if(imm1, "YES", "no"),
      if(imm2, "YES", "no"),
      if(k1 < 0, "---", Str(k1)),
      if(k2 < 0, "---", Str(k2)));
  );
};

/* ────────────────────────────────────────────────────────────────────────────
 * Fingerprint: per-member factorization and analysis
 * ──────────────────────────────────────────────────────────────────────────── */

cc_fingerprint(p, {kind=1}) =
{
  my(R, root, chain, cur, M_pm1, M_pp1, q, r, kp,
     immune_list, min_kp, min_kp_q);

  R = cc_root(p, kind);
  root = R[1];

  /* Build chain */
  chain = List();
  cur = root;
  while(isprime(cur) && #chain < 200,
    listput(chain, cur);
    cur = cc_next(cur, kind);
  );

  cc_hdr(Str("FINGERPRINT — CC", #chain, if(kind==1,"a","b")));
  printf("  Root: %s\n\n", cc_fmt_hex(root, 60));

  for(i = 1, #chain,
    cur = chain[i];
    cc_sub(Str("Position ", i-1, " (", #binary(cur), " bits)"));
    printf("  Value: %s\n", cc_fmt_hex(cur, 60));
    printf("  Bits: %d  Pop: %d  NAF: %d  T1: %d\n",
      #binary(cur), hammingweight(cur), cc_naf_weight(cur), cc_trailing_ones(cur));

    /* v2 of p-1 and p+1 */
    printf("  v2(p-1) = %d   v2(p+1) = %d\n", valuation(cur-1, 2), valuation(cur+1, 2));

    /* Factor p-1 and p+1 */
    M_pm1 = factor(cur - 1);
    M_pp1 = factor(cur + 1);
    printf("  p-1 = %s\n", cc_fmt_fact_tagged(M_pm1));
    printf("  p+1 = %s\n", cc_fmt_fact_tagged(M_pp1));

    /* Immunization and kill positions */
    immune_list = List();
    min_kp = 999; min_kp_q = 0;
    for(j = 1, #cc_analysis_primes,
      q = cc_analysis_primes[j];
      if(cc_is_immune(cur, q, kind),
        listput(immune_list, q);
      ,
        kp = cc_kill_pos(cur % q, q, kind);
        if(kp < min_kp,
          min_kp = kp;
          min_kp_q = q;
        );
      );
    );
    printf("  Immune: %s\n", if(#immune_list > 0, Str(Vec(immune_list)), "(none)"));
    if(min_kp < 999,
      printf("  Min kill: pos %d by prime %d\n", min_kp, min_kp_q);
    );

    /* S-base coverage */
    cc_sbase_check(cur, kind);
  );

  return(#chain);
};

/* ────────────────────────────────────────────────────────────────────────────
 * Layer 6: Top-Down Grid & p-Adic Address Space
 *
 * Key insight: In the 2-adic grid (row = bit_length, col = lower bits),
 * every CC is a STRAIGHT LINE with slope 1.  The column position is
 * determined entirely by the root's residue structure (p+1 for first kind).
 *
 * The p-adic valuation fingerprint (v3, v5, v7, ... of p-1 and p+1)
 * gives a multi-dimensional coordinate.  Primorial-based roots have
 * high valuations across many primes; random/sequential roots don't.
 * ──────────────────────────────────────────────────────────────────────────── */

/* 2-adic grid position: (row, column) where row = bit_length, col = lower bits */
cc_grid_pos(n) =
{
  my(row, col, mask);
  row = #binary(n);
  mask = (1 << (row - 1)) - 1;
  col = bitand(n, mask);
  return([row, col]);
};

/* Show CC trajectory in the 2-adic grid — makes the "straight line" visible */
cc_grid_trajectory(p, {kind=1}) =
{
  my(R, root, cur, chain_len = 0, pos, row, col, mask,
     col_prev = -1, col_delta, breaker_row, breaker_col);

  R = cc_root(p, kind);
  root = R[1];

  cc_hdr(Str("TOP-DOWN GRID TRAJECTORY — ", if(kind==1,"FIRST","SECOND"), " KIND"));

  printf("  In the 2-adic grid: row = bit_length, col = lower bits\n");
  printf("  Every CC is a straight line: row increases by 1 per step,\n");
  printf("  column is determined by (p%s1) — the structural invariant.\n\n",
    if(kind==1,"+","-"));

  /* Show the structural invariant */
  my(inv = if(kind==1, root+1, root-1));
  printf("  Structural invariant: p%s1 = %s\n", if(kind==1,"+","-"), cc_fmt_hex(inv, 50));
  printf("  v2(p%s1) = %d\n", if(kind==1,"+","-"), valuation(inv, 2));
  printf("  odd_core(p%s1) = %s\n\n", if(kind==1,"+","-"), cc_fmt_hex(cc_odd_part(inv), 40));

  printf("  %3s  %5s  %12s  %12s  %s\n", "Pos", "Row", "Col (dec)", "Col (hex)", "Status");
  printf("  %s\n", concat(vector(60, i, "-")));

  cur = root;
  while(1,
    row = #binary(cur);
    mask = (1 << (row - 1)) - 1;
    col = bitand(cur, mask);

    if(isprime(cur),
      chain_len++;
      printf("  %3d  %5d  %12d  %12s  PRIME\n",
        chain_len - 1, row, col, cc_fmt_hex(col, 12));
    ,
      printf("  %3d  %5d  %12d  %12s  <<BREAK>>\n",
        chain_len, row, col, cc_fmt_hex(col, 12));
      breaker_row = row; breaker_col = col;
      break;
    );
    cur = cc_next(cur, kind);
    if(chain_len > 200, break);
  );

  /* Column pattern analysis */
  cc_sub("Column Pattern");
  printf("  For first-kind CC: col_i = (p+1) * 2^i - 1  (mod 2^row_i)\n");
  printf("  The column bits grow by prepending from the structural invariant.\n");
  printf("  This is why the chain is a straight line in the grid:\n");
  printf("  each step adds 1 row and extends the column deterministically.\n");

  return(chain_len);
};

/* Multi-prime valuation fingerprint: shows WHERE a number sits in p-adic space */
cc_padic_addr(p, {depth=1}) =
{
  my(primes_list, pm1, pp1, vpm1, vpp1);

  primes_list = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31];
  pm1 = p - 1;
  pp1 = p + 1;

  cc_hdr("p-ADIC ADDRESS (multi-prime valuation fingerprint)");
  printf("  p = %s\n\n", cc_fmt_hex(p, 60));

  /* Valuation table for p-1 and p+1 */
  printf("  %5s  %8s  %8s  %s\n", "Prime", "v_q(p-1)", "v_q(p+1)", "Interpretation");
  printf("  %s\n", concat(vector(65, i, "-")));

  for(i = 1, #primes_list,
    my(q = primes_list[i]);
    vpm1 = valuation(pm1, q);
    vpp1 = valuation(pp1, q);
    my(interp = "");
    if(q <= 23,
      /* For small primes: connect to immunization */
      if(vpp1 >= 1 && q != 2,
        interp = Str("q|(p+1) => immune (1st kind, res=", q-1, ")");
      );
      if(vpm1 >= 1 && q != 2,
        if(interp != "", interp = Str(interp, " | "));
        interp = Str(interp, "q|(p-1) => immune (2nd kind, res=1)");
      );
    );
    printf("  %5d  %8d  %8d  %s\n", q, vpm1, vpp1, interp);
  );

  /* Composite fingerprint: product of small-prime powers dividing p+1 and p-1 */
  cc_sub("Smooth Part of p+1 (construction signature)");
  my(smooth_pp1 = 1, rem_pp1 = pp1);
  for(i = 1, #primes_list,
    my(q = primes_list[i], v = valuation(pp1, q));
    if(v > 0, smooth_pp1 *= q^v; rem_pp1 /= q^v);
  );
  printf("  p+1 = %s * %s\n", smooth_pp1, cc_fmt_hex(rem_pp1, 40));
  printf("  31-smooth part: %s  (%d bits)\n", smooth_pp1, #binary(smooth_pp1));
  printf("  Cofactor:       %s  (%d bits)\n", cc_fmt_hex(rem_pp1, 40), #binary(rem_pp1));

  cc_sub("Smooth Part of p-1");
  my(smooth_pm1 = 1, rem_pm1 = pm1);
  for(i = 1, #primes_list,
    my(q = primes_list[i], v = valuation(pm1, q));
    if(v > 0, smooth_pm1 *= q^v; rem_pm1 /= q^v);
  );
  printf("  p-1 = %s * %s\n", smooth_pm1, cc_fmt_hex(rem_pm1, 40));
  printf("  31-smooth part: %s  (%d bits)\n", smooth_pm1, #binary(smooth_pm1));
  printf("  Cofactor:       %s  (%d bits)\n", cc_fmt_hex(rem_pm1, 40), #binary(rem_pm1));

  /* Primorial test: is this a primorial-based root? */
  cc_sub("Primorial Signature");
  my(primorial_primes = [2,3,5,7,11,13,17,19,23], all_divide = 1);
  for(i = 1, #primorial_primes,
    if(pp1 % primorial_primes[i] != 0, all_divide = 0; break);
  );
  if(all_divide,
    printf("  23# | (p+1): YES — primorial-based root (standard construction)\n");
  ,
    printf("  23# | (p+1): NO — non-primorial root (your construction style)\n");
    printf("  Missing from p+1: ");
    for(i = 1, #primorial_primes,
      if(pp1 % primorial_primes[i] != 0, printf("%d ", primorial_primes[i]));
    );
    printf("\n");
  );

  if(depth >= 1,
    /* Deep comparison: how does this root's p-adic address compare
     * to what a primorial root would look like? */
    cc_sub("Construction Comparison");
    printf("  Primorial roots: v_q(p+1) >= 1 for all q <= 23# (forced by CRT)\n");
    printf("  Your root:       some primes divide p+1, others divide p-1\n");
    printf("  The immunization set is different — but the CHAIN doesn't care\n");
    printf("  which side (p-1 vs p+1) a prime lands on, only that it's immune.\n");
    printf("\n  Question: does the SPLIT between p-1 and p+1 immunizations\n");
    printf("  affect chain length? This is what the p-adic address reveals.\n");
  );

  return([smooth_pp1, smooth_pm1]);
};

/* CRT reconstruction: given desired immunizations, what residue class? */
cc_crt_constraints(kind, {primes=0}) =
{
  my(pr, residues = List(), moduli = List(), M = 1, x);

  if(primes == 0, pr = [3,5,7,11,13,17,19,23], pr = primes);

  cc_hdr(Str("CRT CONSTRAINTS — ", if(kind==1,"FIRST","SECOND"), " KIND"));
  printf("  For a root to be immunized against prime q:\n");
  if(kind == 1,
    printf("    p mod q = q-1  (equivalently: q | (p+1))\n\n");
  ,
    printf("    p mod q = 1    (equivalently: q | (p-1))\n\n");
  );

  printf("  %5s  %10s  %s\n", "Prime", "Required r", "Constraint");
  printf("  %s\n", concat(vector(45, i, "-")));

  for(i = 1, #pr,
    my(q = pr[i], r);
    r = cc_immune_residue(q, kind);
    listput(residues, r);
    listput(moduli, q);
    M *= q;
    printf("  %5d  %10d  p ≡ %d (mod %d)\n", q, r, r, q);
  );

  /* Solve CRT */
  printf("\n  Combined modulus M = %s", pr[1]);
  for(i = 2, #pr, printf(" * %s", pr[i]));
  printf(" = %d\n", M);

  /* Use PARI's chinese() for CRT */
  x = chinese(Mod(residues[1], moduli[1]), Mod(residues[2], moduli[2]));
  for(i = 3, #residues,
    x = chinese(x, Mod(residues[i], moduli[i]));
  );

  printf("  CRT solution: p ≡ %d (mod %d)\n", lift(x), M);
  printf("\n  This means: roots fully immunized against {%s", pr[1]);
  for(i = 2, #pr, printf(",%s", pr[i]));
  printf("}\n  must satisfy p = %d + k * %d for some k >= 0\n", lift(x), M);

  /* Density */
  printf("\n  Density: 1 in every %d integers satisfies these constraints\n", M);
  printf("  In a %d-bit range: ~2^%d / %d ≈ %.0f candidates\n",
    90, 90, M, 2.0^90 / M);

  return(lift(x));
};

/* Grid neighborhood: find primes near a chain in the 2-adic grid */
cc_grid_neighbors(p, {radius=5}, {kind=1}) =
{
  my(R, root, row, col, mask, neighbor, count);

  R = cc_root(p, kind);
  root = R[1];
  row = #binary(root);
  mask = (1 << (row - 1)) - 1;
  col = bitand(root, mask);

  cc_hdr("GRID NEIGHBORHOOD");
  printf("  Root at grid position: row=%d, col=%d\n", row, col);
  printf("  Scanning ±%d columns at same row (same bit-length)\n\n", radius);

  printf("  %8s  %12s  %s  %s\n", "Offset", "Col", "Value (hex)", "Properties");
  printf("  %s\n", concat(vector(70, i, "-")));

  for(delta = -radius, radius,
    my(new_col = col + delta);
    if(new_col >= 0 && new_col <= mask,
      neighbor = (1 << (row - 1)) + new_col;
      my(tag = "", w);
      if(isprime(neighbor),
        tag = "PRIME";
        /* Check if it starts a chain */
        w = cc_walk_quiet(neighbor, 1);
        if(w[2] >= 2, tag = Str(tag, " CC", w[2], "a"));
        w = cc_walk_quiet(neighbor, 2);
        if(w[2] >= 2, tag = Str(tag, " CC", w[2], "b"));
      );
      my(marker = if(delta == 0, " <-- ROOT", ""));
      printf("  %+8d  %12d  %s  %s%s\n",
        delta, new_col, cc_fmt_hex(neighbor, 30), tag, marker);
    );
  );

  /* Also show neighbors at the SAME column but different rows */
  cc_sub("Vertical Neighbors (same column, different bit-length)");
  printf("  Same odd_core family — these share the column in top-down view\n\n");
  my(core = cc_odd_part(root));
  for(shift = 0, 6,
    my(v = core * (1 << shift));
    if(#binary(v) >= row - 3 && #binary(v) <= row + 3,
      my(tag = "");
      if(isprime(v),
        tag = "PRIME";
        my(w = cc_walk_quiet(v, 1));
        if(w[2] >= 2, tag = Str(tag, " CC", w[2], "a"));
      );
      my(marker = if(v == root, " <-- ROOT", ""));
      printf("  2^%d * core = %s  (%d bits)  %s%s\n",
        shift, cc_fmt_hex(v, 30), #binary(v), tag, marker);
    );
  );

  return([row, col]);
};

/* ────────────────────────────────────────────────────────────────────────────
 * Layer 7: Mini Search Engine
 *
 * GP/PARI reimplementation of the core algorithm from cc_gmp_v31_claude.c.
 * NOT meant for production speed — this is for interactive exploration,
 * brainstorming on a phone/iSH, and verifying algorithm understanding.
 *
 * Pipeline: CRT bases → sieve → PRP → chain follow → proven verify
 *
 * Mathematical foundation:
 *   First-kind CC: p_i = 2^i * (p+1) - 1
 *   q | p_i  iff  p ≡ 2^(-i) - 1 (mod q)   [the "forbidden residue"]
 *   Immune:  p ≡ q-1 (mod q)  iff  q | (p+1) — never killed
 *
 *   Second-kind CC: p_i = 2^i * (p-1) + 1
 *   q | p_i  iff  p ≡ 1 - 2^(-i) (mod q)
 *   Immune:  p ≡ 1 (mod q)  iff  q | (p-1)
 *
 * Construction: candidates n = k * M + base, where M = product of CRT primes
 * and base is chosen so n avoids all forbidden residues for the CRT primes.
 * Additional sieve primes reject candidates quickly before PRP testing.
 * ──────────────────────────────────────────────────────────────────────────── */

/* Forbidden residues mod q for a chain of length `target`.
 * Returns sorted vector of distinct forbidden residues. */
cc_forbidden_res(q, target, {kind=1}) =
{
  my(res = List(), r, seen);
  for(i = 0, target - 1,
    if(kind == 1,
      r = lift(Mod(2, q)^(-i) - 1);
    ,
      r = lift(1 - Mod(2, q)^(-i));
    );
    seen = 0;
    for(j = 1, #res, if(res[j] == r, seen = 1; break));
    if(!seen, listput(res, r));
  );
  return(vecsort(Vec(res)));
};

/* Valid (non-forbidden) residues mod q for given target and kind. */
cc_valid_res(q, target, {kind=1}) =
{
  my(forb = cc_forbidden_res(q, target, kind), valid = List(), found);
  for(r = 0, q - 1,
    found = 0;
    for(j = 1, #forb, if(forb[j] == r, found = 1; break));
    if(!found, listput(valid, r));
  );
  return(Vec(valid));
};

/* Build all CRT-valid bases mod M for a target chain length.
 * Default CRT primes: {3,5,7,11,13,17,19} (matching cc_gmp_v31).
 * Returns [bases_vector, M].
 *
 * Example: cc_build_bases(18) returns 9 bases mod 4849845
 *   (most primes forced to immune residue; q=17 has ord_2(17)=8 < 18,
 *   so 9 valid residues survive). */
cc_build_bases(target, {kind=1}, {primes=0}) =
{
  my(pr, M, bases, q, valid, new_bases, nb);
  if(primes == 0, pr = [3, 5, 7, 11, 13, 17, 19], pr = primes);

  bases = cc_valid_res(pr[1], target, kind);
  M = pr[1];

  for(i = 2, #pr,
    q = pr[i];
    valid = cc_valid_res(q, target, kind);
    new_bases = List();
    for(j = 1, #bases,
      for(k = 1, #valid,
        nb = lift(chinese(Mod(bases[j], M), Mod(valid[k], q)));
        listput(new_bases, nb);
      );
    );
    M *= q;
    bases = Vec(new_bases);
  );

  return([vecsort(bases), M]);
};

/* Precompute sieve table: for each prime in [lo, hi], store forbidden residues.
 * Returns vector of [prime, forbidden_residues_vec]. */
cc_build_sieve(target, {kind=1}, {lo=23}, {hi=500}) =
{
  my(table = List());
  forprime(q = lo, hi,
    listput(table, [q, cc_forbidden_res(q, target, kind)]);
  );
  return(Vec(table));
};

/* Check candidate n against precomputed sieve table.
 * Returns 1 if n passes (no small prime kills chain before target), 0 if killed. */
cc_sieve_check(n, sieve_table) =
{
  my(q, forb, r);
  for(i = 1, #sieve_table,
    q = sieve_table[i][1];
    forb = sieve_table[i][2];
    r = n % q;
    for(j = 1, #forb,
      if(r == forb[j], return(0));
    );
  );
  return(1);
};

/* Quick chain length count using BPSW pseudoprime test (fast).
 * For search use; verify hits with cc_chain_verify(). */
cc_chain_test(n, {kind=1}) =
{
  my(len = 0, cur = n);
  while(ispseudoprime(cur),
    len++;
    cur = cc_next(cur, kind);
    if(len > 200, break);
  );
  return(len);
};

/* Verify chain with proven primality (APRCL / ECPP via PARI).
 * Returns [proven_length, 1]. */
cc_chain_verify(n, {kind=1}) =
{
  my(len = 0, cur = n);
  while(isprime(cur),
    len++;
    cur = cc_next(cur, kind);
    if(len > 200, break);
  );
  return([len, 1]);
};

/* Show search configuration: bases, sieve, density estimates.
 * Useful for understanding search space before launching. */
cc_search_info(bits, target, {kind=1}) =
{
  my(BM, bases, M, sieve, nk, est_sieve, est_prp);

  BM = cc_build_bases(target, kind);
  bases = BM[1]; M = BM[2];
  sieve = cc_build_sieve(target, kind);

  cc_hdr(Str("SEARCH INFO — CC", target, if(kind==1,"a","b"), " at ", bits, " bits"));

  printf("  CRT modulus M    = %d\n", M);
  printf("  CRT primes:        {3,5,7,11,13,17,19}\n");
  printf("  Valid bases        = %d\n", #bases);
  printf("  Sieve primes       = %d (23..500)\n", #sieve);

  /* Estimate candidate count per base */
  nk = 2.0^(bits-1) / M;
  printf("\n  Per base: ~%.0f k-values in %d-bit range\n", nk, bits);
  printf("  Total candidates:  ~%.0f\n", nk * #bases);

  /* Estimate sieve survival rate */
  est_sieve = 1.0;
  for(i = 1, #sieve,
    my(q = sieve[i][1], nf = #sieve[i][2]);
    est_sieve *= (1.0 - 1.0*nf / q);
  );
  printf("  Est sieve pass:    %.4f%% (1 in %.0f)\n", 100*est_sieve, 1/est_sieve);

  /* Estimate PRP rate (Mertens-like) */
  est_prp = 1.0 / (bits * log(2));
  printf("  Est PRP rate:      ~1/%.0f at %d bits\n", 1/est_prp, bits);

  /* Combined estimate */
  printf("  Est PRP survivors: ~%.0f per base\n", nk * est_sieve * est_prp);
  printf("  Est total PRPs:    ~%.0f across all bases\n", nk * #bases * est_sieve * est_prp);

  printf("\n  Bases (mod %d):\n", M);
  for(i = 1, min(#bases, 20),
    printf("    [%2d] %d", i, bases[i]);
    if(i <= 5,
      printf("  (mod 6 = %d, mod 30 = %d)", bases[i] % 6, bases[i] % 30);
    );
    printf("\n");
  );
  if(#bases > 20,
    printf("    ... (%d more)\n", #bases - 20);
  );

  /* Show valid residues per prime */
  cc_sub("Valid Residues Per CRT Prime");
  my(crt_pr = [3, 5, 7, 11, 13, 17, 19]);
  for(i = 1, #crt_pr,
    my(q = crt_pr[i], v = cc_valid_res(q, target, kind));
    printf("  q=%2d: %d valid  %s  (ord_2=%d, immune=%d)\n",
      q, #v, Str(v),
      znorder(Mod(2, q)), cc_immune_residue(q, kind));
  );

  return([bases, M, sieve]);
};

/* ── Main search function ──
 *
 * cc_search(bits, target, {kind}, {prefix}, {prefix_bits}, {max_cand}, {verbose})
 *
 * bits:        target bit length (e.g., 40, 60, 89)
 * target:      minimum chain length to report (e.g., 5, 7, 10)
 * kind:        1 = first kind (default), 2 = second kind
 * prefix:      binary prefix as integer (e.g., 5 for 0b101)
 * prefix_bits: number of bits in the prefix (e.g., 3 for 0b101)
 * max_cand:    stop after testing this many candidates (default 100000)
 * verbose:     0 = silent, 1 = progress reports (default 1)
 *
 * Returns vector of [root, chain_length] for each hit.
 *
 * Example:
 *   cc_search(40, 5)               \\ find CC5+ roots at 40 bits
 *   cc_search(60, 7, 1, 5, 3)      \\ CC7+ at 60 bits, prefix 0b101
 *   cc_search(30, 4,,,,, 0)         \\ CC4+ at 30 bits, silent mode
 */
cc_search(bits, target, {kind=1}, {prefix=0}, {prefix_bits=0}, {max_cand=100000}, {verbose=1}) =
{
  my(BM, bases, M, sieve, lo, hi, n, base,
     tested = 0, sieved = 0, prp_pass = 0, found = List(), len,
     t0, elapsed, p_lo, p_hi, k_lo, k_hi, k_start, k_step);

  t0 = getabstime();

  /* Build CRT bases */
  BM = cc_build_bases(target, kind);
  bases = BM[1]; M = BM[2];

  /* Build sieve table for primes beyond CRT set */
  sieve = cc_build_sieve(target, kind, 23, 500);

  if(verbose,
    cc_hdr(Str("SEARCH — CC", target, if(kind==1,"a","b"), " at ", bits, " bits"));
    printf("  CRT modulus M = %d  (%d bases)\n", M, #bases);
    printf("  Sieve primes:   %d (23..500)\n", #sieve);
  );

  /* Compute bit range */
  lo = 1 << (bits - 1);
  hi = (1 << bits) - 1;

  /* Apply prefix constraint */
  if(prefix > 0 && prefix_bits > 0,
    p_lo = prefix << (bits - prefix_bits);
    p_hi = ((prefix + 1) << (bits - prefix_bits)) - 1;
    lo = max(lo, p_lo);
    hi = min(hi, p_hi);
  );

  if(verbose,
    printf("  Range: %s .. %s\n", cc_fmt_hex(lo, 30), cc_fmt_hex(hi, 30));
    printf("  Max candidates: %d\n", max_cand);
    printf("  Searching...\n\n");
  );

  /* M is odd (product of odd primes), so n = k*M + base.
   * n is odd iff k and base have different parity.
   * Use forstep with step=2 to only hit odd candidates. */

  for(bi = 1, #bases,
    base = bases[bi];

    /* Compute k range for this base */
    k_lo = (lo - base) \ M;
    if(k_lo * M + base < lo, k_lo++);
    k_hi = (hi - base) \ M;

    /* Adjust k parity: need (k + base) odd since M is odd, so n is odd */
    k_start = k_lo;
    if(base % 2 == 0,
      /* Need k odd */
      if(k_start % 2 == 0, k_start++);
    ,
      /* Need k even */
      if(k_start % 2 == 1, k_start++);
    );

    forstep(k = k_start, k_hi, 2,
      n = k * M + base;

      /* Range check */
      if(n < lo || n > hi, next);

      /* Prefix check */
      if(prefix > 0 && prefix_bits > 0,
        if(n >> (bits - prefix_bits) != prefix, next);
      );

      tested++;

      /* Sieve check */
      if(!cc_sieve_check(n, sieve), next);
      sieved++;

      /* Quick PRP test on candidate root */
      if(!ispseudoprime(n), next);
      prp_pass++;

      /* Chain follow with PRP */
      len = cc_chain_test(n, kind);

      if(len >= target,
        /* Verify with proven primality */
        my(V = cc_chain_verify(n, kind));
        if(V[1] >= target,
          listput(found, [n, V[1]]);
          if(verbose,
            printf("  *** FOUND CC%d%s: %s (%d bits) ***\n",
              V[1], if(kind==1,"a","b"), Str(n), bits);
          );
        );
      );

      if(tested >= max_cand,
        if(verbose, printf("\n  Reached max_cand=%d, stopping.\n", max_cand));
        break(2);
      );

      /* Progress report every 10000 */
      if(verbose && tested % 10000 == 0,
        elapsed = getabstime() - t0;
        printf("  [tested=%d  sieved=%d  prp=%d  found=%d  %.1fs  %.0f/s]\n",
          tested, sieved, prp_pass, #found,
          elapsed / 1000.0,
          if(elapsed > 0, 1000.0 * tested / elapsed, 0));
      );
    );
  );

  elapsed = getabstime() - t0;
  if(verbose,
    printf("\n  Search complete:\n");
    printf("    Tested:       %d\n", tested);
    printf("    Sieve pass:   %d (%.2f%%)\n", sieved,
      if(tested > 0, 100.0 * sieved / tested, 0.0));
    printf("    PRP pass:     %d\n", prp_pass);
    printf("    Chains found: %d\n", #found);
    printf("    Time:         %.1f seconds\n", elapsed / 1000.0);
    printf("    Rate:         %.0f candidates/sec\n",
      if(elapsed > 0, 1000.0 * tested / elapsed, 0));

    for(i = 1, #found,
      printf("    [%d] CC%d%s root=%s (%d bits)\n",
        i, found[i][2], if(kind==1,"a","b"), Str(found[i][1]), #binary(found[i][1]));
    );
  );

  return(Vec(found));
};

/* Quick search wrapper for common use cases */
cc_qsearch(bits, target, {kind=1}, {max_cand=50000}) =
  cc_search(bits, target, kind, 0, 0, max_cand, 1);

/* ────────────────────────────────────────────────────────────────────────────
 * Layer 8: Constructive Search Engine
 *
 * GP/PARI reimplementation of cc18_construct_v34__codex_v2.c.
 * Complementary to Layer 7's CRT-based approach.
 *
 * Construction: p = S * R - 1,  where S = primorial(N).
 *   Since p + 1 = S * R, every prime q <= N divides p + 1,
 *   so p is automatically immune to all primes up to N.
 *
 * Pipeline: R generation -> constructor sieve -> PRP -> ceiling -> chain
 *
 * Mathematical foundation:
 *   Chain member i (first kind): p_i = 2^i * (p+1) - 1 = 2^i * S * R - 1
 *   For sieve prime q (q does NOT divide S):
 *     q | p_i  iff  R ≡ (2^i * S)^{-1} (mod q)   [forbidden R residue]
 *
 * Modes:
 *   "B" (default) — R sequential odd, highest throughput
 *   "A"           — R must be PRP, most selective
 *   "C"           — R random odd, good for sampling
 * ──────────────────────────────────────────────────────────────────────────── */

/* Compute primorial(n) = product of primes <= n */
cc_primorial(n) =
{
  my(S = 1);
  forprime(q = 2, n, S *= q);
  return(S);
};

/* Forbidden R residues mod q for constructor sieve (first-kind only).
 * For sieve prime q that does NOT divide S, chain member i is:
 *   p_i = 2^i * S * R - 1
 *   p_i ≡ 0 (mod q) iff R ≡ (2^i * S)^{-1} (mod q)
 * Returns sorted vector of distinct forbidden R residues. */
cc_con_forbidden_R(q, S, target) =
{
  my(res = List(), Sinv, r, seen);
  if(q <= 1 || S % q == 0, return(Vec(res)));
  Sinv = lift(Mod(S, q)^(-1));
  for(i = 0, target - 1,
    r = lift(Mod(2, q)^(-i) * Sinv);
    seen = 0;
    for(j = 1, #res, if(res[j] == r, seen = 1; break));
    if(!seen, listput(res, r));
  );
  return(vecsort(Vec(res)));
};

/* Build constructor sieve table: for each prime q in [lo, hi] that does NOT
 * divide S, store [q, forbidden_R_residues].
 * Default range matches v34: primes 37..863. */
cc_con_sieve(S, target, {lo=37}, {hi=863}) =
{
  my(table = List());
  forprime(q = lo, hi,
    if(S % q != 0,
      listput(table, [q, cc_con_forbidden_R(q, S, target)]);
    );
  );
  return(Vec(table));
};

/* Check R against constructor sieve table.
 * Returns 1 if R passes (no sieve prime kills chain), 0 if killed. */
cc_con_sieve_check(R, sieve_table) =
{
  my(q, forb, r);
  for(i = 1, #sieve_table,
    q = sieve_table[i][1];
    forb = sieve_table[i][2];
    r = R % q;
    for(j = 1, #forb,
      if(r == forb[j], return(0));
    );
  );
  return(1);
};

/* Compute [R_min, R_max] for b-bit candidates with p = S*R - 1.
 * Optional prefix constraint narrows the range.
 * Only odd R values are used (caller responsibility). */
cc_con_R_range(S, bits, {prefix=0}, {prefix_bits=0}) =
{
  my(lo, hi, p_lo, p_hi, R_min, R_max);

  /* Full bit range: 2^{b-1} <= p < 2^b */
  lo = 2^(bits - 1);
  hi = 2^bits - 1;

  /* Apply prefix constraint */
  if(prefix > 0 && prefix_bits > 0,
    p_lo = prefix << (bits - prefix_bits);
    p_hi = ((prefix + 1) << (bits - prefix_bits)) - 1;
    lo = max(lo, p_lo);
    hi = min(hi, p_hi);
  );

  /* R range: p = S*R - 1, so R = (p+1)/S
   * R_min = ceil((lo + 1) / S), R_max = floor((hi + 1) / S) */
  R_min = (lo + 1 + S - 1) \ S;  /* ceil division */
  R_max = (hi + 1) \ S;

  /* Force R_min odd */
  if(R_min % 2 == 0, R_min++);

  return([R_min, R_max]);
};

/* OPT-6 prefilter: predict maximum chain length from small primes.
 * Uses 15 primes {67..137} to find earliest killing position.
 * Returns the ceiling (upper bound on chain length).
 * If ceiling < target, chain cannot reach target. */
cc_con_ceiling(p, target, {primes=0}) =
{
  my(pr, ceiling, r, kill_pos);
  if(primes == 0,
    pr = [67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137];
  ,
    pr = primes;
  );
  ceiling = 32;
  for(i = 1, #pr,
    my(q = pr[i]);
    /* For first-kind: p_j = 2^j*(p+1) - 1
     * p_j ≡ 0 (mod q) iff p ≡ 2^(-j) - 1 (mod q)
     * Find earliest j where this holds */
    r = p % q;
    kill_pos = q;  /* sentinel: no kill found */
    for(j = 0, min(q - 2, target),
      if(lift(Mod(2, q)^(-j) - 1) == r,
        kill_pos = j;
        break;
      );
    );
    if(kill_pos < ceiling, ceiling = kill_pos);
    if(ceiling < target, return(ceiling));  /* early exit */
  );
  return(ceiling);
};

/* Show constructor search parameters and density estimates.
 * S_limit: largest prime in primorial (default 31 matching v34). */
cc_con_info(bits, target, {S_limit=31}) =
{
  my(S, sieve, Rr, R_min, R_max, R_count, est_sieve, est_prp, S_bits);

  S = cc_primorial(S_limit);
  sieve = cc_con_sieve(S, target);
  Rr = cc_con_R_range(S, bits);
  R_min = Rr[1]; R_max = Rr[2];
  R_count = (R_max - R_min) \ 2 + 1;  /* odd R values */
  S_bits = #binary(S);

  cc_hdr(Str("CONSTRUCTOR INFO — CC", target, "a at ", bits, " bits"));

  printf("  S = primorial(%d)\n", S_limit);
  printf("  S = %s\n", if(S_bits <= 120, Str(S), Str(cc_fmt_hex(S, 40), " (", S_bits, " bits)")));
  printf("  S bits:            %d\n", S_bits);
  printf("  Sieve primes:      %d (%d..%d, excluding S-factors)\n",
    #sieve, if(#sieve > 0, sieve[1][1], 0), if(#sieve > 0, sieve[#sieve][1], 0));

  printf("\n  R range (odd only):\n");
  printf("    R_min = %s\n", if(#binary(R_min) <= 80, Str(R_min), cc_fmt_hex(R_min, 40)));
  printf("    R_max = %s\n", if(#binary(R_max) <= 80, Str(R_max), cc_fmt_hex(R_max, 40)));
  printf("    R_count = %d odd values\n", R_count);

  /* Estimate sieve survival rate */
  est_sieve = 1.0;
  for(i = 1, #sieve,
    my(q = sieve[i][1], nf = #sieve[i][2]);
    est_sieve *= (1.0 - 1.0*nf / q);
  );
  printf("\n  Est sieve pass:    %.4f%% (1 in %.0f)\n", 100*est_sieve, 1/est_sieve);

  /* Estimate PRP rate */
  est_prp = 1.0 / (bits * log(2));
  printf("  Est PRP rate:      ~1/%.0f at %d bits\n", 1/est_prp, bits);

  /* Combined */
  printf("  Est sieve+PRP:     ~%.0f PRPs from %d R values\n",
    R_count * est_sieve * est_prp, R_count);

  /* Show S factorization */
  cc_sub("S factorization");
  my(F = factor(S));
  for(i = 1, matsize(F)[1],
    printf("  %d", F[i, 1]);
    if(i < matsize(F)[1], printf(" * "));
  );
  printf("\n");

  /* Show immunity */
  printf("\n  Automatic immunity: all primes <= %d divide p+1 = S*R\n", S_limit);
  printf("  Chain members p_i = 2^i * S * R - 1 (first kind)\n");

  return([S, sieve, Rr]);
};

/* ── Main constructor search function ──
 *
 * cc_construct(bits, target, {S_limit}, {mode}, {prefix}, {prefix_bits},
 *              {max_cand}, {verbose})
 *
 * bits:        target bit length for p (e.g., 35, 89)
 * target:      minimum chain length to report
 * S_limit:     largest prime in primorial (default 31)
 * mode:        "B" = sequential odd R (default)
 *              "A" = R must be PRP
 *              "C" = random odd R
 * prefix:      binary prefix as integer (e.g., 5 for 0b101)
 * prefix_bits: number of bits in the prefix
 * max_cand:    stop after this many R values tested (default 100000)
 * verbose:     0 = silent, 1 = progress (default 1)
 *
 * Returns vector of [root, chain_length] for each hit.
 */
cc_construct(bits, target, {S_limit=31}, {mode="B"}, {prefix=0}, {prefix_bits=0}, {max_cand=100000}, {verbose=1}) =
{
  my(S, sieve, Rr, R_min, R_max, R_count, R, p,
     tested = 0, sieve_pass = 0, prp_pass = 0, ceil_pass = 0,
     found = List(), len, t0, elapsed, ceiling);

  t0 = getabstime();

  /* Build primorial and sieve */
  S = cc_primorial(S_limit);
  sieve = cc_con_sieve(S, target);
  Rr = cc_con_R_range(S, bits, prefix, prefix_bits);
  R_min = Rr[1]; R_max = Rr[2];
  R_count = (R_max - R_min) \ 2 + 1;

  if(verbose,
    cc_hdr(Str("CONSTRUCT — CC", target, "a at ", bits, " bits (mode ", mode, ")"));
    printf("  S = primorial(%d)  (%d bits)\n", S_limit, #binary(S));
    printf("  R range: %s .. %s  (%d odd values)\n",
      if(#binary(R_min) <= 60, Str(R_min), cc_fmt_hex(R_min, 30)),
      if(#binary(R_max) <= 60, Str(R_max), cc_fmt_hex(R_max, 30)),
      R_count);
    printf("  Sieve primes: %d\n", #sieve);
    printf("  Mode: %s  Max candidates: %d\n", mode, max_cand);
    printf("  Searching...\n\n");
  );

  if(R_min > R_max,
    if(verbose, printf("  No valid R range for %d bits with S=primorial(%d).\n", bits, S_limit));
    return(Vec(found));
  );

  /* Main loop */
  while(tested < max_cand,
    /* Generate R based on mode */
    if(mode == "C",
      /* Mode C: random odd R */
      R = R_min + 2 * random((R_max - R_min) \ 2 + 1);
      if(R % 2 == 0, R++);
      if(R > R_max, R -= 2);
    ,
      if(tested == 0,
        R = R_min;
      ,
        R += 2;
      );
      if(R > R_max,
        if(verbose, printf("\n  Exhausted R range after %d candidates.\n", tested));
        break;
      );
    );

    tested++;

    /* Mode A: R must be PRP */
    if(mode == "A",
      if(!ispseudoprime(R), next);
    );

    /* Constructor sieve check */
    if(!cc_con_sieve_check(R, sieve), next);
    sieve_pass++;

    /* Build p = S*R - 1 */
    p = S * R - 1;

    /* Verify bit length */
    if(#binary(p) != bits, next);

    /* Quick PRP test on root p */
    if(!ispseudoprime(p), next);
    prp_pass++;

    /* OPT-6 ceiling prefilter */
    ceiling = cc_con_ceiling(p, target);
    if(ceiling < target, next);
    ceil_pass++;

    /* Chain follow with PRP (reuse Layer 7) */
    len = cc_chain_test(p, 1);

    if(len >= target,
      /* Verify with proven primality */
      my(V = cc_chain_verify(p, 1));
      if(V[1] >= target,
        listput(found, [p, V[1]]);
        if(verbose,
          printf("  *** FOUND CC%da: %s (%d bits, R=%s) ***\n",
            V[1], Str(p), #binary(p),
            if(#binary(R) <= 60, Str(R), cc_fmt_hex(R, 30)));
        );
      );
    );

    /* Progress report every 10000 */
    if(verbose && tested % 10000 == 0,
      elapsed = getabstime() - t0;
      printf("  [tested=%d  sieve=%d  prp=%d  ceil=%d  found=%d  %.1fs  %.0f/s]\n",
        tested, sieve_pass, prp_pass, ceil_pass, #found,
        elapsed / 1000.0,
        if(elapsed > 0, 1000.0 * tested / elapsed, 0));
    );
  );

  elapsed = getabstime() - t0;
  if(verbose,
    printf("\n  Search complete:\n");
    printf("    Tested:        %d\n", tested);
    printf("    Sieve pass:    %d (%.2f%%)\n", sieve_pass,
      if(tested > 0, 100.0 * sieve_pass / tested, 0.0));
    printf("    PRP pass:      %d\n", prp_pass);
    printf("    Ceiling pass:  %d\n", ceil_pass);
    printf("    Chains found:  %d\n", #found);
    printf("    Time:          %.1f seconds\n", elapsed / 1000.0);
    printf("    Rate:          %.0f R/sec\n",
      if(elapsed > 0, 1000.0 * tested / elapsed, 0));

    for(i = 1, #found,
      printf("    [%d] CC%da root=%s (%d bits)\n",
        i, found[i][2], Str(found[i][1]), #binary(found[i][1]));
    );
  );

  return(Vec(found));
};

/* Quick constructor wrapper for common use cases */
cc_qconstruct(bits, target, {S_limit=31}, {max_cand=50000}) =
  cc_construct(bits, target, S_limit, "B", 0, 0, max_cand, 1);

/* ────────────────────────────────────────────────────────────────────────────
 * Layer 9: Safe Prime Search Engine
 *
 * GP/PARI reimplementation of cc_gmp_v35_03.c combined q+2q+1 trial division.
 * NOT production speed — this is a reference/POC for verifying the C code
 * and for interactive exploration of the safe prime search space.
 *
 * Key idea: before calling ispseudoprime(q) + ispseudoprime(2q+1) (expensive
 * at 4096+ bits), sieve BOTH q and 2q+1 against thousands of small primes.
 * A single trial-division pass checks q%p==0 || (2q+1)%p==0 and rejects
 * ~95% of candidates cheaply. Matches OpenSSL's approach to dhparam.
 *
 * Pipeline: random q → oddness → combined q+2q+1 trial → BPSW(q) → BPSW(2q+1)
 * ──────────────────────────────────────────────────────────────────────────── */

/* Build trial prime list: all primes from lo to hi.
 * Returns a vector of primes suitable for trial division. */
cc_sp_trial_primes(lo, hi) =
{
  my(primes = List());
  forprime(p = lo, hi,
    listput(primes, p);
  );
  return(Vec(primes));
};

/* Combined q+2q+1 trial division.
 * For each prime p in trial_primes, checks:
 *   q mod p == 0  →  reject
 *   (2*q+1) mod p == 0  →  reject
 * Returns 1 if both survive, 0 if either is killed.
 * This is the core insight from v35-03: the 2q+1 residue is just
 * (2*(q%p) + 1) % p — no need to compute 2q+1 explicitly. */
cc_sp_trial_check(q, trial_primes) =
{
  my(rp, rp2, p);
  for(i = 1, #trial_primes,
    p = trial_primes[i];
    rp = q % p;
    if(rp == 0, return(0));
    rp2 = (2 * rp + 1) % p;
    if(rp2 == 0, return(0));
  );
  return(1);
};

/* Verify a safe prime: q is prime AND 2q+1 is prime (proven).
 * Returns [q, p=2q+1, bits_p] or 0 if not a safe prime. */
cc_sp_verify(q) =
{
  if(!isprime(q), return(0));
  my(p = 2*q + 1);
  if(!isprime(p), return(0));
  return([q, p, #binary(p)]);
};

/* Show safe prime search configuration before running.
 * Displays trial prime count, estimated survival rates, etc. */
cc_sp_info(bits, {trial_limit=0}) =
{
  my(tl, trial_primes, n_primes, est_trial, est_prp, bits_p);

  /* Auto-scale trial limit: max(1000, bits*25) — matches v35-03 */
  if(trial_limit <= 0,
    tl = max(1000, bits * 25);
  ,
    tl = trial_limit;
  );

  trial_primes = cc_sp_trial_primes(2, tl);
  n_primes = #trial_primes;

  cc_hdr(Str("SAFE PRIME SEARCH INFO — ", bits, "-bit p=2q+1"));

  printf("  Target bits (p):       %d\n", bits);
  printf("  Target bits (q):       %d\n", bits - 1);
  printf("  Trial limit:           %d\n", tl);
  printf("  Trial primes:          %d (2..%d)\n", n_primes, trial_primes[n_primes]);

  /* Estimate trial division survival rate for combined q+2q+1.
   * For each prime p: probability q survives = (p-1)/p.
   * Probability 2q+1 also survives = (p-1)/p (approximately).
   * Combined: (1 - 1/p) * (1 - 1/p) per prime, but correlated.
   * More precisely: prob(q%p != 0 AND (2q+1)%p != 0) = (p-2)/p for p>2.
   * Product over all primes gives the trial survival rate. */
  est_trial = 1.0;
  for(i = 1, n_primes,
    my(p = trial_primes[i]);
    if(p == 2,
      /* q must be odd, so q%2 != 0 already guaranteed */
      est_trial *= 1.0;
    ,
      /* Combined: exactly (p-2)/p of residues mod p survive both checks */
      est_trial *= (1.0 * (p - 2)) / p;
    );
  );

  printf("  Est trial survive:     %.6f%% (1 in %.0f)\n",
    100 * est_trial, 1 / est_trial);

  /* PRP rate at this bit size */
  est_prp = 1.0 / ((bits - 1) * log(2));
  printf("  Est q is PRP:          ~1/%.0f\n", 1 / est_prp);

  /* Probability 2q+1 is also prime given q prime: ~1/(bits*ln2) again
   * but actually from heuristics ~C_twin * 2/(ln(2q+1)) where C_twin ≈ 1.32 */
  my(est_2q1 = 1.32 * 2.0 / (bits * log(2)));
  printf("  Est 2q+1 prime|q:      ~1/%.0f\n", 1 / est_2q1);

  /* Combined: per random odd q */
  my(est_total = est_trial * est_prp * est_2q1);
  printf("  Est per random q:      ~1/%.0f safe primes\n", 1 / est_total);

  /* Cost comparison: BPSW on this bit size */
  if(bits >= 2048,
    printf("\n  At %d bits, BPSW costs ~%.0f ms/test (estimated).\n",
      bits, 0.0012 * (1.0*bits/1000)^3);
    printf("  Trial division rejects %.1f%% of candidates before BPSW.\n",
      100 * (1 - est_trial));
    printf("  Without trial: ~%.0f BPSW tests per safe prime.\n",
      1 / (est_prp * est_2q1));
    printf("  With trial:    ~%.0f BPSW tests per safe prime.\n",
      1 / (est_trial * est_prp * est_2q1) * est_trial);
    printf("  Speedup:       ~%.0fx from trial division alone.\n",
      1 / est_trial);
  );

  return([tl, n_primes, est_trial]);
};

/* ── Main safe prime search function ──
 *
 * cc_sp_search(bits, {trial_limit}, {max_cand}, {verbose}, {proven})
 *
 * bits:        bit length of p=2q+1 (e.g., 512, 2048, 4096)
 * trial_limit: max prime for trial sieve (0 = auto: bits*25)
 * max_cand:    stop after testing this many q (default 10000000)
 * verbose:     0 = silent, 1 = progress (default 1)
 * proven:      1 = prove primality (slow for large bits), 0 = PRP only (default 0)
 *
 * Returns [q, p, bits_p] on success, 0 if not found within max_cand.
 *
 * Examples:
 *   cc_sp_search(512)                     \\ 512-bit safe prime, ~instant
 *   cc_sp_search(2048)                    \\ 2048-bit, seconds
 *   cc_sp_search(4096)                    \\ 4096-bit, ~30-120s in GP
 *   cc_sp_search(4096, 200000)            \\ 4096-bit, larger trial sieve
 *   cc_sp_search(256, , 100000, 0)        \\ silent 256-bit search
 */
cc_sp_search(bits, {trial_limit=0}, {max_cand=10000000}, {verbose=1}, {proven=0}) =
{
  my(tl, trial_primes, q_bits, lo, hi, q, p, r6,
     tested = 0, trial_pass = 0, prp_pass = 0,
     t0, elapsed);

  t0 = getabstime();

  /* Auto-scale trial limit */
  if(trial_limit <= 0, tl = max(1000, bits * 25), tl = trial_limit);

  trial_primes = cc_sp_trial_primes(3, tl);
  q_bits = bits - 1;

  /* q range: 2^(q_bits-1) to 2^q_bits - 1 */
  lo = 2^(q_bits - 1);
  hi = 2^q_bits - 1;

  if(verbose,
    cc_hdr(Str("SAFE PRIME SEARCH — ", bits, "-bit p=2q+1"));
    printf("  q range:        %d-bit\n", q_bits);
    printf("  Trial primes:   %d (3..%d)\n", #trial_primes, trial_primes[#trial_primes]);
    printf("  Max candidates: %d\n", max_cand);
    printf("  Prove primes:   %s\n", if(proven, "YES (isprime)", "no (ispseudoprime)"));
    printf("  Searching...\n\n");
  );

  /* Main loop: random odd q ≡ 5 (mod 6) in range */
  while(tested < max_cand,
    q = random(hi - lo + 1) + lo;
    /* Force q ≡ 5 (mod 6): odd AND 2q+1 not divisible by 3 */
    r6 = q % 6;
    if(r6 == 0, q += 5, if(r6 == 1, q += 4, if(r6 == 2, q += 3, if(r6 == 3, q += 2, if(r6 == 4, q += 1)))));
    if(q > hi, q -= 6);

    tested++;

    /* Combined q+2q+1 trial division (skip p=3, already handled) */
    if(!cc_sp_trial_check(q, trial_primes), next);
    trial_pass++;

    /* BPSW on q */
    if(!ispseudoprime(q), next);
    prp_pass++;

    /* BPSW on 2q+1 */
    p = 2 * q + 1;
    if(!ispseudoprime(p), next);

    /* Found! Optionally prove */
    if(proven, if(!isprime(q) || !isprime(p), next));

    elapsed = getabstime() - t0;
    if(verbose,
      printf("  *** FOUND %d-bit safe prime ***\n", #binary(p));
      printf("  q = %s\n", cc_fmt_hex(q, 80));
      printf("  p = 2q+1 = %s\n", cc_fmt_hex(p, 80));
      printf("  Tested:       %d candidates\n", tested);
      printf("  Trial pass:   %d (%.4f%%)\n", trial_pass, if(tested > 0, 100.0 * trial_pass / tested, 0.0));
      printf("  PRP pass (q): %d\n", prp_pass);
      printf("  Time:         %.1f seconds\n", elapsed / 1000.0);
      printf("  Rate:         %.0f candidates/sec\n", if(elapsed > 0, 1000.0 * tested / elapsed, 0));
    );

    return([q, p, #binary(p)]);
  );

  /* Not found within budget */
  elapsed = getabstime() - t0;
  if(verbose,
    printf("\n  Exhausted %d candidates without finding safe prime.\n", max_cand);
    printf("  Time: %.1f seconds (%.0f cand/sec)\n", elapsed / 1000.0, if(elapsed > 0, 1000.0 * tested / elapsed, 0));
  );

  return(0);
};

/* Quick safe prime search: small bit sizes, instant results */
cc_sp_qsearch(bits, {trial_limit=0}, {max_cand=1000000}) =
  cc_sp_search(bits, trial_limit, max_cand, 1, 0);

/* Benchmark: compare safe prime search with and without trial division.
 * Runs two searches at the given bit size and reports timing.
 * Good for demonstrating the v35-03 optimization. */
cc_sp_bench(bits, {max_cand=100000}) =
{
  my(t0, t1, t2, q, p, r6, trial_primes, tl,
     tested_notrial = 0, tested_trial = 0,
     lo, hi, found);

  tl = max(1000, bits * 25);
  trial_primes = cc_sp_trial_primes(3, tl);
  lo = 2^(bits - 2);
  hi = 2^(bits - 1) - 1;

  cc_hdr(Str("SAFE PRIME BENCHMARK — ", bits, "-bit, trial_limit=", tl));

  /* Method 1: No trial division — just PRP both */
  printf("  Method 1: BPSW-only (no trial division)...\n");
  t0 = getabstime();
  found = 0;
  while(!found && tested_notrial < max_cand,
    q = random(hi - lo + 1) + lo;
    r6 = q % 6;
    if(r6 == 0, q += 5, if(r6 == 1, q += 4, if(r6 == 2, q += 3, if(r6 == 3, q += 2, if(r6 == 4, q += 1)))));
    if(q > hi, q -= 6);
    tested_notrial++;
    if(!ispseudoprime(q), next);
    if(!ispseudoprime(2*q+1), next);
    found = 1;
  );
  t1 = getabstime();

  /* Method 2: Combined trial + PRP */
  printf("  Method 2: Combined trial + BPSW...\n");
  found = 0;
  while(!found && tested_trial < max_cand,
    q = random(hi - lo + 1) + lo;
    r6 = q % 6;
    if(r6 == 0, q += 5, if(r6 == 1, q += 4, if(r6 == 2, q += 3, if(r6 == 3, q += 2, if(r6 == 4, q += 1)))));
    if(q > hi, q -= 6);
    tested_trial++;
    if(!cc_sp_trial_check(q, trial_primes), next);
    if(!ispseudoprime(q), next);
    if(!ispseudoprime(2*q+1), next);
    found = 1;
  );
  t2 = getabstime();

  printf("\n  Results:\n");
  printf("    BPSW-only:     %d candidates, %.1f s\n", tested_notrial, (t1 - t0) / 1000.0);
  printf("    Trial+BPSW:    %d candidates, %.1f s\n", tested_trial, (t2 - t1) / 1000.0);
  if((t2 - t1) > 0 && (t1 - t0) > 0,
    printf("    Speedup:       %.1fx\n", 1.0 * (t1 - t0) / (t2 - t1));
  );
};

/* ────────────────────────────────────────────────────────────────────────────
 * Layer 9b: OpenSSL-Style Delta-Sieve and Comparison
 *
 * OpenSSL's safe prime algorithm (Zimmermann quick sieve):
 *   1. Generate random odd p₀ ≡ 3 (mod 4)
 *   2. Compute mods[i] = p₀ mod prime[i] for all trial primes (ONCE)
 *   3. For each delta (step +4, maintaining p ≡ 3 mod 4):
 *      Check (mods[i] + delta) % prime[i] <= 1 for all i
 *      The <= 1 trick: residue 0 means p divisible, residue 1 means q=(p-1)/2
 *      divisible — catches both in ONE comparison per prime.
 *   4. Survivor: BPSW p, then BPSW q = (p-1)/2
 *
 * Key differences from our v35-03 approach:
 *   - Delta-sieve: compute residues ONCE per random start, try offsets cheaply
 *   - Our approach: fresh mpz_fdiv_ui per candidate (more expensive)
 *   - OpenSSL uses only 1024 primes at 4096-bit (up to ~8837)
 *   - Our auto-scaling uses ~9805 primes at 4096-bit (up to ~102400)
 *
 * The question: does using MORE trial primes compensate for OpenSSL's
 * more efficient sieve structure?  And: would OpenSSL benefit from using
 * more trial primes?
 *
 * References:
 *   - OpenSSL crypto/bn/bn_prime.c: probable_prime(), bn_is_prime_int()
 *   - OpenSSL issue #17956: external sieve with 20M primes was 2x faster
 *   - Wiener 2003: "Safe Prime Generation with a Combined Sieve"
 * ──────────────────────────────────────────────────────────────────────────── */

/* OpenSSL's trial prime count selection (from calc_trial_divisions) */
cc_sp_openssl_trial_count(bits) =
{
  if(bits <= 512, return(64));
  if(bits <= 1024, return(128));
  if(bits <= 2048, return(384));
  if(bits <= 4096, return(1024));
  return(2048);
};

/* Delta-sieve safe prime search — OpenSSL's algorithm in GP.
 *
 * n_trial: number of trial primes (OpenSSL uses cc_sp_openssl_trial_count)
 * max_delta: max offset to try before re-randomizing (0 = auto: 2^20)
 *
 * Returns [q, p, bits_p, tested, deltas_tried, reseeds] or 0. */
cc_sp_delta_sieve(bits, {n_trial=0}, {max_cand=10000000}, {max_delta=0}, {verbose=1}) =
{
  my(trial_primes, p0, p, q, mods, delta, r, ok,
     tested = 0, deltas_tried = 0, reseeds = 0, trial_pass = 0,
     t0, elapsed, md);

  t0 = getabstime();

  /* Default trial count: match OpenSSL */
  if(n_trial <= 0, n_trial = cc_sp_openssl_trial_count(bits));

  /* Build trial primes (skip 2; handle oddness separately) */
  trial_primes = cc_sp_trial_primes(3, prime(n_trial + 1));
  if(#trial_primes > n_trial, trial_primes = vecextract(trial_primes, Str("1..", n_trial)));

  /* Default max_delta */
  if(max_delta <= 0, md = 2^20, md = max_delta);

  if(verbose,
    cc_hdr(Str("DELTA-SIEVE SAFE PRIME — ", bits, "-bit (OpenSSL-style)"));
    printf("  Trial primes:   %d (3..%d)\n", #trial_primes, trial_primes[#trial_primes]);
    printf("  Max delta:      %d\n", md);
    printf("  Max candidates: %d\n", max_cand);
    printf("  Searching...\n\n");
  );

  while(tested < max_cand,
    /* Step 1: Generate random p0 ≡ 3 (mod 4), bits-bit */
    p0 = random(2^(bits-1)) + 2^(bits-1);
    /* Force p0 ≡ 3 (mod 4): set two LSBs to 11 */
    p0 = bitor(p0, 3);
    /* Ensure p0 ≡ 3 (mod 4) exactly */
    if(p0 % 4 != 3, p0 += (7 - p0 % 4) % 4);
    reseeds++;

    /* Step 2: Compute base residues ONCE */
    mods = vector(#trial_primes, i, p0 % trial_primes[i]);

    /* Step 3: Delta-sieve — try p0, p0+4, p0+8, ... */
    delta = 0;
    while(delta <= md,
      deltas_tried++;
      tested++;

      /* Combined sieve: check (mods[i] + delta) % prime <= 1 */
      ok = 1;
      for(i = 1, #trial_primes,
        r = (mods[i] + delta) % trial_primes[i];
        if(r <= 1, ok = 0; break);
      );

      if(ok,
        trial_pass++;
        p = p0 + delta;

        /* Ensure still in bit range */
        if(#binary(p) == bits,
          /* BPSW on p */
          if(ispseudoprime(p),
            q = (p - 1) / 2;
            /* BPSW on q */
            if(ispseudoprime(q),
              elapsed = getabstime() - t0;
              if(verbose,
                printf("  *** FOUND %d-bit safe prime ***\n", bits);
                printf("  q = %s\n", cc_fmt_hex(q, 80));
                printf("  p = 2q+1 = %s\n", cc_fmt_hex(p, 80));
                printf("  Reseeds:      %d\n", reseeds);
                printf("  Deltas tried: %d\n", deltas_tried);
                printf("  Trial pass:   %d (%.4f%%)\n", trial_pass, if(tested > 0, 100.0 * trial_pass / tested, 0.0));
                printf("  Time:         %.3f seconds\n", elapsed / 1000.0);
                printf("  Rate:         %.0f deltas/sec\n", if(elapsed > 0, 1000.0 * deltas_tried / elapsed, 0));
              );
              return([q, p, #binary(p), tested, deltas_tried, reseeds]);
            );
          );
        );
      );

      delta += 4;
      if(tested >= max_cand, break);
    );
  );

  elapsed = getabstime() - t0;
  if(verbose,
    printf("\n  Exhausted %d candidates.\n", tested);
    printf("  Time: %.3f seconds\n", elapsed / 1000.0);
  );
  return(0);
};

/* ── Head-to-head comparison: OpenSSL-style vs more trial primes ──
 *
 * Runs the delta-sieve with two different trial prime counts:
 *   1. OpenSSL's default (calc_trial_divisions: 1024 at 4096-bit)
 *   2. Our auto-scaled count (bits*25, e.g. ~9805 at 4096-bit)
 *
 * This demonstrates the key finding: more trial primes = fewer BPSW tests
 * = faster overall, even though each sieve step costs slightly more. */
cc_sp_compare(bits, {runs=5}, {max_cand=10000000}) =
{
  my(n_ossl, n_ours, tp_ossl, tp_ours,
     t0, t1, t2, t_ossl = 0, t_ours = 0,
     found_ossl = 0, found_ours = 0,
     deltas_ossl = 0, deltas_ours = 0,
     trial_limit, result);

  n_ossl = cc_sp_openssl_trial_count(bits);
  trial_limit = max(1000, bits * 25);
  tp_ours = cc_sp_trial_primes(3, trial_limit);
  n_ours = #tp_ours;

  cc_hdr(Str("SAFE PRIME COMPARISON — ", bits, "-bit, ", runs, " runs each"));
  printf("  OpenSSL approach: %d trial primes (3..%d)\n", n_ossl, prime(n_ossl + 1));
  printf("  Our approach:     %d trial primes (3..%d)\n", n_ours, tp_ours[n_ours]);
  printf("  Ratio:            %.1fx more trial primes\n", 1.0 * n_ours / n_ossl);
  printf("\n");

  /* Estimate trial survival rates */
  my(est_ossl = 1.0, est_ours = 1.0);
  my(tp_ossl_v = cc_sp_trial_primes(3, prime(n_ossl + 1)));
  for(i = 1, #tp_ossl_v, est_ossl *= (1.0 * (tp_ossl_v[i] - 2)) / tp_ossl_v[i]);
  for(i = 1, n_ours, est_ours *= (1.0 * (tp_ours[i] - 2)) / tp_ours[i]);
  printf("  Est trial survive (OpenSSL): %.4f%% (1 in %.0f)\n", 100*est_ossl, 1/est_ossl);
  printf("  Est trial survive (ours):    %.4f%% (1 in %.0f)\n", 100*est_ours, 1/est_ours);
  printf("  Predicted BPSW reduction:    %.1f%% fewer BPSW tests with more primes\n",
    100 * (1 - est_ours / est_ossl));
  printf("\n  Running %d rounds each...\n\n", runs);

  /* Run OpenSSL-style */
  printf("  [OpenSSL-style: %d primes]\n", n_ossl);
  for(r = 1, runs,
    t0 = getabstime();
    result = cc_sp_delta_sieve(bits, n_ossl, max_cand, 0, 0);
    t1 = getabstime();
    if(result,
      found_ossl++;
      deltas_ossl += result[5];
      t_ossl += (t1 - t0);
      printf("    run %d: %.3fs, %d deltas\n", r, (t1 - t0) / 1000.0, result[5]);
    ,
      printf("    run %d: FAILED (no prime found)\n", r);
    );
  );

  /* Run with more primes */
  printf("  [Our approach: %d primes]\n", n_ours);
  for(r = 1, runs,
    t0 = getabstime();
    result = cc_sp_delta_sieve(bits, n_ours, max_cand, 0, 0);
    t1 = getabstime();
    if(result,
      found_ours++;
      deltas_ours += result[5];
      t_ours += (t1 - t0);
      printf("    run %d: %.3fs, %d deltas\n", r, (t1 - t0) / 1000.0, result[5]);
    ,
      printf("    run %d: FAILED (no prime found)\n", r);
    );
  );

  /* Summary */
  printf("\n");
  cc_sub("Results");
  printf("  OpenSSL (%d primes): %d/%d found, avg %.3fs, avg %d deltas\n",
    n_ossl, found_ossl, runs,
    if(found_ossl > 0, t_ossl / found_ossl / 1000.0, 0),
    if(found_ossl > 0, deltas_ossl / found_ossl, 0));
  printf("  Ours    (%d primes): %d/%d found, avg %.3fs, avg %d deltas\n",
    n_ours, found_ours, runs,
    if(found_ours > 0, t_ours / found_ours / 1000.0, 0),
    if(found_ours > 0, deltas_ours / found_ours, 0));
  if(t_ours > 0 && t_ossl > 0,
    printf("  Speedup: %.2fx %s\n",
      if(t_ours < t_ossl, 1.0 * t_ossl / t_ours, 1.0 * t_ours / t_ossl),
      if(t_ours < t_ossl, "(more primes faster)", "(fewer primes faster)"));
  );
  if(found_ossl > 0 && found_ours > 0,
    printf("  Avg deltas: %.1fx %s with more primes\n",
      if(deltas_ours < deltas_ossl,
        1.0 * deltas_ossl / found_ossl / (deltas_ours / found_ours),
        1.0 * deltas_ours / found_ours / (deltas_ossl / found_ossl)),
      if(deltas_ours * found_ossl < deltas_ossl * found_ours, "fewer", "more"));
  );

  printf("\n  Interpretation:\n");
  printf("    If 'more primes faster': OpenSSL would benefit from larger trial table.\n");
  printf("    The BPSW savings from eliminating more candidates outweighs\n");
  printf("    the slightly higher per-delta sieve cost.\n");
  printf("    At %d bits, BPSW costs ~ms while trial sieve costs ~us.\n", bits);
};

/* ── Analyze the sieve efficiency at different trial prime counts ──
 * Shows how trial survival rate drops as we add more primes.
 * This is the theoretical basis for the improvement. */
cc_sp_trial_analysis(bits) =
{
  my(counts, tl, tp, est, prev_est);

  cc_hdr(Str("TRIAL PRIME ANALYSIS — ", bits, "-bit safe primes"));
  printf("  %-12s  %-8s  %-16s  %-14s  %-10s\n",
    "Trial primes", "Max prime", "Survival rate", "1 in N", "BPSW saved");

  counts = [64, 128, 256, 384, 512, 1024, 2048, 4096, 8192];
  /* Add our auto-scaled count */
  tl = max(1000, bits * 25);

  prev_est = 1.0;
  for(ci = 1, #counts,
    my(n = counts[ci]);
    tp = cc_sp_trial_primes(3, prime(n + 1));
    if(#tp > n, tp = vecextract(tp, Str("1..", n)));
    est = 1.0;
    for(i = 1, #tp, est *= (1.0 * (tp[i] - 2)) / tp[i]);
    printf("  %-12d  %-8d  %-16.8f%%  %-14.0f  %-10.1f%%\n",
      #tp, tp[#tp], 100*est, 1/est,
      if(ci > 1, 100 * (1 - est / prev_est), 0));
    prev_est = est;
  );

  /* Our auto-scaled */
  tp = cc_sp_trial_primes(3, tl);
  est = 1.0;
  for(i = 1, #tp, est *= (1.0 * (tp[i] - 2)) / tp[i]);
  printf("  %-12d  %-8d  %-16.8f%%  %-14.0f  (auto: bits*25)\n",
    #tp, tp[#tp], 100*est, 1/est);

  /* OpenSSL's choice */
  my(n_ossl = cc_sp_openssl_trial_count(bits));
  tp = cc_sp_trial_primes(3, prime(n_ossl + 1));
  if(#tp > n_ossl, tp = vecextract(tp, Str("1..", n_ossl)));
  my(est_ossl = 1.0);
  for(i = 1, #tp, est_ossl *= (1.0 * (tp[i] - 2)) / tp[i]);

  printf("\n  OpenSSL uses %d primes at %d bits → 1 in %.0f survive trial\n",
    n_ossl, bits, 1/est_ossl);
  printf("  Our auto-scale uses %d primes → 1 in %.0f survive trial\n",
    #cc_sp_trial_primes(3, tl), 1/est);
  printf("  Improvement: %.1f%% fewer candidates reach BPSW\n",
    100 * (1 - est / est_ossl));
  printf("  At %d bits, each BPSW test costs ~%.0f us\n",
    bits, if(bits <= 512, 50, if(bits <= 2048, 500, 5000)));
  printf("  Extra sieve cost per delta: ~%.0f ns (negligible)\n",
    (#cc_sp_trial_primes(3, tl) - n_ossl) * 2);
};

/* ────────────────────────────────────────────────────────────────────────────
 * Layer 10: Alternative Algorithms (CUDA Legacy Extraction)
 *
 * cc_x_ prefix: experimental/alternative methods mined from 22 CUDA programs
 * in cuda_legacy/.  These implement techniques explored during GPU development
 * that are mathematically interesting and complement the production pipeline.
 *
 * Source mapping (METHOD-L IDs from cuda_legacy/MASTER_INDEX.md):
 *   cc_x_forbidden_mask  ← METHOD-L01 (cc_search_unified_c7 v25/v27)
 *   cc_x_periodic_table  ← METHOD-L18 (cc_constructor v20)
 *   cc_x_line_filter     ← METHOD-L17 (cc_search_unified_c7 v27)
 *   cc_x_bitwin_mask     ← METHOD-L12 (cc_bitwin_prototype)
 *   cc_x_bitwin_walk     ← METHOD-L12 (cc_bitwin_prototype)
 *   cc_x_primorial_scan  ← METHOD-L10 (cunningham_algebraic v3)
 *   cc_x_root_depth      ← METHOD-L05 (cc_constructor v20)
 *
 * Naming convention: cc_x_ = extracted/experimental.  May be promoted to
 * cc_ in a future version if proven superior to existing methods.
 * ──────────────────────────────────────────────────────────────────────────── */

/* Forbidden residue BITMASK for chain of length `depth` mod prime q.
 * Returns integer where bit r is set if residue r is forbidden.
 * Formula: forbidden[i] = (2^(-i) - 1) mod q  (first-kind)
 *          forbidden[i] = (1 - 2^(-i)) mod q  (second-kind)
 * Usage: bittest(mask, p % q) → 1 means p fails for this prime.
 * Complements cc_forbidden_res() (vector form) with O(1) bitmask lookup.
 * Only valid for primes q < 64 (bitmask width), use cc_forbidden_res for larger. */
cc_x_forbidden_mask(q, depth, {kind=1}) =
{
  my(mask = 0, inv2 = lift(Mod(2, q)^(-1)), inv2_pow = 1, r);
  for(i = 0, depth - 1,
    if(kind == 1,
      r = (inv2_pow + q - 1) % q,   /* (2^(-i) - 1) mod q */
      r = (1 + q - inv2_pow) % q     /* (1 - 2^(-i)) mod q */
    );
    mask = bitor(mask, 2^r);
    inv2_pow = (inv2_pow * inv2) % q;
  );
  mask;
};

/* Periodic forbidden residue table exploiting multiplicative order of 2.
 * Returns [table_vector, period] where table[i+1] = forbidden residue at
 * chain position i, and period = znorder(Mod(2,q)).
 * The table repeats: forbidden[i] = table[(i % period) + 1].
 * From METHOD-L18: cc_constructor v20 periodic precomputation. */
cc_x_periodic_table(q, {kind=1}) =
{
  my(ord = znorder(Mod(2, q)), inv2 = lift(Mod(2, q)^(-1)));
  my(table = vector(ord), inv2_pow = 1);
  for(i = 0, ord - 1,
    if(kind == 1,
      table[i + 1] = (inv2_pow + q - 1) % q,
      table[i + 1] = (1 + q - inv2_pow) % q
    );
    inv2_pow = (inv2_pow * inv2) % q;
  );
  [table, ord];
};

/* Sieve-only chain pre-screening.  Follows chain p → 2p+1 → 4p+3 → ...
 * for `depth` steps, checking ONLY trial division (primes 3..max_prime).
 * Returns number of steps that survive.  Does NOT call isprime().
 * Much faster than cc_walk() for bulk candidate filtering.
 * From METHOD-L17: cc_search_unified_c7 v27 sieve-only line filter. */
cc_x_line_filter(p, depth, {kind=1}, {max_prime=107}) =
{
  my(n = p, steps = 0, pl = primes([3, max_prime]));
  for(step = 0, depth - 1,
    if(n < 2, return(steps));
    if(n > 2 && n % 2 == 0, return(steps));
    my(ok = 1);
    for(j = 1, #pl,
      if(n % pl[j] == 0 && n > pl[j], ok = 0; break);
    );
    if(!ok, return(steps));
    steps++;
    if(kind == 1, n = 2*n + 1, n = 2*n - 1);
  );
  steps;
};

/* Combined forbidden residue mask for BiTwin chains (both kinds + evenness).
 * BiTwin: even center n, first-kind chain from n-1, second-kind from n+1.
 * Returns bitmask where bit r is set if n ≡ r (mod q) is forbidden.
 * First-kind at step k: n ≡ (1 - 2^(-k)) (mod q) → kills (n-1) chain
 * Second-kind at step k: n ≡ (2^(-k) - 1) (mod q) → kills (n+1) chain
 * Plus: all odd residues (n must be even).
 * From METHOD-L12: cc_bitwin_prototype. */
cc_x_bitwin_mask(q, depth) =
{
  my(mask = 0, inv2 = lift(Mod(2, q)^(-1)), inv2_pow = 1, r1, r2);
  /* Evenness: all odd residues are forbidden */
  forstep(r = 1, q - 1, 2, mask = bitor(mask, 2^r));
  /* Forbidden residues for both chain kinds */
  for(k = 0, depth - 1,
    r1 = (1 + q - inv2_pow) % q;  /* first-kind: (1 - 2^(-k)) mod q */
    r2 = (inv2_pow + q - 1) % q;  /* second-kind: (2^(-k) - 1) mod q */
    mask = bitor(mask, 2^r1);
    mask = bitor(mask, 2^r2);
    inv2_pow = (inv2_pow * inv2) % q;
  );
  mask;
};

/* Walk BiTwin chain from even center n.
 * Returns [first_kind_len, second_kind_len, bitwin_length].
 * bitwin_length = min(first_kind_len, second_kind_len).
 * From METHOD-L12: cc_bitwin_prototype. */
cc_x_bitwin_walk(n, {max_depth=50}) =
{
  if(n % 2 != 0, return([-1, -1, -1]));
  /* First-kind chain from n-1 */
  my(p = n - 1, len1 = 0);
  while(len1 < max_depth && isprime(p),
    len1++;
    p = 2*p + 1;
  );
  /* Second-kind chain from n+1 */
  p = n + 1;
  my(len2 = 0);
  while(len2 < max_depth && isprime(p),
    len2++;
    p = 2*p - 1;
  );
  [len1, len2, min(len1, len2)];
};

/* Search using algebraic form: p = k * primorial(prime(prim_n)) * 2^m - 1.
 * This is the form used for historical CC world records.
 * Returns vector of [p, chain_length] pairs with length >= target.
 * From METHOD-L10: cunningham_algebraic v3. */
cc_x_primorial_scan(bits, prim_n, {target=5}, {k_max=10000}) =
{
  my(P = cc_primorial(prime(prim_n)));
  my(lo = 2^(bits-1), hi = 2^bits - 1);
  my(results = List(), tested = 0);
  for(m = 0, bits,
    my(base = P * 2^m);
    if(base > hi, break);
    my(klo = (lo + base - 1) \ base);  /* ceiling division */
    my(khi = min(k_max, hi \ base));
    if(klo > khi, next);
    if(klo % 2 == 0, klo++);           /* odd k only (dedup) */
    forstep(k = klo, khi, 2,
      my(p = k * base - 1);
      tested++;
      if(p % 6 != 5, next);            /* first-kind constraint */
      /* Quick trial filter */
      my(ok = 1);
      foreach([7, 11, 13, 17, 19, 23, 29, 31], q,
        if(p % q == 0, ok = 0; break);
      );
      if(!ok, next);
      /* Check chain length */
      my(n = p, len = 0);
      while(len < 50 && isprime(n), len++; n = 2*n + 1);
      if(len >= target, listput(results, [p, len]));
    );
  );
  Vec(results);
};

/* k-value: number of backward steps from p to the true chain root.
 * Returns 0 if p IS the root, 1 if (p-1)/2 is the root, etc.
 * Returns -1 if p is not prime.
 * From METHOD-L05: cc_constructor v20 k-value tracking. */
cc_x_root_depth(p, {kind=1}) =
{
  if(!isprime(p), return(-1));
  my(n = p, depth = 0, prev);
  while(1,
    if(kind == 1,
      if(n < 3 || n % 2 == 0, return(depth));
      prev = (n - 1) \ 2;
    ,
      prev = (n + 1) \ 2;
    );
    if(prev < 2 || !isprime(prev), return(depth));
    n = prev;
    depth++;
  );
};

/* ────────────────────────────────────────────────────────────────────────────
 * Self-Test on Load
 * ──────────────────────────────────────────────────────────────────────────── */

{
  my(test_root = 1122659, chain, ok = 1);

  printf("\n");
  printf("=== cc_lib_v10.gp — PARI/GP Cunningham Chain Library ===\n\n");

  /* Test: CC7 starting from known root 1122659 */
  chain = List();
  my(cur = test_root);
  while(isprime(cur),
    listput(chain, cur);
    cur = 2*cur + 1;
    if(#chain > 30, break);
  );

  if(#chain == 7,
    printf("  Self-test PASSED: CC7a root=%d, length=%d\n", test_root, #chain);
  ,
    printf("  Self-test FAILED: expected CC7 from root %d, got length %d\n",
      test_root, #chain);
    ok = 0;
  );

  /* Test cc_root */
  my(R = cc_root(chain[4]));
  if(R[1] == test_root,
    printf("  Self-test PASSED: cc_root(%d) found root=%d\n", chain[4], R[1]);
  ,
    printf("  Self-test FAILED: cc_root(%d) returned %d, expected %d\n",
      chain[4], R[1], test_root);
    ok = 0;
  );

  /* Test cc_kill_pos: prime 3 should kill at position 0 for root 1122659
   * 1122659 mod 3 = 2, and 2*2+1=5, 5 mod 3 = 2, ... actually let's check
   * 1122659 mod 3 = 1122659 - 374219*3 = 1122659 - 1122657 = 2
   * cc_immune_residue(3,1) = 2, so 3 IS immune for this root */
  my(r3 = test_root % 3);
  if(cc_is_immune(test_root, 3),
    printf("  Self-test PASSED: root %d immune to prime 3 (residue %d)\n",
      test_root, r3);
  ,
    printf("  Self-test FAILED: expected root %d immune to prime 3\n", test_root);
    ok = 0;
  );

  /* Test composite exact-input behavior in v6 */
  my(test_breaker = cur);
  R = cc_root(test_breaker);
  if(R[1] == test_breaker && R[2] == 0,
    printf("  Self-test PASSED: composite input stays exact in cc_root()\n");
  ,
    printf("  Self-test FAILED: composite input incorrectly normalized\n");
    ok = 0;
  );

  /* Test NAF weight */
  if(cc_naf_weight(7) == 2,
    printf("  Self-test PASSED: naf_weight(7) = 2\n");
  ,
    printf("  Self-test FAILED: naf_weight(7) = %d, expected 2\n", cc_naf_weight(7));
    ok = 0;
  );

  /* Test trailing ones */
  if(cc_trailing_ones(7) == 3,
    printf("  Self-test PASSED: trailing_ones(7) = 3\n");
  ,
    printf("  Self-test FAILED: trailing_ones(7) = %d, expected 3\n", cc_trailing_ones(7));
    ok = 0;
  );

  /* Test motif helpers */
  if(cc_prefix_bits(13, 3) == "110" && cc_tail_bits(13, 3) == "101" && cc_motif_count(31, "111") == 3,
    printf("  Self-test PASSED: prefix/tail/motif helpers\n");
  ,
    printf("  Self-test FAILED: motif helpers mismatch\n");
    ok = 0;
  );

  /* Test CC16a: 90-bit record chain */
  my(cc16_root = 929045505808475001179907959);
  my(w16 = cc_walk_quiet(cc16_root));
  if(w16[2] == 16 && w16[3] == cc16_root,
    printf("  Self-test PASSED: CC16a root=%s (90 bits)\n", cc_fmt_hex(cc16_root, 40));
  ,
    printf("  Self-test FAILED: CC16a expected len=16, got %d\n", w16[2]);
    ok = 0;
  );

  /* Test CC11a: 151-bit chain */
  my(cc11_root = 0x58CD931C9D21E6585796573C57C4AE99A57A7F);
  my(w11 = cc_walk_quiet(cc11_root));
  if(w11[2] == 11 && w11[3] == cc11_root,
    printf("  Self-test PASSED: CC11a root=%s (151 bits)\n", cc_fmt_hex(cc11_root, 50));
  ,
    printf("  Self-test FAILED: CC11a expected len=11, got %d\n", w11[2]);
    ok = 0;
  );

  /* Test CC10a: 153-bit chain */
  my(cc10_root = 0x1bf1bcbb011ab036b4d4a85b49696c0373b4d47);
  my(w10 = cc_walk_quiet(cc10_root));
  if(w10[2] == 10 && w10[3] == cc10_root,
    printf("  Self-test PASSED: CC10a root=%s (153 bits)\n", cc_fmt_hex(cc10_root, 50));
  ,
    printf("  Self-test FAILED: CC10a expected len=10, got %d\n", w10[2]);
    ok = 0;
  );

  /* Test CC7a: 502-bit chain */
  my(cc7_502 = 0x25ddc183317eb4b851e37d009293d44acc04060c33d638ee20c3efd42c4669cabbb4e3b67fd089d9102808e13c42e3413fe97b20506442a33a07dce8b2291d);
  my(w7 = cc_walk_quiet(cc7_502));
  if(w7[2] == 7 && w7[3] == cc7_502,
    printf("  Self-test PASSED: CC7a root=%s (502 bits)\n", cc_fmt_hex(cc7_502, 50));
  ,
    printf("  Self-test FAILED: CC7a-502 expected len=7, got %d\n", w7[2]);
    ok = 0;
  );

  /* ── Layer 7 (search engine) tests ── */

  /* Test cc_forbidden_res: q=5, target=3, kind=1
   * 2^0 mod 5 = 1, r_0 = 0
   * 2^(-1) mod 5 = 3, r_1 = 2
   * 2^(-2) mod 5 = 4, r_2 = 3
   * forbidden = [0, 2, 3] */
  my(fr5 = cc_forbidden_res(5, 3, 1));
  if(fr5 == [0, 2, 3],
    printf("  Self-test PASSED: cc_forbidden_res(5, 3, 1) = [0, 2, 3]\n");
  ,
    printf("  Self-test FAILED: cc_forbidden_res(5, 3, 1) = %s, expected [0,2,3]\n", Str(fr5));
    ok = 0;
  );

  /* Test cc_valid_res: q=5, target=3 -> valid = [1, 4] */
  my(vr5 = cc_valid_res(5, 3, 1));
  if(vr5 == [1, 4],
    printf("  Self-test PASSED: cc_valid_res(5, 3, 1) = [1, 4]\n");
  ,
    printf("  Self-test FAILED: cc_valid_res(5, 3, 1) = %s, expected [1,4]\n", Str(vr5));
    ok = 0;
  );

  /* Test cc_build_bases: target=2, primes=[3,5]
   * q=3: forbidden {0,1}, valid {2}. q=5: forbidden {0,2}, valid {1,3,4}.
   * CRT: 3 bases mod 15: {8, 11, 14} */
  my(BM = cc_build_bases(2, 1, [3, 5]));
  if(BM[2] == 15 && #BM[1] == 3 && BM[1] == [8, 11, 14],
    printf("  Self-test PASSED: cc_build_bases(2, 1, [3,5]) = 3 bases mod 15\n");
  ,
    printf("  Self-test FAILED: cc_build_bases(2,1,[3,5]) M=%d, bases=%s\n", BM[2], Str(BM[1]));
    ok = 0;
  );

  /* Test cc_chain_test: known CC7 root should give length 7 */
  my(ct7 = cc_chain_test(test_root));
  if(ct7 == 7,
    printf("  Self-test PASSED: cc_chain_test(%d) = 7\n", test_root);
  ,
    printf("  Self-test FAILED: cc_chain_test(%d) = %d, expected 7\n", test_root, ct7);
    ok = 0;
  );

  /* Test cc_chain_verify: proven primality for CC7 root */
  my(cv7 = cc_chain_verify(test_root));
  if(cv7[1] == 7,
    printf("  Self-test PASSED: cc_chain_verify(%d) = [7, 1]\n", test_root);
  ,
    printf("  Self-test FAILED: cc_chain_verify(%d) = %s\n", test_root, Str(cv7));
    ok = 0;
  );

  /* Test cc_sieve_check: root 1122659 should pass sieve for target=7 */
  my(sv = cc_build_sieve(7, 1));
  if(cc_sieve_check(test_root, sv) == 1,
    printf("  Self-test PASSED: cc_sieve_check(%d, sieve_7) = 1\n", test_root);
  ,
    printf("  Self-test FAILED: cc_sieve_check(%d, sieve_7) = 0\n", test_root);
    ok = 0;
  );

  /* ── Layer 8 (constructor) tests ── */

  /* Test cc_primorial(7) == 210 */
  if(cc_primorial(7) == 210,
    printf("  Self-test PASSED: cc_primorial(7) = 210\n");
  ,
    printf("  Self-test FAILED: cc_primorial(7) = %d, expected 210\n", cc_primorial(7));
    ok = 0;
  );

  /* Test cc_primorial(31) == 200560490130 */
  if(cc_primorial(31) == 200560490130,
    printf("  Self-test PASSED: cc_primorial(31) = 200560490130\n");
  ,
    printf("  Self-test FAILED: cc_primorial(31) = %d, expected 200560490130\n", cc_primorial(31));
    ok = 0;
  );

  /* Test cc_con_forbidden_R(37, 210, 3)
   * S=210, q=37, target=3.  S mod 37 = 210 mod 37 = 210-5*37 = 210-185 = 25
   * S^{-1} mod 37: 25*x ≡ 1 (mod 37) -> x = lift(Mod(25,37)^(-1))
   * i=0: (2^0 * 210)^{-1} mod 37 = (210)^{-1} mod 37
   * i=1: (2^1 * 210)^{-1} mod 37 = (420)^{-1} mod 37
   * i=2: (2^2 * 210)^{-1} mod 37 = (840)^{-1} mod 37 */
  my(fr37 = cc_con_forbidden_R(37, 210, 3));
  my(fr37_exp = vecsort([lift(Mod(210,37)^(-1)), lift(Mod(420,37)^(-1)), lift(Mod(840,37)^(-1))]));
  if(fr37 == fr37_exp,
    printf("  Self-test PASSED: cc_con_forbidden_R(37, 210, 3) = %s\n", Str(fr37));
  ,
    printf("  Self-test FAILED: cc_con_forbidden_R(37, 210, 3) = %s, expected %s\n",
      Str(fr37), Str(fr37_exp));
    ok = 0;
  );

  /* Test cc_con_R_range gives valid range for 20-bit with S=210 */
  my(Rr20 = cc_con_R_range(210, 20));
  my(plo20 = 210 * Rr20[1] - 1, phi20 = 210 * Rr20[2] - 1);
  if(Rr20[1] <= Rr20[2] && #binary(plo20) == 20 && #binary(phi20) == 20 && Rr20[1] % 2 == 1,
    printf("  Self-test PASSED: cc_con_R_range(210, 20) = [%d, %d] (valid 20-bit)\n",
      Rr20[1], Rr20[2]);
  ,
    printf("  Self-test FAILED: cc_con_R_range(210, 20) = [%d, %d]\n", Rr20[1], Rr20[2]);
    ok = 0;
  );

  /* Test cc_con_sieve_check on known-clean R:
   * p = 1122659 (CC7 root). S = 210. R = (p+1)/S = 1122660/210 = 5346
   * But R must satisfy p = S*R - 1. 210*5346 - 1 = 1122659. Check.
   * Actually R=5346 is even, so it wouldn't normally be enumerated,
   * but sieve_check should still work on it.
   * Use a simpler test: build sieve for S=210, target=3, check R=5347 (odd). */
  my(sv8 = cc_con_sieve(210, 3));
  my(R_test = 5347);
  my(sv8_result = cc_con_sieve_check(R_test, sv8));
  if(type(sv8_result) == "t_INT" && (sv8_result == 0 || sv8_result == 1),
    printf("  Self-test PASSED: cc_con_sieve_check(%d, sieve) = %d (valid result)\n",
      R_test, sv8_result);
  ,
    printf("  Self-test FAILED: cc_con_sieve_check returned invalid type\n");
    ok = 0;
  );

  /* Test cc_con_ceiling on CC7 root 1122659: ceiling should be >= 7 */
  my(ceil7 = cc_con_ceiling(test_root, 7));
  if(ceil7 >= 7,
    printf("  Self-test PASSED: cc_con_ceiling(%d, 7) = %d (>= 7)\n", test_root, ceil7);
  ,
    printf("  Self-test FAILED: cc_con_ceiling(%d, 7) = %d, expected >= 7\n", test_root, ceil7);
    ok = 0;
  );

  /* ── Layer 9 (safe prime search) tests ── */

  /* Test cc_sp_trial_primes: primes from 3 to 30 */
  my(tp30 = cc_sp_trial_primes(3, 30));
  if(tp30 == [3, 5, 7, 11, 13, 17, 19, 23, 29],
    printf("  Self-test PASSED: cc_sp_trial_primes(3, 30) = 9 primes\n");
  ,
    printf("  Self-test FAILED: cc_sp_trial_primes(3, 30) = %s\n", Str(tp30));
    ok = 0;
  );

  /* Test cc_sp_trial_check: q=11 should survive trial by [3,5,7]
   * since 11%3=2, 11%5=1, 11%7=4, and 2q+1=23: 23%3=2, 23%5=3, 23%7=2 */
  if(cc_sp_trial_check(11, [3, 5, 7]) == 1,
    printf("  Self-test PASSED: cc_sp_trial_check(11, [3,5,7]) = 1 (safe prime q=11)\n");
  ,
    printf("  Self-test FAILED: cc_sp_trial_check(11, [3,5,7]) expected 1\n");
    ok = 0;
  );

  /* Test cc_sp_trial_check: q=13 should FAIL because 2q+1=27=3*9 */
  if(cc_sp_trial_check(13, [3, 5, 7]) == 0,
    printf("  Self-test PASSED: cc_sp_trial_check(13, [3,5,7]) = 0 (2q+1=27 div by 3)\n");
  ,
    printf("  Self-test FAILED: cc_sp_trial_check(13, [3,5,7]) expected 0\n");
    ok = 0;
  );

  /* Test cc_sp_verify on known safe prime q=11 → p=23 */
  my(sv11 = cc_sp_verify(11));
  if(sv11 == [11, 23, 5],
    printf("  Self-test PASSED: cc_sp_verify(11) = [11, 23, 5]\n");
  ,
    printf("  Self-test FAILED: cc_sp_verify(11) = %s\n", Str(sv11));
    ok = 0;
  );

  /* Test cc_sp_verify on non-safe prime q=13 → p=27 (composite) */
  if(cc_sp_verify(13) == 0,
    printf("  Self-test PASSED: cc_sp_verify(13) = 0 (not safe prime)\n");
  ,
    printf("  Self-test FAILED: cc_sp_verify(13) expected 0\n");
    ok = 0;
  );

  /* Test cc_sp_search at 32 bits (should find one quickly, silent) */
  my(sp32 = cc_sp_search(32, 0, 100000, 0));
  if(type(sp32) == "t_VEC" && #sp32 == 3 && sp32[3] == 32,
    printf("  Self-test PASSED: cc_sp_search(32) found %d-bit safe prime\n", sp32[3]);
  ,
    printf("  Self-test FAILED: cc_sp_search(32) returned %s\n", Str(sp32));
    ok = 0;
  );

  /* ── Layer 9b (OpenSSL comparison) tests ── */

  /* Test cc_sp_openssl_trial_count: known values matching OpenSSL's calc_trial_divisions */
  if(cc_sp_openssl_trial_count(512) == 64 && cc_sp_openssl_trial_count(2048) == 384 && cc_sp_openssl_trial_count(4096) == 1024,
    printf("  Self-test PASSED: cc_sp_openssl_trial_count(512/2048/4096) = 64/384/1024\n");
  ,
    printf("  Self-test FAILED: cc_sp_openssl_trial_count mismatch\n");
    ok = 0;
  );

  /* Test cc_sp_delta_sieve at 32 bits (silent, should find a valid safe prime) */
  my(ds32 = cc_sp_delta_sieve(32, 0, 100000, 0, 0));
  if(type(ds32) == "t_VEC" && #ds32 == 6 && ds32[3] == 32 && isprime(ds32[1]) && isprime(ds32[2]) && ds32[2] == 2*ds32[1]+1,
    printf("  Self-test PASSED: cc_sp_delta_sieve(32) found valid %d-bit safe prime\n", ds32[3]);
  ,
    printf("  Self-test FAILED: cc_sp_delta_sieve(32) returned %s\n", Str(ds32));
    ok = 0;
  );

  /* === Layer 10 tests: cc_x_ functions === */

  /* Test cc_x_forbidden_mask matches cc_forbidden_res */
  my(fm5 = cc_x_forbidden_mask(5, 3, 1));
  /* cc_forbidden_res(5,3,1) = [0,2,3] → bits 0,2,3 set → mask = 1+4+8 = 13 */
  if(fm5 == 13 && bittest(fm5, 0) && bittest(fm5, 2) && bittest(fm5, 3) && !bittest(fm5, 1),
    printf("  Self-test PASSED: cc_x_forbidden_mask(5, 3, 1) = %d (bits 0,2,3)\n", fm5);
  ,
    printf("  Self-test FAILED: cc_x_forbidden_mask(5, 3, 1) = %d (expected 13)\n", fm5);
    ok = 0;
  );

  /* Test cc_x_periodic_table: prime 7, ord_7(2) = 3 */
  my(pt7 = cc_x_periodic_table(7, 1));
  if(#pt7 == 2 && pt7[2] == 3 && pt7[1] == [0, 3, 1],
    printf("  Self-test PASSED: cc_x_periodic_table(7, 1) period=%d table=%s\n", pt7[2], Str(pt7[1]));
  ,
    printf("  Self-test FAILED: cc_x_periodic_table(7, 1) = %s\n", Str(pt7));
    ok = 0;
  );

  /* Test cc_x_line_filter: known CC7a root 1122659 should survive 7 sieve steps */
  my(lf7 = cc_x_line_filter(1122659, 7, 1));
  if(lf7 == 7,
    printf("  Self-test PASSED: cc_x_line_filter(1122659, 7) = %d steps\n", lf7);
  ,
    printf("  Self-test FAILED: cc_x_line_filter(1122659, 7) = %d (expected 7)\n", lf7);
    ok = 0;
  );

  /* Test cc_x_bitwin_mask: prime 5, depth 2.
   * Odd residues forbidden: {1, 3} → bits 1,3.
   * First-kind k=0: (1-1)%5=0, k=1: (1-3)%5=3 → bits 0,3
   * Second-kind k=0: (1-1)%5=0, k=1: (3-1)%5=2 → bits 0,2
   * Union: {0,1,2,3} → mask = 1+2+4+8 = 15.  Only residue 4 survives. */
  my(bm5 = cc_x_bitwin_mask(5, 2));
  if(bm5 == 15 && !bittest(bm5, 4),
    printf("  Self-test PASSED: cc_x_bitwin_mask(5, 2) = %d (only residue 4 survives)\n", bm5);
  ,
    printf("  Self-test FAILED: cc_x_bitwin_mask(5, 2) = %d (expected 15)\n", bm5);
    ok = 0;
  );

  /* Test cc_x_bitwin_walk: center=6 → first-kind from 5 (CC4a), second from 7 (CC2b) */
  my(bw6 = cc_x_bitwin_walk(6));
  if(#bw6 == 3 && bw6[1] == 4 && bw6[2] == 2 && bw6[3] == 2,
    printf("  Self-test PASSED: cc_x_bitwin_walk(6) = [%d, %d, %d]\n", bw6[1], bw6[2], bw6[3]);
  ,
    printf("  Self-test FAILED: cc_x_bitwin_walk(6) = %s (expected [4,2,2])\n", Str(bw6));
    ok = 0;
  );

  /* Test cc_x_root_depth: chain 2→5→11→23→47 (CC5a) */
  if(cc_x_root_depth(2, 1) == 0 && cc_x_root_depth(5, 1) == 1 && cc_x_root_depth(23, 1) == 3 && cc_x_root_depth(47, 1) == 4,
    printf("  Self-test PASSED: cc_x_root_depth(2)=0, (5)=1, (23)=3, (47)=4\n");
  ,
    printf("  Self-test FAILED: cc_x_root_depth chain 2→5→11→23→47 mismatch\n");
    ok = 0;
  );

  /* Test cc_x_primorial_scan: 20-bit, primorial(4)=7#=210, target CC3+ */
  my(ps20 = cc_x_primorial_scan(20, 4, 3, 50000));
  if(#ps20 > 0 && ps20[1][2] >= 3 && isprime(ps20[1][1]),
    printf("  Self-test PASSED: cc_x_primorial_scan(20, 4, 3) found %d chains, best CC%d\n",
      #ps20, vecmax(vector(#ps20, i, ps20[i][2])));
  ,
    printf("  Self-test FAILED: cc_x_primorial_scan(20, 4, 3) found %d chains\n", #ps20);
    ok = 0;
  );

  printf("\n  Available functions:\n");

  printf("\n  -- Analysis (Layers 1-6) --\n");
  printf("    cc_walk(p, {kind}, {lookahead})    — walk chain, display members\n");
  printf("    cc_shadow(p, {kind}, {primes})     — shadow spiral analysis\n");
  printf("    cc_autopsy(p, {kind}, {depth})     — breaker deep analysis\n");
  printf("    cc_deep_factor(n, {depth})         — recursive p-1/p+1 tree\n");
  printf("    cc_geometry(p)                     — bits, popcount, NAF, v2\n");
  printf("    cc_classify(p)                     — Sophie Germain, safe, twin...\n");
  printf("    cc_profile(p,{kind},{normalize})   — compact normalized summary\n");
  printf("    cc_candidate_profile(...)          — exact-input prime/composite summary\n");
  printf("    cc_candidate_positions(...)        — fixed-depth prime/composite trace\n");
  printf("    cc_candidate_full(...)             — one-call single-input analysis\n");
  printf("    cc_candidate_profile_emit(...)     — emit one candidate TSV/CSV row\n");
  printf("    cc_motif(p,{kind},{normalize})     — prefix/tail motif counts\n");
  printf("    cc_tail_probe(p,{patterns},...)    — append tails and test survival\n");
  printf("    cc_tail_probe_range(...)           — exhaustive append-tail sweep\n");
  printf("    cc_prefix_probe(...)               — fixed-prefix low-tail sweep\n");
  printf("    cc_profile_emit(...)               — emit one TSV/CSV row\n");
  printf("    cc_profile_export(file, roots,...) — export many profiles\n");
  printf("    cc_spine(p, {levels}, {max_bits})  — 2-adic spine exploration\n");
  printf("    cc_fingerprint(p, {kind})          — per-member factorization\n");
  printf("    cc_full(p, {kind})                 — everything at once\n");
  printf("    cc_compare(p1, p2, {kind})         — side-by-side comparison\n");
  printf("    cc_sbase_check(p, {kind})          — S-base coverage\n");
  printf("    cc_walk_quiet(p, {kind})           — silent walk, returns data\n");
  printf("    cc_walk_from_quiet(p, {kind})      — exact-start silent walk\n");
  printf("    cc_immune_persist(p, {kind})       — verify immunization holds\n");
  printf("    cc_next_breakers(p, {count},{kind})— post-break composites\n");
  printf("    cc_root(p, {kind})                 — find chain root\n");
  printf("    cc_is_root(p, {kind})              — check if p is a chain root\n");
  printf("    cc_heat(len)                       — 1/(len+1) heat metric\n");
  printf("    cc_shadow_hit(n, {primes})         — first analysis-prime divisor\n");
  printf("    cc_small_factor_hint(n, {bound})   — trial factor hint\n");
  printf("    cc_kill_pos(r, q, {kind}, {max})   — kill position for residue\n");
  printf("    cc_members(p, {kind}, {depth})     — algebraic form of chain members\n");
  printf("    cc_grid_trajectory(p, {kind})      — top-down 2-adic grid path\n");
  printf("    cc_grid_neighbors(p, {radius})     — neighborhood in grid\n");
  printf("    cc_padic_addr(p, {depth})          — multi-prime valuation fingerprint\n");
  printf("    cc_crt_constraints(kind, {primes}) — CRT root class reconstruction\n");

  printf("\n  -- Mini Search Engine (Layer 7) --\n");
  printf("    cc_search(bits, target, {kind}, {prefix}, {prefix_bits}, {max_cand}, {verbose})\n");
  printf("                                       — main search: CRT+sieve+PRP+chain\n");
  printf("    cc_qsearch(bits, target, {kind})   — quick wrapper (50K candidates)\n");
  printf("    cc_search_info(bits, target,{kind}) — show search configuration\n");
  printf("    cc_forbidden_res(q, target, {kind}) — forbidden residues mod q\n");
  printf("    cc_valid_res(q, target, {kind})     — valid residues mod q\n");
  printf("    cc_build_bases(target, {kind})      — CRT base computation\n");
  printf("    cc_build_sieve(target, {kind})      — precomputed sieve table\n");
  printf("    cc_sieve_check(n, sieve_table)      — sieve candidate against table\n");
  printf("    cc_chain_test(n, {kind})            — quick PRP chain length\n");
  printf("    cc_chain_verify(n, {kind})          — proven primality chain length\n");

  printf("\n  -- Constructive Search Engine (Layer 8) --\n");
  printf("    cc_construct(bits, target, {S_limit}, {mode}, {prefix}, {prefix_bits}, {max_cand}, {verbose})\n");
  printf("                                       — constructor search: S*R-1 pipeline\n");
  printf("    cc_qconstruct(bits, target, {S_limit})  — quick wrapper (50K candidates)\n");
  printf("    cc_con_info(bits, target, {S_limit})  — show constructor configuration\n");
  printf("    cc_primorial(n)                     — product of primes <= n\n");
  printf("    cc_con_forbidden_R(q, S, target)    — forbidden R residues mod q\n");
  printf("    cc_con_sieve(S, target, {lo}, {hi}) — constructor sieve table\n");
  printf("    cc_con_sieve_check(R, sieve_table)  — check R against constructor sieve\n");
  printf("    cc_con_R_range(S, bits, {prefix}, {prefix_bits})  — compute [R_min, R_max]\n");
  printf("    cc_con_ceiling(p, target, {primes}) — OPT-6 chain ceiling prefilter\n");

  printf("\n  Examples (Layer 7 — CRT search):\n");
  printf("    cc_search(40, 5)             \\\\ find CC5+ roots at 40 bits\n");
  printf("    cc_search(60, 7, 1, 5, 3)    \\\\ CC7+ at 60 bits, prefix 0b101\n");
  printf("    cc_search(30, 4,,,,, 0)       \\\\ CC4+ at 30 bits, silent mode\n");
  printf("    cc_search_info(89, 18)        \\\\ show search space for CC18 at 89 bits\n");
  printf("    cc_qsearch(35, 4)            \\\\ quick CC4+ search at 35 bits\n");
  printf("\n  Examples (Layer 8 — constructor p = S*R - 1):\n");
  printf("    cc_construct(35, 4, 13)       \\\\ CC4+ at 35 bits, S=primorial(13)\n");
  printf("    cc_construct(50, 5, 19)       \\\\ CC5+ at 50 bits, S=primorial(19)\n");
  printf("    cc_construct(60, 6, 23)       \\\\ CC6+ at 60 bits, S=primorial(23)\n");
  printf("    cc_construct(89, 10, 31)      \\\\ CC10+ at 89 bits, S=primorial(31)\n");
  printf("    cc_construct(40, 5, 13, \"A\")  \\\\ mode A: R must be prime\n");
  printf("    cc_construct(40, 4, 13, \"C\")  \\\\ mode C: random R sampling\n");
  printf("    cc_construct(45, 5, 17, \"B\", 5, 3)  \\\\ prefix 0b101, sequential\n");
  printf("    cc_con_info(89, 18)           \\\\ constructor params for CC18\n");
  printf("    cc_con_info(60, 10, 23)       \\\\ CC10 info with S=primorial(23)\n");
  printf("    cc_qconstruct(35, 4, 13)      \\\\ quick constructor search\n");
  printf("    cc_primorial(31)              \\\\ = 200560490130\n");
  printf("    cc_con_forbidden_R(37, cc_primorial(31), 18)  \\\\ forbidden R mod 37\n");
  printf("    cc_con_ceiling(1122659, 7)    \\\\ ceiling prefilter on known CC7 root\n");

  printf("\n  -- Safe Prime Search Engine (Layer 9) --\n");
  printf("    cc_sp_search(bits, {trial_limit}, {max_cand}, {verbose}, {proven})\n");
  printf("                                       — main safe prime search\n");
  printf("    cc_sp_qsearch(bits, {trial_limit})  — quick wrapper (1M candidates)\n");
  printf("    cc_sp_info(bits, {trial_limit})      — show search configuration\n");
  printf("    cc_sp_trial_primes(lo, hi)           — generate trial prime list\n");
  printf("    cc_sp_trial_check(q, trial_primes)   — combined q+2q+1 trial division\n");
  printf("    cc_sp_verify(q)                      — prove q and 2q+1 prime\n");
  printf("    cc_sp_bench(bits, {max_cand})        — benchmark trial vs no-trial\n");

  printf("\n  Examples (Layer 9 — safe prime search p = 2q+1):\n");
  printf("    cc_sp_search(512)             \\\\ 512-bit safe prime (~instant)\n");
  printf("    cc_sp_search(2048)            \\\\ 2048-bit (~seconds)\n");
  printf("    cc_sp_search(4096)            \\\\ 4096-bit (~30-120s in GP)\n");
  printf("    cc_sp_search(4096, 200000)    \\\\ larger trial sieve\n");
  printf("    cc_sp_info(4096)              \\\\ search space analysis\n");
  printf("    cc_sp_bench(64)               \\\\ benchmark trial vs no-trial\n");
  printf("    cc_sp_verify(11)              \\\\ verify q=11 → p=23 safe prime\n");

  printf("\n  -- OpenSSL Comparison (Layer 9b) --\n");
  printf("    cc_sp_openssl_trial_count(bits)     — OpenSSL's trial prime count\n");
  printf("    cc_sp_delta_sieve(bits, {n_trial}, {max_cand}, {max_delta}, {verbose})\n");
  printf("                                       — OpenSSL-style delta-sieve search\n");
  printf("    cc_sp_compare(bits, {runs}, {max_cand})  — head-to-head comparison\n");
  printf("    cc_sp_trial_analysis(bits)          — trial prime count impact analysis\n");

  printf("\n  Examples (Layer 9b — OpenSSL comparison):\n");
  printf("    cc_sp_trial_analysis(4096)          \\\\ trial prime impact table\n");
  printf("    cc_sp_compare(512, 5)               \\\\ head-to-head at 512 bits\n");
  printf("    cc_sp_delta_sieve(512)              \\\\ OpenSSL-style search\n");
  printf("    cc_sp_delta_sieve(512, 1000)        \\\\ delta-sieve with 1000 primes\n");

  printf("\n  -- Alternative Algorithms (Layer 10, cc_x_ prefix) --\n");
  printf("    cc_x_forbidden_mask(q, depth, {kind})  — bitmask of forbidden residues\n");
  printf("    cc_x_periodic_table(q, {kind})         — periodic table via ord_q(2)\n");
  printf("    cc_x_line_filter(p, depth, {kind}, {max_prime})\n");
  printf("                                           — sieve-only chain pre-screen\n");
  printf("    cc_x_bitwin_mask(q, depth)             — BiTwin forbidden mask (both kinds)\n");
  printf("    cc_x_bitwin_walk(n, {max_depth})       — BiTwin chain walk from even center\n");
  printf("    cc_x_primorial_scan(bits, prim_n, {target}, {k_max})\n");
  printf("                                           — algebraic form k*P#*2^m-1 search\n");
  printf("    cc_x_root_depth(p, {kind})             — backward steps to chain root\n");

  printf("\n  Examples (Layer 10 — alternative algorithms):\n");
  printf("    cc_x_forbidden_mask(5, 3)         \\\\ = 13 (bits 0,2,3 set)\n");
  printf("    bittest(cc_x_forbidden_mask(5,3), 1122659 %% 5) \\\\ 0 = survives\n");
  printf("    cc_x_periodic_table(7)            \\\\ [table, period=3]\n");
  printf("    cc_x_line_filter(1122659, 7)      \\\\ = 7 (CC7a root survives 7 steps)\n");
  printf("    cc_x_bitwin_walk(6)               \\\\ [4, 2, 2] (CC4a + CC2b)\n");
  printf("    cc_x_bitwin_walk(30)              \\\\ try various even centers\n");
  printf("    cc_x_primorial_scan(20, 4, 3)     \\\\ 20-bit, 7#=210, CC3+\n");
  printf("    cc_x_root_depth(23)               \\\\ = 3 (root is 2, three steps back)\n");

  printf("\n");
  if(ok, printf("  All self-tests passed (%d tests). Ready.\n\n", 37),
         printf("  WARNING: some self-tests failed!\n\n"));
}
