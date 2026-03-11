\\ Copyright (c) 2026 Nenad Mićić <nenad@micic.be>
\\ SPDX-License-Identifier: Apache-2.0
\\
\\ validate_chains.gp — Validate Cunningham Chain (first kind) files
\\ Usage: gp -q -s 512m validate_chains.gp

{
dir = "./";  \\ directory containing CC_8.txt .. CC_16.txt
files = [8,9,10,11,12,13,14,15,16];

total_ok = 0;
total_fail = 0;
total_longer = 0;

parse_hex(s) =
{
  my(v = 0);
  for(i = 1, #s,
    my(c = Vecsmall(s)[i], d = 0);
    if(c >= 48 && c <= 57, d = c - 48,
       c >= 97 && c <= 102, d = c - 87,
       c >= 65 && c <= 70, d = c - 55,
       next);
    v = v * 16 + d;
  );
  v;
};

cc_len(p) =
{
  my(q = p, n = 0);
  while(ispseudoprime(q), n++; q = 2*q + 1);
  n;
};

for(fi = 1, #files,
  L = files[fi];
  fname = Str(dir, "all_CC", L, ".txt");
  printf("=== Validating CC%d ===\n", L);

  ok = 0; fail = 0; longer = 0; count = 0;

  fp = fileopen(fname, "r");
  while(1,
    s = iferr(filereadstr(fp), E, "");
    if(s == "" || type(s) != "t_STR", break);

    \\ strip 0x prefix
    if(#s >= 2 && (s[1] == "0") && (s[2] == "x" || s[2] == "X"),
      s = s[3..#s]
    );
    if(#s == 0, next);

    p = parse_hex(s);
    if(p <= 1, next);
    count++;

    clen = cc_len(p);

    if(clen >= L,
      ok++;
      if(clen > L, longer++;
        printf("  #%d: chain=%d (>%d!) p=0x%s\n", count, clen, L, s);
      );
    ,
      fail++;
      printf("  FAIL #%d: chain=%d (expected >=%d) p=0x%s\n", count, clen, L, s);
    );

    if(count % 10000 == 0, printf("  ... %d checked so far\n", count));
  );
  fileclose(fp);

  printf("  Results: %d checked, %d OK, %d LONGER, %d FAIL\n\n", count, ok, longer, fail);
  total_ok += ok;
  total_fail += fail;
  total_longer += longer;
);

printf("=== TOTAL: %d OK, %d longer-than-expected, %d FAIL ===\n", total_ok, total_longer, total_fail);
if(total_fail == 0, printf("ALL CHAINS VALID!\n"), printf("*** %d FAILURES DETECTED ***\n", total_fail));
}
quit
