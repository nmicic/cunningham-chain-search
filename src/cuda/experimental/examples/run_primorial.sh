#!/usr/bin/env bash
# Example one-shot launcher — PRIMORIAL-QUOTIENT fingerprint pool (3000 entries).
#
# Run from the repo root after `make -j`:
#     bash examples/run_primorial.sh
#
# Re-attach to the running search:
#     tmux attach -t cc_v16_primorial
#
# Stop the search:
#     tmux kill-session -t cc_v16_primorial
#
# Structural pool: every D of the form X#/{Y_1,...,Y_k} for X ∈ {41,43,47,53,
# 59,61,67,71}, k ∈ {0,1,2,3}, with {2,3,5} always kept. Captures exotic
# primorial-quotient patterns (e.g. 41#/31, 71#/31) that the historical-
# frequency pools (exp_d, top10) filter out as singleton-frequency shapes.
#
# This script runs in `--q-band-mode exhaustive`: at startup the engine
# classifies each V by the bit-width of its Q-range. Variants whose range
# fits in ≤ 40 bits get sequential deterministic coverage (exhausts in
# seconds–hours, no birthday collisions); the rest stay on random
# sampling. Banner reports the split, e.g.
#   q-band-mode: exhaustive  threshold=2^40 Q-range  (sequential V's=1250  random V's=1750)
# `--dedup-p0` (default on) collapses repeat-p_0 emits (from random
# resampling **and** multi-V framings of the same chain) to HIT-DUP
# lines.
set -uo pipefail

POOL_NAME="primorial"
POOL_FILE="pools/known_fingerprints_primorial.txt"
SESSION="cc_v16_${POOL_NAME}"
BIN="${BIN:-./bin/cc20_first_kind_immune_v16}"

if tmux has-session -t "$SESSION" 2>/dev/null; then
    echo "ERROR: tmux session '$SESSION' already exists." >&2
    echo "  Attach:  tmux attach -t $SESSION" >&2
    echo "  Kill:    tmux kill-session -t $SESSION" >&2
    exit 1
fi

[ -x "$BIN" ] || { echo "ERROR: $BIN missing — build with 'make -j'" >&2; exit 1; }
[ -f "$POOL_FILE" ] || { echo "ERROR: pool file '$POOL_FILE' missing — run from repo root" >&2; exit 1; }

# Auto-tune prove-threads: keep ~4 cores for K1/host work
N_CPU="$(nproc)"
PROVE_THREADS=$(( N_CPU > 4 ? N_CPU - 4 : 1 ))
[ "$PROVE_THREADS" -gt 24 ] && PROVE_THREADS=24

TS="$(date -u +%Y%m%dT%H%M%SZ)"
OUTDIR="run/cc_v16_${POOL_NAME}/run_${TS}_b107"
mkdir -p "$OUTDIR"

tmux new-session -d -s "$SESSION" "stdbuf -oL $BIN --q-iter --seed-name cc_v16_${POOL_NAME} --seed-pool-file ${POOL_FILE} --immune-prime 5 --exp-start 0 --target 21 --depth 21 --sieve-max 503 --bits-min 90 --bits-max 107 --cpu-prove --prove-threads ${PROVE_THREADS} --prove-queue-depth 16 --streams 4 --q-band-mode exhaustive --exhaustive-max-q-bits 40 --dedup-p0 --min-report-len 10 --time 0 --report 30 --gpu 0 --log ${OUTDIR}/hit.log --checkpoint ${OUTDIR}/checkpoint.txt 2> ${OUTDIR}/stderr.log | tee ${OUTDIR}/stdout.log"

echo "Started tmux session: $SESSION  (pool=$POOL_NAME, pt=$PROVE_THREADS, band=[90,107], target=21)"
echo "Outdir:  $OUTDIR"
echo "Attach:  tmux attach -t $SESSION"
echo "Tail:    tail -F $OUTDIR/stdout.log"
