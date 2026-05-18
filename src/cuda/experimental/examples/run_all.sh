#!/usr/bin/env bash
# Example one-shot launcher — FULL fingerprint pool (1,372 entries, deduplicated D_V values).
#
# Run from the repo root after `make -j`:
#     bash examples/run_all.sh
#
# Re-attach to the running search:
#     tmux attach -t cc_v16_all
#
# Stop the search:
#     tmux kill-session -t cc_v16_all
#
# Broadest D_V coverage (1,372 shapes). Better diversity but spreads GPU time
# thinner than top10/top100; use when the small pools have been worked hard.
set -uo pipefail

POOL_NAME="all"
POOL_FILE="pools/known_fingerprints_all.txt"
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

N_CPU="$(nproc)"
PROVE_THREADS=$(( N_CPU > 4 ? N_CPU - 4 : 1 ))
[ "$PROVE_THREADS" -gt 24 ] && PROVE_THREADS=24

TS="$(date -u +%Y%m%dT%H%M%SZ)"
OUTDIR="run/cc_v16_${POOL_NAME}/run_${TS}_b107"
mkdir -p "$OUTDIR"

tmux new-session -d -s "$SESSION" "stdbuf -oL $BIN --q-iter --seed-name cc_v16_${POOL_NAME} --seed-pool-file ${POOL_FILE} --immune-prime 5 --exp-start 0 --target 21 --depth 21 --sieve-max 503 --bits-min 90 --bits-max 107 --cpu-prove --prove-threads ${PROVE_THREADS} --prove-queue-depth 16 --streams 4 --q-order random --dedup-p0 --min-report-len 10 --time 0 --report 30 --gpu 0 --log ${OUTDIR}/hit.log --checkpoint ${OUTDIR}/checkpoint.txt 2> ${OUTDIR}/stderr.log | tee ${OUTDIR}/stdout.log"

echo "Started tmux session: $SESSION  (pool=$POOL_NAME, pt=$PROVE_THREADS, band=[90,107])"
echo "Outdir:  $OUTDIR"
echo "Attach:  tmux attach -t $SESSION"
echo "Tail:    tail -F $OUTDIR/stdout.log"
