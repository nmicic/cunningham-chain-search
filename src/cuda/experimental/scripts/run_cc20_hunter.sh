#!/usr/bin/env bash
# Copyright (c) 2026 Nenad Mićić <nenad@micic.be>
# SPDX-License-Identifier: Apache-2.0
#
# =============================================================================
# CC20/CC21-scale search preset for cc20_first_kind_immune_v16
#
# Strategy: q-iter against the exponent-aware 4248-fingerprint pool,
# bit-band [90, 107] hard-enforced (target-21 top < 2^127, u128-safe).
# CC10+ live HIT to stderr + hit.log. Continuous loop; --time 0 = run forever.
#
# Run from the repo root after `make -j`:
#     bash scripts/run_cc20_hunter.sh
#
# Edit the knobs below to match your hardware.
# =============================================================================

set -uo pipefail

CMD="${CMD:-./bin/cc20_first_kind_immune_v16}"

# --- Tunable knobs (edit here) -----------------------------------------------
POOL_FILE="${POOL_FILE:-pools/known_fingerprints_exp_d.txt}"
SEED_NAME="${SEED_NAME:-cc20_hunter}"
TARGET="${TARGET:-21}"              # chain length we report on
DEPTH="${DEPTH:-21}"                # GPU sieve depth
SIEVE_MAX="${SIEVE_MAX:-503}"       # sieve prime ceiling
# Auto-tune prove-threads: nproc - 4 (keep cores for K1/host), capped at 24
N_CPU="$(nproc)"
DEFAULT_PT=$(( N_CPU > 4 ? N_CPU - 4 : 1 ))
[ "$DEFAULT_PT" -gt 24 ] && DEFAULT_PT=24
PROVE_THREADS="${PROVE_THREADS:-$DEFAULT_PT}"
MIN_REPORT_LEN="${MIN_REPORT_LEN:-10}"
REPORT_INTERVAL="${REPORT_INTERVAL:-30}"
IMMUNE_PRIME="${IMMUNE_PRIME:-5}"   # smallest q forced into immunity
EXP_START="${EXP_START:-0}"
BITS_MIN="${BITS_MIN:-90}"          # hard bit-band floor
BITS_MAX="${BITS_MAX:-107}"         # hard bit-band ceiling (target-21 top < 2^127)
GPU="${GPU:-0}"
TIME_BUDGET="${TIME_BUDGET:-0}"     # 0 = run forever

# --- Output paths ------------------------------------------------------------
RUN_TS="$(date -u +%Y%m%dT%H%M%SZ)"
OUTDIR="run/${SEED_NAME}_${RUN_TS}"
mkdir -p "$OUTDIR"
HIT_LOG="${OUTDIR}/hit.log"
STDOUT_LOG="${OUTDIR}/stdout.log"
STDERR_LOG="${OUTDIR}/stderr.log"
CHECKPOINT="${OUTDIR}/checkpoint.txt"

# --- Verify binary -----------------------------------------------------------
if [ ! -x "$CMD" ]; then
    echo "ERROR: $CMD not found or not executable" >&2
    echo "Build with: make -j" >&2
    exit 1
fi
[ -f "$POOL_FILE" ] || { echo "ERROR: pool file '$POOL_FILE' missing — run from repo root" >&2; exit 1; }

# --- Signal handling ---------------------------------------------------------
trap "echo 'Caught signal, exiting...'; exit 0" INT TERM QUIT

echo "=== Cunningham-chain hunter: target=$TARGET depth=$DEPTH sieve=$SIEVE_MAX pt=$PROVE_THREADS ==="
echo "Pool:    $POOL_FILE"
echo "Bits:    [$BITS_MIN, $BITS_MAX]"
echo "Outdir:  $OUTDIR"
echo "Tail:    tail -F $STDOUT_LOG"
echo ""

# Notes on Q sampling:
#   --q-order random          : exp_d's Q-ranges are astronomical, random
#                               is the only practical sampling mode.
#   --q-band-mode fixed       : default; honours --q-order globally.
#                               For mixed-range pools (e.g. primorial),
#                               flip to `exhaustive` — see examples/run_primorial.sh.
#   --dedup-p0                : default ON; collapses repeat-p_0 emits
#                               to compact `HIT-DUP` lines in hit.log.
#                               Pass --no-dedup-p0 if you want every emit verbose.
stdbuf -oL "$CMD" \
    --q-iter \
    --seed-name "$SEED_NAME" \
    --seed-pool-file "$POOL_FILE" \
    --immune-prime "$IMMUNE_PRIME" \
    --exp-start "$EXP_START" \
    --target "$TARGET" \
    --depth "$DEPTH" \
    --sieve-max "$SIEVE_MAX" \
    --bits-min "$BITS_MIN" \
    --bits-max "$BITS_MAX" \
    --q-order random \
    --q-band-mode fixed \
    --dedup-p0 \
    --prove-queue-depth 16 \
    --streams 4 \
    --cpu-prove --prove-threads "$PROVE_THREADS" \
    --min-report-len "$MIN_REPORT_LEN" \
    --report "$REPORT_INTERVAL" \
    --time "$TIME_BUDGET" \
    --gpu "$GPU" \
    --log "$HIT_LOG" \
    --checkpoint "$CHECKPOINT" \
    2> "$STDERR_LOG" | tee "$STDOUT_LOG"
