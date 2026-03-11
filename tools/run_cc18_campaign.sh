#!/bin/bash
# Copyright (c) 2026 Nenad Mićić <nenad@micic.be>
# SPDX-License-Identifier: Apache-2.0
#
# =============================================================================
# CC18a Campaign Runner (originally for cc_gmp_v31_claude, now cc_gmp_v33_03)
#
# Strategy: Random-chunk sweep across bit sizes 89-99.
# Each bit size runs for 15 minutes, then rotates.
# Full cycle: 11 bit sizes × 15 min = ~2.75 hours, then repeats.
#
# Hardware: AMD Ryzen 9 7950X3D + RTX 5090 (CPU-only for now)
# =============================================================================

CMD=./cc_v31_claude
NAME=$(basename "$CMD")

# --- Tunable parameters ---
BITS_MIN=89
BITS_MAX=99
RUN_SECONDS=900          # 15 minutes per bit size
THREADS=30               # 30 of 32 threads (leaves 2 for OS)
PIN_BASE=0               # Pin starting from core 0
LOG_THRESHOLD=7           # Log CC7+ chains
REPORT_INTERVAL=5         # Progress every 5s
LINE_DEPTH=17             # Full CC18 line-sieve
CHUNK_TILES=500           # Tiles per random chunk (default)

# --- CCD0-only mode (uncomment to test V-Cache-only performance) ---
# THREADS=16
# PIN_BASE=0
# Note: 7950X3D CCD0 has 96 MB V-Cache. Test if 16 threads on CCD0
# beats 30 threads across both CCDs.

# --- Signal handling ---
trap "echo 'Caught signal, exiting...'; exit 0" INT TERM QUIT

# --- Verify binary exists ---
if [ ! -x "$CMD" ]; then
    echo "ERROR: $CMD not found or not executable"
    echo "Build with: gcc -O3 -march=native -flto -o cc_v31_claude src/cpu/cc_gmp_v33_03.c -lgmp -lpthread -lm"
    exit 1
fi

echo "=== CC18a Campaign: ${BITS_MIN}-${BITS_MAX} bit, ${THREADS} threads, ${RUN_SECONDS}s per size ==="
echo "Binary: $CMD"
echo ""

while true; do
    for j in $(seq "$BITS_MIN" "$BITS_MAX"); do

        outdir="${NAME}_$(date +%Y)/$(date +%m)/$(date +%d)/$(date +%H)"
        mkdir -p "$outdir"

        ts="$(date +%s)"
        report="${outdir}/cc18a_${j}bit_${ts}.txt"
        teeout="${outdir}/cc18a_${j}bit_console_${ts}.txt"

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting ${j}-bit run (${RUN_SECONDS}s)..."

        timeout --foreground -s 15 -k 5 "$RUN_SECONDS" \
            stdbuf -oL -eL "$CMD" \
                --target 18 \
                --bits "$j" \
                --threads "$THREADS" \
                --pin --pin-base "$PIN_BASE" \
                --log "$LOG_THRESHOLD" \
                --line-depth "$LINE_DEPTH" \
                --chunk-tiles "$CHUNK_TILES" \
                --report "$REPORT_INTERVAL" \
                --continuous \
                --output "$report" \
            |& tee -ai "$teeout"

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] ${j}-bit run complete."
        sleep 0.5
    done

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] === Full cycle complete, restarting ==="
done
