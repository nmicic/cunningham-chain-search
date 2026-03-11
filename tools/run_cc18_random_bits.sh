#!/bin/bash
# Copyright (c) 2026 Nenad Mićić <nenad@micic.be>
# SPDX-License-Identifier: Apache-2.0
#
# CC18 First-Kind Random Bit-Size Search — cc18_filter_cuda_CpC_v13
#
# Picks random bit sizes from 86-107, runs 5-minute sessions continuously.
# Usage: ./run_cc18_random_bits.sh <gpu_id>
# Example: ./run_cc18_random_bits.sh 0

CMD=./cc18_filter
NAME=$(basename "$CMD")
GPU=$1

trap "exit" TERM
trap "exit" QUIT
trap "exit" INT

while [ 1 ]; do
  j=$((RANDOM % 22 + 86))
#j=88
  outdir="${NAME}_$(date +%Y)/$(date +%m)/$(date +%d)/$(date +%H)"
  mkdir -p "$outdir"

  ts="$(date +%s)"
  report="${outdir}/gpu${GPU}_report_bits${j}_${ts}.txt"
  teeout="${outdir}/gpu${GPU}_tee_bits${j}_${ts}.txt"

  timeout --foreground -s 15 -k 5 305 \
    stdbuf -oL -eL "$CMD" \
      --time 300 --target 18 --bits "$j" --depth 18 \
      --log 10 --report 5 --output "$report" \
      --continuous --batch 524288 --blocks-per-sm 16 \
      --gpu "$GPU" --prove-threads 32 \
      |& tee -ai "$teeout"

  sleep 0.1
done
