#!/bin/bash
# Copyright (c) 2026 Nenad Mićić <nenad@micic.be>
# SPDX-License-Identifier: Apache-2.0
#
# CC17 Second-Kind Sequential Search — cc18_filter_cuda_CpC_v13b
#
# Loops through bit sizes 72-76, sequential prefix sweep.
# Usage: ./run_cc17b_sequential.sh <gpu_id> <prefix>
# Example: ./run_cc17b_sequential.sh 0 0b10

CMD=./cc18_filterb
NAME=$(basename "$CMD")
GPU=$1
PREFIX=$2

trap "exit" TERM
trap "exit" QUIT
trap "exit" INT

while [ 1 ]; do
  for j in 72 73 74 75 76; do
    outdir="${NAME}_$(date +%Y)/$(date +%m)/$(date +%d)/$(date +%H)"
    mkdir -p "$outdir"

    ts="$(date +%s)"
    report="${outdir}/gpu${GPU}_report_bits${j}_${ts}.txt"
    teeout="${outdir}/gpu${GPU}_tee_bits${j}_${ts}.txt"

    stdbuf -oL -eL ./cc18_filterb \
      --target 17 --bits "$j" --depth 17 \
      --log 12 --report 5 --output "$report" \
      --prefix-mode sequential --prefix "$PREFIX" \
      --batch 524288 --blocks-per-sm 16 \
      --gpu "$GPU" --prove-threads 60 \
      |& tee -ai "$teeout"

    sleep 0.1
  done
done
