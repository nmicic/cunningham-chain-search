#!/usr/bin/env bash
# Fetch the CC10+ root snapshot CSV used by the fingerprint pool builders
# and the external-prefix validator.
#
# Run from the repo root:
#     bash data/fetch_data.sh
#
# Source: https://github.com/nmicic/cunningham-chain-data
# Release: "Second aggregate first-kind snapshot (2026-03-19)"
set -euo pipefail

URL="https://github.com/nmicic/cunningham-chain-data/releases/download/v2026-03-19-snapshot/cc10plus_roots_snapshot_2026-03-19.csv.gz"
GZ="data/cc10plus_roots_snapshot_2026-03-19.csv.gz"
CSV="data/cc10plus_roots_snapshot_2026-03-19.csv"

if [ -f "$CSV" ]; then
    echo "Already present: $CSV"
    exit 0
fi

mkdir -p data

if command -v wget >/dev/null 2>&1; then
    wget -O "$GZ" "$URL"
elif command -v curl >/dev/null 2>&1; then
    curl -L -o "$GZ" "$URL"
else
    echo "ERROR: neither wget nor curl is installed" >&2
    exit 1
fi

gunzip "$GZ"
echo "Fetched: $CSV"
ls -lh "$CSV"
