# `data/` — CC10+ root snapshot

This directory holds the **CC10+ root snapshot CSV** that the fingerprint
pool builders read. It is not bundled with the source — fetch it from
the dedicated data repo before running any builder or validator.

## Fetch

```bash
# from the repo root
bash data/fetch_data.sh
```

The script fetches `cc10plus_roots_snapshot_2026-03-19.csv.gz` (~5 MB
gzipped, ~52 MB uncompressed) from the release page and gunzips it in
place. Re-running is a no-op if the CSV is already present.

## Source

[`github.com/nmicic/cunningham-chain-data`](https://github.com/nmicic/cunningham-chain-data)
— release **"Second aggregate first-kind snapshot (2026-03-19)"**.

Mirror the data repo if you need to host it elsewhere; the engine itself
doesn't care about the URL, only the file at `data/cc10plus_roots_snapshot_2026-03-19.csv`.

## What this file contains

One row per known Cunningham chain of the first kind with length ≥ 10,
columns:

| column | meaning |
|---|---|
| `cc` | chain length |
| `root_hex` | `p_0` as hex |
| `digits` | decimal digit count of `p_0` |
| `bits` | bit length of `p_0` |

Aggregated from the historical record corpus through 2026-03-19.

## Scripts that read this file

| script | what it does |
|---|---|
| `scripts/build_known_fingerprints.py` | Builds `pools/known_fingerprints.txt` (top-100) |
| `scripts/build_known_fingerprints_all.py` | Builds `pools/known_fingerprints_all.txt` (all squarefree shapes) |
| `scripts/build_known_fingerprints_exp.py` | Builds `pools/known_fingerprints_exp_d.txt` (exponent-aware, 4248 variants) |
| `scripts/verify_fingerprints.py` | Re-derives the pool from scratch and cross-checks distribution |
| `tools/snapshot_replay.py` | External-prefix structural validator — checks pool covers historical chains |

All of them default to reading `./data/cc10plus_roots_snapshot_2026-03-19.csv`.
Override with `--csv PATH` if you keep the data elsewhere.
