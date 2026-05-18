# Example launchers — one-command CC21 hunt

Five self-contained scripts to start a v16 search on this machine. Each one:

1. Picks a different fingerprint pool (the only knob worth varying for a first run).
2. Starts the engine in a detached `tmux` session so it survives logout / SSH drops.
3. Writes logs to `run/cc_v16_<pool>/run_<UTC-ts>_b107/`.
4. Auto-tunes `--prove-threads` for your CPU (≈ nproc − 4, capped at 24).

| Script | Pool file | Entries | When to use |
|---|---|---|---|
| `run10.sh` | `pools/known_fingerprints_top10.txt` | 10 | Highest per-GPU-second hit rate — narrowest structural search. **Good first run.** |
| `run100.sh` | `pools/known_fingerprints.txt` | 100 | Balanced: more shapes, slightly slower per-GPU-second hit rate. |
| `run_all.sh` | `pools/known_fingerprints_all.txt` | 1,372 | Deduplicated historical `D_V` values, including singleton-frequency shapes the `exp_d` pool drops. |
| `run_exp.sh` | `pools/known_fingerprints_exp_d.txt` | 4248 | Largest historical-frequency pool, ~96.78 % coverage of historical CC10+. |
| `run_primorial.sh` | `pools/known_fingerprints_primorial.txt` | 3000 | **Structural** pool: every `X#/{Y_1,..,Y_k}` (X up to 71#, drop 0–3 small primes). Captures exotic primorial-quotient patterns like `41#/31` and `71#/31` that historical pools miss. |

## Quick start

```bash
# from this folder's parent, after `make -j`
bash examples/run10.sh
```

Then either:

```bash
tmux attach -t cc_v16_top10                          # watch it live
tail -F run/cc_v16_top10/run_*_b107/stdout.log       # tail the periodic report
tail -F run/cc_v16_top10/run_*_b107/hit.log          # tail CC10+ hits
```

## Stop or restart

```bash
tmux kill-session -t cc_v16_top10        # stop
bash examples/run10.sh                   # restart (refuses if session is alive)
```

## What's hard-coded

These scripts target **CC21-scale searches** with conservative bit safety:

- `--target 21 --depth 21`
- `--bits-min 90 --bits-max 107` → chain top `p_20 ≤ 2^127` (tight u128 edge, validated at startup)
- `--q-order random` → two GPUs on the same pool diverge structurally
- `--immune-prime 5` → keep `D_V` variants that exclude 7 (76 % of CC10+ shapes)
- `--sieve-max 503`, `--prove-queue-depth 16`, `--streams 4` → empirical optima

Want a different bit band, target length, or hardware tuning? Copy the script
and edit. See `../CONFIGURATION.md` for the full CLI reference.
