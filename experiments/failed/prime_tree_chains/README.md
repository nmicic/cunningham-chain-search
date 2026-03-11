# prime_tree_chains index

**Latest file**: `00_prime_tree_chains__20260117-181212.cu`

## Functional spec
Prime Tree Chains traces yellow chains along tree edges with transitions of the form p → 2p or 2p+1.

## Problem / capability
Describe sequences that climb tree edges while staying prime, providing compile/run instructions for exploration.

## Key methods / implementation notes
- Right-edge transitions p → 2p+1 are Cunningham chains of the second kind (CC2)
- Left edge p → 2p is modeled but usually dies immediately for odd primes
- Examples given for lengths
- Binary tree view of prime progressions

## Build / run
- **Compile**: `nvcc -O3 -arch=sm_86 prime_tree_chains.cu -o prime_tree_chains`
- **Usage tip**: `./prime_tree_chains --exp 30 --min-chain 5` to scan for chains.

## Notes
CC2 behavior appears through right-edge chains, so this is Cunningham-related visualization tooling.

## Cunningham relation
- **Status**: CC-related visualization (contains CC2 paths)
- **Why**: Tracks both `2p` and `2p+1` edges; right-edge chains are CC2.
