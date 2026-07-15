# Running the Yannakakis executor at larger scales (server notes)

Notes for scaling `exec.cpp` (SMS join tree → COUNT(*) execution) beyond the
32 GB laptop runs recorded in README.md. Numbers below extrapolate from the
Apple M-series measurements; expect ~1.5–2× more ns/value per core on typical
x86 servers, and re-measure before trusting projections.

## Build & smoke test

The repo's original target is Ubuntu + GCC 14 (the two Apple-clang
portability fixes in `src/` are harmless there). No cmake needed:

```
g++ -std=c++2c -O3 -march=native -DNDEBUG -w compare/exec.cpp \
    src/QueryGraph.cpp src/UnionFind.cpp -o compare/exec
compare/exec verify        # ~seconds; must print 128 OK, no FAIL
```

## Invocation

```
compare/exec <shape> <n> <rows>           one point, exactly n tables
compare/exec <shape> 0 <rows> [maxN]      geometric sweep 1000..maxN
EXEC_MEM_GIB=200 compare/exec ...         raise the database size guard
```

- `rows` = rows per table (0 = tiny default ~5.5). Values are uint8, domain
  D=4, so a table is `rows × arity` bytes.
- Default sweep cap is the 1 s planning budget size (chain/clique/star 16.2M,
  tree 4.4M); pass `maxN` to go beyond — nothing else limits n until memory.
- The guard refuses databases > `EXEC_MEM_GIB` (default 20) GiB of *values*;
  set it to ~50–60% of server RAM (values are only part of the peak, see
  below).
- CSV goes to stdout (`type,n,rows,values,gen_s,plan_s,exec_s,ns_per_value,
  count,overflow`), flushed per point — safe to `nohup ... > out.csv &` and
  tail. Each point regenerates data from a seed derived from n; runs are
  reproducible and independent.

## Memory model (peak, bytes)

Let n = tables, R = rows/table, I = incidences (chain/star/tree: 2(n−1);
clique: n). D = 4.

| what | size | notes |
|---|---|---|
| values | I·(R+1.5) | uint8; the dominant term for R ≳ 100 |
| message arena | 8·D·n = 32n | uint64 count per (table, key value) |
| hypergraph + SMS + exec CSRs | ~40–60n | evtx, vedge, gamma, coff/clist, sattr, ... |
| query graph construction | ~100–150n *transient* | `genGraph`'s `vector<vector<unsigned>>` + QueryGraph, freed before execution — but sets the peak at small R |

Rule of thumb peak: **chain/star/tree ≈ 2nR + ~150n; clique ≈ nR + ~150n.**

Worked examples:
- chain n=100M, R=100 → ~20 GB values + ~15 GB overhead ≈ 35 GB (64 GB box)
- chain n=500M, R=100 → ~100 GB + ~75 GB ≈ 175 GB (256 GB box)
- chain n=16.2M, R=5000 → ~160 GB + ~2.5 GB (row-heavy, overhead negligible)
- clique n=1B, R=100 → ~100 GB + ~150 GB ≈ 250 GB — the ~150n overhead
  dominates; do the direct-construction fix below first

## Time model (single-threaded, M-series; scale up for x86)

- datagen: ~0.5–1 ns/value (power-of-two D fast path: 8 values per rng draw)
- SMS plan: ~17 ns/table (chain/star/clique); random trees ~190 ns/table
  (cache-hostile MCS ordering)
- execute: ~1.0 ns/value chain, ~1.65 star, ~3.4 clique; tree ~2.7 and slowly
  growing with n (cache misses on parent/child access)

E.g. chain n=500M, R=100 (100 GB of values): gen ~100 s, plan ~9 s,
exec ~200 s per trial. Everything is one core; run the four shapes as
separate processes in parallel if RAM allows (they don't share anything).

## Hard ceilings in the current code

- **Incidences must fit uint32** (`Hypergraph::eoff`, `voff` in sms.hpp):
  chain/star/tree n ≤ ~2.1B, clique n ≤ ~4.2B. Memory binds first anyway.
- **`genGraph`'s adjacency lists** (`vector<vector<unsigned>>`, ~56 B/table
  transient) are the silly bottleneck past n ≈ 100M for chain/star/tree.
  Fix: build the Hypergraph directly in `buildHypergraph` (exec.cpp) instead
  of going through QueryGraph — for chain, table i has attrs {i−1, i}
  (clipped at the ends); for star, table 0 has attrs {0..n−2} and leaf i has
  {i−1}; clique already constructs directly (`intersectionHypergraph`).
  The executor doesn't care about attribute numbering, only that each query
  edge is one shared attr. Trees still need `QueryGraph::randomTree`.
- **Message arena is 32 B/table** (uint64 × D=4). At n=2B that's 64 GB. If it
  binds, switch counts to uint32 mod 2^31−1 (halves it; keep the two-tree
  verification, it works under any modulus).
- Counts are reported mod 2^61−1 — do NOT "simplify" to plain uint64
  wraparound; huge counts are ≡ 0 mod 2^64 (their 2-adic valuation exceeds
  64), which silently blinds the two-tree cross-check.

## Practical tips

- Pin to a core (`taskset -c N`) and use 3+ points per size if you want
  paper-grade numbers; the current driver does one trial per point.
- Transparent hugepages help the big streaming arrays; `madvise` mode plus
  default glibc malloc is fine (all large allocations are single vectors).
- On multi-socket boxes, `numactl --interleave=all` avoids one node's memory
  filling up (single-threaded but > single-node-sized data).
- Memory-check a config cheaply: the guard prints the values size and refuses
  before allocating, so a dry run with small `EXEC_MEM_GIB` tells you the
  values footprint without OOMing anything.
- `verify` covers n ≤ 20K. After code changes for scale, also re-run one
  medium point (e.g. `exec chain 1000000 100`) and check the count matches
  the committed CSVs before burning hours on big runs (same n + rows + code
  ⇒ same data ⇒ same count).
