# SMS vs. Adaptive LinDP: join tree construction within a time budget

This directory compares **SMS** (sort + maximal cardinality search), a
linear-time join tree construction algorithm, against the **adaptive LinDP**
join order optimizer implemented in this repository, using the experimental
setup of Fig. 15 of the BTW 2025 paper: for each query shape (chain, clique,
star, tree), *how large a query can each algorithm optimize within a given
time budget?*

![comparison](fig15_comparison.png)

## Max relations within a 1 s budget (single-threaded, Apple M-series)

| shape  | SMS (ours) | adaptive LinDP (full) | adaptive DP kernel | old DP kernel | vs LinDP | vs DP kernel |
|--------|-----------:|----------------------:|-------------------:|--------------:|---------:|-------------:|
| chain  | 16.2M      | 793                   | 19.9K              | 2.5K          | ~20,000× | ~800×        |
| clique | 16.2M      | 296                   | 12.2K              | 2.1K          | ~55,000× | ~1,300×      |
| star   | 16.2M      | 4.4K                  | 8M                 | 43.7K         | ~3,700×  | 2.0×         |
| tree   | 4.4M       | 1.8K                  | 3.1M               | 59.2K         | ~2,500×  | 1.4×         |

At the sizes the baseline needs ~1 s for, SMS takes: 2K chain ≈ 0.03 ms,
2K clique ≈ 0.02 ms, 3M star ≈ 105 ms, 1.5M tree ≈ 210 ms.

## The algorithms

**SMS** takes a *hypergraph* (one hyperedge per relation, one vertex per join
attribute) plus relation sizes, and produces a *join tree* for Yannakakis-style
semi-join processing:

1. Sort relations by size (LSD radix sort — linear time; the digit width
   scales with the input so tiny queries don't pay a fixed histogram cost).
2. Run Tarjan & Yannakakis' maximum cardinality search, which is essentially
   Prim's algorithm over the hypergraph, breaking ties toward larger
   relations. Runs in time linear in the hypergraph size.

The parent map produced by MCS *is* the join tree. A validator checks the
join tree property (for every attribute, the relations containing it form a
connected subtree); it passes on all acyclic inputs tested (n = 3 … 20K, all
shapes, 3 size-randomizations each).

**Adaptive LinDP** (this repo) takes a *query graph* (relations, binary join
predicates with selectivities) and produces an optimal *binary join operator
tree* restricted to linearizations: IKKBZ computes one linearization per root,
and an adaptive DP finds the best tree per linearization, sharing work across
linearizations. The "adaptive DP kernel" rows above run a *single* DP pass
over one fixed linearization — that is what Fig. 15 of the paper plots as
"new" / "old" (the full optimizer allocates an n×n table and cannot reach the
figure's star n = 3M).

## Setup notes

- **Same machine, same harness.** Both sides use this repo's query generators
  (`genTestData` in `src/main.cpp`), its repeat-until-stable timing harness
  (`src/Benchmark.cpp`), and its size-growth schedule (grow n geometrically
  until optimization exceeds ~1 s; 3 trials per size with re-randomized
  weights/sizes). The published Fig. 15 numbers come from different hardware,
  so only same-machine ratios are meaningful; the reproduced kernel curves
  match the published shapes.
- **Hypergraph encodings.** Chain/star/tree: one attribute per query graph
  edge, shared by its two relations. Clique: encoded as the intersection
  query R1(x), …, Rn(x) — every pair joined on the same variable, which is
  what a clique query graph means here; the hypergraph is n unary hyperedges
  (α-acyclic, size O(n)). Under the alternative reading (a distinct variable
  per pair) the hypergraph is cyclic and no join tree exists at all, while
  the baseline's input — the query graph — is identical either way.
- **What is timed.** For LinDP: everything after query graph construction.
  For SMS: everything after hypergraph construction, including all working
  memory allocation, the radix sort, and MCS. Relation sizes are uniform in
  [1, 10^6] per trial.
- **Outputs differ.** SMS produces a join tree for semi-join processing;
  LinDP produces a cost-optimal linearized operator tree under its cost
  model. Both solve "plan this join query", but SMS does no cost-based
  join-order search — its guarantee is structural (a valid join tree,
  size-prioritized) and it is the appropriate object for Yannakakis-style
  execution.

## Executing the plans (Yannakakis COUNT(*))

`exec.cpp` checks that the join trees SMS produces at the *largest* sizes are
actually executable: it materializes real relations (domain [0, D=4) per
attribute, D "diagonal" rows per table so the join is never empty, plus 0–3
random rows), builds the join tree from the actual table sizes, and computes
`COUNT(*)` of the full join in **one bottom-up Yannakakis pass** (semi-join +
COUNT pushdown fused; COUNT(*) needs no top-down pass). Counts are reported
mod 2^61−1 — the true counts are astronomical by design.

Measured at the largest sizes SMS plans within 1 s (Apple M-series,
single-threaded; `values` = Σ arity·rows). Default tiny tables (~5.5 rows
each) and 100-rows-per-table (`exec <shape> 0 100`):

| shape  | n (tables) | rows  | plan (SMS) | execute COUNT(*) | ns/value |
|--------|-----------:|------:|-----------:|-----------------:|---------:|
| chain  | 16.2M      | 89.1M | 0.30 s     | 0.52 s           | 2.9      |
| clique | 16.2M      | 89.1M | 0.23 s     | 0.70 s           | 7.9      |
| star   | 16.2M      | 89.1M | 0.31 s     | 0.60 s           | 3.2      |
| tree   | 4.4M       | 24.2M | 0.85 s     | 1.12 s           | 23.2     |
| chain  | 16.2M      | 1.64B | 0.28 s     | 3.29 s           | 1.00     |
| clique | 16.2M      | 1.64B | 0.16 s     | 5.62 s           | 3.42     |
| star   | 16.2M      | 1.64B | 0.25 s     | 5.43 s           | 1.65     |
| tree   | 4.4M       | 447M  | 0.83 s     | 2.41 s           | 2.70     |

Execution is linear-ish across the whole sweep (n = 1K … 16.2M): at 100
rows/table ns/value is flat at 1.0 (chain), 1.65 (star), 3.4 (clique), and
creeps 1.1 → 2.7 for random trees (parent/child accesses fall out of cache;
bigger constant, still near-linear) — i.e. ~1.6 *billion* rows counted in
3–6 s, roughly 300–500M rows/s single-threaded. With tiny tables per-table
overhead dominates (~3–8 ns/value; trees ~23). The executor iterates rows
outer / children inner so wide tables (a star center has one attribute per
child) stream memory sequentially instead of doing one strided column scan
per child.

Correctness: 120 brute-force comparisons (all shapes, n = 2…8, random D,
counts match exact join enumeration) plus a two-tree check at n = 1K/20K
(same database, two *different* SMS trees from reshuffled declared sizes must
give identical counts — they do, mod 2^61−1) — run `exec verify`. Every
executed tree is also validated with the join tree property checker first.

## Files

- `sms.hpp` — C++ port of SMS
  (from `Linear_time_semi_join_planning/sms/{radix_sort,mcs}.py`), the
  hypergraph encodings, the join tree validator, and the query generators;
  shared by `bench.cpp` and `exec.cpp`.
- `bench.cpp` — benchmark driver: SMS plus wrappers around this
  repo's LinDP / DP-kernel implementations.
  Usage: `bench <algo> <shape> [n]` with algo ∈ {sms, sms-verify,
  lindp-new-{none,basic,sorted}, lindp-old-sorted, dp-new, dp-old}.
- `exec.cpp` — Yannakakis COUNT(*) executor for SMS join trees (see above).
  Usage: `exec verify` or `exec <shape> [n] [rows] [maxN]` (n=0 sweeps;
  rows = rows per table, default ~5.5); sweeps write
  `results/exec_<shape>.csv` / `results/exec100_<shape>.csv`
  (`type,n,rows,values,gen_s,plan_s,exec_s,ns_per_value,count,overflow`).
- `SERVER-NOTES.md` — how to run the executor at larger scales on a server:
  memory/time models, hard ceilings, invocation, tips.
- `run_all.sh`, `run_dp.sh` — resumable drivers for the full suite.
- `results/*.csv` — raw measurements (`algo,type,n,m,trial,time`).
- `plot.py` — builds `fig15_comparison.png` and the summary table.

Build (no cmake needed):

```
c++ -std=c++2c -O3 -march=native -DNDEBUG -w compare/bench.cpp \
    src/Benchmark.cpp src/BumpAlloc.cpp src/DP.cpp src/LinDP.cpp \
    src/QueryGraph.cpp src/IKKBZ.cpp src/MDQ.cpp src/UnionFind.cpp \
    -o compare/bench
c++ -std=c++2c -O3 -march=native -DNDEBUG -w compare/exec.cpp \
    src/QueryGraph.cpp src/UnionFind.cpp -o compare/exec
compare/run_all.sh && compare/run_dp.sh && python3 compare/plot.py
```

(Two one-line portability fixes to `src/` — `unsigned __int128` and a missing
`<algorithm>` include — make the repo build with Apple clang.)
