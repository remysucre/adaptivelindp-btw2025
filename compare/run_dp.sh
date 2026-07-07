#!/bin/bash
# Fig-15 single-pass DP kernels (eff/mcm); resumable via .done markers.
cd "$(dirname "$0")/.."
mkdir -p compare/results
for algo in dp-new dp-old; do
  for t in chain clique star tree; do
    out="compare/results/${algo}_${t}.csv"
    if [ -f "$out.done" ]; then continue; fi
    echo "=== $algo $t $(date +%H:%M:%S) ==="
    ./compare/bench "$algo" "$t" > "$out" && touch "$out.done"
  done
done
echo "DP-DONE $(date +%H:%M:%S)"
