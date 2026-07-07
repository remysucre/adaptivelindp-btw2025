#!/bin/bash
# Sequential benchmark suite; resumable (skips runs whose .done marker exists).
cd "$(dirname "$0")/.."
mkdir -p compare/results
for algo in lindp-new-sorted lindp-new-none lindp-new-basic lindp-old-sorted sms; do
  for t in chain clique star tree; do
    out="compare/results/${algo}_${t}.csv"
    if [ -f "$out.done" ]; then continue; fi
    echo "=== $algo $t $(date +%H:%M:%S) ==="
    ./compare/bench "$algo" "$t" > "$out" && touch "$out.done"
  done
done
echo "ALL-DONE $(date +%H:%M:%S)"
