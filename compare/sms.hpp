#pragma once
// SMS (sort + Tarjan-Yannakakis maximum cardinality search) join tree
// construction, shared between bench.cpp (timing) and exec.cpp (execution).
// C++ port of sms/{radix_sort,mcs}.py from ../Linear_time_semi_join_planning,
// non-active-sorting variant.
#include "../src/QueryGraph.hpp"
#include "../src/Random.hpp"
#include <algorithm>
#include <bit>
#include <cstdint>
#include <cstring>
#include <numeric>
#include <string_view>
#include <vector>

//===========================================================================
// Hypergraph encoding of a join query: one hyperedge per table, one vertex
// per join attribute. Every query graph edge (u,v) becomes an attribute
// shared by exactly tables u and v.
//===========================================================================
struct Hypergraph {
    uint32_t m = 0;                  // #hyperedges (tables)
    uint32_t nV = 0;                 // #vertices (attributes)
    std::vector<uint32_t> eoff;      // m+1 offsets into evtx
    std::vector<uint32_t> evtx;      // incidence lists (vertex ids), per edge
    std::vector<uint32_t> sizes;     // table cardinalities
};

inline Hypergraph fromQueryGraph(const fastoptim::QueryGraph& q) {
    Hypergraph h;
    h.m = static_cast<uint32_t>(q.size());
    h.eoff.assign(h.m + 1, 0);
    // degree = number of incident query graph edges
    for (unsigned i = 0; i < q.size(); i++) {
        for (auto& e : q.getEdges(i)) {
            h.eoff[i + 1]++;
            h.eoff[e.target + 1]++;
        }
    }
    for (uint32_t i = 0; i < h.m; i++)
        h.eoff[i + 1] += h.eoff[i];
    h.evtx.resize(h.eoff[h.m]);
    std::vector<uint32_t> cur(h.eoff.begin(), h.eoff.end() - 1);
    uint32_t vid = 0;
    for (unsigned i = 0; i < q.size(); i++) {
        for (auto& e : q.getEdges(i)) {
            h.evtx[cur[i]++] = vid;
            h.evtx[cur[e.target]++] = vid;
            vid++;
        }
    }
    h.nV = vid;
    h.sizes.assign(h.m, 1);
    return h;
}

// Clique query graphs arise from the intersection query R1(x),...,Rn(x):
// every pair joined, all on the same variable. The hypergraph is n unary
// hyperedges sharing one vertex (alpha-acyclic, size O(n)) — unlike the
// pairwise-distinct-variables reading, which is cyclic and has no join tree.
inline Hypergraph intersectionHypergraph(uint32_t n) {
    Hypergraph h;
    h.m = n;
    h.nV = 1;
    h.eoff.resize(n + 1);
    for (uint32_t i = 0; i <= n; i++)
        h.eoff[i] = i;
    h.evtx.assign(n, 0);
    h.sizes.assign(n, 1);
    return h;
}

inline void randomizeSizes(Hypergraph& h, fastoptim::Random& rng) {
    for (auto& s : h.sizes)
        s = 1 + static_cast<uint32_t>(rng.nextRange(1'000'000));
}

//===========================================================================
// SMS: sort-MCS
//===========================================================================
struct SMSResult {
    std::vector<uint32_t> gamma;      // parent of each hyperedge (UINT32_MAX for root)
    std::vector<uint32_t> res;        // final numbering (1-based)
    std::vector<uint32_t> labelOrder; // edges in labeling order (parents first)
    uint32_t root = 0;
    uint64_t checksum = 0;
};

inline void radixSortBySize(const std::vector<uint32_t>& sizes, std::vector<uint32_t>& order) {
    // LSD radix sort of edge ids by 32-bit size key. The digit width scales
    // with m (#buckets ~ m, clamped to [2^8, 2^16]) and passes cover only the
    // bits where keys differ, so tiny inputs don't pay a fixed 2^16-bucket
    // prefix scan (~30us) per call.
    const size_t m = order.size();
    if (m < 2)
        return;
    std::vector<uint32_t> tmp(m);
    uint32_t lo = ~0u, hi = 0;
    for (size_t i = 0; i < m; i++) {
        lo &= sizes[order[i]];
        hi |= sizes[order[i]];
    }
    const int bits = std::bit_width(lo ^ hi); // highest differing bit
    if (!bits)
        return; // all keys equal
    const int digit = std::clamp(static_cast<int>(std::bit_width(m)), 8, 16);
    const size_t B = size_t(1) << digit;
    std::vector<uint32_t> count(B);
    for (int shift = 0; shift < bits; shift += digit) {
        memset(count.data(), 0, B * sizeof(uint32_t));
        for (size_t i = 0; i < m; i++)
            count[(sizes[order[i]] >> shift) & (B - 1)]++;
        uint32_t sum = 0;
        for (size_t d = 0; d < B; d++) {
            uint32_t c = count[d];
            count[d] = sum;
            sum += c;
        }
        for (size_t i = 0; i < m; i++)
            tmp[count[(sizes[order[i]] >> shift) & (B - 1)]++] = order[i];
        std::swap(order, tmp);
    }
}

inline SMSResult smsRun(const Hypergraph& h) {
    const uint32_t m = h.m;
    const uint32_t nV = h.nV;
    SMSResult out;

    // --- step 1: sort hyperedges by size (ascending), linear time ---
    std::vector<uint32_t> order(m);
    std::iota(order.begin(), order.end(), 0u);
    radixSortBySize(h.sizes, order);

    // --- build v2e with each vertex's edge list ascending by size ---
    // (mcs.py preprocessing(): v2e[v] sorted so reverse iteration visits
    //  heavier edges first)
    std::vector<uint32_t> voff(static_cast<size_t>(nV) + 1, 0);
    for (uint32_t v : h.evtx)
        voff[v + 1]++;
    for (uint32_t v = 0; v < nV; v++)
        voff[v + 1] += voff[v];
    std::vector<uint32_t> vedge(h.evtx.size());
    {
        std::vector<uint32_t> cur(voff.begin(), voff.end() - 1);
        for (uint32_t oi = 0; oi < m; oi++) {
            uint32_t e = order[oi];
            for (uint32_t k = h.eoff[e]; k < h.eoff[e + 1]; k++)
                vedge[cur[h.evtx[k]]++] = e;
        }
    }

    // --- step 2: maximum cardinality search (Tarjan & Yannakakis) ---
    std::vector<uint32_t> vSelected(m, 0);
    std::vector<uint8_t> eOut(m, 0);
    std::vector<uint8_t> vOut(nV, 0);
    std::vector<uint32_t> gamma(m, UINT32_MAX);
    std::vector<uint32_t> labelOrder;
    labelOrder.reserve(m);
    // buckets[c] holds edges with c marked vertices; lazily deleted
    std::vector<std::vector<uint32_t>> buckets(1);
    buckets[0] = order; // ascending by size; pop from back -> largest first

    int64_t j = 0;
    uint32_t e = buckets[0].back();
    buckets[0].pop_back();
    out.root = e;
    bool first = true;
    while (j >= 0) {
        if (!first) {
            e = buckets[j].back();
            buckets[j].pop_back();
        }
        first = false;
        eOut[e] = 1;
        labelOrder.push_back(e);
        for (uint32_t k = h.eoff[e]; k < h.eoff[e + 1]; k++) {
            uint32_t v = h.evtx[k];
            if (vOut[v])
                continue;
            for (uint32_t p = voff[v + 1]; p-- > voff[v];) { // heaviest first
                uint32_t s = vedge[p];
                if (eOut[s])
                    continue;
                gamma[s] = e;
                uint32_t c = ++vSelected[s];
                if (c < h.eoff[s + 1] - h.eoff[s]) { // still maximal
                    if (c >= buckets.size())
                        buckets.resize(c + 1);
                    buckets[c].push_back(s);
                    if (static_cast<int64_t>(c) > j)
                        j = c;
                }
            }
            vOut[v] = 1;
        }
        while (j >= 0) {
            if (buckets[j].empty()) {
                j--;
                continue;
            }
            if (eOut[buckets[j].back()])
                buckets[j].pop_back();
            else
                break;
        }
    }

    // MCS numbering (the actual impl has no second, per-sibling sort phase)
    std::vector<uint32_t> res(m);
    for (uint32_t i = 0; i < static_cast<uint32_t>(labelOrder.size()); i++)
        res[labelOrder[i]] = i + 1;

    uint64_t chk = 0;
    for (uint32_t i = 0; i < m; i++)
        chk = chk * 0x9e3779b97f4a7c15ull + gamma[i] + res[i];
    out.checksum = chk;
    out.gamma = std::move(gamma);
    out.res = std::move(res);
    out.labelOrder = std::move(labelOrder);
    return out;
}

//===========================================================================
// Join tree validity check (not timed): for every attribute, the tables
// containing it must form a connected subtree of gamma.
//===========================================================================
inline bool verifyJoinTree(const Hypergraph& h, const SMSResult& r) {
    // build v2e (unsorted is fine)
    std::vector<uint32_t> voff(static_cast<size_t>(h.nV) + 1, 0);
    for (uint32_t v : h.evtx)
        voff[v + 1]++;
    for (uint32_t v = 0; v < h.nV; v++)
        voff[v + 1] += voff[v];
    std::vector<uint32_t> vedge(h.evtx.size());
    {
        std::vector<uint32_t> cur(voff.begin(), voff.end() - 1);
        for (uint32_t e = 0; e < h.m; e++)
            for (uint32_t k = h.eoff[e]; k < h.eoff[e + 1]; k++)
                vedge[cur[h.evtx[k]]++] = e;
    }
    // check gamma forms a tree rooted at root (no cycles, all reach root)
    {
        std::vector<uint8_t> state(h.m, 0); // 0 unknown, 1 ok
        state[r.root] = 1;
        std::vector<uint32_t> path;
        for (uint32_t e = 0; e < h.m; e++) {
            uint32_t x = e;
            path.clear();
            while (!state[x]) {
                path.push_back(x);
                if (r.gamma[x] == UINT32_MAX)
                    return false; // non-root without parent
                x = r.gamma[x];
                if (path.size() > h.m)
                    return false; // cycle
            }
            for (uint32_t y : path)
                state[y] = 1;
        }
    }
    // connectedness: within the edge set of a vertex, exactly one member has
    // its parent outside the set (or is the root)
    std::vector<uint8_t> mark(h.m, 0);
    for (uint32_t v = 0; v < h.nV; v++) {
        uint32_t lo = voff[v], hi = voff[v + 1];
        for (uint32_t k = lo; k < hi; k++)
            mark[vedge[k]] = 1;
        uint32_t tops = 0;
        for (uint32_t k = lo; k < hi; k++) {
            uint32_t e = vedge[k];
            if (e == r.root || !mark[r.gamma[e]])
                tops++;
        }
        for (uint32_t k = lo; k < hi; k++)
            mark[vedge[k]] = 0;
        if (hi > lo && tops != 1)
            return false;
    }
    return true;
}

//===========================================================================
// Query generation, identical to genTestData() in src/main.cpp
//===========================================================================
inline fastoptim::QueryGraph genGraph(std::string_view type, uint64_t n) {
    using fastoptim::QueryGraph;
    if (type == "clique") {
        std::vector<std::vector<unsigned>> adj(n);
        for (unsigned i = 0; i + 1 < n; i++)
            for (unsigned j = i + 1; j < n; j++)
                adj[i].push_back(j);
        return QueryGraph::fromAdjList(adj);
    } else if (type == "chain") {
        std::vector<std::vector<unsigned>> adj(n);
        for (unsigned i = 0; i + 1 < n; i++)
            adj[i].push_back(i + 1);
        return QueryGraph::fromAdjList(adj);
    } else if (type == "star") {
        std::vector<std::vector<unsigned>> adj(n);
        for (unsigned i = 1; i < n; i++)
            adj[0].push_back(i);
        return QueryGraph::fromAdjList(adj);
    } else { // tree: chainLen=0.5, density=0, seed as in main.cpp
        fastoptim::Random rng(1024);
        return QueryGraph::randomTree(static_cast<unsigned>(n), 0.5, 0.0, rng);
    }
}
