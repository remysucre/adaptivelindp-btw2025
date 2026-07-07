// Benchmark driver comparing adaptive LinDP (baseline, this repo) against
// SMS (sort + Tarjan-Yannakakis maximum cardinality search), the linear-time
// join tree construction algorithm from ../Linear_time_semi_join_planning.
//
// Methodology follows doBenchmark() in src/main.cpp / Fig. 15 of the paper:
// for each query graph type (chain, clique, star, tree), grow n geometrically
// and measure per-optimization time until it exceeds ~1s. The resulting
// time-vs-n curve is inverted downstream into "largest query optimizable
// within a time budget".
//
// Usage: bench <algo> <type>
//   algo: lindp-new-none | lindp-new-basic | lindp-new-sorted |
//         lindp-old-sorted | sms | sms-verify
//   type: chain | clique | star | tree
// Writes CSV (algo,type,n,m,trial,time) to stdout.
#include "../src/Benchmark.hpp"
#include "../src/BumpAlloc.hpp"
#include "../src/DP.hpp"
#include "../src/IKKBZ.hpp"
#include "../src/LinDP.hpp"
#include "../src/Plan.hpp"
#include "../src/QueryGraph.hpp"
#include "../src/Random.hpp"
#include <cinttypes>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <numeric>
#include <string_view>
#include <vector>

using namespace std;
using namespace fastoptim;

//===========================================================================
// Baseline: LinDP, copied from testLinDP<false> in src/main.cpp
//===========================================================================
static void runLinDP(const QueryGraph& qg, LinDPOptions opt) {
    struct PlanC : Plan {
        float cardinality = 0.0f;
        uint64_t bs = 0;
    };
    BumpAlloc<PlanC> alloc;
    vector<Plan*> basePlans;
    basePlans.reserve(qg.size());
    for (unsigned i = 0; i < qg.size(); i++) {
        auto* p = &alloc.emplace();
        basePlans.push_back(p);
        p->cardinality = qg.nodes[i].cardinality;
        p->cost = 0.0f;
        p->bs = 1ull << (i % 64);
    }
    auto combine = [&](Plan* targetr, Plan* leftr, Plan* rightr, float sel) -> Plan* {
        auto* target = static_cast<PlanC*>(targetr);
        auto* left = static_cast<PlanC*>(leftr);
        auto* right = static_cast<PlanC*>(rightr);
        if (!target) {
            target = alloc.allocate();
            target->cardinality = left->cardinality * sel * right->cardinality;
        }
        target->cost = left->cost + right->cost + target->cardinality;
        target->bs = left->bs | right->bs;
        return target;
    };
    [[maybe_unused]] auto* best = LinDP::optimize(qg, basePlans, combine, opt);
}

//===========================================================================
// Fig. 15's actual subjects: one DP recomputation over the identity
// linearization; copied from mcm()/eff() in src/main.cpp with prep=false.
//===========================================================================
static void runSingleDP(const QueryGraph& g, bool newAlgo) {
    IncrementalDP dp(g);
    vector<unsigned> nodes(g.size());
    for (unsigned i = 0; i < g.size(); i++)
        nodes[i] = g.size() - 1 - i;
    auto cb = [&](unsigned, unsigned, unsigned, float) {};
    if (newAlgo)
        dp.updateNew(nodes, cb);
    else
        dp.updateOld(nodes, cb);
}

//===========================================================================
// Hypergraph encoding of a join query: one hyperedge per table, one vertex
// per join attribute. Every query graph edge (u,v) becomes an attribute
// shared by exactly tables u and v.
//===========================================================================
struct Hypergraph {
    uint32_t m = 0;            // #hyperedges (tables)
    uint32_t nV = 0;           // #vertices (attributes)
    vector<uint32_t> eoff;     // m+1 offsets into evtx
    vector<uint32_t> evtx;     // incidence lists (vertex ids), per edge
    vector<uint32_t> sizes;    // table cardinalities
};

static Hypergraph fromQueryGraph(const QueryGraph& q) {
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
    vector<uint32_t> cur(h.eoff.begin(), h.eoff.end() - 1);
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
static Hypergraph intersectionHypergraph(uint32_t n) {
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

static void randomizeSizes(Hypergraph& h, Random& rng) {
    for (auto& s : h.sizes)
        s = 1 + static_cast<uint32_t>(rng.nextRange(1'000'000));
}

//===========================================================================
// SMS: sort-MCS-sort, C++ port of sms/mcs.py + sms/sms.py
// (../Linear_time_semi_join_planning), non-active-sorting variant.
//===========================================================================
struct SMSResult {
    vector<uint32_t> gamma;        // parent of each hyperedge (UINT32_MAX for root)
    vector<uint32_t> res;          // final numbering (1-based)
    uint32_t root = 0;
    uint64_t checksum = 0;
};

static void radixSortBySize(const vector<uint32_t>& sizes, vector<uint32_t>& order) {
    // LSD radix sort of edge ids by 32-bit size key, two 16-bit passes
    const size_t m = order.size();
    vector<uint32_t> tmp(m);
    static constexpr size_t B = 1 << 16;
    vector<uint32_t> count(B);
    for (int pass = 0; pass < 2; pass++) {
        const int shift = pass * 16;
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
        swap(order, tmp);
    }
}

static SMSResult smsRun(const Hypergraph& h) {
    const uint32_t m = h.m;
    const uint32_t nV = h.nV;
    SMSResult out;

    // --- step 1: sort hyperedges by size (ascending), linear time ---
    vector<uint32_t> order(m);
    iota(order.begin(), order.end(), 0u);
    radixSortBySize(h.sizes, order);

    // --- build v2e with each vertex's edge list ascending by size ---
    // (mcs.py preprocessing(): v2e[v] sorted so reverse iteration visits
    //  heavier edges first)
    vector<uint32_t> voff(static_cast<size_t>(nV) + 1, 0);
    for (uint32_t v : h.evtx)
        voff[v + 1]++;
    for (uint32_t v = 0; v < nV; v++)
        voff[v + 1] += voff[v];
    vector<uint32_t> vedge(h.evtx.size());
    {
        vector<uint32_t> cur(voff.begin(), voff.end() - 1);
        for (uint32_t oi = 0; oi < m; oi++) {
            uint32_t e = order[oi];
            for (uint32_t k = h.eoff[e]; k < h.eoff[e + 1]; k++)
                vedge[cur[h.evtx[k]]++] = e;
        }
    }

    // --- step 2: maximum cardinality search (Tarjan & Yannakakis) ---
    vector<uint32_t> vSelected(m, 0);
    vector<uint8_t> eOut(m, 0);
    vector<uint8_t> vOut(nV, 0);
    vector<uint32_t> gamma(m, UINT32_MAX);
    vector<uint32_t> labelOrder;
    labelOrder.reserve(m);
    // buckets[c] holds edges with c marked vertices; lazily deleted
    vector<vector<uint32_t>> buckets(1);
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
    vector<uint32_t> res(m);
    for (uint32_t i = 0; i < static_cast<uint32_t>(labelOrder.size()); i++)
        res[labelOrder[i]] = i + 1;

    uint64_t chk = 0;
    for (uint32_t i = 0; i < m; i++)
        chk = chk * 0x9e3779b97f4a7c15ull + gamma[i] + res[i];
    out.checksum = chk;
    out.gamma = std::move(gamma);
    out.res = std::move(res);
    return out;
}

//===========================================================================
// Join tree validity check (not timed): for every attribute, the tables
// containing it must form a connected subtree of gamma.
//===========================================================================
static bool verifyJoinTree(const Hypergraph& h, const SMSResult& r) {
    // build v2e (unsorted is fine)
    vector<uint32_t> voff(static_cast<size_t>(h.nV) + 1, 0);
    for (uint32_t v : h.evtx)
        voff[v + 1]++;
    for (uint32_t v = 0; v < h.nV; v++)
        voff[v + 1] += voff[v];
    vector<uint32_t> vedge(h.evtx.size());
    {
        vector<uint32_t> cur(voff.begin(), voff.end() - 1);
        for (uint32_t e = 0; e < h.m; e++)
            for (uint32_t k = h.eoff[e]; k < h.eoff[e + 1]; k++)
                vedge[cur[h.evtx[k]]++] = e;
    }
    // check gamma forms a tree rooted at root (no cycles, all reach root)
    {
        vector<uint8_t> state(h.m, 0); // 0 unknown, 1 ok
        state[r.root] = 1;
        vector<uint32_t> path;
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
    vector<uint8_t> mark(h.m, 0);
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
static QueryGraph genGraph(string_view type, uint64_t n) {
    if (type == "clique") {
        vector<vector<unsigned>> adj(n);
        for (unsigned i = 0; i + 1 < n; i++)
            for (unsigned j = i + 1; j < n; j++)
                adj[i].push_back(j);
        return QueryGraph::fromAdjList(adj);
    } else if (type == "chain") {
        vector<vector<unsigned>> adj(n);
        for (unsigned i = 0; i + 1 < n; i++)
            adj[i].push_back(i + 1);
        return QueryGraph::fromAdjList(adj);
    } else if (type == "star") {
        vector<vector<unsigned>> adj(n);
        for (unsigned i = 1; i < n; i++)
            adj[0].push_back(i);
        return QueryGraph::fromAdjList(adj);
    } else { // tree: chainLen=0.5, density=0, seed as in main.cpp
        Random rng(1024);
        return QueryGraph::randomTree(static_cast<unsigned>(n), 0.5, 0.0, rng);
    }
}

//===========================================================================
int main(int argc, char* argv[]) {
    if (argc < 3) {
        fprintf(stderr, "usage: %s <algo> <type>\n", argv[0]);
        return 1;
    }
    string_view algo = argv[1];
    string_view type = argv[2];

    const bool isSMS = algo.starts_with("sms");
    const bool verify = (algo == "sms-verify");
    const bool isSingleDP = algo.starts_with("dp-");

    LinDPOptions opt;
    if (algo == "lindp-new-none")
        opt = {.ikkbz = {.transfer = LinerizationTransfer::None}, .dp = {.newAlgo = true}};
    else if (algo == "lindp-new-basic")
        opt = {.ikkbz = {.transfer = LinerizationTransfer::Basic}, .dp = {.newAlgo = true}};
    else if (algo == "lindp-new-sorted")
        opt = {.ikkbz = {.transfer = LinerizationTransfer::Sorted}, .dp = {.newAlgo = true}};
    else if (algo == "lindp-old-sorted")
        opt = {.ikkbz = {.transfer = LinerizationTransfer::Sorted}, .dp = {.newAlgo = false}};
    else if (!isSMS && algo != "dp-new" && algo != "dp-old") {
        fprintf(stderr, "unknown algo\n");
        return 1;
    }

    if (verify) {
        // Join tree validity of SMS output on acyclic inputs
        for (uint64_t n : {3ull, 5ull, 17ull, 100ull, 1024ull, 20000ull}) {
            Hypergraph h = (type == "clique")
                ? intersectionHypergraph(static_cast<uint32_t>(n))
                : fromQueryGraph(genGraph(type, n));
            for (uint64_t trial = 0; trial < 3; trial++) {
                Random rng(n * 131 + trial);
                randomizeSizes(h, rng);
                auto r = smsRun(h);
                bool ok = verifyJoinTree(h, r);
                printf("verify %s n=%" PRIu64 " trial=%" PRIu64 " : %s\n",
                       string(type).c_str(), n, trial, ok ? "OK" : "FAIL");
                if (!ok)
                    return 1;
            }
        }
        return 0;
    }

    auto nextN = [](uint64_t n) -> uint64_t {
        if (n < 16)
            return n + 1;
        if (n < 100000)
            return static_cast<uint64_t>(ceil(n * 1.0625));
        return static_cast<uint64_t>(ceil(n * 1.125));
    };

    // optional argv[3]: measure exactly this n and stop
    uint64_t onlyN = (argc > 3) ? strtoull(argv[3], nullptr, 10) : 0;

    printf("algo,type,n,m,trial,time\n");
    for (uint64_t n = onlyN ? onlyN : 3;; n = nextN(n)) {
        // resource caps
        if (!isSMS && type == "clique" && n * (n - 1) > 1'200'000'000ull)
            break;
        if (n > 150'000'000ull)
            break;
        if (algo == "dp-old" && n > 60000)
            break; // n^2 byte DP table
        if (algo == "dp-new" && n > 20'000'000ull)
            break;
        if (!isSMS && !isSingleDP && n >= 65535)
            break; // LinDP uses 16-bit versions + n^2 table

        QueryGraph q;
        Hypergraph h;
        if (isSMS) {
            h = (type == "clique") ? intersectionHypergraph(static_cast<uint32_t>(n))
                                   : fromQueryGraph(genGraph(type, n));
        } else {
            q = genGraph(type, n);
        }
        chrono::duration<double> avg{};
        size_t totTrials = 3;
        volatile uint64_t sink = 0;
        for (size_t trial = 0; trial < totTrials; trial++) {
            Random rng(n * 131 + trial);
            chrono::duration<double> dur;
            if (isSMS) {
                randomizeSizes(h, rng);
                dur = Benchmark::measure([&] { sink += smsRun(h).checksum; });
            } else if (isSingleDP) {
                q.randomizeWeights(rng);
                dur = Benchmark::measure([&] { runSingleDP(q, algo == "dp-new"); });
            } else {
                q.randomizeWeights(rng);
                dur = Benchmark::measure([&] { runLinDP(q, opt); });
            }
            uint64_t m = isSMS ? h.evtx.size() : q.edges.size();
            printf("%s,%s,%" PRIu64 ",%" PRIu64 ",%zu,%.9f\n",
                   string(algo).c_str(), string(type).c_str(), n, m, trial, dur.count());
            fflush(stdout);
            avg += dur / totTrials;
        }
        if (onlyN || avg >= chrono::duration<double>(1.05))
            break;
    }
    return 0;
}
