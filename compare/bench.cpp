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
#include "sms.hpp"
#include <algorithm>
#include <bit>
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
