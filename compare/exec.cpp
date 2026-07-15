// Executes the join trees produced by SMS (sort + MCS): materializes actual
// relations for a query, builds the join tree from the *actual* table sizes,
// and computes COUNT(*) of the full join with a single bottom-up Yannakakis
// pass (semi-join + COUNT pushdown fused; for COUNT(*) no top-down pass or
// second bottom-up pass is needed). Runtime is linear in database size +
// hypergraph size, up to hashing — and since attribute domains are small the
// per-child "hash maps" are dense arrays, so it is genuinely linear.
//
// Data generation: each attribute ranges over domain [0, D). Every table gets
// D "diagonal" rows (v,v,...,v) — so the join is never empty and counts are
// non-trivial — plus uniformly random rows up to the requested rows/table
// (default D+0..3). Counts are computed modulo
// the Mersenne prime 2^61-1 (reduction is two shifts; counts below 2^61 are
// exact, and an overflow flag reports when reduction happened — with millions
// of tables the true count is astronomical by design). Mod 2^64 wraparound
// would collapse huge counts to 0 (their power-of-two factors exceed 2^64),
// making cross-checks between different trees vacuous; an odd modulus keeps
// them meaningful.
//
// Correctness checks:
//   * small n (all shapes): count equals brute-force join enumeration
//   * large n: two *different* join trees (from different declared sizes)
//     must yield the same count on the same database
//   * every tree is validated with verifyJoinTree before execution
//
// Usage:
//   exec verify                    run correctness checks
//   exec <type> [n] [rows] [maxN]  sweep n geometrically up to maxN (default:
//                                  largest size SMS plans within 1s — chain/
//                                  clique/star 16.2M, tree 4.4M; see
//                                  ../README.md), or measure exactly n if
//                                  given (n=0: sweep);
//                                  rows = rows per table (default: D+0..3)
//   type: chain | clique | star | tree
// Env: EXEC_MEM_GIB=<gib> raises the database size guard (default 20 GiB).
// Sweep writes CSV (type,n,rows,values,gen_s,plan_s,exec_s,ns_per_value,
// count,ovf). See SERVER-NOTES.md for running at larger scales.
#include "../src/QueryGraph.hpp"
#include "../src/Random.hpp"
#include "sms.hpp"
#include <chrono>
#include <cinttypes>
#include <cmath>
#include <cstdio>
#include <string_view>
#include <vector>

using namespace std;
using namespace fastoptim;

//===========================================================================
// Database: one table per hyperedge, row-major values, domain [0, D) per
// attribute (uint8_t, so D <= 255).
//===========================================================================
struct Database {
    vector<uint64_t> valoff; // m+1: offset of table e's values in vals
    vector<uint32_t> nrows;  // m
    vector<uint8_t> vals;    // table e, row r, attr position k:
                             //   vals[valoff[e] + r*arity(e) + k]
    uint64_t totalRows = 0;
};

// baseRows = 0 means the default D + 0..3 rows per table
static Database genDatabase(const Hypergraph& h, uint32_t D, Random& rng,
                            uint32_t baseRows = 0) {
    const uint32_t rbase = max(baseRows, D); // >= D for the diagonal rows
    Database db;
    db.nrows.resize(h.m);
    db.valoff.assign(static_cast<size_t>(h.m) + 1, 0);
    for (uint32_t e = 0; e < h.m; e++) {
        uint32_t arity = h.eoff[e + 1] - h.eoff[e];
        uint32_t rows = rbase + static_cast<uint32_t>(rng.nextRange(4));
        db.nrows[e] = rows;
        db.valoff[e + 1] = db.valoff[e] + uint64_t(rows) * arity;
        db.totalRows += rows;
    }
    uint64_t memGib = 20;
    if (const char* env = getenv("EXEC_MEM_GIB"))
        memGib = strtoull(env, nullptr, 10);
    if (db.valoff[h.m] > memGib << 30) {
        fprintf(stderr,
                "database would be %.1f GiB, refusing (set EXEC_MEM_GIB to raise)\n",
                static_cast<double>(db.valoff[h.m]) / (1ull << 30));
        exit(1);
    }
    db.vals.resize(db.valoff[h.m]);
    for (uint32_t e = 0; e < h.m; e++) {
        uint32_t arity = h.eoff[e + 1] - h.eoff[e];
        uint8_t* base = db.vals.data() + db.valoff[e];
        for (uint32_t v = 0; v < D; v++) // diagonal rows: join never empty
            for (uint32_t k = 0; k < arity; k++)
                base[uint64_t(v) * arity + k] = static_cast<uint8_t>(v);
        uint8_t* p = base + uint64_t(D) * arity;
        uint64_t nrand = uint64_t(db.nrows[e] - D) * arity;
        if ((D & (D - 1)) == 0) { // power of two: 8 uniform values per rng draw
            uint64_t buf = 0;
            unsigned have = 0;
            for (uint64_t idx = 0; idx < nrand; idx++) {
                if (!have) {
                    buf = rng.next();
                    have = 8;
                }
                p[idx] = static_cast<uint8_t>(buf) & static_cast<uint8_t>(D - 1);
                buf >>= 8;
                have--;
            }
        } else {
            for (uint64_t idx = 0; idx < nrand; idx++)
                p[idx] = static_cast<uint8_t>(rng.nextRange(D));
        }
    }
    return db;
}

//===========================================================================
// COUNT(*) via one bottom-up pass over the join tree. For each non-root
// table c, the "message" to its parent maps each value combination of
// key(c) = attrs(c) ∩ attrs(parent(c)) to the COUNT(*) of the join of c's
// subtree, grouped by key(c). By the join tree property, key(c) separates
// c's subtree from the rest, so multiplying messages into parent rows is
// exact. Messages are dense arrays of size D^|key| (all our shapes have
// |key| = 1). Arithmetic is mod 2^61-1 with reduction tracking.
//===========================================================================
static constexpr uint64_t MODP = (1ull << 61) - 1;

struct ExecResult {
    uint64_t count = 0; // mod 2^61-1
    bool overflow = false;
    bool ok = true; // structural checks passed
};

static ExecResult yannakakisCount(const Hypergraph& h, const Database& db,
                                  const SMSResult& r, uint32_t D) {
    const uint32_t m = h.m;
    ExecResult out;

    // children lists (CSR over gamma)
    vector<uint32_t> coff(static_cast<size_t>(m) + 1, 0);
    for (uint32_t e = 0; e < m; e++)
        if (e != r.root)
            coff[r.gamma[e] + 1]++;
    for (uint32_t e = 0; e < m; e++)
        coff[e + 1] += coff[e];
    vector<uint32_t> clist(m ? m - 1 : 0);
    {
        vector<uint32_t> cur(coff.begin(), coff.end() - 1);
        for (uint32_t e = 0; e < m; e++)
            if (e != r.root)
                clist[cur[r.gamma[e]]++] = e;
    }

    // key attrs per non-root edge c: positions (within c's attr list) of the
    // attrs shared with parent(c). Linear: stamp each parent's attrs once.
    vector<uint32_t> sstart(m, 0), slen(m, 0);
    vector<uint32_t> sattr; // concatenated position lists
    sattr.reserve(h.evtx.size());
    {
        vector<uint32_t> stamp(h.nV, UINT32_MAX);
        for (uint32_t f = 0; f < m; f++) {
            if (coff[f] == coff[f + 1])
                continue;
            for (uint32_t k = h.eoff[f]; k < h.eoff[f + 1]; k++)
                stamp[h.evtx[k]] = f;
            for (uint32_t ci = coff[f]; ci < coff[f + 1]; ci++) {
                uint32_t c = clist[ci];
                sstart[c] = static_cast<uint32_t>(sattr.size());
                for (uint32_t k = h.eoff[c]; k < h.eoff[c + 1]; k++)
                    if (stamp[h.evtx[k]] == f)
                        sattr.push_back(k - h.eoff[c]);
                slen[c] = static_cast<uint32_t>(sattr.size()) - sstart[c];
            }
        }
    }

    // message arena: table c gets D^slen[c] uint64 slots
    vector<uint64_t> mapoff(static_cast<size_t>(m) + 1, 0);
    for (uint32_t e = 0; e < m; e++) {
        uint64_t sz = 1;
        if (e == r.root)
            sz = 0;
        else
            for (uint32_t j = 0; j < slen[e]; j++) {
                sz *= D;
                if (sz > (1u << 22)) { // never hit by our shapes (|key| = 1)
                    fprintf(stderr, "join key domain too large for dense messages\n");
                    out.ok = false;
                    return out;
                }
            }
        mapoff[e + 1] = mapoff[e] + sz;
    }
    vector<uint64_t> arena(mapoff[m], 0);

    bool ovf = false;
    auto mul = [&](uint64_t a, uint64_t b) { // a, b < 2^61
        unsigned __int128 z = static_cast<unsigned __int128>(a) * b;
        ovf |= z >= MODP;
        uint64_t res = static_cast<uint64_t>(z & MODP) + static_cast<uint64_t>(z >> 61);
        res = (res & MODP) + (res >> 61);
        return res >= MODP ? res - MODP : res;
    };
    auto add = [&](uint64_t a, uint64_t b) { // a, b < 2^61
        uint64_t res = a + b;
        ovf |= res >= MODP;
        return res >= MODP ? res - MODP : res;
    };

    // bottom-up: reverse MCS label order visits children before parents.
    // Rows are the outer loop: each row's values are scanned once,
    // sequentially — for wide tables (a star center / clique root has one
    // attr per child) the child-outer alternative would do one strided,
    // cache-missing column scan per child.
    vector<uint32_t> pos(h.nV, UINT32_MAX); // attr -> position in current table
    vector<uint32_t> cpp;                   // child key attrs as parent positions
    vector<uint64_t> cppoff;                // per-child offsets into cpp
    vector<const uint64_t*> cmaps;          // per-child message arrays
    uint64_t total = 0;
    for (uint32_t i = m; i-- > 0;) {
        uint32_t f = r.labelOrder[i];
        uint32_t arity = h.eoff[f + 1] - h.eoff[f];
        uint32_t rows = db.nrows[f];
        const uint8_t* base = db.vals.data() + db.valoff[f];
        for (uint32_t k = 0; k < arity; k++)
            pos[h.evtx[h.eoff[f] + k]] = k;
        uint32_t nch = coff[f + 1] - coff[f];
        cpp.clear();
        cppoff.assign(1, 0);
        cmaps.clear();
        for (uint32_t ci = coff[f]; ci < coff[f + 1]; ci++) {
            uint32_t c = clist[ci];
            for (uint32_t j = 0; j < slen[c]; j++)
                cpp.push_back(pos[h.evtx[h.eoff[c] + sattr[sstart[c] + j]]]);
            cppoff.push_back(cpp.size());
            cmaps.push_back(arena.data() + mapoff[c]);
        }
        uint64_t* fmap = (f == r.root) ? nullptr : arena.data() + mapoff[f];
        for (uint32_t rw = 0; rw < rows; rw++) {
            const uint8_t* rowv = base + uint64_t(rw) * arity;
            uint64_t cval = 1;
            for (uint32_t ci = 0; ci < nch; ci++) {
                uint64_t key = 0;
                for (uint64_t t = cppoff[ci]; t < cppoff[ci + 1]; t++)
                    key = key * D + rowv[cpp[t]];
                cval = mul(cval, cmaps[ci][key]);
            }
            if (f == r.root) {
                total = add(total, cval);
            } else {
                uint64_t key = 0;
                for (uint32_t j = 0; j < slen[f]; j++)
                    key = key * D + rowv[sattr[sstart[f] + j]];
                fmap[key] = add(fmap[key], cval);
            }
        }
        for (uint32_t k = 0; k < arity; k++)
            pos[h.evtx[h.eoff[f] + k]] = UINT32_MAX;
    }
    out.count = total;
    out.overflow = ovf;
    return out;
}

//===========================================================================
// Brute force COUNT(*): backtracking join enumeration (mod 2^61-1). Small n only.
//===========================================================================
static uint64_t bruteCount(const Hypergraph& h, const Database& db) {
    vector<int> assign(h.nV, -1);
    uint64_t count = 0;
    vector<uint32_t> newly;
    auto rec = [&](auto&& self, uint32_t e) -> void {
        if (e == h.m) {
            count++;
            return;
        }
        uint32_t arity = h.eoff[e + 1] - h.eoff[e];
        const uint8_t* base = db.vals.data() + db.valoff[e];
        for (uint32_t rw = 0; rw < db.nrows[e]; rw++) {
            size_t mark = newly.size();
            bool ok = true;
            for (uint32_t k = 0; k < arity; k++) {
                uint32_t v = h.evtx[h.eoff[e] + k];
                int val = base[uint64_t(rw) * arity + k];
                if (assign[v] < 0) {
                    assign[v] = val;
                    newly.push_back(v);
                } else if (assign[v] != val) {
                    ok = false;
                    break;
                }
            }
            if (ok)
                self(self, e + 1);
            while (newly.size() > mark) {
                assign[newly.back()] = -1;
                newly.pop_back();
            }
        }
    };
    rec(rec, 0);
    return count;
}

//===========================================================================
static Hypergraph buildHypergraph(string_view type, uint64_t n) {
    return (type == "clique") ? intersectionHypergraph(static_cast<uint32_t>(n))
                              : fromQueryGraph(genGraph(type, n));
}

static double seconds(auto fn) {
    auto t0 = chrono::steady_clock::now();
    fn();
    return chrono::duration<double>(chrono::steady_clock::now() - t0).count();
}

static int runVerify() {
    int failures = 0;
    // brute-force comparison at small n
    for (string_view type : {"chain", "star", "tree", "clique"}) {
        for (uint64_t n : {2ull, 3ull, 4ull, 5ull, 6ull, 8ull}) {
            Hypergraph h = buildHypergraph(type, n);
            for (uint64_t trial = 0; trial < 5; trial++) {
                Random rng(n * 977 + trial);
                uint32_t D = 2 + static_cast<uint32_t>(rng.nextRange(3)); // 2..4
                uint32_t R = (n <= 4 && trial >= 3) ? 12 : 0; // wider tables
                Database db = genDatabase(h, D, rng, R);
                h.sizes.assign(db.nrows.begin(), db.nrows.end());
                SMSResult r = smsRun(h);
                bool treeOk = verifyJoinTree(h, r);
                ExecResult y = yannakakisCount(h, db, r, D);
                uint64_t ref = bruteCount(h, db);
                bool ok = treeOk && y.ok && !y.overflow && y.count == ref;
                printf("verify %-6s n=%" PRIu64 " D=%u R=%u trial=%" PRIu64
                       " yann=%" PRIu64 " brute=%" PRIu64 " : %s\n",
                       string(type).c_str(), n, D, R, trial, y.count, ref,
                       ok ? "OK" : "FAIL");
                failures += !ok;
            }
        }
    }
    // two-tree consistency at larger n: same database, different join trees
    for (string_view type : {"chain", "star", "tree", "clique"}) {
        for (uint64_t n : {1000ull, 20000ull}) {
            Hypergraph h = buildHypergraph(type, n);
            Random rng(n * 31 + 7);
            uint32_t D = 4;
            uint32_t R = (n == 20000) ? 100 : 0; // wide tables on one size
            Database db = genDatabase(h, D, rng, R);
            h.sizes.assign(db.nrows.begin(), db.nrows.end());
            SMSResult r1 = smsRun(h);
            randomizeSizes(h, rng); // declared sizes reshuffled -> different tree
            SMSResult r2 = smsRun(h);
            ExecResult y1 = yannakakisCount(h, db, r1, D);
            ExecResult y2 = yannakakisCount(h, db, r2, D);
            bool ok = verifyJoinTree(h, r1) && verifyJoinTree(h, r2) && y1.ok &&
                      y2.ok && r1.checksum != r2.checksum && y1.count == y2.count;
            printf("two-tree %-6s n=%" PRIu64 " count1=%" PRIu64 " count2=%" PRIu64
                   "%s : %s\n",
                   string(type).c_str(), n, y1.count, y2.count,
                   y1.overflow ? " (mod 2^61-1)" : "", ok ? "OK" : "FAIL");
            failures += !ok;
        }
    }
    return failures ? 1 : 0;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        fprintf(stderr, "usage: %s verify | %s <type> [n]\n", argv[0], argv[0]);
        return 1;
    }
    string_view type = argv[1];
    if (type == "verify")
        return runVerify();

    // default: largest n SMS plans within the 1s budget (compare/README.md)
    uint64_t maxN = (type == "tree") ? 4'400'000ull : 16'200'000ull;
    const uint64_t onlyN = (argc > 2) ? strtoull(argv[2], nullptr, 10) : 0;
    const uint32_t rowsPerTable = (argc > 3) ? static_cast<uint32_t>(strtoul(argv[3], nullptr, 10)) : 0;
    if (argc > 4)
        maxN = strtoull(argv[4], nullptr, 10);
    const uint32_t D = 4;

    printf("type,n,rows,values,gen_s,plan_s,exec_s,ns_per_value,count,overflow\n");
    for (uint64_t n = onlyN ? onlyN : 1000;; n = min(n * 2, maxN)) {
        Hypergraph h;
        Database db;
        Random rng(n * 131);
        double genS = seconds([&] {
            h = buildHypergraph(type, n);
            db = genDatabase(h, D, rng, rowsPerTable);
            h.sizes.assign(db.nrows.begin(), db.nrows.end());
        });
        SMSResult r;
        double planS = seconds([&] { r = smsRun(h); });
        if (!verifyJoinTree(h, r)) {
            fprintf(stderr, "join tree INVALID at n=%" PRIu64 "\n", n);
            return 1;
        }
        ExecResult y;
        double execS = seconds([&] { y = yannakakisCount(h, db, r, D); });
        if (!y.ok)
            return 1;
        uint64_t values = db.valoff[h.m];
        printf("%s,%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%.6f,%.6f,%.6f,%.2f,%" PRIu64 ",%d\n",
               string(type).c_str(), n, db.totalRows, values, genS, planS, execS,
               execS * 1e9 / static_cast<double>(values), y.count, y.overflow ? 1 : 0);
        fflush(stdout);
        if (onlyN || n >= maxN)
            break;
    }
    return 0;
}
