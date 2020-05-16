#include "omp.h"

#include "blackholes.h"

#include <algorithm.h>
#include <graph.h>

#include <atomic>
#include <iostream>
#include <limits>
#include <set>

using namespace std;

void ForceSkipPrecise(const vector<size_t>& tsOrder, const vector<char>& special, std::vector<size_t>& pos) {
    size_t skip = 0;
    for (size_t i = 1; i < pos.size(); i++) {
        size_t p = pos[i];
        size_t v = tsOrder[p];
        if (special[v]) {
            skip = i;
            break;
        }
    }

    for (size_t i = skip + 1; i < pos.size(); i++) {
        pos[i] = tsOrder.size() - pos.size() + i;
    }
}

void ForceSkipFast(const vector<size_t>& tsOrder, const vector<char>& special, std::vector<size_t>& pos) {
    size_t skip = 0;
    for (size_t i = 1; i < pos.size(); i++) {
        size_t p = pos[i];
        size_t v = tsOrder[p];
        if (special[v]) {
            skip = i;
            break;
        }
    }

    size_t skipPos = pos[skip];
    pos[0] = skipPos - 1;
    for (size_t i = 1; i < pos.size(); i++) {
        pos[i] = tsOrder.size() - pos.size() + i;
    }
}

void ForceSkip(const vector<size_t>& tsOrder, const vector<char>& special, std::vector<size_t>& pos, ESkipMode skipMode = ESkipMode::Precise) {
    switch(skipMode) {
        case ESkipMode::Precise:
            ForceSkipPrecise(tsOrder, special, pos);
            break;
        case ESkipMode::Fast:
            ForceSkipFast(tsOrder, special, pos);
            break;
        default:
            ForceSkipPrecise(tsOrder, special, pos);
    }
}

void TopsortSolver(const TGraph& graph, size_t ompThreads, ESkipMode skipMode, bool printDebugInfo) {
    TGraph graphRev = GetReverseGraph(graph);
    TGraph graphUndir = GetUndirGraph(graph);
    vector<size_t> roots = GetRoots(graphRev);
    size_t maxid = 0;
    size_t maxsize = 0;
    // do in parallel!!
    for (size_t i = 0; i < roots.size(); i++) {
        TClosure closure = GetClosure(graph, roots[i]);
        if (closure.size() > maxsize) {
            maxid = i;
            maxsize = closure.size();
        }
    }
    swap(roots[0], roots[maxid]);
    Print("roots", roots, printDebugInfo);
    vector<char> special(graph.size());
    vector<size_t> tsOrder = TopSort(graph, roots, special);
    Print("topsort", tsOrder, printDebugInfo);
    if (ompThreads > 1) {
        BruteForceParallel(graph, graphUndir, tsOrder, special, ompThreads, skipMode);
    } else {
        BruteForce(graph, graphUndir, tsOrder, special, skipMode);
    }
}

void BruteForce(const TGraph& graph, const TGraph& graphUndir, const vector<size_t>& tsOrder) {
    (void)graphUndir;
    const size_t totalVertex = static_cast<size_t>(tsOrder.size());
    vector<size_t> pos;
    pos.reserve(totalVertex);
    for (size_t sampleSize = 1; sampleSize <= totalVertex; sampleSize++) {
        pos.resize(sampleSize);
        for (size_t i = 0; i < pos.size(); ++i) {
            pos[i] = i;
        }
        while (true) {
            set<size_t> bh = GetBlackHole(graph, tsOrder, pos);
            FoundBlackHole(bh);
            if (!BruteNext(pos, totalVertex)) {
                break;
            }
        }
    }
}

bool CheckOutdegree(const TGraph& graph, const set<size_t>& bh) {
    for (size_t v : bh) {
        for (size_t to : graph[v]) {
            if (bh.find(to) == bh.end()) {
                return false;
            }
        }
    }
    return true;
}

bool ShouldSkip(const vector<size_t>& tsOrder,
                const vector<char>& special,
                const vector<size_t>& pos)
{
    for (size_t i = 1; i < pos.size(); i++) {
        size_t p = pos[i];
        size_t v = tsOrder[p];
        if (special[v]) {
            return true;
        }
    }
    return false;
}

void BruteForce(const TGraph& graph, const TGraph& graphUndir, const vector<size_t>& tsOrder, const vector<char>& special, ESkipMode skipMode, bool printDebugInfo) {
    const size_t totalVertex = tsOrder.size();
    const size_t maxSampleSize = totalVertex;
    vector<size_t> pos;
    pos.reserve(totalVertex);
    for (size_t sampleSize = 1; sampleSize <= min(maxSampleSize, totalVertex); sampleSize++) {
        size_t sampleSizeBlackholesFound = 0;
        pos.resize(sampleSize);
        for (size_t i = 0; i < pos.size(); ++i) {
            pos[i] = i;
        }
        bool stop = false;
        while (!stop) {
            set<size_t> bh = GetBlackHole(graph, tsOrder, pos);
            if (bh.size() > 0 && CheckConnectivity(graphUndir, bh)) {
                sampleSizeBlackholesFound++;
                FoundBlackHole(bh);
            } else {
                FilteredCandidate();
            }
            while (true) {
                stop = !BruteNext(pos, totalVertex);
                if (stop) {
                    break;
                }
                bool skip = ShouldSkip(tsOrder, special, pos);
                if (skip) {
                    ForceSkip(tsOrder, special, pos, skipMode);
                    FilteredCandidate();
                }
                if (!skip) {
                    break;
                }
            }
        }

        if (sampleSizeBlackholesFound == 0) {
            break;
        }

        if (printDebugInfo) {
            cout << "sampleSize done " << sampleSize << endl;
        }
    }
}

bool BruteNext(vector<size_t>& pos, size_t vertexCount) {
    size_t sampleSize = pos.size();
    bool res = true;
    for (size_t i = sampleSize; i > 0; --i) {
        size_t& cur = pos[i - 1];
        size_t distFromEnd = sampleSize - i;
        if (cur == vertexCount - 1 - distFromEnd) {
            if (i == 1) {
                res = false;
            }
            continue;
        } else {
            cur++;
            for (size_t j = i; j < sampleSize; j++) {
                pos[j] = pos[j - 1] + 1;
            }
            break;
        }
    }
    return res;
}

size_t ProcessSampleRange(size_t start,
                          size_t end,
                          const TGraph& graph,
                          const TGraph& graphUndir,
                          const vector<size_t>& tsOrder,
                          const vector<char>& special,
                          ESkipMode skipMode,
                          bool)
{
    std::vector<size_t> pos;
    pos.reserve(end);
    for (size_t sampleSize = start; sampleSize < end; ++sampleSize) {
        size_t totalVertex = tsOrder.size();
        size_t sampleSizeBlackholesFound = 0;
        pos.resize(sampleSize);
        for (size_t i = 0; i < pos.size(); ++i) {
            pos[i] = i;
        }
        bool stop = false;
        while (!stop) {
            set<size_t> bh = GetBlackHole(graph, tsOrder, pos);
            if (bh.size() > 0 && CheckConnectivity(graphUndir, bh)) {
                sampleSizeBlackholesFound++;
                FoundBlackHole(bh);
            } else {
                FilteredCandidate();
            }
            while (true) {
                stop = !BruteNext(pos, totalVertex);
                if (stop) {
                    break;
                }
                bool skip = ShouldSkip(tsOrder, special, pos);

                if (skip) {
                    FilteredCandidate();
                    ForceSkip(tsOrder, special, pos, skipMode);
                }
                if (!skip) {
                    break;
                }
            }
        }

        if (sampleSizeBlackholesFound == 0) {
            return sampleSize;
        }
    }
    return std::numeric_limits<size_t>::max();
}

void BruteForceParallel(const TGraph& graph,
                        const TGraph& graphUndir,
                        const vector<size_t>& tsOrder,
                        const vector<char>& special,
                        size_t ompThreads,
                        ESkipMode skipMode,
                        bool printDebugInfo)
{
    std::atomic<size_t> minimalEmptySize(std::numeric_limits<size_t>::max());
    #pragma omp parallel for num_threads(ompThreads) schedule(dynamic)
    for (size_t i = 0; i < tsOrder.size(); i++) {
        size_t start = i;
        size_t end = min(tsOrder.size(), i + 1);
        if (start < end && start < tsOrder.size()) {
            size_t ss = ProcessSampleRange(start + 1, end + 1, graph, graphUndir, tsOrder, special, skipMode, printDebugInfo);
            #pragma omp critical
            {
                size_t mes = minimalEmptySize.load();
                minimalEmptySize.store(min(ss, mes));
            }
        }
    }
    (void)printDebugInfo;
}

void FoundBlackHole(const set<size_t>& bh, bool printDebugInfo) {
    static std::atomic<size_t> bhCount(0);
    size_t oldValue = bhCount.fetch_add(1);
    cout << "bh_count " << oldValue + 1 << "\n";
    if ((oldValue & (1 << 10)) == 0) {
        cout.flush();
    }
    if (printDebugInfo) {
        cout << "BH: ";
        for (size_t x : bh) {
            cout << x << " ";
        }
        cout << endl;
    }
}

void FilteredCandidate() {
    static std::atomic<size_t> filtered(0);
    size_t oldValue = filtered.fetch_add(1);
    //printf("bh_count %zu\n", oldValue + 1);
    cout << "bh_filtered " << oldValue + 1 << "\n";
    if ((oldValue & (1 << 10)) == 0) {
        cout.flush();
    }
}

set<size_t> GetBlackHole(const TGraph& graph, const vector<size_t>& tsOrder, const vector<size_t>& pos) {
    set<size_t> blackhole;
    TUsed used;
    used.assign(graph.size(), 0);
    for (size_t p : pos) {
        size_t v = tsOrder[p];
        TClosure closure = GetClosure(graph, v, used);
        for (size_t x : closure) {
            blackhole.insert(x);
        }
    }
    for (size_t p : pos) {
        size_t x = tsOrder[p];
        if (used[x] >= 2) {
            return {}; // means this blackhole has already been detected;
        }
    }
    return blackhole;
}

void Print(const string& tag, const set<size_t>& v, bool printDebugInfo) {
    if (!printDebugInfo) {
        return;
    }
    cout << tag << ": ";
    for (size_t x : v) {
        cout << x << " ";
    }
    cout << endl;
}

void Print(const string& tag, const vector<size_t>& v, bool printDebugInfo) {
    if (!printDebugInfo) {
        return;
    }
    cout << tag << ": ";
    for (size_t x : v) {
        cout << x << " ";
    }
    cout << endl;
}

void Print(const string& tag, const vector<vector<size_t>>& vv, bool printDebugInfo) {
    if (!printDebugInfo) {
        return;
    }
    cout << tag << " BEGIN" << endl;
    for (size_t i = 0; i < vv.size(); i++) {
        const auto& v = vv[i];
        cout << i << ": ";
        for (size_t x : v) {
            cout << x << " ";
        }
        cout << endl;
    }
    cout << tag + " END" << endl;
}
