#include <iostream>
#include <queue>
#include <chrono>

#include <algorithm.h>
#include <blackholes.h>
#include <common.h>
#include <graph.h>
#include <io.h>


using namespace std;

using TCandidates = set<size_t>;

bool printDebugInfo = false;

set<size_t> GetBlackHoleChinese(const TGraph& graph, const vector<size_t>& tsOrder, const vector<size_t>& pos) {
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

void BruteForce(const TGraph& graph, const TGraph& graphUndir, const vector<size_t>& tsOrder, size_t bhSize) {
    size_t filtered = 0;
    const size_t totalVertex = static_cast<size_t>(tsOrder.size());
    vector<size_t> pos;
    pos.reserve(totalVertex);
    for (size_t sampleSize = 1; sampleSize <= totalVertex; sampleSize++) {
        pos.resize(sampleSize);
        for (size_t i = 0; i < pos.size(); ++i) {
            pos[i] = i;
        }
        while (true) {
            set<size_t> bh = GetBlackHoleChinese(graph, tsOrder, pos);
            if (bh.size() == 0) {
                filtered++;
            }
            if (bh.size() == bhSize && CheckConnectivity(graphUndir, bh)) {
                FoundBlackHole(bh, printDebugInfo);
            }
            if (!BruteNext(pos, totalVertex)) {
                break;
            }
        }
        if (printDebugInfo) {
            cout << "bh_filtered " << filtered << endl;
        }
    }
}

bool HasCandidates(const TClosure& closure, const TCandidates& candidates) {
    for (int v : closure) {
        if (candidates.find(v) == candidates.end()) {
            return false;
        }
    }
    return true;
}

void ChineseSolver(const TGraph& graph, size_t maxBHSize, bool printDebugInfo = false) {
    TGraph graphRev = GetReverseGraph(graph);
    TGraph graphUndir = GetUndirGraph(graph);

    if (maxBHSize == 0) {
        maxBHSize = graph.size();
    } else {
        maxBHSize = min(maxBHSize, graph.size());
    }

    vector<TCandidates> candidates(maxBHSize + 1);
    for (size_t i = 1; i <= maxBHSize; i++) {
        //cerr << "Iter: " << i << endl;
        for (size_t v = 0; v < graph.size(); v++) {
            if (printDebugInfo) {
                cout << graph[v].size() << endl;
            }
            if (graph[v].size() < i) {
                candidates[i].insert(v);
            }
        }
        Print("P", candidates[i], printDebugInfo);

        TCandidates toRemove;
        for (int v : candidates[i]) {
            if (toRemove.find(v) != toRemove.end()) {
                continue;
            }
            if (candidates[i - 1].find(v) != candidates[i - 1].end()) {
                continue;
            }
            bool remove = false;
            for (int s : graph[v]) {
                if (candidates[i].find(s) == candidates[i].end() || toRemove.find(s) != toRemove.end()) {
                    remove = true;
                    break;
                }
            }
            if (remove) {
                TClosure revClosure = GetClosure(graphRev, v);
                toRemove.insert(revClosure.begin(), revClosure.end());
            }
        }

        for (int v : toRemove) {
            auto it = candidates[i].find(v);
            if (it != candidates[i].end()) {
                candidates[i].erase(it);
            }
        }
        toRemove.clear();
        Print("C", candidates[i], printDebugInfo);
// ---------------------------------------------
        for (int v : candidates[i]) {
            if (toRemove.find(v) != toRemove.end()) {
                continue;
            }
            TClosure closure = GetClosure(graph, v);
            size_t cs = closure.size();
            if (printDebugInfo) {
                cout << "cs " << cs << endl;
            }
            if (cs >= i) {
                TClosure revClosure = GetClosure(graphRev, v);
                toRemove.insert(revClosure.begin(), revClosure.end());
                if (cs == i) {
                    FoundBlackHole(closure, printDebugInfo);
                }
            }
        }

        for (int v : toRemove) {
            auto it = candidates[i].find(v);
            if (it != candidates[i].end()) {
                candidates[i].erase(it);
            }
        }

        toRemove.clear();
        Print("F", candidates[i], printDebugInfo);

        vector<size_t> order;
        order.insert(order.end(), candidates[i].begin(), candidates[i].end());
        BruteForce(graph, graphUndir, order, i);
    }
}

int main(int argc, char** argv) {
    if (argc < 4) {
        cout << "use ./chinese_cli inputFileName binaryInput maxBHSize printDebugInfo" << endl;
        return 0;
    }
    const string inputFileName(argv[1]);
    const int binaryInput(atoi(argv[2]));
    size_t maxBHSize = atoi(argv[3]);
    if (argc > 4) {
        printDebugInfo = atoi(argv[4]);
    }

    using namespace std::chrono;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    {

        TGraph graph;
        if (binaryInput) {
            graph = ReadGraphFromBinaryFile(inputFileName);
        } else {
            graph = ReadGraphFromTextFile(inputFileName);
        }

        TGraph graphRev = GetReverseGraph(graph);

        size_t compCnt = 0;
        TComponents comp2Vertex = GetStrongConnectivityComponents(graph, graphRev, compCnt);
        TGraph graphCond = BuildCondensedGraph(compCnt, comp2Vertex, graph);

        ChineseSolver(graphCond, maxBHSize, printDebugInfo);
    }
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
    cout << "It took me: " << time_span.count() << " seconds" << endl;
    return 0;
}
