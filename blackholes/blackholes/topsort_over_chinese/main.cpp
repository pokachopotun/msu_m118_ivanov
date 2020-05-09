#include <omp.h>
#include <iostream>
#include <queue>
#include <chrono>

#include <algorithm.h>
#include <blackholes.h>
#include <common.h>
#include <graph.h>
#include <io.h>

#include <iostream>
#include <set>
#include <vector>

using namespace std;

using TCandidates = set<size_t>;

void ChineseSolveSingleSize(size_t i,
                            const TGraph& graph,
                            size_t ompThreads = 1,
                            bool printDebugInfo = false)
{
    TGraph graphRev = GetReverseGraph(graph);
    TGraph graphUndir = GetUndirGraph(graph);
    TCandidates curCandidates;
    //cerr << "Iter: " << i << endl;
    for (size_t v = 0; v < graph.size(); v++) {
        if (printDebugInfo) {
            cout << graph[v].size() << endl;
        }
        if (graph[v].size() < i) {
            curCandidates.insert(v);
        }
    }
    Print("P", curCandidates, printDebugInfo);

    TCandidates toRemove;
    for (int v : curCandidates) {
        if (toRemove.find(v) != toRemove.end()) {
            continue;
        }
        bool remove = false;
        for (int s : graph[v]) {
            if (curCandidates.find(s) == curCandidates.end() || toRemove.find(s) != toRemove.end()) {
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
        auto it = curCandidates.find(v);
        if (it != curCandidates.end()) {
            curCandidates.erase(it);
        }
    }
    toRemove.clear();
    Print("C", curCandidates, printDebugInfo);
// ---------------------------------------------
    for (int v : curCandidates) {
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
        auto it = curCandidates.find(v);
        if (it != curCandidates.end()) {
            curCandidates.erase(it);
        }
    }

    toRemove.clear();
    Print("F", curCandidates, printDebugInfo);

    {
        TUsed used(graph.size(), 1);
        for (size_t x : curCandidates) {
            used[x] = 0;
        }
        for (size_t v : curCandidates) {
            if (used[v]) {
                continue;
            }
            TClosure closure = SimpleClosureBFS(graphUndir, used, v);
            TGraph wccGraph = CarveComponent(graph, closure);
            TopsortSolver(wccGraph, ompThreads, printDebugInfo);
        }
    }
}

int main(int argc, char** argv) {
    if (argc < 4) {
        cout << "use ./topsort_over_chinese_cli inputFileName binaryInput maxBHSize useCondensation ompThreads printDebugInfo" << endl;
        return 0;
    }
    const string inputFileName(argv[1]);
    const int binaryInput(atoi(argv[2]));
    size_t maxBHSize = atoi(argv[3]);
    bool useCondensation = atoi(argv[4]);
    size_t ompThreads = 1;
    bool printDebugInfo = false;
    if (argc > 5) {
        ompThreads = atoi(argv[5]);
    }
    if (argc > 6) {
        printDebugInfo = atoi(argv[6]);
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

        if (useCondensation) {
            TGraph graphRev = GetReverseGraph(graph);
            size_t compCnt = 0;
            TComponents comp2Vertex = GetStrongConnectivityComponents(graph, graphRev, compCnt);
            TGraph graphCond = BuildCondensedGraph(compCnt, comp2Vertex, graph);
            ChineseSolveSingleSize(maxBHSize, graphCond, ompThreads, printDebugInfo);
        } else {
            ChineseSolveSingleSize(maxBHSize, graph, ompThreads, printDebugInfo);
        }
    }
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
    cout << "It took me: " << time_span.count() << " seconds" << endl;
    return 0;
}
