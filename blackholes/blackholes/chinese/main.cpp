#include <omp.h>
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

void BruteForce(const TGraph& graph, const TGraph& graphUndir, const vector<size_t>& tsOrder, size_t bhSize, bool printDebugInfo = false) {
    const size_t totalVertex = static_cast<size_t>(tsOrder.size());
    vector<size_t> pos;
    pos.reserve(bhSize);
    for (size_t sampleSize = bhSize; sampleSize <= min(bhSize, totalVertex); sampleSize++) {
        pos.resize(sampleSize);
        for (size_t i = 0; i < pos.size(); ++i) {
            pos[i] = i;
        }
        while (true) {
            set<size_t> bh;
            for (size_t p : pos) {
                bh.insert(tsOrder[p]);
            }
            if (bh.size() == bhSize && CheckConnectivity(graphUndir, bh) && CheckOutdegree(graph, bh)) {
                FoundBlackHole(bh, printDebugInfo);
            } else {
                FilteredCandidate();
            }
            if (!BruteNext(pos, totalVertex)) {
                break;
            }
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

void ChineseSolveSingleSize(size_t i,
                            const TGraph& graph,
                            const TGraph graphRev,
                            const TGraph& graphUndir,
                            const TCandidates& prevCandidates,
                            TCandidates& curCandidates,
                            bool printDebugInfo = false)
{
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
        if (prevCandidates.find(v) != prevCandidates.end()) {
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

    vector<size_t> order;
    order.insert(order.end(), curCandidates.begin(), curCandidates.end());
    BruteForce(graph, graphUndir, order, i);
}

void ChineseSolveAllSizesAsync(const TGraph& graph,
                               const TGraph graphRev,
                               const TGraph& graphUndir,
                               vector<TCandidates>& candidates,
                               size_t maxBHSize,
                               size_t numThreads,
                               bool printDebugInfo = false)
{
    TCandidates prevDummy;
    #pragma omp parallel for num_threads(numThreads)
    for (size_t i = 1; i <= maxBHSize; i++) {
        ChineseSolveSingleSize(i, graph, graphRev, graphUndir, prevDummy, candidates[i], printDebugInfo);
    }
}

void ChineseSolveAllSizes(const TGraph& graph,
                          const TGraph graphRev,
                          const TGraph& graphUndir,
                          vector<TCandidates>& candidates,
                          size_t maxBHSize,
                          bool printDebugInfo = false)
{
    for (size_t i = 1; i <= maxBHSize; i++) {
        ChineseSolveSingleSize(i, graph, graphRev, graphUndir, candidates[i - 1], candidates[i], printDebugInfo);
    }
}

void ChineseSolver(const TGraph& graph,
                   size_t maxBHSize,
                   size_t numThreads = 1,
                   bool printDebugInfo = false)
{
    TGraph graphRev = GetReverseGraph(graph);
    TGraph graphUndir = GetUndirGraph(graph);

    if (maxBHSize == 0) {
        maxBHSize = graph.size();
    } else {
        maxBHSize = min(maxBHSize, graph.size());
    }

    vector<TCandidates> candidates(maxBHSize + 1);
    if (numThreads > 1) {
        ChineseSolveAllSizesAsync(graph, graphRev, graphUndir, candidates, maxBHSize, numThreads, printDebugInfo);
    } else {
        ChineseSolveAllSizes(graph, graphRev, graphUndir, candidates, maxBHSize, printDebugInfo);
    }
}

int main(int argc, char** argv) {
    if (argc < 4) {
        cout << "use ./chinese_cli inputFileName binaryInput maxBHSize useCondensation ompThreads printDebugInfo" << endl;
        return 0;
    }
    const string inputFileName(argv[1]);
    const int binaryInput(atoi(argv[2]));
    size_t maxBHSize = atoi(argv[3]);
    bool useCondensation = atoi(argv[4]);
    size_t ompThreads = 1;
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
            ChineseSolver(graphCond, maxBHSize, ompThreads, printDebugInfo);
        } else {
            ChineseSolver(graph, maxBHSize, ompThreads, printDebugInfo);
        }
    }
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
    cout << "It took me: " << time_span.count() << " seconds" << endl;
    return 0;
}
