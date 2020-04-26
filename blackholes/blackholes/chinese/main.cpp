#include <iostream>
//#include <fstream>
//#include <string>
//#include <vector>
//#include <set>
//#include <queue>
//#include <stack>
//#include <algorithm>
#include <chrono>

#include <algorithm.h>
#include <blackholes.h>
#include <common.h>
#include <graph.h>
#include <io.h>


using namespace std;

//using TGraph = vector<vector<int>>;
//using TClosure = set<int>;
using TCandidates = set<int>;
//using TUsed = vector<char>;

//int bhCount = 0;
int printDebugInfo = 0;

void BruteForce(const TGraph& graph, const TGraph& graphUndir, const vector<size_t>& tsOrder, size_t bhSize) {
    size_t filtered = 0;
    const size_t totalVertex = static_cast<size_t>(tsOrder.size());
    vector<size_t> pos;
    pos.reserve(totalVertex);
    for (size_t sampleSize = 1; sampleSize <= totalVertex; sampleSize++) {
        //cerr << "BF Sample size: " << sampleSize << endl;
        pos.resize(sampleSize);
        for (size_t i = 0; i < pos.size(); ++i) {
            pos[i] = i;
        }
        while (true) {
            set<size_t> bh = GetBlackHole(graph, tsOrder, pos);
            if (bh.size() == 0) {
                filtered++;
            }
            if (bh.size() == bhSize && CheckConnectivity(graphUndir, bh)) {
                FoundBlackHole(bh.size());
            }
            if (!BruteNext(pos, totalVertex)) {
                break;
            }
        }
        cout << "filtered " << filtered << endl;
    }
}
/*
    static void Closure(const TGraph& g, int v, TUsed& used) {
        ClosureSize(g, v, used);
        return;
    }

    static void Closure(const TGraph& g, int v) {
        ClosureSize(g,v);
        return;
    } */
/*
    static int ClosureSize(const TGraph& g, int v, TUsed& used) {
        return BFS(g, used, v);
    } */
/*
    static int ClosureSize(const TGraph& g, int v) {
        TUsed used;
        used.resize(g.size());
        return BFS(g, used, v);
    } */

int DegreeOut(const TClosure& closure, const TGraph& gd) {
    long long dout = 0;
    for (int v : closure) {
        for (int to : gd[v]) {
            bool vin = closure.find(v) != closure.end();
            bool toin = closure.find(to) != closure.end();
            dout += vin && !toin;
        }
    }
    return dout;
}

bool HasCandidates(const TClosure& closure, const TCandidates& candidates) {
    for (int v : closure) {
        if (candidates.find(v) == candidates.end()) {
            return false;
        }
    }
    return true;
}
/*
    static int BFS(const TGraph& g, TUsed& used, int v) {
        queue<int> q;
        int closureSize = 0;
        if (!used[v]) {
            q.push(v);
            used[v] = 1;
            closureSize++;
        }
        while(!q.empty()) {
            v = q.front();
            q.pop();
            for (int to : g[v]) {
                if (!used[to]) {
                    used[to] = 1;
                    closureSize++;
                    q.push(to);
                }
            }
        }
        return closureSize;
    } */

void Print(const string& heading, const set<int>& closure) {
    if (!printDebugInfo) {
        return;
    }
    cout << heading << ": ";
    for (int v : closure) {
        cout << v << " ";
    }
    cout << endl;
}

int main(int argc, char** argv) {
    const string inputFileName(argv[1]);
    size_t maxBHSize = atoi(argv[2]);
    if (argc >= 4) {
        printDebugInfo = atoi(argv[3]);
    }

    using namespace std::chrono;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    {

        TGraph graph = ReadGraphFromBinaryFile(inputFileName);
        TGraph graphRev = GetReverseGraph(graph);
        TGraph graphUndir = GetUndirGraph(graph);

        // cerr << "Graph built" << endl;

        if (maxBHSize == 0) {
            maxBHSize = graph.size();
        } else {
            maxBHSize = min(maxBHSize, graph.size());
        }

        vector<TCandidates> candidates(maxBHSize + 1);
        for (size_t i = 1; i <= maxBHSize; i++) {
            cerr << "Iter: " << i << endl;
            for (size_t v = 0; v < graph.size(); v++) {
                if (printDebugInfo) {
                    cout << graph[v].size() << endl;
                }
                if (graph[v].size() < i) {
                    candidates[i].insert(v);
                }
            }
            Print("P", candidates[i]);

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
               /*
                TClosure closure = GetClosure(graph, v);
                Print("closure", closure);
                if (!HasCandidates(closure, candidates[i])) {
                    TClosure revClosure = GetClosure(graphRev, v);
                    toRemove.insert(revClosure.begin(), revClosure.end());
                } */
            }

            for (int v : toRemove) {
                auto it = candidates[i].find(v);
                if (it != candidates[i].end()) {
                    candidates[i].erase(it);
                }
            }
            toRemove.clear();
            Print("C", candidates[i]);
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
                        FoundBlackHole(cs);
                        // closure v is a blackhole
                        // Print("BH_C", closure);
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
            Print("F", candidates[i]);

            vector<size_t> order;
            order.insert(order.end(), candidates[i].begin(), candidates[i].end());
            BruteForce(graph, graphUndir, order, i);
    // ---------------------------------------------
/*            {

                TUsed used;
                used.assign(graph.size(), 1);
                for (int v : candidates[i]) {
                    used[v] = 0;
                }
                for (int v : candidates[i]) {
                    if (used[v]) {
                        continue;
                    }
                    TClosure closure = GetClosure(graphUndir, v, used);
                    size_t cs = closure.size();
                    int dout = DegreeOut(closure, graph);
                    if (cs == i && dout == 0) {
                        bhCount++;
                        FoundBlackHole(cs);
                        // restricted closure of v is blackhole
                        // Print("BH_F", closure);
                    }
                }
            }

    // ---------------------------------------------
            cerr << "Going undir" << endl;
            {
                TUsed used;
                used.assign(graph.size(), 1);
                for (int v : candidates[i]) {
                    used[v] = 0;
                }
                for (int v : candidates[i]) {
                    if (used[v]) {
                        continue;
                    }
                    TClosure closure = GetClosure(graphUndir, v, used);
                    closure.push_back(v);
                    size_t cs = closure.size();
                    int dout = DegreeOut(closure, graph, graphRev);
                    if (cs == i && dout == 0) {
                        // restricted closure of v is blackhole
                        Print("BH", closure);
                    }
                }
            } */
        }
    }

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
    cout << "It took me: " << time_span.count() << " seconds" << endl;

    return 0;
}
