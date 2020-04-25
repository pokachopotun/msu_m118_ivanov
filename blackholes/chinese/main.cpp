#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <queue>
#include <stack>
#include <algorithm>
#include <chrono>


using namespace std;

using TGraph = vector<vector<int>>;
using TClosure = set<int>;
using TCandidates = set<int>;
using TUsed = vector<char>;

int bhCount = 0;
int filtered = 0;

int printDebugInfo = 0;

class Solution {
public:

    static bool CheckConnectivity(const TGraph& graphUndir, const set<int>& bh) {
        TUsed used;
        used.assign(graphUndir.size(), 1);
        for (int v : bh) {
            used[v] = 0;
        }
        int v = *bh.begin();
        TClosure closure = GetClosure(graphUndir, v, used);
        return closure.size() == bh.size();
    }

    static set<int> GetBlackHole(const TGraph& graph, const vector<int>& tsOrder, const vector<int>& pos) {
        set<int> blackhole;
        TUsed used;
        used.assign(graph.size(), 0);
        for (int p : pos) {
            int v = tsOrder[p];
            TClosure closure = GetClosure(graph, v, used);
            for (int x : closure) {
                blackhole.insert(x);
            }
        }
        for (int p : pos) {
            int x = tsOrder[p];
            if (used[x] >= 2) {
                return {}; // means this blackhole has already been detected;
            }
        }
        return blackhole;
    }

    static bool BruteNext(vector<int>& pos, int vertexCount) {
        int sampleSize = pos.size();
        bool res = true;
        for (int i = sampleSize - 1; i >= 0; --i) {
            int& cur = pos[i];
            int distFromEnd = sampleSize - 1 - i;
            if (cur == vertexCount - 1 - distFromEnd) {
                if (i == 0) {
                    res = false;
                }
                continue;
            } else {
                cur++;
                for (int j = i + 1; j < sampleSize; j++) {
                    pos[j] = pos[j - 1] + 1;
                }
                break;
            }
        }
        return res;
    }

    static void BruteForce(const TGraph& graph, const TGraph& graphUndir, const vector<int>& tsOrder, int bhSize) {
        const int totalVertex = static_cast<int>(tsOrder.size());
        vector<int> pos;
        pos.reserve(totalVertex);
        for (int sampleSize = 1; sampleSize <= totalVertex; sampleSize++) {
            //cerr << "BF Sample size: " << sampleSize << endl;
            pos.resize(sampleSize);
            for (int i = 0; i < pos.size(); ++i) {
                pos[i] = i;
            }
            while (true) {
                set<int> bh = GetBlackHole(graph, tsOrder, pos);
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
            cout << "bhcount " << bhCount << endl;
            cout << "filtered " << filtered << endl;
        }
    }

    static void FoundBlackHole(size_t size) {
        bhCount++;
        cout << "bhcount " << bhCount << endl;
    }

    static TGraph ReadGraphFromFile(const string& filename) {
        fstream file(filename, ios::in | ios::binary);
        int n;
        long long m;
        // read header
        file.read((char*)(&n), sizeof(int));
        file.read((char*)(&m), sizeof(long long));

        TGraph graph(n);

        vector<set<int>> tmp(n);
        for(long long i = 0; i < m; i++) {
            int x, y;
            file.read((char*)(&x), sizeof(int));
            file.read((char*)(&y), sizeof(int));
            if (x == y) {
                continue;
            }
            tmp[x].insert(y);
        }
        for(int i = 0; i < tmp.size(); i++) {
            const auto& t = tmp[i];
            for (int x : t) {
                graph[i].push_back(x);
            }
        }
        return graph;
    }

    static TGraph GetReverseGraph(const TGraph& graph) {
        TGraph graphRev(graph.size());
        for (int v = 0; v < graph.size(); v++) {
            for (int to : graph[v]) {
                graphRev[to].push_back(v);
            }
        }
        return graphRev;
    }

    static TGraph GetUndirectedGraph(const TGraph& graph) {
        TGraph graphUndir(graph.size());
        for (int v = 0; v < graph.size(); v++) {
            for (int to : graph[v]) {
                graphUndir[to].push_back(v);
                graphUndir[v].push_back(to);
            }
        }
        return graphUndir;
    }

    static void Closure(const TGraph& g, int v, TUsed& used) {
        ClosureSize(g, v, used);
        return;
    }

    static void Closure(const TGraph& g, int v) {
        ClosureSize(g,v);
        return;
    }

    static int ClosureSize(const TGraph& g, int v, TUsed& used) {
        return BFS(g, used, v);
    }

    static int ClosureSize(const TGraph& g, int v) {
        TUsed used;
        used.resize(g.size());
        return BFS(g, used, v);
    }

    static TClosure GetClosure(const TGraph& g, int v, TUsed& used) {
        return ClosureBFS(g, used, v);
    }

    static TClosure GetClosure(const TGraph& g, int v) {
        TUsed used;
        used.resize(g.size());
        return ClosureBFS(g, used, v);
    }

    static int DegreeOut(const TClosure& closure, const TGraph& gd) {
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

    static bool HasCandidates(const TClosure& closure, const TCandidates& candidates) {
        for (int v : closure) {
            if (candidates.find(v) == candidates.end()) {
                return false;
            }
        }
        return true;
    }

private:
    static void dfs_non_recursive(const TGraph& g, vector<char>& used, vector<int>& order, int v) { // shitty code. write recursive or refactor
        vector<int> pos(g.size());
        stack<int> st;
        st.push(v);
        while (!st.empty()) {
            v = st.top();
            used[v] = 1;
            if (pos[v] < g[v].size()) {
                int to = g[v][pos[v]];
                pos[v]++;
                if (used[to] > 0) {
                    continue;
                }
                st.push(to);
                continue;
            }
            st.pop();
            order.push_back(v);
        }
    }

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
    }

    static TClosure ClosureBFS(const TGraph& g, TUsed& used, int v) {
        queue<int> q;
        TClosure closure;
        int closureSize = 0;
        if (!used[v]) {
            q.push(v);
            closure.insert(v);
        }
        used[v] += 1;
        while(!q.empty()) {
            v = q.front();
            q.pop();
            for (int to : g[v]) {
                if (!used[to]) {
                    used[to] += 1;
                    q.push(to);
                    closure.insert(to);
                } else {
                    used[to] += 1;
                }
            }
        }
        return closure;
    }

};

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
    int maxBHSize = atoi(argv[2]);
    if (argc >= 4) {
        printDebugInfo = atoi(argv[3]);
    }

    using namespace std::chrono;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    {

        TGraph graph = Solution::ReadGraphFromFile(inputFileName);
        TGraph graphRev = Solution::GetReverseGraph(graph);
        TGraph graphUndir = Solution::GetUndirectedGraph(graph);

        // cerr << "Graph built" << endl;

        if (maxBHSize == -1) {
            maxBHSize = graph.size();
        } else {
            maxBHSize = min(maxBHSize, static_cast<int>(graph.size()));
        }

        vector<TCandidates> candidates(maxBHSize + 1);
        for (int i = 1; i <= maxBHSize; i++) {
            cerr << "Iter: " << i << endl;
            for (int v = 0; v < graph.size(); v++) {
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
                    TClosure revClosure = Solution::GetClosure(graphRev, v);
                    toRemove.insert(revClosure.begin(), revClosure.end());
                }
               /* 
                TClosure closure = Solution::GetClosure(graph, v);
                Print("closure", closure);
                if (!Solution::HasCandidates(closure, candidates[i])) {
                    TClosure revClosure = Solution::GetClosure(graphRev, v);
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
                TClosure closure = Solution::GetClosure(graph, v);
                size_t cs = closure.size();
                if (printDebugInfo) {
                    cout << "cs " << cs << endl;
                }
                if (cs >= i) {
                    TClosure revClosure = Solution::GetClosure(graphRev, v);
                    toRemove.insert(revClosure.begin(), revClosure.end());
                    if (cs == i) {
                        Solution::FoundBlackHole(cs);
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

            vector<int> order;
            order.insert(order.end(), candidates[i].begin(), candidates[i].end());
            Solution::BruteForce(graph, graphUndir, order, i);
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
                    TClosure closure = Solution::GetClosure(graphUndir, v, used);
                    size_t cs = closure.size();
                    int dout = Solution::DegreeOut(closure, graph);
                    if (cs == i && dout == 0) {
                        bhCount++;
                        Solution::FoundBlackHole(cs);
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
                    TClosure closure = Solution::GetClosure(graphUndir, v, used);
                    closure.push_back(v);
                    size_t cs = closure.size();
                    int dout = Solution::DegreeOut(closure, graph, graphRev);
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
