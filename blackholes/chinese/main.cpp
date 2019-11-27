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
 
class Solution {
public:

    static TGraph ReadGraphFromFile(const string& filename) {
        ifstream fin(filename);
        int n, m;
        fin >> n >> m;
        TGraph graph(n);
        for(int i = 0; i < m; i++) {
            int x, y;
            fin >> x >> y;
            graph[x].push_back(y);
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
            used[v] = 1;
            closure.insert(v);
        }
        while(!q.empty()) {
            v = q.front();
            q.pop();
            for (int to : g[v]) {
                if (!used[to]) {
                    used[to] = 1;
                    q.push(to);
                    closure.insert(to);
                }
            }
        }
        return closure;
    }

};

void Print(const string& heading, const set<int>& closure) {
    cout << heading << ": ";
    for (int v : closure) {
        cout << v << " ";
    }
    cout << endl;
}

int main(int argc, char** argv) { 
    const string inputFileName(argv[1]);
    int maxBHSize = atoi(argv[2]);

    using namespace std::chrono;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    {

        TGraph graph = Solution::ReadGraphFromFile(inputFileName);
        TGraph graphRev = Solution::GetReverseGraph(graph);
        TGraph graphUndir = Solution::GetUndirectedGraph(graph);

        cerr << "Graph built" << endl;

        if (maxBHSize == -1) {
            maxBHSize = graph.size();
        } else {
            maxBHSize = min(maxBHSize, static_cast<int>(graph.size()));
        }
 
        vector<TCandidates> candidates(maxBHSize + 1);
        for (int i = 1; i <= maxBHSize; i++) {
            cerr << "Iter: " << i << endl;
            for (int v = 0; v < graph.size(); v++) {
                if (graph[v].size() < i) {
                    candidates[i].insert(v);
                }
            }
            Print("P", candidates[i]);
            
            TCandidates toRemove;
            for (int v : candidates[i]) {
                if (candidates[i - 1].find(v) != candidates[i - 1].end()) {
                    continue;
                }
                TClosure closure = Solution::GetClosure(graph, v);
                if (!Solution::HasCandidates(closure, candidates[i])) {
                    TClosure revClosure = Solution::GetClosure(graphRev, v);
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
            Print("C", candidates[i]);
    // ---------------------------------------------
            for (int v : candidates[i]) {
                TClosure closure = Solution::GetClosure(graph, v);
                size_t cs = closure.size();
                if (cs >= i) {
                    TClosure revClosure = Solution::GetClosure(graphRev, v);
                    toRemove.insert(revClosure.begin(), revClosure.end());
                    if (cs == i) {
                        // closure v is a blackhole
                        Print("BH_C", closure);
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
    // ---------------------------------------------
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
                    size_t cs = closure.size();
                    int dout = Solution::DegreeOut(closure, graph);
                    if (cs == i && dout == 0) {
                        // restricted closure of v is blackhole
                        Print("BH_F", closure);                    
                    }
                }
            }
/*
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
