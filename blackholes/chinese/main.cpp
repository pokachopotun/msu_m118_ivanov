#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <queue>
#include <stack>


using namespace std;

using TGraph = vector<vector<int>>;
using TClosure = vector<int>;
using TCandidates = set<int>;
using TUsed = vector<char>;
 
class Solution {
public:

    static TGraph ReadGraphFromFile(const string& filename) {
        ifstream fin(filename, std::ifstream::in);
        int n, m;
        TGraph graph(n);
        fin >> n >> m;
        for(int i = 0; i < m; i++) {
            int x, y;
            cin >> x >> y;
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

    static int DegreeOut(const TClosure& closure, const TGraph& gd, const TGraph& gr) {
        long long dout = 0;
        for (int v : closure) {
           dout += static_cast<long long>(gd[v].size()) - static_cast<long long>(gr[v].size());
        }
        return dout;
    } 

    static bool HasCandidates(const TClosure& closure, const TCandidates& candidates) {
        for (int v : closure) {
            if (candidates.find(v) != candidates.end()) {
                return true;
            }
        }
        return false;
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
            closure.push_back(v);
        }
        while(!q.empty()) {
            v = q.front();
            q.pop();
            for (int to : g[v]) {
                if (!used[to]) {
                    used[to] = 1;
                    q.push(to);
                    closure.push_back(v);
                }
            }
        }
        return closure;
    }

};

int main(int argc, char** argv) { 
    const string inputFileName(argv[1]);
    const int maxBHSize = atoi(argv[2]);

    {

        int maxBHSize;

        TGraph graph = Solution::ReadGraphFromFile(inputFileName);
        TGraph graphRev = Solution::GetReverseGraph(graph);
        TGraph graphUndir = Solution::GetUndirectedGraph(graph);
        
        vector<TCandidates> candidates(maxBHSize + 1);
        for (int i = 1; i <= maxBHSize; i++) {
            for (int v = 0; v < graph.size(); v++) {
                if (graph[v].size() < i) {
                    candidates[i].insert(v);
                }
            }

            TCandidates toRemove;
            for (int v : candidates[i]) {
                if (candidates[i - 1].find(v) != candidates[i - 1].end()) {
                    continue;
                }
                TClosure closure = Solution::GetClosure(graph, v);
                if (Solution::HasCandidates(closure, candidates[i])) {
                    TClosure revClosure = Solution::GetClosure(graphRev, v);
                    for (int x : revClosure) {
                        toRemove.insert(x);
                    }
                } 
            }

            for (int v : toRemove) {
                auto it = candidates[i].find(v);
                if (it != candidates[i].end()) {
                    candidates[i].erase(v);
                }
            }
            toRemove.clear();
    // ---------------------------------------------
            for (int v : candidates[i]) {
                int cs = Solution::ClosureSize(graph, v);
                if (cs >= i) {
                    if (cs == i) {
                        // closure v is a blackhole
                    }
                    TClosure closure = Solution::GetClosure(graphRev, v);
                    toRemove.insert(closure.begin(), closure.end());
                }
            }

            for (int v : toRemove) {
                auto it = candidates[i].find(v);
                if (it != candidates[i].end()) {
                    candidates[i].erase(v);
                }
            }
            toRemove.clear();
    // ---------------------------------------------
            TUsed used;
            used.resize(graph.size());
            for (int v = 0; v < candidates[i].size(); v++) {
                if (used[v]) {
                    continue;
                }
                TClosure closure = Solution::GetClosure(graphUndir, v, used);
                size_t cs = closure.size();
                int dout = Solution::DegreeOut(closure, graph, graphRev);
                if (cs == i && dout == 0) {
                    // restricted closure of v is blackhole
                }
            }
        }
    }

    return 0;
}
