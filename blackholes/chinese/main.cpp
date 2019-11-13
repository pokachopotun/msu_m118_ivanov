#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <queue>
#include <stack>


using namespace std;

using Graph = vector<vector<int>>;

using Candidates = set<int>;
using Potential = set<int>;
using Final = set<int>;

class Solution {
public:

    static Graph ReadGraphFromFile(const string& filename) {
        ifstream fin(filename, std::ifstream::in);
        int n, m;
        Graph graph(n);
        fin >> n >> m;
        for(int i = 0; i < m; i++) {
            int x, y;
            cin >> x >> y;
            graph[x].push_back(y);
        }
        return graph;
    }

    static Graph GetReverseGraph(const Graph& graph) {
        Graph graphRev(graph.size());
        for (int v = 0; v < graph.size(); v++) {
            for (int to : graph[v]) {
                graphRev[to].push_back(v);
            }
        }
        return graphRev;
    }

    static Graph GetUndirectedGraph(const Graph& graph) {
        Graph graphUndir(graph.size());
        for (int v = 0; v < graph.size(); v++) {
            for (int to : graph[v]) {
                graphUndir[to].push_back(v);
                graphUndir[v].push_back(to);
            }
        }
        return graphUndir;
    }


    static bool checkSuccessors(const Graph& g, const Candidates& cand, int v) {
        vector<char> used(g.size());
        queue<int> q;
        q.push(v);
        used[v] = 1;
        while(!q.empty()) {
            v = q.front();
            q.pop();
            if (cand.find(v) == cand.end()) {
               return false;
            } 
            for (int to : g[v]) {
                if (!used[to]) {
                    used[to] = 1;
                    q.push(to);
                }
            }
        }
        return true;
    }

    static int closureSize(const Graph& g, int v) {
        vector<char> used(g.size());
        queue<int> q;
        q.push(v);
        used[v] = 1;
        int cnt = 0;
        while(!q.empty()) {
            v = q.front();
            q.pop();
            cnt++;
            for (int to : g[v]) {
                if (!used[to]) {
                    used[to] = 1;
                    q.push(to);
                }
            }
        }
        return cnt;
    }

    static void removePredecessors(const Graph& g, Candidates& toRemove, int v) {
        vector<char> used(g.size());
        queue<int> q;
        q.push(v);
        used[v] = 1;
        while(!q.empty()) {
            v = q.front();
            q.pop();
            toRemove.insert(v);
            for (int to : g[v]) {
                if (!used[to]) {
                    used[to] = 1;
                    q.push(to);
                }
            }
        }
    }

    static int restrictedClosureSize(const Graph& g, const Candidates& cand, int v, vector<char>& used) {
        queue<int> q;
        q.push(v);
        used[v] = 1;
        int cnt = 0;
        while(!q.empty()) {
            v = q.front();
            q.pop();
            cnt++;
            for (int to : g[v]) {
                if (cand.find(to) != cand.end()) {
                    used[to] = 1;
                }
                if (!used[to]) {
                    used[to] = 1;
                    q.push(to);
                }
            }
        }
        return cnt;
    }
        
    static int degreeOut(const Graph& g, const Candidates& cand, int v, vector<char>& used, const Graph& gd, const Graph& gr) {
        queue<int> q;
        q.push(v);
        used[v] = 1;
        int dout = 0;
        while(!q.empty()) {
            v = q.front();
            q.pop();
            dout += gd[v].size() - gr[v].size();
            for (int to : g[v]) {
                if (cand.find(to) != cand.end()) {
                    used[to] = 1;
                }
                if (!used[to]) {
                    used[to] = 1;
                    q.push(to);
                }
            }
        }
        return dout;
    } 

private:
    static void dfs_non_recursive(const Graph& g, vector<char>& used, vector<int>& order, int v) { // shitty code. write recursive or refactor
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
};

int main(int argc, char** argv) { 
    const string inputFileName(argv[1]);
    const int maxBHSize = atoi(argv[2]);

    {

        int maxBHSize;

        Graph graph = Solution::ReadGraphFromFile(inputFileName);
        Graph graphRev = Solution::GetReverseGraph(graph);
        Graph graphUndir = Solution::GetUndirectedGraph(graph);
        
        vector<Candidates> candidates(maxBHSize + 1);
        for (int i = 1; i <= maxBHSize; i++) {
            for (int v = 0; v < graph.size(); v++) {
                if (graph[v].size() < i) {
                    candidates[i].insert(v);
                }
            }

            Candidates toRemove;
            for (int v : candidates[i]) {
                if (candidates[i - 1].find(v) != candidates[i - 1].end()) {
                    continue;
                }
                if (!Solution::checkSuccessors(graph, candidates[i], v)) {
                    toRemove.insert(v);
                    Solution::removePredecessors(graph, toRemove, v);
                } 
            }

            for (int v : toRemove) {
                auto it = candidates[i].find(v);
                if (it != candidates[i].end()) {
                    candidates[i].erase(it);
                }
            }
            toRemove.clear();
    // ---------------------------------------------
            for (int v : candidates[i]) {
                int cs = Solution::closureSize(graph, v);
                if (cs >= i) {
                    if (cs == i) {
                        // closure v is a blackhole
                    }
                    toRemove.insert(v);
                    Solution::removePredecessors(graphRev, toRemove, v);
                }
            }

            for (int v : toRemove) {
                auto it = candidates[i].find(v);
                if (it != candidates[i].end()) {
                    candidates[i].erase(it);
                }
            }
            toRemove.clear();

            vector<char> used;
            for (int v : candidates[i]) {
                if (used[v]) {
                    continue;
                }
                int cs = Solution::restrictedClosureSize(graphUndir, candidates[i], v, used);
                int dout = Solution::degreeOut(graphUndir, candidates[i], v, used, graph, graphRev);
                if (cs == i && dout == 0) {
                    // restricted closure of v is blackhole
                }
            }
        }
    }

    return 0;
}
