#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <chrono>
#include <fstream>
#include <stack>
#include <queue>

using namespace std;

using TUsed = vector<char>;
using TGraph = vector<vector<int>>;
using TComponents = vector<int>;
using TSccMembers = vector<vector<int>>;
using TClosure = set<int>;

class Solution {
public:
    static TGraph BuildCondensedGraph(int compCnt, const TComponents& comp, const TGraph& graph) {
        TGraph condensed(compCnt);
        for (int v = 0; v < graph.size(); v++) {
            for (int to : graph[v]) {
                if(comp[v] == comp[to]) {
                    continue;
                }
                condensed[comp[v]].push_back(comp[to]);
            }
        }
        for (auto& cond : condensed) {
            sort(cond.begin(), cond.end());
            auto it = std::unique(cond.begin(), cond.end());
            cond.resize(std::distance(cond.begin(), it));
        }
        return condensed;
    }

    static TComponents GetStrongConnectivityComponents(const TGraph& graph, const TGraph& graphRev) {
        TComponents scc(graph.size());
        vector<int> order = TopSort(graph);
        vector<char> used(graphRev.size());
        int cnt = 0;
        for (int v : order) {
            if (used[v]) {
                continue;
            }
            markSingleSCC(graphRev, used, scc, v, cnt);
            cnt++;
        }
        return scc;
    }

    static vector<int> TopSort(const TGraph& g){
        vector<char> used(g.size());
        vector<int> order;
        for (int i = 0; i < g.size(); i++) {
            if(used[i]) {
                continue;
            }
            dfs(g, used, order, i);
        }
        reverse(order.begin(), order.end());
        return order;
    }

    static vector<int> TopSort(const TGraph& g, const vector<int>& roots, vector<char>& special) {
        vector<char> used(g.size());
        vector<int> order;
        vector<int> count(g.size());
        for (int v : roots) {
            dfs_recursive(g, used, order, v, count, special);
        }
        return order;
    }

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

    static int GetCompCnt(const TComponents& comp) {
        int compCnt = 0;
        for (int c : comp) {
            compCnt = max(c, compCnt);
        }
        return compCnt + 1; 
    }

    static TSccMembers GetSccMembers(const TComponents& comp) {
        TSccMembers memb(GetCompCnt(comp));
        for (int v = 0; v < comp.size(); v++) {
            memb[comp[v]].push_back(v);
        }
        return memb;
    } 

    static int dfs_recursive(const TGraph& g, vector<char>& used, vector<int>& order, int v, vector<int>& count, vector<char>& special) {
        used[v] = 1;
        for (int to : g[v]) {
            if (!used[to]) { 
                count[v] += dfs_recursive(g, used, order, to, count, special);
            }
        }
        count[v] += 1;
        order.push_back(v);
        special[v] = char(order.size() == count[v]);
        return count[v];
    }
        
    static void dfs(const TGraph& g, vector<char>& used, vector<int>& order, int v) {
        vector<int> pos(g.size());
        stack<int> st;
        st.push(v);
        while (!st.empty()) {
            v = st.top();
            used[v] = 1;
            if (pos[v] < g[v].size()) {
                int to = g[v][pos[v]];
                pos[v]++;
                if (used[to]) {
                    continue;
                }
                st.push(to);
                continue;
            }
            st.pop();
            order.push_back(v);
        }
    }

    static void markSingleSCC(const TGraph& g, vector<char>& used, TComponents& comp, int v, int compId) {
        int depth = 0;
        vector<int> pos(g.size());
        stack<int> st;
        st.push(v);
        while (!st.empty()){
            v = st.top();
            used[v] = 1;
            if (pos[v] < g[v].size()) {
                int to = g[v][pos[v]];
                pos[v]++;
                if (used[to]) {
                    continue;
                }
                st.push(to);
                continue;
            }
            st.pop();
            comp[v] = compId;
        }
    }
    static void Print(const string& tag, const set<int>& v) {
        cout << tag << ": ";
        for (int x : v) {
            cout << x << " ";
        }
        cout << endl;
    }

    static void Print(const string& tag, const vector<int>& v) {
        cout << tag << ": ";
        for (int x : v) {
            cout << x << " ";
        }
        cout << endl;
    }

    static void Print(const string& tag, const vector<vector<int>>& vv) {
        cout << tag << " BEGIN" << endl;
        for (int i = 0; i < vv.size(); i++) {
            const auto& v = vv[i]; 
            cout << i << ": ";
            for (int x : v) {
                cout << x << " ";
            }
            cout << endl;
        }
        cout << tag + " END" << endl;
    }

    static TClosure GetClosure(const TGraph& g, int v, TUsed& used) {
        return ClosureBFS(g, used, v);
    }

    static TClosure GetClosure(const TGraph& g, int v) {
        TUsed used;
        used.resize(g.size());
        return ClosureBFS(g, used, v);
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

    static set<int> GetBlackHole(const TGraph& graph, const vector<int>& tsOrder, const vector<int>& pos) {
        set<int> blackhole;
        TUsed used;
        used.resize(graph.size());
        for (int p : pos) {
            int v = tsOrder[p];
            TClosure closure = GetClosure(graph, v, used);
            for (int x : closure) {
                blackhole.insert(x);
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
                

    static void BruteForce(const TGraph& graph, const vector<int>& tsOrder) {
        const int totalVertex = static_cast<int>(tsOrder.size());
        vector<int> pos;
        pos.reserve(totalVertex);
        for (int sampleSize = 1; sampleSize <= totalVertex; sampleSize++) {
            pos.resize(sampleSize);
            for (int i = 0; i < pos.size(); ++i) {
                pos[i] = i;
            }
            while (true) {
                set<int> bh = GetBlackHole(graph, tsOrder, pos);
                Print("BH", bh);
                if (!BruteNext(pos, totalVertex)) {
                    break;
                }
            }
        }
    }

    static void BruteForce(const TGraph& graph, const vector<int>& tsOrder, const vector<char>& special) {
        const int totalVertex = static_cast<int>(tsOrder.size());
        vector<int> pos;
        pos.reserve(totalVertex);
        for (int sampleSize = 1; sampleSize <= totalVertex; sampleSize++) {
            pos.resize(sampleSize);
            for (int i = 0; i < pos.size(); ++i) {
                pos[i] = i;
            }
            bool stop = false;
            while (!stop) {
                set<int> bh = GetBlackHole(graph, tsOrder, pos);
                Print("BH", bh);
                while (true) {
                    stop = !BruteNext(pos, totalVertex);
                    if (stop) {
                        break;
                    }
                    bool skip = false;
                    for (int i = 1; i < pos.size(); i++) {
                        int p = pos[i];
                        int v = tsOrder[p];
                        if (special[v]) {
                            skip = true;
                            break;
                        }
                    }
                    if (!skip) {
                        break;
                    }
                }
            }
        }
    }
    
    static vector<int> GetRoots(const TGraph& graphRev) {
        vector<int> roots;
        for (int v = 0; v < graphRev.size(); v++) {
            if (graphRev[v].size() == 0) {
                roots.push_back(v);
            }
        }
        return roots;
    }
};

int main(int argc, char * argv[]) {
    if (argc < 3){
            cout << "use ./main input.txt optimize";
    }
    
    const string inputFilePath(argv[1]);
    const int optimize(atoi(argv[2]));

    // #ifdef _DEBUG
        //ofstream fout(outputFilePath, std::ifstream::out);
    // #endif

    using namespace std::chrono;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    {

        TGraph graph = Solution::ReadGraphFromFile(inputFilePath);
        TGraph graphRev = Solution::GetReverseGraph(graph);

        TComponents comp2Vertex = Solution::GetStrongConnectivityComponents(graph, graphRev);
        TSccMembers compMembers = Solution::GetSccMembers(comp2Vertex);
        Solution::Print("components", compMembers);
        int compCnt = Solution::GetCompCnt(comp2Vertex); TGraph graphCond = Solution::BuildCondensedGraph(compCnt, comp2Vertex, graph); 
        if (optimize) {
            TGraph graphCondRev = Solution::GetReverseGraph(graphCond);
            vector<int> roots = Solution::GetRoots(graphCondRev);
            vector<char> special(graphCond.size());
            vector<int> tsOrder = Solution::TopSort(graphCond, roots, special);
            Solution::Print("topsort", tsOrder);
            Solution::BruteForce(graphCond, tsOrder, special);
        } else {
            vector<int> tsOrder = Solution::TopSort(graphCond);
            Solution::Print("topsort", tsOrder);
            Solution::BruteForce(graphCond, tsOrder);
        }
            
        
/*
        int total = 0;
        fout << inputFilePath << endl;
        fout << "Nodes " << n << "Edges " << m << endl;
        fout << "BlackHoles cnt: " << endl; 
        for( int i = 0 ; i < blackHoles.size(); i++ ){
                if( blackHoles[i].size() > 0 ){
                        total += blackHoles[i].size();
                        fout <<" Size: " << i << " cnt "<< blackHoles[i].size() << endl;
                }
        }
        fout << "Total: " << total << endl;
*/
    }
    high_resolution_clock::time_point t2 = high_resolution_clock::now();

    duration<double> time_span = duration_cast<duration<double>>( t2 - t1 );
    cout << "It took me: "<< time_span.count() <<  " seconds" << endl;
    return 0;
}
