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
#include <map>

using namespace std;

using TUsed = vector<char>;
using TGraph = vector<vector<int>>;
using TComponents = vector<int>;
using TSccMembers = vector<vector<int>>;
using TClosure = set<int>;

int bhCount = 0;
int filtered = 0;
int skipped = 0;

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
        for (int i = 0; i < roots.size(); i++) {
            int v = roots[i];
            dfs_recursive(g, used, order, v, count, special);
        }
        return order;
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

 /*   static TGraph ReadGraphFromFile(const string& filename) {
        ifstream fin(filename);
        int n, m;
        fin >> n >> m;
        TGraph graph(n);

        vector<set<int>> tmp(n);
        for(int i = 0; i < m; i++) {
            int x, y;
            fin >> x >> y;
            tmp[x].insert(y);
        }
        for(int i = 0; i < tmp.size(); i++) {
            const auto& t = tmp[i];
            for (int x : t) {
                graph[i].push_back(x);
            }
        }
        return graph;
    } */

    static TGraph GetReverseGraph(const TGraph& graph) {
        TGraph graphRev(graph.size());
        for (int v = 0; v < graph.size(); v++) {
            for (int to : graph[v]) {
                graphRev[to].push_back(v);
            }
        }
        return graphRev;
    }

    static TGraph GetUndirGraph(const TGraph& graph) {
        TGraph graphUndir(graph.size());
        for (int v = 0; v < graph.size(); v++) {
            for (int to : graph[v]) {
                graphUndir[to].push_back(v);
                graphUndir[v].push_back(to);
            }
        }
        return graphUndir;
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

    static void PrintStat(const string& tag, const vector<vector<int>>& vv) {
        map<int, int> sizes;
        for (int i = 0; i < vv.size(); i++) {
            const auto& v = vv[i]; 
            sizes[v.size()]++;    
        }
        cout << tag << " BEGIN" << endl;
        for (const auto& p : sizes) {
            cout << "Size: " << p.first << " cnt " << p.second << endl;
        }
        cout << tag + " END" << endl;
    }

    static void PrintStat(const string& tag, const vector<TGraph>& v) {
        map<int, int> sizes;
        for (const auto& x : v) {
            sizes[x.size()]++;    
        }
        cout << tag << " BEGIN" << endl;
        for (const auto& p : sizes) {
            cout << "Size: " << p.first << " cnt " << p.second << endl;
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

    static void FoundBlackHole(size_t size) {
        bhCount++;
        cout << "bhcount " << bhCount << endl;
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
                return {}; // means this blackhole has already been deteected;
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
                

    static void BruteForce(const TGraph& graph, const TGraph& graphUndir, const vector<int>& tsOrder) {
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
                FoundBlackHole(bh.size());
                // cout << "BH count " << bhCount << endl;
                // Print("BH", bh);
                if (!BruteNext(pos, totalVertex)) {
                    break;
                }
            }
        }
    }

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

    static void BruteForce(const TGraph& graph, const TGraph& graphUndir, const vector<int>& tsOrder, const vector<char>& special) {
        const int totalVertex = static_cast<int>(tsOrder.size());
        const int maxSampleSize = totalVertex;
        vector<int> pos;
        pos.reserve(totalVertex);
        for (int sampleSize = 1; sampleSize <= min(maxSampleSize, totalVertex); sampleSize++) {
            int sampleSizeCount = 0;
            // cerr << "BF Sample size: " << sampleSize << endl;
            pos.resize(sampleSize);
            for (int i = 0; i < pos.size(); ++i) {
                pos[i] = i;
            }
            bool stop = false;
            while (!stop) {
                set<int> bh = GetBlackHole(graph, tsOrder, pos);
                if (bh.size() == 0) {
                    filtered++;
                }
                if (bh.size() > 0 && CheckConnectivity(graphUndir, bh)) {
                    sampleSizeCount++;
                    FoundBlackHole(bh.size());
                    // cout << "BH count " << bhCount << endl;
                    // Print("BH", bh);
                }
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
                            skipped += 1;
                            break;
                        }
                    }
                    if (!skip) {
                        break;
                    }
                }
            }

            if (sampleSizeCount == 0) {
                break;
            }

            cout << "bh_count " << bhCount << endl;
            cout << "bh_filtered " << filtered << endl;
            cout << "bh_skipped " << skipped << endl;
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

    static TGraph CarveComponent(const TGraph& graph, const TClosure& component) {
        TGraph weak;
        weak.resize(component.size());    
        map<int, int> mapping;
        int cnt = 0;
        for (int v : component) {
            mapping[v] = cnt;
            cnt++;
        }

        for (const auto& p : mapping) {
            int old = p.first;
            int neww = p.second;
            for (int v : graph[old]) {
                if (component.find(v) != component.end()) {
                    weak[neww].push_back(mapping[v]);
                }
            }
        }
        return weak;
    }

        

    static vector<TGraph> GetWeakComponents(const TGraph& graph, const TGraph& graphUndir) {
        TUsed used;
        used.resize(graph.size());
        vector<TGraph> res;
        for (int v = 0; v < graph.size(); ++v) {
            if (used[v]) {
                continue;
            }
            TClosure component = GetClosure(graphUndir, v, used);
            res.push_back(CarveComponent(graph, component));
        }
        return res;
    }


    static void Solve(const TGraph& graphCond) {
        TGraph graphCondRev = Solution::GetReverseGraph(graphCond);
        TGraph graphCondUndir = Solution::GetUndirGraph(graphCond);
        vector<int> roots = Solution::GetRoots(graphCondRev);
        int maxid = 0;
        int maxsize = 0;
        for (int i = 0; i < roots.size(); i++) {
            TClosure closure = GetClosure(graphCond, roots[i]);
            if (closure.size() > maxsize) {
                maxid = i;
                maxsize = closure.size();
            }
        }
        swap(roots[0], roots[maxid]);
        // Solution::Print("roots", roots);
        vector<char> special(graphCond.size());
        vector<int> tsOrder = Solution::TopSort(graphCond, roots, special);
//           Solution::Print("topsort", tsOrder);
        Solution::BruteForce(graphCond, graphCondUndir, tsOrder, special);
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
        // Solution::Print("components", compMembers);
        Solution::PrintStat("SCC component sizes", compMembers);
        int compCnt = Solution::GetCompCnt(comp2Vertex); 
        TGraph graphCond = Solution::BuildCondensedGraph(compCnt, comp2Vertex, graph);
        if (optimize) {
            TGraph graphCondUndir = Solution::GetUndirGraph(graphCond);
            vector<TGraph> weaks = Solution::GetWeakComponents(graphCond, graphCondUndir);
            //Solution::PrintStat("Weakly connected component sizes", weaks);
            //if (false) 
            for (const TGraph& weak : weaks) {
                Solution::Solve(weak);
            }
        } else {
            TGraph graphCondUndir = Solution::GetUndirGraph(graphCond);
            vector<int> tsOrder = Solution::TopSort(graphCond);
  //          Solution::Print("topsort", tsOrder);
            Solution::BruteForce(graphCond, graphCondUndir, tsOrder);
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
