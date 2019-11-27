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

using namespace std;

using TGraph = vector<vector<int>>;
using TComponents = vector<int>;
using TSccMembers = vector<vector<int>>;

class Solution {
public:
    static TGraph BuildCondensedGraph(int compCnt, const TComponents& comp, const TGraph& graph) {
        TGraph condensed(compCnt, vector<int>());
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
        
private:
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
        while(!st.empty()){
            v = st.top();
            used[v] = 1;
            if(pos[v] < g[v].size()){
                int to = g[v][pos[v]];
                pos[v]++;
                if(used[to]){
                    continue;
                }
                st.push( to );
                continue;
            }
            st.pop();
            comp[v] = compId;
        }
    }
};

int main(int argc, char * argv[]) {
        if( argc < 3){
                cout << "use ./main input.txt output.txt";
        }
        
        const string inputFilePath(argv[1]);
        const string outputFilePath(argv[2]);

// #ifdef _DEBUG
        ofstream fout(outputFilePath, std::ifstream::out);
// #endif

        using namespace std::chrono;
        high_resolution_clock::time_point t1 = high_resolution_clock::now();
        {

            TGraph graph = Solution::ReadGraphFromFile(inputFilePath);
            TGraph graphRev = Solution::GetReverseGraph(graph);

            TComponents comp2Vertex = Solution::GetStrongConnectivityComponents(graph, graphRev);
            TSccMembers compMembers = Solution::GetSccMembers(comp2Vertex);
            int compCnt = Solution::GetCompCnt(comp2Vertex);
            TGraph graphCond = Solution::BuildCondensedGraph(compCnt, comp2Vertex, graph);

            vector<int> tsOrder = Solution::TopSort(graphCond);
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
