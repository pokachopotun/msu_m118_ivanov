#include "algorithm.h"

#include <graph.h>

#include <algorithm>
#include <stack>

using namespace std;

size_t dfs_recursive_with_special(const TGraph& g, vector<char>& used, vector<size_t>& order, size_t v, vector<size_t>& count, vector<char>& special) {
    used[v] = 1;
    for (size_t to : g[v]) {
        if (!used[to]) {
            count[v] += dfs_recursive_with_special(g, used, order, to, count, special);
        }
    }
    count[v] += 1;
    order.push_back(v);
    special[v] = char(order.size() == count[v]);
    return count[v];
}

vector<size_t> TopSort(const TGraph& g, const vector<size_t>& roots, vector<char>& special) {
    vector<char> used(g.size());
    vector<size_t> order;
    vector<size_t> count(g.size());
    for (size_t v : roots) {
        dfs_recursive_with_special(g, used, order, v, count, special);
    }
    return order;
}

void dfs_non_recursive(const TGraph& g, vector<char>& used, vector<size_t>& order, size_t v) {
    vector<size_t> pos(g.size());
    stack<size_t> st;
    st.push(v);
    while (!st.empty()) {
        v = st.top();
        used[v] = 1;
        if (pos[v] < g[v].size()) {
            size_t to = g[v][pos[v]];
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

vector<size_t> TopSort(const TGraph& g){
    vector<char> used(g.size());
    vector<size_t> order;
    for (size_t i = 0; i < g.size(); i++) {
        if(used[i]) {
            continue;
        }
        dfs_non_recursive(g, used, order, i);
    }
    reverse(order.begin(), order.end());
    return order;
}

TComponents GetStrongConnectivityComponents(const TGraph& graph, const TGraph& graphRev, size_t& componentCount) {
    auto&& markSingleSCC = [](const TGraph& g, vector<char>& used, TComponents& comp, size_t v, size_t compId) {
        vector<size_t> pos(g.size());
        stack<size_t> st;
        st.push(v);
        while (!st.empty()){
            v = st.top();
            used[v] = 1;
            if (pos[v] < g[v].size()) {
                size_t to = g[v][pos[v]];
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
    };

    TComponents scc(graph.size());
    vector<size_t> order = TopSort(graph);
    vector<char> used(graphRev.size());
    componentCount = 0;
    for (size_t v : order) {
        if (used[v]) {
            continue;
        }
        markSingleSCC(graphRev, used, scc, v, componentCount);
        componentCount++;
    }
    return scc;
}

bool CheckConnectivity(const TGraph& graphUndir, const set<size_t>& bh) {
    TUsed used;
    used.assign(graphUndir.size(), 1);
    for (size_t v : bh) {
        used[v] = 0;
    }
    size_t v = *bh.begin();
    TClosure closure = GetClosure(graphUndir, v, used);
    return closure.size() == bh.size();
}
