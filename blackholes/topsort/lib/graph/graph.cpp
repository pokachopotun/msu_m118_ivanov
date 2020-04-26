#include "graph.h"

#include <algorithm>
#include <map>
#include <queue>

using namespace std;

TGraph GetReverseGraph(const TGraph& graph) {
    TGraph graphRev(graph.size());
    for (size_t v = 0; v < graph.size(); v++) {
        for (size_t to : graph[v]) {
            graphRev[to].push_back(v);
        }
    }
    return graphRev;
}

TGraph GetUndirGraph(const TGraph& graph) {
    TGraph graphUndir(graph.size());
    for (size_t v = 0; v < graph.size(); v++) {
        for (size_t to : graph[v]) {
            graphUndir[to].push_back(v);
            graphUndir[v].push_back(to);
        }
    }
    return graphUndir;
}

TClosure GetClosure(const TGraph& g, size_t v, TUsed& used) {
    return ClosureBFS(g, used, v);
}

TClosure GetClosure(const TGraph& g, size_t v) {
    TUsed used;
    used.resize(g.size());
    return ClosureBFS(g, used, v);
}

TClosure ClosureBFS(const TGraph& g, TUsed& used, size_t v) {
    queue<size_t> q;
    TClosure closure;
    if (!used[v]) {
        q.push(v);
        closure.insert(v);
    }
    used[v] += 1;
    while(!q.empty()) {
        v = q.front();
        q.pop();
        for (size_t to : g[v]) {
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

TGraph CarveComponent(const TGraph& graph, const TClosure& component) {
    TGraph weak;
    weak.resize(component.size());
    map<size_t, size_t> mapping;
    size_t cnt = 0;
    for (size_t v : component) {
        mapping[v] = cnt;
        cnt++;
    }

    for (const auto& p : mapping) {
        size_t old = p.first;
        size_t neww = p.second;
        for (size_t v : graph[old]) {
            if (component.find(v) != component.end()) {
                weak[neww].push_back(mapping[v]);
            }
        }
    }
    return weak;
}

vector<TGraph> GetWeakComponents(const TGraph& graph, const TGraph& graphUndir) {
    TUsed used;
    used.resize(graph.size());
    vector<TGraph> res;
    for (size_t v = 0; v < graph.size(); ++v) {
        if (used[v]) {
            continue;
        }
        TClosure component = GetClosure(graphUndir, v, used);
        res.push_back(CarveComponent(graph, component));
    }
    return res;
}

vector<size_t> GetRoots(const TGraph& graphRev) {
    vector<size_t> roots;
    for (size_t v = 0; v < graphRev.size(); v++) {
        if (graphRev[v].size() == 0) {
            roots.push_back(v);
        }
    }
    return roots;
}

TGraph BuildCondensedGraph(size_t compCnt, const TComponents& comp, const TGraph& graph) {
    TGraph condensed(compCnt);
    for (size_t v = 0; v < graph.size(); v++) {
        for (size_t to : graph[v]) {
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

size_t GetCompCnt(const TComponents& comp) {
    size_t compCnt = 0;
    for (size_t c : comp) {
        compCnt = max(c, compCnt);
    }
    return compCnt + 1;
}

TSccMembers GetSccMembers(const TComponents& comp) {
    TSccMembers memb(GetCompCnt(comp));
    for (size_t v = 0; v < comp.size(); v++) {
        memb[comp[v]].push_back(v);
    }
    return memb;
}
