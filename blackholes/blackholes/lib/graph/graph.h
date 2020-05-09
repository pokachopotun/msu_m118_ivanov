#pragma once
#include <common.h>

#include <vector>

TGraph GetReverseGraph(const TGraph& graph);
TGraph GetUndirGraph(const TGraph& graph);
TClosure GetClosure(const TGraph& g, size_t v, TUsed& used);
TClosure GetClosure(const TGraph& g, size_t v);
TClosure ClosureBFS(const TGraph& g, TUsed& used, size_t v);
std::vector<TGraph> GetWeakComponents(const TGraph& graph, const TGraph& graphUndir);
std::vector<size_t> GetRoots(const TGraph& graphRev);
TGraph BuildCondensedGraph(size_t compCnt, const TComponents& comp, const TGraph& graph);
size_t GetCompCnt(const TComponents& comp);
TSccMembers GetSccMembers(const TComponents& comp);

TClosure SimpleClosureBFS(const TGraph& g, TUsed& used, size_t s);

TGraph CarveComponent(const TGraph& graph, const TClosure& component);
