#pragma once

#include <common.h>

#include <vector>
#include <set>

std::vector<size_t> TopSort(const TGraph& g, const std::vector<size_t>& roots, std::vector<char>& special);
std::vector<size_t> TopSort(const TGraph& g);
TComponents GetStrongConnectivityComponents(const TGraph& graph, const TGraph& graphRev, size_t& componentCount);
bool CheckConnectivity(const TGraph& graphUndir, const std::set<size_t>& bh);
