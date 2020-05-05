#pragma once
#include <cstddef>
#include <vector>
#include <set>

using TUsed = std::vector<int>;
using TGraph = std::vector<std::vector<size_t>>;
using TComponents = std::vector<size_t>;
using TSccMembers = std::vector<std::vector<size_t>>;
using TClosure = std::set<size_t>;
