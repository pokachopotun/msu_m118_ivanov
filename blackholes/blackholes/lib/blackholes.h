#pragma once
#include <common.h>

#include <cstddef>

#include <vector>
#include <set>
#include <string>

void TopsortSolver(const TGraph& graphCond, size_t ompThreads = 1, bool printDebugInfo = false);

bool CheckOutdegree(const TGraph& graph, const std::set<size_t>& bh);

void BruteForce(const TGraph& graph, const TGraph& graphUndir, const std::vector<size_t>& tsOrder);
void BruteForce(const TGraph& graph, const TGraph& graphUndir, const std::vector<size_t>& tsOrder, const std::vector<char>& special, bool printDebugInfo = false);
void BruteForceParallel(const TGraph& graph, const TGraph& graphUndir, const std::vector<size_t>& tsOrder, const std::vector<char>& special, size_t ompThreads = 1, bool printDebugInfo = false);
bool BruteNext(std::vector<size_t>& pos, size_t vertexCount);

void FoundBlackHole(const std::set<size_t>& bh, bool printDebugInfo = false);
void FilteredCandidate();
std::set<size_t> GetBlackHole(const TGraph& graph, const std::vector<size_t>& tsOrder, const std::vector<size_t>& pos);

void Print(const std::string& tag, const std::set<size_t>& v, bool printDebugInfo);
void Print(const std::string& tag, const std::vector<size_t>& v, bool printDebugInfo);
void Print(const std::string& tag, const std::vector<std::vector<size_t>>& vv, bool printDebugInfo);

