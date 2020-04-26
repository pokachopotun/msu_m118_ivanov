#pragma once

#include <common.h>

#include <cstddef>

#include <vector>
#include <set>
#include <string>

void Solve(const TGraph& graphCond);

void BruteForce(const TGraph& graph, const TGraph& graphUndir, const std::vector<size_t>& tsOrder);
void BruteForce(const TGraph& graph, const TGraph& graphUndir, const std::vector<size_t>& tsOrder, const std::vector<char>& special);
bool BruteNext(std::vector<size_t>& pos, size_t vertexCount);

void FoundBlackHole(size_t size);
std::set<size_t> GetBlackHole(const TGraph& graph, const std::vector<size_t>& tsOrder, const std::vector<size_t>& pos);

void Print(const std::string& tag, const std::set<size_t>& v);
void Print(const std::string& tag, const std::vector<size_t>& v);
void Print(const std::string& tag, const std::vector<std::vector<size_t>>& vv);
