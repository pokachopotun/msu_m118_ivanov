#include "blackholes.h"

#include <algorithm.h>
#include <graph.h>

#include <iostream>

using namespace std;

void Solve(const TGraph& graphCond) {
    TGraph graphCondRev = GetReverseGraph(graphCond);
    TGraph graphCondUndir = GetUndirGraph(graphCond);
    vector<size_t> roots = GetRoots(graphCondRev);
    size_t maxid = 0;
    size_t maxsize = 0;
    for (size_t i = 0; i < roots.size(); i++) {
        TClosure closure = GetClosure(graphCond, roots[i]);
        if (closure.size() > maxsize) {
            maxid = i;
            maxsize = closure.size();
        }
    }
    swap(roots[0], roots[maxid]);
    // Solution::Print("roots", roots);
    vector<char> special(graphCond.size());
    vector<size_t> tsOrder = TopSort(graphCond, roots, special);
//           Solution::Print("topsort", tsOrder);
    BruteForce(graphCond, graphCondUndir, tsOrder, special);
}

void BruteForce(const TGraph& graph, const TGraph& graphUndir, const vector<size_t>& tsOrder) {
    (void)graphUndir;
    const size_t totalVertex = static_cast<size_t>(tsOrder.size());
    vector<size_t> pos;
    pos.reserve(totalVertex);
    for (size_t sampleSize = 1; sampleSize <= totalVertex; sampleSize++) {
        //cerr << "BF Sample size: " << sampleSize << endl;
        pos.resize(sampleSize);
        for (size_t i = 0; i < pos.size(); ++i) {
            pos[i] = i;
        }
        while (true) {
            set<size_t> bh = GetBlackHole(graph, tsOrder, pos);
            FoundBlackHole(bh.size());
            // cout << "BH count " << bhCount << endl;
            // Print("BH", bh);
            if (!BruteNext(pos, totalVertex)) {
                break;
            }
        }
    }
}

void BruteForce(const TGraph& graph, const TGraph& graphUndir, const vector<size_t>& tsOrder, const vector<char>& special) {
    size_t filtered = 0;
    size_t skipped = 0;
    const size_t totalVertex = static_cast<size_t>(tsOrder.size());
    const size_t maxSampleSize = totalVertex;
    vector<size_t> pos;
    pos.reserve(totalVertex);
    for (size_t sampleSize = 1; sampleSize <= min(maxSampleSize, totalVertex); sampleSize++) {
        size_t sampleSizeCount = 0;
        // cerr << "BF Sample size: " << sampleSize << endl;
        pos.resize(sampleSize);
        for (size_t i = 0; i < pos.size(); ++i) {
            pos[i] = i;
        }
        bool stop = false;
        while (!stop) {
            set<size_t> bh = GetBlackHole(graph, tsOrder, pos);
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
                for (size_t i = 1; i < pos.size(); i++) {
                    size_t p = pos[i];
                    size_t v = tsOrder[p];
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

        //cout << "bh_count " << bhCount << endl;
        cout << "bh_filtered " << filtered << endl;
        cout << "bh_skipped " << skipped << endl;
    }
}

bool BruteNext(vector<size_t>& pos, size_t vertexCount) {
    size_t sampleSize = pos.size();
    bool res = true;
    for (size_t i = sampleSize; i > 0; --i) {
        size_t& cur = pos[i - 1];
        size_t distFromEnd = sampleSize - i;
        if (cur == vertexCount - 1 - distFromEnd) {
            if (i == 1) {
                res = false;
            }
            continue;
        } else {
            cur++;
            for (size_t j = i; j < sampleSize + 1; j++) {
                pos[j] = pos[j - 1] + 1;
            }
            break;
        }
    }
    return res;
}

void FoundBlackHole(size_t size) {
    (void)size;
    static size_t bhCount = 0;
    bhCount++;
    cout << "bh_count" << bhCount << endl;
}

set<size_t> GetBlackHole(const TGraph& graph,
                         const vector<size_t>& tsOrder, const vector<size_t>& pos) {
    set<size_t> blackhole;
    TUsed used;
    used.assign(graph.size(), 0);
    for (size_t p : pos) {
        size_t v = tsOrder[p];
        TClosure closure = GetClosure(graph, v, used);
        for (size_t x : closure) {
            blackhole.insert(x);
        }
    }
    for (size_t p : pos) {
        size_t x = tsOrder[p];
        if (used[x] >= 2) {
            return {}; // means this blackhole has already been deteected;
        }
    }
    return blackhole;
}

void Print(const string& tag, const set<size_t>& v) {
    cout << tag << ": ";
    for (size_t x : v) {
        cout << x << " ";
    }
    cout << endl;
}

void Print(const string& tag, const vector<size_t>& v) {
    cout << tag << ": ";
    for (size_t x : v) {
        cout << x << " ";
    }
    cout << endl;
}

void Print(const string& tag, const vector<vector<size_t>>& vv) {
    cout << tag << " BEGIN" << endl;
    for (size_t i = 0; i < vv.size(); i++) {
        const auto& v = vv[i];
        cout << i << ": ";
        for (size_t x : v) {
            cout << x << " ";
        }
        cout << endl;
    }
    cout << tag + " END" << endl;
}
