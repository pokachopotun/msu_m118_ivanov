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
    static TGraph ReadGraphFromFile(const string& filename) {
        fstream file(filename, ios::in | ios::binary);
        int n;
        long long m;
        // read header
        file.read((char*)(&n), sizeof(int));
        file.read((char*)(&m), sizeof(long long));

        TGraph graph(n);

        vector<set<int>> tmp(n);
        int e = 0;
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
            e += t.size();
            for (int x : t) {
                graph[i].push_back(x);
            }
        }
        cout << filename << " " << n << " " << e << endl;
        return graph;
    }
};

int main(int argc, char * argv[]) {
    if (argc < 2){
        cout << "use ./main input.txt";
    }

    const string inputFilePath(argv[1]);

    TGraph graph = Solution::ReadGraphFromFile(inputFilePath);

    return 0;
}
