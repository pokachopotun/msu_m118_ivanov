#include "io.h"

using namespace std;

TGraph ReadGraphFromBinaryStream(std::istream& file) {
    int n;
    long long m;
    file.read((char*)(&n), sizeof(int));
    file.read((char*)(&m), sizeof(long long));
    TGraph graph(n);
    vector<set<int>> tmp(n);
    for(long long i = 0; i < m; i++) {
        int x, y;
        file.read((char*)(&x), sizeof(int));
        file.read((char*)(&y), sizeof(int));
        if (x == y) { // omit self loops
            continue;
        }
        tmp[x].insert(y);
    }
    for(size_t i = 0; i < tmp.size(); i++) {
        const auto& t = tmp[i];
        for (int x : t) {
            graph[i].push_back(static_cast<size_t>(x));
        }
    }
    return graph;
}

TGraph ReadGraphFromBinaryFile(const std::string& filename) {
    ifstream file(filename, ios::in | ios::binary);
    return ReadGraphFromBinaryStream(file);
}

TGraph ReadGraphFromTextStream(std::istream& stream) {
    int n = 0;
    long long m = 0;
    stream >> n >> m;
    TGraph graph(n);

    vector<set<int>> tmp(n);
    for(int i = 0; i < m; i++) {
        int x, y;
        stream >> x >> y;
        tmp[x].insert(y);
    }
    for(size_t i = 0; i < tmp.size(); i++) {
        const auto& t = tmp[i];
        for (int x : t) {
            graph[i].push_back(static_cast<size_t>(x));
        }
    }
    return graph;
}

TGraph ReadGraphFromTextFile(const string& filename) {
    ifstream stream(filename);
    return ReadGraphFromTextStream(stream);
}
