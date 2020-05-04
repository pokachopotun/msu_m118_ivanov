#include <algorithm.h>
#include <common.h>
#include <blackholes.h>
#include <graph.h>
#include <io.h>

#include <iostream>
#include <string>
#include <chrono>

using namespace std;

int main(int argc, char * argv[]) {
    if (argc < 4){
        cout << "use ./topsort_cli input.txt binaryInput divideAndConquer ompThreads" << endl;
        return 0;
    }

    const string inputFilePath(argv[1]);
    const int binaryInput(atoi(argv[2]));
    const int divideAndConquer(atoi(argv[3]));
    const size_t ompThreads(atoi(argv[4]));

    using namespace std::chrono;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    {
        TGraph graph;
        if (binaryInput) {
            graph = ReadGraphFromBinaryFile(inputFilePath);
        } else {
            graph = ReadGraphFromTextFile(inputFilePath);
        }
        TGraph graphRev = GetReverseGraph(graph);

        size_t compCnt = 0;
        TComponents comp2Vertex = GetStrongConnectivityComponents(graph, graphRev, compCnt);
        TGraph graphCond = BuildCondensedGraph(compCnt, comp2Vertex, graph);

        if (divideAndConquer) {
            TGraph graphCondUndir = GetUndirGraph(graphCond);
            vector<TGraph> weaks = GetWeakComponents(graphCond, graphCondUndir);
            for (const TGraph& weak : weaks) {
                TopsortSolver(weak, ompThreads);
            }
        } else {
            TopsortSolver(graphCond, ompThreads);
        }
    }
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
    cout << "It took me: "<< time_span.count() <<  " seconds" << endl;
    return 0;
}
