#include <iostream>
#include <chrono>

#include <algorithm.h>
#include <blackholes.h>
#include <common.h>
#include <graph.h>
#include <io.h>

using namespace std;

bool printDebugInfo = false;

int main(int argc, char** argv) {
    if (argc < 4) {
        cout << "use ./condensation_cli inputFileName binaryInput printDebugInfo" << endl;
        return 0;
    }
    const string inputFileName(argv[1]);
    const int binaryInput(atoi(argv[2]));
    if (argc > 3) {
        printDebugInfo = atoi(argv[3]);
    }

    using namespace std::chrono;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    {
        TGraph graph;
        if (binaryInput) {
            graph = ReadGraphFromBinaryFile(inputFileName);
        } else {
            graph = ReadGraphFromTextFile(inputFileName);
        }

        TGraph graphRev = GetReverseGraph(graph);
        size_t compCnt = 0;
        TComponents comp2Vertex = GetStrongConnectivityComponents(graph, graphRev, compCnt);
        TGraph graphCond = BuildCondensedGraph(compCnt, comp2Vertex, graph);

        auto&& countEdges = [](const TGraph& graph) -> size_t {
            size_t edges = 0;
            for (const vector<size_t>& e : graph) {
                edges += e.size();
            }
            return edges;
        };

        cout << "orig " << graph.size() << " " << countEdges(graph) << endl;
        cout << "cond " << graphCond.size() << " " << countEdges(graphCond) << endl;
    }
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
    cout << "It took me: " << time_span.count() << " seconds" << endl;
    return 0;
}
