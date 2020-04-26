#include <algorithm.h>
#include <common.h>
#include <blackholes.h>
#include <graph.h>
#include <io.h>

#include <iostream>
#include <string>
//#include <vector>
//#include <set>
//#include <algorithm>
//#include <cmath>
//#include <iomanip>
#include <chrono>
//#include <fstream>
//#include <stack>
//#include <queue>
//#include <map>

using namespace std;

int main(int argc, char * argv[]) {
    if (argc < 3){
        cout << "use ./main input.txt optimize";
    }

    const string inputFilePath(argv[1]);
    const int optimize(atoi(argv[2]));

    // #ifdef _DEBUG
        //ofstream fout(outputFilePath, std::ifstream::out);
    // #endif

    using namespace std::chrono;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    {
        TGraph graph = ReadGraphFromBinaryFile(inputFilePath);
        TGraph graphRev = GetReverseGraph(graph);

        TComponents comp2Vertex = GetStrongConnectivityComponents(graph, graphRev);
        TSccMembers compMembers = GetSccMembers(comp2Vertex);
        // Print("components", compMembers);
        //PrintStat("SCC component sizes", compMembers);
        size_t compCnt = GetCompCnt(comp2Vertex);
        TGraph graphCond = BuildCondensedGraph(compCnt, comp2Vertex, graph);
        if (optimize) {
            TGraph graphCondUndir = GetUndirGraph(graphCond);
            vector<TGraph> weaks = GetWeakComponents(graphCond, graphCondUndir);
            //PrintStat("Weakly connected component sizes", weaks);
            //if (false)
            for (const TGraph& weak : weaks) {
                Solve(weak);
            }
        } else {
            TGraph graphCondUndir = GetUndirGraph(graphCond);
            vector<size_t> tsOrder = TopSort(graphCond);
  //          Print("topsort", tsOrder);
            BruteForce(graphCond, graphCondUndir, tsOrder);
        }
/*
        int total = 0;
        fout << inputFilePath << endl;
        fout << "Nodes " << n << "Edges " << m << endl;
        fout << "BlackHoles cnt: " << endl;
        for( int i = 0 ; i < blackHoles.size(); i++ ){
                if( blackHoles[i].size() > 0 ){
                        total += blackHoles[i].size();
                        fout <<" Size: " << i << " cnt "<< blackHoles[i].size() << endl;
                }
        }
        fout << "Total: " << total << endl;
*/
    }
    high_resolution_clock::time_point t2 = high_resolution_clock::now();

    duration<double> time_span = duration_cast<duration<double>>( t2 - t1 );
    cout << "It took me: "<< time_span.count() <<  " seconds" << endl;
    return 0;
}
