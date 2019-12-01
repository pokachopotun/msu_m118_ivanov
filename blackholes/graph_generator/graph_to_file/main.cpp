#include <iostream>
#include <fstream>
#include <string>

using namespace std;

int main(int argc, char** argv) {

    if (argc < 3) {
        cout << "Use ./graph_to_file file.bin file.graph" << endl;
        return 0;
    }

    const string file_name(argv[1]);
    const string output_file_name(argv[2]);

    // open file
    fstream file(file_name, ios::in | ios::binary);
    ofstream out(output_file_name, ios::out);

    int vertices_count = 0;
    long long edges_count = 0;

    // read header
    file.read((char*)(&vertices_count), sizeof(int));
    file.read((char*)(&edges_count), sizeof(long long));

    // print graph
//    cout << "Graph has " << vertices_count << " vertices" << endl;
 //   cout << "Graph has " << edges_count << " edges" << endl;

    out << vertices_count << " " << edges_count << endl;

    for(long long i = 0; i < edges_count; i++){
            int src_id = 0, dst_id = 0;
            float weight = 0;

            file.read((char*)(&src_id), sizeof(int));
            file.read((char*)(&dst_id), sizeof(int));
//            file.read((char*)(&weight), sizeof(float)); // remove it for unweighed graph

  //          cout << src_id << " " << dst_id << " | " << weight << endl;
            out << src_id << " " << dst_id << endl;
    }

    file.close();
    out.close();
}
