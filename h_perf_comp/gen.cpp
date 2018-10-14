#include <matrix.h>
#include <iostream>
#include <fstream>
#include <random>
#include <string>

#include <cstdlib>

using namespace std;

int main(int argc, char* argv[]){
	if(argc < 4){
		cout << "use ./gen N M K outfile.dat" << endl;
		return 0;
	}
	const size_t N = atoi(argv[1]);
	const size_t M = atoi(argv[2]);
	const size_t K = atoi(argv[3]);
	// const string outputFilename(argv[4]);

	Matrix m;
	m.Generate(N, M, K);
	m.PrintToTxt();
	cout << "_row: "<< endl;
	m.PrintBuffers();
	return 0;
}
