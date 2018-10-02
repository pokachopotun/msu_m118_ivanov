#include <matrix.h>
#include <iostream>
#include <fstream>
#include <random>
#include <string>

#include <cstdlib>

using namespace std;

int main(int argc, char* argv[]){
	if(argc < 6){
		cout << "use ./fromtxt f N M outfile.dat txt.in " << endl;
		return 0;
	}
	const char type = argv[1][0];
	const size_t N = atoi(argv[2]);
	const size_t M = atoi(argv[3]);
	const string txtFileName(argv[4]);
	const string outputFilename(argv[5]);

	Matrix m(N, M, type, nullptr);
	m.ReadInTxt(txtFileName);
	m.WriteToFile(outputFilename);
	return 0;
}