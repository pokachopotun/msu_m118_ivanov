#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;


int main(int argc, char* argv[]){
	if(argc < 4){
		cout << "use ./mul A.in B.in C.out" << endl;
		return 0;
	}

	const double eps = 1e-6;

	int rows[2], cols[2];
	int mSize[2];
	double * buf[2];
	for(int i = 0; i < 2; i++){
		std::ifstream ifs(argv[i + 1], std::ios::binary);
		ifs.read(reinterpret_cast<char *>( &rows[i] ), sizeof(int)/sizeof(char));
		ifs.read(reinterpret_cast<char *>( &cols[i] ), sizeof(int)/sizeof(char));
	
		mSize[i] = cols[i] * rows[i];
		int bufSize = mSize[i] * sizeof(double) / sizeof(char);
		buf[i] = new double[mSize[i]];

		char * charBuf = reinterpret_cast<char *>(buf[i]);
		ifs.read(charBuf, bufSize);
	}

	double * ansBuf = new double[rows[0] * cols[1]];
	for(int i = 0; i < rows[0] * cols[1]; i++){
		ansBuf[i] = 0.0;
	}

	for(int i=0 ;i < rows[0]; i++){
		for(int j =0; j < cols[1]; j++){
			for(int k = 0; k < rows[1]; k++){
				ansBuf[i * cols[1] + j] += buf[0][i * cols[0] + k] * buf[1][k * cols[1] + j];
			}
		}
	}
	{
		std::ofstream ofs(argv[3], std::ios::binary);
		ofs.write(reinterpret_cast<char *>( &rows[0] ), sizeof(int)/sizeof(char));
		ofs.write(reinterpret_cast<char *>( &cols[1] ), sizeof(int)/sizeof(char));
	
		int mSize = cols[1] * rows[0];
		int bufSize = mSize * sizeof(double) / sizeof(char);

		char * charBuf = reinterpret_cast<char *>(ansBuf);
		ofs.write(charBuf, bufSize);
	}
	return 0;
}
