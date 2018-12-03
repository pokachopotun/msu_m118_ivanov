#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;


int main(int argc, char* argv[]){
	if(argc < 3){
		cout << "use ./compare A.in B.in" << endl;
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

	bool res = true;
	res = res && rows[0] == rows[1];
	res = res && cols[0] == cols[1];
	for(int i = 0; i < mSize[0]; i++){
		res = res && fabs( buf[0][i] - buf[1][i] ) < eps;
	}
	
	cout << boolalpha << res << endl;
	return 0;
}
