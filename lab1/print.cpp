#include <iostream>
#include <fstream>
#include <string>

using namespace std;

int main(int argc, char * argv[]){
	if(argc < 2){
		cout << "use ./print A.dat" << endl;
		return 0;
	}

	string inputFileName(argv[1]);
	ifstream ifs(argv[1], ios::binary);
	
	char type;
	size_t N;
	size_t M;

	ifs.read(reinterpret_cast<char *>( &type ), sizeof(char));
	ifs.read(reinterpret_cast<char *>( &N ), sizeof(size_t)/sizeof(char));
	ifs.read(reinterpret_cast<char *>( &M ), sizeof(size_t)/sizeof(char));

	const size_t matrixSize = N * M;
	size_t bufferSize = 0;
	if( type == 'f' ){
		bufferSize = matrixSize * sizeof(float)/sizeof(char);
	}else{
		bufferSize = matrixSize * sizeof(double)/sizeof(char);
	}

	char * charBuf = new char[bufferSize];
	float * floatBuf = reinterpret_cast<float *> (charBuf);
	double * doubleBuf = reinterpret_cast<double *> (charBuf);
		
	ifs.read(charBuf, bufferSize);

	if(type == 'f'){
		for(size_t i = 0 ; i < N; i++){
			for(size_t j = 0; j < M; j++){
				size_t pos = i * M + j;
				cout << floatBuf[pos] << " ";
			}	
			cout << endl;
		}
	}else{
		for(size_t i=0 ; i < N; i++){
			for(size_t j = 0; j < M; j++){
				size_t pos = i * M + j;
				cout << doubleBuf[pos] << " ";
			}	
			cout << endl;
		}
	}

	delete[] charBuf;
	return 0;
}
