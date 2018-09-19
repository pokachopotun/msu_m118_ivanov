#include <iostream>
#include <fstream>
#include <random>
#include <string>

#include <cstdlib>

using namespace std;

int main(int argc, char* argv[]){
	if(argc < 5){
		cout << "use ./gen f N M outfile.dat" << endl;
		return 0;
	}
	const char type = argv[1][0];
	const size_t N = atoi(argv[2]);
	const size_t M = atoi(argv[3]);
	const string outputFilename(argv[4]);
	
	const size_t matrixSize = N * M;
	size_t bufferSize = 0;
	if( type == 'f' ){
		bufferSize = matrixSize * sizeof(float)/sizeof(char);
	}else{
		bufferSize = matrixSize * sizeof(double)/sizeof(char);
	}
	// cout << bufferSize << endl;
	ofstream ofs(outputFilename, ios::binary);
	
	char * charbuf = new char[bufferSize];

	double * doubleBuf = reinterpret_cast< double * >(charbuf);
	float * floatBuf = reinterpret_cast< float * >(charbuf);

	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<float> dis(0, 1);
	if(type == 'f'){
		for(size_t i = 0 ; i < N; i++){
			for(size_t j = 0; j < M; j++){
				size_t pos = i * M + j;
				// cout << i << " " << j << " " << pos << endl;
				floatBuf[pos] = static_cast<float>( dis(gen) );
			}	
		}
	}else{
		for(size_t i=0 ; i < N; i++){
			for(size_t j = 0; j < M; j++){
				size_t pos = i * M + j;
				doubleBuf[pos] = static_cast<double>( dis(gen) );
			}	
		}
	}
	// cout << "dumping" << endl;
	ofs.write(&type, sizeof(char));
	ofs.write(reinterpret_cast< const char* >( &N ), sizeof(size_t)/sizeof(char));
	ofs.write(reinterpret_cast< const char* >( &M ), sizeof(size_t)/sizeof(char));
	ofs.write(charbuf, bufferSize);
	delete[] charbuf;
	return 0;
}
