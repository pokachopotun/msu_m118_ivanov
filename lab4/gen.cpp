#include <iostream>
#include <fstream>
#include <random>
#include <string>
#include <random>

#include <cstdlib>

using namespace std;

int main(int argc, char* argv[]){
	if(argc < 4){
		cout << "use ./gen N M outfile.dat" << endl;
		return 0;
	}

	const int rows = atoi(argv[1]);
	const int cols = atoi(argv[2]);
	const string outputFilename(argv[3]);

	double * buf = new double[rows * cols];
	char * charBuf = reinterpret_cast<char*>(buf);

	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<float> dis(0, 1);
	static_cast<double>( dis(gen) );
	for(int i = 0; i < rows * cols; i++){
		buf[i] = static_cast<double>( dis(gen) );
	}

	std::ofstream ofs(outputFilename, std::ios::binary);
	ofs.write(reinterpret_cast< const char* >( &rows ), sizeof(int)/sizeof(char));
	ofs.write(reinterpret_cast< const char* >( &cols ), sizeof(int)/sizeof(char));
	ofs.write(charBuf, rows * cols * sizeof(double)/sizeof(char));

	cout << "Generation done" << endl;
	return 0;
}
