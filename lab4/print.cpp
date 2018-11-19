#include <iostream>
#include <fstream>


using namespace std;

int main(int argc, char * argv[]){

	if(argc < 2){
		cout << "use ./print A.dat" << endl;
		return 0;
	}
	const string inputFileName(argv[1]);

	std::ifstream ifs(inputFileName, std::ios::binary);

	int rows, cols;
	ifs.read(reinterpret_cast<char *>( &rows ), sizeof(int)/sizeof(char));
	ifs.read(reinterpret_cast<char *>( &cols ), sizeof(int)/sizeof(char));
	int mSize = cols * rows;
	int bufSize = mSize * sizeof(double) / sizeof(char);
	double * buf = new double[mSize];

	char * charBuf = reinterpret_cast<char *>(buf);
	ifs.read(charBuf, bufSize);
	
	cout << rows << " " << cols << endl;
	for(int i =0 ;i < rows; i++){
		for(int j = 0; j < cols; j++){
			cout << buf[i * cols + j] << " ";
		}
		cout << endl;
	}
	return 0;
}