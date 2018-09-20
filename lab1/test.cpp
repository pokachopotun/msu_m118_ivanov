#include <matrix.h>

#include <iostream>
#include <fstream>


using namespace std;

int main(int argc, char * argv[]){

	if(argc < 2){
		cout << "use ./print A.dat" << endl;
		return 0;
	}
	string inputFileName(argv[1]);
	Matrix m(inputFileName);
	cout << m;
	return 0;
}