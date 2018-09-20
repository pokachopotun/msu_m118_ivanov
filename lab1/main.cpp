#include <matrix.h>

#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char * argv[]) {
	
	if(argc < 4){
		cout << "use ./main A.in B.in C.out" << endl;
		return 0;
	}
	cout << "main.cpp" << endl;
	string inputFileA(argv[1]);
	string inputFileB(argv[2]);
	string outputFileName(argv[3]);

	Matrix a(inputFileA);
	Matrix b(inputFileB);
	std::cout << "Init done" << std::endl;
	Matrix c = Multiply(a, b, "ijk" );
	cout << "Multiply done " << endl;
	c.WriteToFile(outputFileName);
	a.clean();
	b.clean();
	c.clean();
	return 0;
}
