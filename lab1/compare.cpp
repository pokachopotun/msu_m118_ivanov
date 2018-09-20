#include <iostream>
#include <fstream>
#include <matrix.h>

using namespace std;


int main(int argc, char* argv[]){
	if(argc < 3){
		cout << "use ./compare A.in B.in" << endl;
		return 0;
	}

	cout << "compare.cpp" << endl;

	const string inFileNameA(argv[1]);
	const string inFileNameB(argv[2]);

	Matrix a(inFileNameA);
	Matrix b(inFileNameB);

	cout << boolalpha << Compare(a, b) << endl;

	a.clean();
	b.clean();
	return 0;
}
