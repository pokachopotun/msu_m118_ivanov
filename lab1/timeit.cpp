#include <matrix.h>

#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

int main(int argc, char * argv[]) {
	
	if(argc < 6){
		cout << "use ./timeit f 3 3 ijk stats.out" << endl;
		return 0;
	}
	cout << "timeit.cpp" << endl;
	
	char type(atoi(argv[1]));
	size_t rows(atoi(argv[2]));
	size_t cols(atoi(argv[3]));
	string mode(argv[4]);
	string outputFileName(argv[5]);

	Matrix a(rows, cols, type);
	Matrix b(rows, cols, type);
	std::cout << "Init done" << std::endl;
	vector<double> times;
	for(int i=0 ; i < 10; i++) {
		a.MakeRandom();
		b.MakeRandom();
		Matrix c = Multiply( a, b, mode );
		times.push_back(MyTimer::GetSeconds());
		cout << "Multiply done in " << times.back() << " seconds" << endl;
		c.clean();
	}
	//c.WriteToFile(outputFileName);
	
	a.clean();
	b.clean();
	// c.clean();
	return 0;
}
