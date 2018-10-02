#include <matrix.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>

using namespace std;

int main(int argc, char * argv[]) {
	
	if(argc < 5){
		cout << "use ./timeit f 3 3 ijk" << endl;
		return 0;
	}
	map<string, int> mp;
	mp["ijk"] = 1;
	mp["ikj"] = 2;
	mp["jik"] = 3;
	mp["jki"] = 4;
	mp["kij"] = 5;
	mp["kji"] = 6;


	const int reps = 10;

	
	char type(atoi(argv[1]));
	size_t rows(atoi(argv[2]));
	size_t cols(atoi(argv[3]));
	string mode(argv[4]);


	Matrix a(rows, cols, type);
	Matrix b(rows, cols, type);
	// std::cout << "Init done" << std:

	vector<double> times;
	for(int i=0 ; i < reps; i++) {
		a.MakeRandom();
		b.MakeRandom();
		Matrix c = Multiply( a, b, mode );
		times.push_back(MyTimer::GetSeconds());
		
		c.clean();
	}
	cout << mp[mode] << " " << accumulate( times.begin(), times.end(), 0.0)/times.size() << endl;
	//c.WriteToFile(outputFileName);
	a.clean();
	b.clean();
	// c.clean();
	return 0;
}
