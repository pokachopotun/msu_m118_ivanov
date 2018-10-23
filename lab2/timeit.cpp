#include <papi.h>
#include <iostream>
#include <memory>
#include <testsized.h>
#include <chrono>

using namespace std;
typedef std::chrono::high_resolution_clock Clock;

const int numEvents = 1;
const int reps = 1;


void timeit(int mSize, int bSize, const string& mode){
	double time = 0;
	for(int i=0 ;i < reps; i++){
		unique_ptr< TestSized > a( new TestSized(mSize) );
		auto t1 = Clock::now();
		a->Multiply(bSize, mode);
		auto t2 = Clock::now();
		double s = std::chrono::duration_cast<chrono::milliseconds>(t2 - t1).count();
		time += s;
	}
	// printf("Size %d Block %d Mode %s Time %6.1f ms\n", mSize, bSize, mode.c_str(), time);
	printf("Size %d Block %d Mode %s Time %6.1f ms\n", mSize, bSize, mode.c_str(), time);
}

int main(int argc, char * argv[]) {
	
	if(argc < 4){
		cout << "use ./timeit mSize bSize ijk" << endl;
		return 0;
	}

	const int mSize = atoi(argv[1]);
	const int bSize = atoi(argv[2]);
	const string mode(argv[3]);

	timeit(mSize, bSize, mode);

	return 0;
}
