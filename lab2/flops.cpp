#include <papi.h>
#include <iostream>
#include <memory>
#include <testsized.h>
#include <chrono>

using namespace std;
typedef std::chrono::high_resolution_clock Clock;

const int reps = 1;

int main(int argc, char * argv[]) {
	
	if(argc < 4){
		cout << "use ./timeit mSize bSize ijk" << endl;
		return 0;
	}

	const int mSize = atoi(argv[1]);
	const int bSize = atoi(argv[2]);
	const string mode(argv[3]);
	float rtime, ptime;
	long long flpops;
	float mflops;
	unique_ptr< TestSized > a(new TestSized(mSize));
	PAPI_flops(&rtime, &ptime, &flpops, &mflops);
	a->Multiply(bSize,mode);	
	PAPI_flops(&rtime, &ptime, &flpops, &mflops);
	printf("Size %d Block %d Mode %s rtime %f ptime %f flpops %lld mflops %f \n", mSize, bSize, mode.c_str(), rtime, ptime, flpops, mflops);
	return 0;
}
