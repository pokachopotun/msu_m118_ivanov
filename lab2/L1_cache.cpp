#include <papi.h>
#include <iostream>
#include <memory>
#include <testsized.h>
#include <chrono>

using namespace std;
typedef std::chrono::high_resolution_clock Clock;

const int numEvents = 1;

int main(int argc, char * argv[]) {
	
	if(argc < 4){
		cout << "use ./timeit mSize bSize ijk" << endl;
		return 0;
	}

	const int mSize = atoi(argv[1]);
	const int bSize = atoi(argv[2]);
	const string mode(argv[3]);
	int Events[numEvents] = {PAPI_L1_DCM};
	long long values[numEvents];
	unique_ptr< TestSized > a(new TestSized(mSize));
	if(PAPI_start_counters(Events, numEvents)!=PAPI_OK){
		return 1;
	}
	a->Multiply(bSize,mode);	
	if(PAPI_stop_counters(values, numEvents) != PAPI_OK){
		return 2;
	}
	printf("Size %d Block %d Mode %s L1_DCM %lld \n", mSize, bSize, mode.c_str(), values[0]);
	return 0;
}
