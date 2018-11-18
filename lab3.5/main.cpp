#include <pthread.h>

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <string>
#include <chrono>

using namespace std;

int gFirst, gLast, nThreads;

struct threadParams_t{
	int rank;
	double time;
	vector<int>* numbers;
};

void* threadFunc(void * args){

	auto time = high_resolution_clock::now()

	threadParams_t* threadParams = (threadParams_t * )args;

	int first = gFirst;
	int last = gLast;

	int rank = threadParams->rank;
	int size = nThreads;
	
	int total = last + 1 - first;
	int part = total / size + int( total % size != 0 );
	first += rank * part;
	last = min( last + 1, first + part);

	if(first > gLast){
		time = high_resolution_clock::now() - time;
	}else{
	
		cerr << rank << " " << first << " " << last << endl;
		int n = sqrt(last);
		vector<bool> prime(n + 1, true);
		prime[0] = prime[1] = false;
		for(int i = 2; i <= n; ++i){
			if(prime[i]){
				if(static_cast<long long>(i) * i <= n){
					for(int j = i * i; j <= n; j += i){
						prime[j] = false;
					}
				}
			}
		}
		vector<bool> primeLocal(last - first, true);
		if(first == 1){
			primeLocal[0] = false;
		}
		for(int i = 1; i <= n; i++){
			if(!prime[i]){
				continue;
			}
			int s = ( first / i ) * i;
			if( s < first ){
				s += i;
			}
			for(int j = s; j < last; j += i ){
				primeLocal[j - first] = false;
			}
			if(i == s){
				primeLocal[i - first] = true;
			}
		}
		for(int i = 0; i < primeLocal.size(); i++){
			if(primeLocal[i]){
				threadParams->numbers->push_back( i + first );
			}
		}
		time = high_resolution_clock::now() - time;
	}
	threadParams->time = static_cast<double>( time );
	return 0;
}

int main(int argc, char* argv[]){
	if(argc < 4){
		int n = sqrt(1);
		cout << 
		cout << "Use ./main nThreads first last output.txt" << endl;
		return 0;
	}

	nThreads = atoi(argv[1]);
	gFirst = atoi(argv[2]);
	gLast = atoi(argv[3]);
	const string outputFilename(argv[4]);

	pthread_t* threadIds = new pthread_t[nThreads];
	threadParams_t* threadParams = new threadParams_t[nThreads];
	for(int i = 0; i < nThreads; i++){
		threadParams[i].rank = i;
		threadParams[i].numbers = new vector<int>();
	}

	int* cnt = (int*)calloc(nThreads, sizeof(int));
	int gCnt = 0;

	for(int i =0 ; i < nThreads; i++){
		pthread_create(&threadIds[i], NULL, threadFunc, &threadParams[i]);
	}

	for (int i = 0; i < nThreads; i++){
       pthread_join(threadIds[i], NULL);
	}
	int totalNums = 0;
	ofstream fout(outputFilename.c_str());
	for(int i =0; i < nThreads; i++){
		totalNums += threadParams[i].numbers->size();
		for(int j = 0; j < threadParams[i].numbers->size(); j++){
			fout << threadParams[i].numbers->operator[](j) << " ";
		}
		
	}
	fout << endl;
	printf("Nodes %d Range %d %d Nums %d maxTime %f sumTime %f \n", size, gFirst, gLast, totalNums, localTime, sumTime );
	return 0;

//----------------------------

	// ofstream fout(outputFilename.c_str());
	// for(int i = 0; i < lastNum; i++){
	// 	fout << numbers[i] << " ";
	// }
	// double sumTime = localTime;
	// for(int sender = 1; sender < size; sender++){
	// 	MPI_Status status;
	// 	int recvCnt;
	// 	double recvTime;
	// 	MPI_Recv(&recvTime, 1, MPI_DOUBLE, sender, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	// 	MPI_Recv(numbers, part, MPI_INTEGER, sender, 0, MPI_COMM_WORLD, &status);
	// 	MPI_Get_count(&status, MPI_INTEGER, &recvCnt);
	// 	localTime = max(localTime, recvTime);
	// 	sumTime += recvTime;
	// 	lastNum += recvCnt;
	// 	for( int i =0; i < recvCnt; i++){
	// 		fout << numbers[i] << " ";
	// 	}
	// }
	// fout << endl;
	

	return 0;
}
