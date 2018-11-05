#include <mpi.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <string>

using namespace std;

int main(int argc, char* argv[]){
	if(argc < 4){
		cout << "Use main.cpp first last output.txt" << endl;
		return 0;
	}
	MPI_Init(&argc, &argv);

	double localTime = MPI_Wtime();

	int first = atoi(argv[1]);
	int last = atoi(argv[2]);
	const string outputFilename(argv[3]);
	int rank, size;
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	int total = last + 1 - first;
	int part = total / size + int( total % size != 0 );
	first += rank * part;
	last = min( last + 1, first + part);
	
	// cerr << rank << " " << first << " " << last << endl;
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
	int * numbers = (int * ) calloc(part, sizeof(int));
	int lastNum = 0;
	for(int i = 0; i < primeLocal.size(); i++){
		if(primeLocal[i]){
			numbers[lastNum++] = i + first;
		}
	}

	localTime = MPI_Wtime() - localTime;

	if(rank > 0){
		MPI_Send(&localTime, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		MPI_Send(numbers, lastNum, MPI_INTEGER, 0, 0, MPI_COMM_WORLD);
	}
	if(rank == 0){
		ofstream fout(outputFilename.c_str());
		for(int i = 0; i < lastNum; i++){
			fout << numbers[i] << " ";
		}
		double sumTime = localTime;
		for(int sender = 1; sender < size; sender++){
			MPI_Status status;
			int recvCnt;
			double recvTime;
			MPI_Recv(&recvTime, 1, MPI_DOUBLE, sender, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(numbers, part, MPI_INTEGER, sender, 0, MPI_COMM_WORLD, &status);
			MPI_Get_count(&status, MPI_INTEGER, &recvCnt);
			localTime = max(localTime, recvTime);
			sumTime += recvTime;
			lastNum += recvCnt;
			for( int i =0; i < recvCnt; i++){
				fout << numbers[i] << " ";
			}
		}
		fout << endl;
		printf("Nums %d maxTime %f sumTime %f \n", lastNum, localTime, sumTime );
	}
	MPI_Finalize();
	return 0;
}
