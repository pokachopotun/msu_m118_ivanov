#include <cstdlib>
#include <mpi.h>
#include <iostream>
#include <string>
#include <cmath>
using namespace std;

double frand() {
	return static_cast<double>(rand()) / static_cast <double> (RAND_MAX);
}

int main( int argc, char* argv[] ) {
	
	MPI_Init(&argc, &argv);	
	if( argc < 4 ) {
		cout << "Use ./main numIter exchangeStep systemId" << endl;
		MPI_Finalize();
		return 0;
	}
	const int numIter = atoi( argv[1] );
	const int exStep = atoi( argv[2] );
	const int systemId = atoi( argv[3] );

	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	string data;
	pair<double, string> L1[2];
	
	switch( systemId ) {
		case 1: 
			data = "a";
			// 1.0 1.0
			L1[0] = make_pair<double, string>( 1.0 , "ab" );
			L1[1] = make_pair<double, string>( 1.0, "bc" );
			break;
		case 2:
			data = "a";
			L1[0] = make_pair<double, string>( 1.0 , "aa" ); // 0.001;
			break;
		case 3:	
			data = "a";
				// 0.01 0.01
			L1[0] = make_pair<double, string>( 1.0, "bc");
			L1[1] = make_pair<double, string>( 1.0 , "a");
			break;
		default:
			cout << " systemId should be in range 1:3" << endl;
			MPI_Finalize();
			return 0;
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	double time = MPI_Wtime();

	if( rank != size / 2 ){
		data = "";
	}

	for( int iter = 0; iter < numIter; iter++ ) 
	{
		if( iter % exStep == 0 ){

			int nextRank = (rank + 1 + size ) % size;
			int prevRank = (rank - 1 + size ) % size;
			MPI_Request req[20];
			{
				int my = data.size();
				int left = 0;
				int right = 0;
				MPI_Isend(&my, 1, MPI_INTEGER, nextRank, 0, MPI_COMM_WORLD,&req[0]);
				MPI_Irecv(&left, 1, MPI_INTEGER, prevRank, 0, MPI_COMM_WORLD,&req[1]);
				MPI_Isend(&my, 1, MPI_INTEGER, prevRank, 0, MPI_COMM_WORLD,&req[2]);
				MPI_Irecv(&right, 1, MPI_INTEGER, nextRank, 0, MPI_COMM_WORLD,&req[3]);

				MPI_Waitall(4, req, MPI_STATUSES_IGNORE);
				
				// left right
				int difLeft = ( my - left ) / 2;
				int difRight = ( right - my ) / 2;
				string sendBuf;
				string recvBuf;
				int cnt = 0;
				if( difLeft > 0 && prevRank < rank ){
					sendBuf = data.substr(0, difLeft);
					data = data.substr( difLeft, int( data.size() ) - difLeft );
					MPI_Isend( &sendBuf[0], difLeft, MPI_CHAR, prevRank, 0, MPI_COMM_WORLD, &req[4 + cnt]);
					cnt++;
				}
				if( difRight > 0 && nextRank > rank ){
					data.resize( my + difRight );
					MPI_Irecv( &data[my], difRight, MPI_CHAR, nextRank, 0, MPI_COMM_WORLD, &req[4 + cnt]); cnt++;
				}
				MPI_Waitall(cnt, &req[4], MPI_STATUSES_IGNORE);
			}
			
			{
				int my = data.size();
				int left = 0;
				int right = 0;
				MPI_Isend(&my, 1, MPI_INTEGER, nextRank, 0, MPI_COMM_WORLD,&req[6]);
				MPI_Irecv(&left, 1, MPI_INTEGER, prevRank, 0, MPI_COMM_WORLD,&req[7]);
				MPI_Isend(&my, 1, MPI_INTEGER, prevRank, 0, MPI_COMM_WORLD,&req[8]);
				MPI_Irecv(&right, 1, MPI_INTEGER, nextRank, 0, MPI_COMM_WORLD,&req[9]);

				MPI_Waitall(4, &req[6], MPI_STATUSES_IGNORE);
				
				// left right
				int difLeft = ( left - my ) / 2;
				int difRight = ( my - right ) / 2;
				string sendBuf;
				string recvBuf;
				int cnt = 0;
				if( difRight > 0 && nextRank > rank ){
					sendBuf = data.substr( int( data.size() ) - difRight, difRight);
					data.resize( int( data.size() ) - difRight );
					MPI_Isend( &sendBuf[0], difRight, MPI_CHAR, nextRank, 0, MPI_COMM_WORLD, &req[10 + cnt]); cnt++;
				}
				if( difLeft > 0 && prevRank < rank ){
					recvBuf.resize( difLeft );
					MPI_Irecv( &recvBuf[0], difLeft, MPI_CHAR, prevRank, 0, MPI_COMM_WORLD, &req[10 + cnt]); cnt++;
				}
				MPI_Waitall(cnt, &req[10], MPI_STATUSES_IGNORE);
				data = recvBuf + data;
			}
		} // end data exchange 
		
		string newData;

		for( int i =0 ; i < data.size(); i++ ){
			char c = data[i] - 'a';
			double r = frand();
			if( r <= L1[c].first ) {
				newData += L1[c].second;
			} else {
				newData += data[i];
			}
		}
		data = newData;

		for ( int i = 0; i < size; i++ ){
			if ( i == rank ) {
			//	printf("data: iter %d rank %d load %d\n", iter, rank, int(data.size()) );
				cout << "debug: rank " << rank << " iter " << iter << " " << data << endl;
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	time = MPI_Wtime() - time;
	if( rank == 0 ) {
		printf("final: System L%d nProc %d nIter %d stepExchange %d time %f\n", systemId, size, numIter, exStep, time );
	}
	MPI_Finalize();
	return 0;
}
