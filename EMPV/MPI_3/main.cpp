#include <cstdlib>
#include <mpi.h>
#include <iostream>
#include <string>
using namespace std;


void setxy( char * data, int gx, int gy, int fieldSize ){
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	const int localLines = fieldSize / size;

	int localSize = localLines * fieldSize;
	int gpos = gx * fieldSize + gy;
	int trank = gpos / localSize;
	if( rank == trank ){
		int pos = gpos % localSize;
		data[fieldSize + pos] = 1;
	//	if(rank == 9)
		for ( int j = 0; j < 3; j++ ) {
		for( int i = 0; i < fieldSize; i++ ) {
			cout << int(data[j * fieldSize + i]) << " ";
		}
		cout << endl;
		}
		cout << "set " << gx << " " << gy << " at " << rank << " " << pos << endl; 
	}
}

int main( int argc, char* argv[] ){
	
	MPI_Init(&argc, &argv);	
	if( argc < 5 ) {
		cout << "Use ./main fieldSize numIter writeFile pattern" << endl;
		MPI_Finalize();
		return 0;
	}
	const int fieldSize = atoi( argv[1] );
	const int numIter = atoi( argv[2] );
	const string writeFile( argv[3] );
	const string pattern(argv[4]);

	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	if ( fieldSize % size != 0 ) {
		cout << "Make sure fieldSize is evenly divided by P" << endl;
		MPI_Finalize();
		return 0;
	}

	const int localLines = fieldSize / size;

	int localSize = localLines * fieldSize;
	char * data = (char * )calloc( localSize + 2 * fieldSize, sizeof(char) );
	char * newdata = (char * )calloc( localSize + 2 * fieldSize, sizeof(char) );
	for( int i = 0; i < localSize + 2 * fieldSize; i++ ) 
		data[i] = 0;	
	
	//init	
	if( pattern == "glyder"){
		setxy(data, 0, 0, fieldSize);
		setxy(data, 2, 2, fieldSize);
		setxy(data, 3, 3, fieldSize);
		setxy(data, 4, 3, fieldSize);
		setxy(data, 4, 2, fieldSize);
		setxy(data, 4, 1, fieldSize);
		setxy(data, 9, 9, fieldSize);
		setxy(data, 9, 8, fieldSize);
		setxy(data, 9, 7, fieldSize);	
	}else{
		for( int i = 0; i < localSize + 2 * fieldSize; i++ ){
			data[i] = static_cast<char>(rand() % 2);
		}
	}


	int nextRank = ( rank + 1 ) % size;
	int prevRank = ( rank - 1 + size ) % size;

	int prevLine = 0;
	int nextLine = fieldSize + fieldSize * localLines;
	int myFirst = fieldSize;
	int myLast = nextLine - fieldSize;
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	double time = MPI_Wtime();

	for( int iter = 0; iter < numIter; iter++ ){
//		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Request reqs[4];
		MPI_Isend(&data[myFirst], fieldSize, MPI_CHAR, prevRank, 0, MPI_COMM_WORLD,&reqs[0]);
		MPI_Isend(&data[myLast], fieldSize, MPI_CHAR, nextRank, 0, MPI_COMM_WORLD, &reqs[1]);
		MPI_Irecv(&data[prevLine], fieldSize, MPI_CHAR, prevRank, 0, MPI_COMM_WORLD, &reqs[2]);
		MPI_Irecv(&data[nextLine], fieldSize, MPI_CHAR, nextRank, 0, MPI_COMM_WORLD, &reqs[3]);
		MPI_Waitall(4, reqs, MPI_STATUSES_IGNORE);

		for( int row = 1 ; row < localLines + 1; row++ )
		{
			for( int  col = 0; col < fieldSize; col++ )
			{
				int alive = 0;
				int dead = 0;
				int val = data[ row * fieldSize + col ];
				for( int i = -1; i <= 1; i++ )
				{
					for( int j = -1; j <= 1; j++ )
					{
						if( i == 0 && i == j ) continue;

						int x = row + i; 
						int y = ( col + j + fieldSize ) % fieldSize ; 
						int pos = fieldSize * x + y;
						int val = data[pos];
						if( val == 0 ) 
							dead++;
						else 
							alive++;
					}
				}
				char res = 0;
				if( val == 0 ){ 
					if ( alive == 3 ) res = 1;
				} else {
					if ( alive == 2 || alive == 3 ) res = 1; else res = 0;
				}
				newdata[row * fieldSize + col] = res;
			}
		}
		
		swap( data, newdata );
	}
	MPI_Barrier(MPI_COMM_WORLD);
	time = MPI_Wtime() - time;
	if( rank == 0 ) {
		printf("fieldSize %d nProc %d nIter %d time %f\n", fieldSize, size, numIter, time );
	}
	MPI_File f;
	if ( writeFile.size() > 1 ) {
		MPI_File_open( MPI_COMM_WORLD, writeFile.data(), MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &f );
		MPI_File_set_view( f, 0, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);
		MPI_File_write_at_all( f, localSize * rank, &data[myFirst], localSize, MPI_CHAR,  MPI_STATUS_IGNORE);
		cout << rank << " local size " << localSize << endl;
	//	cout << "rank " << rank << ": ";
	/*	for( int i = 0; i < localLines; i++ ){
			for( int j = 0; j < fieldSize; j++ ) {
				cout << int(data[myFirst + i * fieldSize + j]);
			}
			cout << endl;
		} */
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_File_close(&f);
	delete[] data;
	delete[] newdata;
	MPI_Finalize();
	return 0;
}
