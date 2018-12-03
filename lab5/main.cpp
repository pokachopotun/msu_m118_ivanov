#include <mpi.h>
#include <iostream>

using namespace std;

int getPowerOfTwo(int a){
	int ans = 0;
	while(a){
		ans++;
		a >>= 1;
	}
	return ans;
}

int main(int argc, char* argv[]){
	if( argc < 4 ){
		cout << "use ./main A.in B.in C.out" << endl;
		return 0;
	}

	MPI_Init( &argc, &argv );

	const string AinputFileName(argv[1]);
	const string BinputFileName(argv[2]);
	const string CoutputFileName(argv[3]);

	int mCols, mRows;

	int worldRank;
	int worldSize;
	MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
	MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

	int pRow, pCol, pSlice, pSection;
	int nProc = 1 << ( getPowerOfTwo(worldSize) / 3 );

	{
		int sliceSize = nProc * nProc;
		pSlice = worldRank / sliceSize;
		int inSliceId = worldRank % sliceSize;
		pRow = inSliceId / nProc;
		pCol = inSliceId % nProc;
		pSection = pRow * nProc + pCol;
	}	

	MPI_Comm commRow, commCol, commSlice, commSection;

	MPI_Comm_split(MPI_COMM_WORLD, pSlice, worldRank, &commSlice);
	MPI_Comm_split(commSlice, pRow, worldRank, &commRow);
	MPI_Comm_split(commSlice, pCol, worldRank, &commCol);
	MPI_Comm_split(MPI_COMM_WORLD, pSection, worldRank, &commSection);


	MPI_File AinputFile, BinputFile, CoutputFile;

	if( pSlice == 0 ){
		MPI_File_open(commSlice, AinputFileName.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &AinputFile);
		MPI_File_open(commSlice, BinputFileName.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &BinputFile);
		MPI_File_open(commSlice, CoutputFileName.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &CoutputFile);

		MPI_File_set_view(AinputFile, 0, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);
		MPI_File_set_view(BinputFile, 0, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);

		MPI_File_read_all(AinputFile, &mRows, 1, MPI_INTEGER, MPI_STATUS_IGNORE);
		MPI_File_read_all(AinputFile, &mCols, 1, MPI_INTEGER, MPI_STATUS_IGNORE);
	}
	MPI_Bcast(&mRows, 1, MPI_INTEGER, 0, commSection);
	MPI_Bcast(&mCols, 1, MPI_INTEGER, 0, commSection);

	int localRows = mRows / nProc;
	int localCols = mCols / nProc;
	int localSize = localRows * localCols;

	// read in data
	double * A = new double[localSize];
	double * B = new double[localSize];
	double * C = new double[localSize];

	if( pSlice == 0 )
	{ // read in A
		int mShift = mCols * pRow * localRows + pCol * localCols + 1;
		MPI_File_set_view(AinputFile, 8 * mShift, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);
		for(int row = 0; row < localRows; row++ ){
			MPI_File_read_at_all(AinputFile, 8 * ( row * mCols ), &A[row * localCols], localCols, MPI_DOUBLE, MPI_STATUS_IGNORE);
		}
	}

	if( pSlice == 0 )
	{ // read in B
		int mShift = mCols * pRow * localRows + pCol * localCols + 1;
		MPI_File_set_view(BinputFile, 8 * mShift, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);
		for(int row = 0; row < localRows; row++ ){
			MPI_File_read_at_all(BinputFile, 8 * ( row * mCols ), &B[row * localCols], localCols, MPI_DOUBLE, MPI_STATUS_IGNORE);
		}
	}

	double time = MPI_Wtime();

	{ // init A
		if( pSlice == 0 && pCol > 0 ){
			MPI_Send(A, localSize, MPI_DOUBLE, pCol, 0, commSection);
		}

		if( pSlice > 0 && pSlice == pCol ) {
			MPI_Recv(A, localSize, MPI_DOUBLE, 0, 0, commSection, MPI_STATUS_IGNORE);
		}
		MPI_Bcast(A, localSize, MPI_DOUBLE, pSlice, commRow);
	}

	{ // init B
		if( pSlice == 0 && pRow > 0 ){
			MPI_Send(B, localSize, MPI_DOUBLE, pRow, 0, commSection);
		}

		if( pSlice > 0 && pSlice == pRow ) {
			MPI_Recv(B, localSize, MPI_DOUBLE, 0, 0, commSection, MPI_STATUS_IGNORE);
		}

		MPI_Bcast(B, localSize, MPI_DOUBLE, pSlice, commCol);
	}

	{ // init C with zeroes
		for( int i = 0; i < localRows; i++ ){
			for( int j = 0; j < localCols; j++ ){
				int posC = localCols * i + j;
				C[posC] = 0.0;
			}
		}
	}

	// multiply local
	for( int i = 0; i < localRows; i++ ){
		for( int j = 0; j < localCols; j++ ){
			for( int k = 0; k < localCols; k++ ){
				int posA = localCols * i + k;
				int posB = localCols * k + j;
				int posC = localCols * i + j;
				C[posC] += A[posA] * B[posB];
			}
		}
	}

	// reduce C to slice 0
	if( pSlice == 0 ){
		MPI_Reduce(MPI_IN_PLACE, C, localSize, MPI_DOUBLE, MPI_SUM, 0, commSection);
	}else{
		MPI_Reduce(C, 0, localSize, MPI_DOUBLE, MPI_SUM, 0, commSection);
	}
	// printf("%d \n", pSlice);
	if( pSlice == 0 )
	{ // write C to file
		MPI_File_set_view(CoutputFile, 0, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);
		if( worldRank == 0 ){
			MPI_File_write_at(CoutputFile, 0, &mRows, 1, MPI_INTEGER, MPI_STATUS_IGNORE);
			MPI_File_write_at(CoutputFile, 4, &mCols, 1, MPI_INTEGER, MPI_STATUS_IGNORE);
		}
		int mShift = mCols * pRow * localRows + pCol * localCols + 1;
		MPI_File_set_view(CoutputFile, 8 * mShift, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);
		for(int row = 0; row < localRows; row++ ){
			MPI_File_write_at_all(CoutputFile, 8 * ( row * mCols ), &C[row * localCols], localCols, MPI_DOUBLE, MPI_STATUS_IGNORE);
		}
	}

	time = MPI_Wtime() - time;

	if( worldRank == 0 ){
		printf("Nodes %d rows %d cols %d Time %f\n", worldSize, mRows, mCols, time);
	}

	MPI_Finalize();
	return 0;
}