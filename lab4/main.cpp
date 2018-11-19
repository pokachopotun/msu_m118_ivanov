#include <mpi.h>

#include <iostream>
#include <string>
#include <vector>


using namespace std;

int main(int argc, char* argv[]){

	if(argc < 4){
		cout << "Use ./main A.in B.in C.out" << endl;
		return 0;
	}
	MPI_Init(&argc, &argv);

	const string AinputFileName(argv[1]);
	const string BinputFileName(argv[2]);
	const string CoutputFileName(argv[3]);

	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	int w = 0, h = 0;

	MPI_File AinputFile, BinputFile, CoutputFile;
	MPI_File_open(MPI_COMM_WORLD, AinputFileName.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &AinputFile);
	MPI_File_open(MPI_COMM_WORLD, BinputFileName.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &BinputFile);
	MPI_File_open(MPI_COMM_WORLD, CoutputFileName.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &CoutputFile);

	MPI_File_set_view(AinputFile, 0, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);

	MPI_File_read_all(AinputFile, &h, 1, MPI_INTEGER, MPI_STATUS_IGNORE);
	MPI_File_read_all(AinputFile, &w, 1, MPI_INTEGER, MPI_STATUS_IGNORE);


	printf("rank %d h %d w %d \n", rank, h, w);
	if(true){
		int hLocal = h / size; 
		int wLocal = w;
		int bwLocal = w / size;
		int localSize = hLocal * wLocal;

		double* A = new double[localSize];
		double* B = new double[bwLocal];
		double* C = new double[bwLocal];

		MPI_File_set_view(AinputFile, 8, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);
		MPI_File_read_at_all(AinputFile, 8 * localSize * rank, A, localSize, MPI_DOUBLE, MPI_STATUS_IGNORE);

		// double * tmp = new double[h * w];
		// for(int i =0 ;i < h; i++){
		// 	for(int j = 0; j < w; j++){
		// 		cout << tmp[i * w + j] << " ";
		// 	}
		// 	cout << endl;
		// }

		// return 0;

		for(int i = 0; i < bwLocal; i++){
			C[i]  = 0.0;
		}

		int rowStartLocal = rank * hLocal;
		int rowEndLocal = rowStartLocal + hLocal;

		int curRank = rank;
		for( int iter =0 ; iter < size; iter++, curRank = (curRank + 1 + size) % size ){

			MPI_File_set_view(BinputFile, 8, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);
			MPI_File_read_at_all(BinputFile, 8 * curRank * bwLocal, B, bwLocal, MPI_DOUBLE, MPI_STATUS_IGNORE);

			if(rank == 0){
				for(int i =0 ;i < bwLocal; i++){
					cout << B[i] << " ";
				}
				cout << endl;
			}

			int colStartLocal = curRank * bwLocal;
			int colEndLocal = colStartLocal + bwLocal;

			for(int row = 0; row < hLocal; row++){
				for(int col = colStartLocal; col < colEndLocal; col++ ){
					int aIndex = row * w + col;
					int cIndex = row;
					int bIndex = col - colStartLocal;
					C[cIndex] += B[bIndex] * A[aIndex];
				}
			}
		}

		MPI_File_set_view(CoutputFile, 0, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);
		if(rank == 0){
			int outRows = h;
			int outCols = 1;
			MPI_File_write_at(CoutputFile, 0, &outRows, 1, MPI_INTEGER, MPI_STATUS_IGNORE);
			MPI_File_write_at(CoutputFile, 4, &outCols, 1, MPI_INTEGER, MPI_STATUS_IGNORE);
		}
		MPI_File_set_view(CoutputFile, 8, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);
		MPI_File_write_at_all(CoutputFile, 8 * rank * bwLocal, C, bwLocal, MPI_DOUBLE, MPI_STATUS_IGNORE);
	}else{

	}

	MPI_Finalize();

	return 0;
}