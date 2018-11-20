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

	double time = MPI_Wtime();

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


	// printf("rank %d h %d w %d \n", rank, h, w);
	if( w <= h ){
		// printf("there\n");
		int hLocal = h / size; 
		int wLocal = w;
		int chLocal = hLocal;
		int localSize = hLocal * wLocal;


		double* A = new double[localSize];
		double* B = new double[wLocal];
		double* C = new double[chLocal];

		MPI_File_set_view(AinputFile, 8, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);
		MPI_File_read_at_all(AinputFile, 8 * localSize * rank, A, localSize, MPI_DOUBLE, MPI_STATUS_IGNORE);

		MPI_File_set_view(BinputFile, 8, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);
		MPI_File_read_at_all(BinputFile, 0, B, w, MPI_DOUBLE, MPI_STATUS_IGNORE);
		// for(int i =0 ; i < hLocal; i++){
		// 	for(int j =0 ; j < wLocal; j++){
		// 		cout << A[i * wLocal + j] <<  " ";
		// 	}
		// 	cout << endl;
		// }
		// cout << endl;
		// for(int i =0 ; i < wLocal; i++){
		// 	cout << B[i] <<  " ";
		// }
		// cout << endl;
		// for(int i = 0; i < hLocal; i++){
		// 	C[i]  = 0.0;
		// }
		MPI_Barrier(MPI_COMM_WORLD);
		time = MPI_Wtime();

		for(int row = 0; row < hLocal; row++){
			for(int col = 0; col < w; col++ ){
				int aIndex = row * w + col;
				int cIndex = row;
				int bIndex = col;
				C[cIndex] += B[bIndex] * A[aIndex];
			}
		}

		MPI_Barrier(MPI_COMM_WORLD);
		time = MPI_Wtime() - time;

		MPI_File_set_view(CoutputFile, 0, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);
		if(rank == 0){
			int outRows = h;
			int outCols = 1;
			MPI_File_write_at(CoutputFile, 0, &outRows, 1, MPI_INTEGER, MPI_STATUS_IGNORE);
			MPI_File_write_at(CoutputFile, 4, &outCols, 1, MPI_INTEGER, MPI_STATUS_IGNORE);
		}
		// MPI_Barrier(MPI_COMM_WORLD);
		// MPI_File_set_view(CoutputFile, 8, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);
		// printf("%d \n", bwLocal);
		MPI_File_write_at_all(CoutputFile, 8 + 8 * rank * chLocal, C, chLocal, MPI_DOUBLE, MPI_STATUS_IGNORE);
	}else{
		int hLocal = h; 
		int wLocal = w / size;
		int bwLocal = w / size;
		int localSize = hLocal * wLocal;

		double* A = new double[localSize];
		double* B = new double[bwLocal];
		double* C = new double[w];
		// double* C1 = new double[bwLocal];
		MPI_File_set_view(AinputFile, 8, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);
		for(int row = 0; row < h; row++ ){
			MPI_File_read_at_all(AinputFile, 8 * ( row * w + wLocal * rank ), &A[wLocal * row], wLocal, MPI_DOUBLE, MPI_STATUS_IGNORE);
		}
		MPI_File_set_view(BinputFile, 8, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);
		MPI_File_read_at_all(BinputFile, 8 * rank * bwLocal, B, bwLocal, MPI_DOUBLE, MPI_STATUS_IGNORE);


		for(int i = 0; i < w; i++){
			C[i]  = 0.0;
		}

		MPI_Barrier(MPI_COMM_WORLD);
		time = MPI_Wtime();

		for(int row = 0; row < h; row++){
			for(int col = 0; col < wLocal; col++ ){
				int aIndex = row * wLocal + col;
				int cIndex = row;
				int bIndex = col;
				C[cIndex] += B[bIndex] * A[aIndex];
			}
		}
		
		if(rank == 0){
			MPI_Reduce(MPI_IN_PLACE, C, w, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		}else{
			MPI_Reduce(C, 0, w, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		}

		time = MPI_Wtime() - time;
		
		MPI_File_set_view(CoutputFile, 0, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);
		if(rank == 0){
			int outRows = h;
			int outCols = 1;
			MPI_File_write_at(CoutputFile, 0, &outRows, 1, MPI_INTEGER, MPI_STATUS_IGNORE);
			MPI_File_write_at(CoutputFile, 4, &outCols, 1, MPI_INTEGER, MPI_STATUS_IGNORE);
		}
		
		if( rank == 0 ){
			MPI_File_write_at(CoutputFile, 8, C, w, MPI_DOUBLE, MPI_STATUS_IGNORE);
		}
	}

	if( rank == 0 )
		printf("Nodes %d rows %d cols %d Time %f\n", size, h, w, time);

	MPI_Finalize();

	return 0;
}