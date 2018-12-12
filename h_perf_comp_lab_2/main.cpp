#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <mpi.h>

using namespace std;

template< typename T > 
class MyVec{
public:
	MyVec(size_t _size = 0)
				 : size(_size), data( nullptr ){
				 	Assign(size, 0.0);
				 }
	size_t Size()
	{
		return size;
	}
	const size_t Size() const
	{
		return size;
	}
	void Resize(size_t newSize)
	{
		if( data != nullptr )
			delete[] data;
		size = newSize;
		if(size > 0){
			data = new T[newSize];
			for( int i = 0; i < Size(); i++){
				(*this)[i] = 0.0;
 			}
		}
		else
			data = nullptr;
	}
	void Assign(size_t newSize, T val)
	{
		Resize(newSize);
		for( int i = 0; i < Size(); i++){
			(*this)[i] = val;
 		}
	}
	T& operator[](size_t pos)
	{
		return data[pos];
	}
	const T& operator[](size_t pos) const
	{
		return data[pos];
	}
	~MyVec(){
		if( data != nullptr )
			delete[] data;
	}
private:
	size_t size;
	T* data;
};

void axpby(MyVec<double>& X, const MyVec<double>& Y, double a, double b){
	for(int i=0 ; i < X.Size(); i++){
		X[i] = a * X[i] + b * Y[i];
	}
}

double dot(const MyVec<double>& X, const MyVec<double>& Y){

	double ans = 0;
	for(int i=0 ; i < X.Size(); i++){
		ans += X[i] * Y[i];
	}
	MPI_Allreduce(MPI_IN_PLACE, &ans, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	// cout << "sadf " << ans << endl;
	return ans;
}

void CopyVector(const MyVec<double>& X, MyVec<double>& Y){
	for(int i=0 ; i < X.Size(); i++){
		Y[i] = X[i];
	}
}

void SetVector(MyVec<double>& X, const double c){
	for(int i=0 ; i < X.Size(); i++){
		X[i] = c;
	}
}

class Matrix{
public:
	Matrix(int N, int M, int K, int Px, int Py, int Pz){
		Generate(N, M, K, Px, Py, Pz);
	}
	// create inverse matrix
	Matrix(const Matrix& rhs){
		_nRows = rhs._nRows;
		_nCols = rhs._nCols;
		_nRowsLocal = rhs._nRowsLocal;
		
		_rowLocal.Resize(_nRowsLocal + 1);
		_col.Resize(_nRowsLocal);
		_val.Resize(_nRowsLocal);

		for(int i =0 ; i <= _nRowsLocal; i++){
			_rowLocal[i] = i;
		}

		for(int i =0 ; i < _nRowsLocal; i++){
			_col[i] = i;
		}

		for(int i=0 ;i < rhs._nRowsLocal; i++){
			for(int j = rhs._rowLocal[i]; j < rhs._rowLocal[i + 1]; j++){
				if( rhs._col[j] == i){
					_val[i] = 1.0 / rhs._val[j];
					break;
				}
			}
		}
		_nonZerocnt = (int)_rowLocal.Size() - 1;
	}

	double MakeVal(double col, double row){
		return sin(col + row + 1.0);
	}
	void Generate(int Nx, int Ny, int Nz, int Px, int Py, int Pz){
		using namespace std;

		int commRank, commSize;
		MPI_Comm_rank(MPI_COMM_WORLD, &commRank);
		MPI_Comm_size(MPI_COMM_WORLD, &commSize);

		int NxLocal = Nx / Px, NyLocal = Ny/Py, NzLocal = Nz/Pz;
		int N = Nx * Ny * Nz;
		_ownedRows.Assign( N, -1);


		_owners.Assign( N, -1);
		for(int k = 0; k < Nz; k++){
			for(int j = 0; j < Ny; j++){
				for(int i = 0; i < Nx; i++){
					int sid = k * (Nx * Ny) + j * Nx + i;
					int pid = k / NzLocal * (Px * Py) + j / NyLocal * Px + i / NxLocal;
					_owners[sid] = pid;
				}
			}
		}

		cout << commRank <<" owners: ";
		for( int i =0 ; i < N; i++){
			cout << _owners[i] << " ";
		}
		cout << endl;
		exit(0);
		int px, py, pz;
		{
			int slice = Px * Py;
			pz = commRank / slice;
			int pxpy = commRank % slice;
			py = pxpy / Py;
			px = pxpy % Py;
		}
		// cout << commRank << ": p " << px << " " << py << " " << pz << endl;
		_nRowsLocal = NxLocal * NyLocal * NzLocal;
		// cout << "Generation: " << endl;
		_nonZerocnt  = 0;
		_nCols = _nRows = Nx * Ny * Nz;
		_rowLocal.Resize(_nRowsLocal + 1);

		int kRange[]{ pz * NzLocal, (pz+ 1) * NzLocal };
		int yRange[]{ py * NyLocal, (py+ 1) * NyLocal };
		int xRange[]{ px * NxLocal, (px+ 1) * NxLocal };
		for(int k = kRange[0]; k < kRange[1]; k++){
			for(int j = yRange[0]; j < yRange[1]; j++){
				for(int i = xRange[0]; i < xRange[1]; i++){
					int kLocal = k - kRange[0];
					int jLocal = j - yRange[0];
					int iLocal = i - xRange[0];
					int idGlobal = k * (Nx * Ny) + j * Nx + i;

					int idLocal = kLocal * (NxLocal * NyLocal) + jLocal * NxLocal + iLocal;
					// cout << "global " << idGlobal << " local " << idLocal << endl;
					_ownedRows[idGlobal] = idLocal;
					// _L2G[idLocal] = idGlobal;
					// cout << id << endl;
					int cnt_cur = 
					+ int(k > 0) 
					+ int(j > 0)
					+ int(i > 0)
					+ 1
					+ int(i < Nx - 1)
					+ int(j < Ny - 1)
					+ int(k < Nz - 1);
					_rowLocal[idLocal] = _nonZerocnt;
					// _G2L[idGlobal] = idLocal;
					_nonZerocnt += cnt_cur;
				}
			}
		}
		// for(int i = 0; i < _ownedRows.Size(); i++){
		// 	cout << _ownedRows[i] << " ";
		// }
		// cout << endl;

		_rowLocal[_nRowsLocal] = _nonZerocnt;
		_col.Resize(_nonZerocnt);
		_val.Resize(_nonZerocnt);

		for(int k = kRange[0]; k < kRange[1]; k++){
			for(int j = yRange[0]; j < yRange[1]; j++){
				for(int i = xRange[0]; i < xRange[1]; i++){
					int row = k * (Nx * Ny) + j * Nx + i;
					int kLocal = k - kRange[0];
					int jLocal = j - yRange[0];
					int iLocal = i - xRange[0];
					int rowLocal = kLocal * (NxLocal * NyLocal) + jLocal * NxLocal + iLocal;
					// cout << row << endl;
					int rowStartPos = _rowLocal[rowLocal];
					int cnt = 0;
					if(k > 0){
						int c = row - Nx * Ny;
						_col[rowStartPos + cnt] = c;
						_val[rowStartPos + cnt] = MakeVal(c, row);
						cnt++;
					}
					if(j > 0){
						int c = row - Nx;
						_col[rowStartPos + cnt] = c;
						_val[rowStartPos + cnt] = MakeVal(c, row);
						cnt++;
					}
					if(i > 0){
						int c = row - 1;
						_col[rowStartPos + cnt] = c;
						_val[rowStartPos + cnt] = MakeVal(c, row);
						cnt++;
					}
					_col[rowStartPos + cnt] = row;
					_val[rowStartPos + cnt] = MakeVal(i, row);
					cnt++;
					if(i < Nx - 1){
						int c = row + 1;
						_col[rowStartPos + cnt] = c;
						_val[rowStartPos + cnt] = MakeVal(c, row);
						cnt++;
					}
					if(j < Ny - 1){
						int c = row + Nx;
						_col[rowStartPos + cnt] = c;
						_val[rowStartPos + cnt] = MakeVal(c, row);
						cnt++;

					}
					if(k < Nz - 1){
						int c = row + Nx * Ny;
						_col[rowStartPos + cnt] = c;
						_val[rowStartPos + cnt] = MakeVal(c, row);
						cnt++;
					}
					int pos = _rowLocal[rowLocal];	
					int pp = 0;
					double sum = 0;
					for(;pos < _rowLocal[rowLocal + 1]; pos++){
						if(_col[pos] == row){
							pp = pos;
						}else{
							sum += fabs(_val[pos]);
						}
					}
					_val[pp] = 1.1 * sum;
				}
			}
		}
	}



	void PrintToTxt(){
		using namespace std;
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		cout << rank << ": ";
		for(int i =0 ;i < _nRowsLocal; i++){
			int pos = _rowLocal[i];
			for(int j =0 ; j < _nCols; j++){
				if(pos < _rowLocal[i + 1] && _col[pos] == j){
					// cout << " 1";
					cout << _val[pos];
					pos++;
				}else{
					cout << " 0.00";
				}
				cout << " ";
			}
			cout << endl;
		}
	}

	void PrintOwning(){
		using namespace std;
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		cout << rank << ": ";
		for( int i =0 ; i < _ownedRows.Size(); i++){
			cout << _ownedRows[i] << " ";
		}
		cout << endl;
	}

	void PrintBuffers(){
		using namespace std;
		cout << "Row " << endl;
		for(int i =0 ;i < _nRows + 1; i++){
			cout << _rowLocal[i] << " ";
		}
		cout << endl;
		cout << "val " << endl;
		for(int i =0 ;i < _nonZerocnt; i++){
			cout << _val[i] << " ";
		}
		cout << endl;
		for(int i =0 ;i < _nRows; i++){
			for(int j = _rowLocal[i] ; j < _rowLocal[i + 1]; j++){
				cout << _col[j]<< " ";
			}
			cout << endl;
		}
	}

	size_t MatrixSize() const;

	size_t Cols() const;
	
	size_t Rows() const;

	friend std::ostream& operator<<(std::ostream& ofs, const Matrix& m);

	virtual ~Matrix(){
	}

	void WriteToFile(const std::string& outFileName);

	friend Matrix Multiply(const Matrix& lhs, const Matrix& rhs, const std::string& mode);
	friend bool Compare(const Matrix& lhs, const Matrix& rhs);
	
	void SpMV(const MyVec<double>& X, MyVec<double>& Y){
		SetVector(Y, 0.0);
		for(int i =0; i < _nRowsLocal; i++){
			for(int j = _rowLocal[i]; j < _rowLocal[i + 1]; j++){
				Y[i] += X[_col[j]] * _val[j];
			}
		}
	}  
	void Clean();

	int GetNonZeroCnt(){
		return _nonZerocnt;
	}
private:
	size_t _nRows, _nCols, _nRowsLocal;
	int _nonZerocnt;
	MyVec<double> _val;
	MyVec<int> _rowLocal;
	MyVec<int> _col;
	MyVec<int> _ownedRows;
	MyVec<int> _owners;
};


int Solver_BiCGSTAB( Matrix& A,  const int Nx, const int Ny,const int Nz,
								const int Px, const int Py, const int Pz,
	double tol = 1e-8, int maxit = 1000, bool info = false){

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int N = Nx * Ny * Nz;
	double mineps = 1e-15;
	double RhoMin = 1e-60;
	double Rhoi_1=1.0, alphai = 1.0, wi = 1.0, betai_1 = 1.0, 
	Rhoi_2=1.0, alphai_1 = 1.0, wi_1 = 1.0;
	Matrix DD(A);
	// DD.PrintToTxt();
	// exit(0);
	MyVec<double> BB(N), XX(N), PP(N), 
		PP2(N), RR(N), RR2(N), TT(N), VV(N), SS(N), SS2(N);

	for(int i=0 ;i < BB.Size(); i++){
		BB[i] = sin(static_cast<double>(i));
	}

	int nit;

	CopyVector(BB, RR);
	CopyVector(BB, RR2);

	double initres = sqrt(dot(RR,RR));
	double eps = max(mineps, tol*initres);
	//cout << "eps " << eps << endl; 
	double res = initres;

	int I;
	for(I=0; I<maxit; I++){
		if(info) 
			printf("It %d: res = %e tol=%e\n", I,res, res/initres);
		if(res<eps){
			// cout << "71" << endl;
			break;
		}
		if(res>initres/mineps){
			// cout << "75" << endl;
			return -1;
		}
		if(I==0) 
			Rhoi_1 = initres*initres;
		else 
			Rhoi_1 = dot(RR2,RR);
		if(fabs(Rhoi_1)<RhoMin){
			// cout << "83" << endl;
			return -1;
		}
		if(I==0){
			CopyVector(RR, PP);
			// cout << "RR" << dot(RR, RR) << endl;
			// cout << "PP" << dot(PP, PP) << endl;
		}
		else{
			betai_1=(Rhoi_1*alphai_1)/(Rhoi_2*wi_1);
			axpby(PP, RR, betai_1, 1.0); // p=r+betai_1*(p-w1*v)
			axpby(PP, VV, 1.0, -wi_1*betai_1);
		}
		// cout << "PP" << dot(PP, PP) << endl;
		// cout << "RR" << dot(RR, RR) << endl;
		DD.SpMV(PP,PP2);
		A.SpMV(PP2,VV);

		// DD.PrintToTxt();
		// cout << "PP2 " << dot(PP2, PP2) << endl;
		
		alphai = dot(RR2,VV);
		
		if(fabs(alphai)<RhoMin){
			return -3;
		}
		alphai = Rhoi_1/alphai;
		// SS=RR; // s=r-alphai*v
		CopyVector(RR, SS);
		axpby(SS, VV, 1.0, -alphai);
		DD.SpMV(SS, SS2);
		A.SpMV(SS2, TT);
		
		wi = dot(TT, TT);
		
		if(fabs(wi)<RhoMin){
			// cout << "112" << endl;
			return -4;
		}
		wi = dot(TT, SS) / wi;
		if(fabs(wi)<RhoMin) {
			// cout << "117" << endl;
			return -5;
		}
		// x=x+alphai*p2+wi*s2
		axpby(XX, PP2,1.0,alphai);
		axpby(XX, SS2,1.0,wi);
		
		//RR=SS; // r=s-wi*t
		CopyVector(SS, RR);
		axpby(RR, TT, 1.0, -wi);
		alphai_1=alphai;
		Rhoi_2=Rhoi_1;
		wi_1=wi;
		res = sqrt(dot(RR,RR));
	}
	if(info && rank == 0 ) 
		printf("Solver_BiCGSTAB: outres: %g\n",res);
	return I;
}

// void test_SpMV( const int Nx,const int My,const int Kz, const int niter ){

// 	using namespace std;
// 	int N = Nx * My * Kz;
// 	Matrix A(Nx, Ny, Nz, );
// 	vector<double> X(N), Y(N);
// 	for(int i =0 ;i < X.Size(); i++){
// 		Y[i] = X[i] = sin(static_cast<double>(i));
// 	}
// 	double totalFLOP = niter * ( N + A.GetNonZeroCnt() * 2 ) * 1e-9;
// 	printf("test SpMV N = %d Nx=%d My=%d Kz = %d niter = %d \n", N, Nx, My, Kz, niter);
// 	vector<double> time(4, 0);
// 	for(int t= 0; t < 4; t++){
// 		omp_set_nRowsum_threads(t + 1);
// 		time[t] = omp_get_wtime();

// 		for(int i =0 ; i < niter; i++){
// 			A.SpMV(X, Y);
// 		}

// 		time[t] = omp_get_wtime() - time[t];

// 		printf("SpMV threads = %d time = %6.3fs GFLOPS = %6.3f  speedup = %fX\n", t + 1, time[t], totalFLOP / time[t], time[0]/time[t]);
// 	}

// }

// void test_axpby( const int Nx,const int My,const int Kz, const int niter ){
	
	

// 	using namespace std;
	
// 	int N = Nx * My * Kz;

// 	vector<double> X(N), Y(N);
// 	for(int i =0 ;i < X.Size(); i++){
// 		Y[i] = X[i] = sin(static_cast<double>(i));
// 	}

// 	double totalFLOP = static_cast<double>( niter ) * static_cast<double>( N * 3 ) * 1e-9;
	
// 	printf("test axpby N = %d Nx=%d My=%d Kz = %d niter = %d \n", N, Nx, My, Kz, niter);
// 	vector<double> time(4, 0);
// 	for(int t= 0; t < 4; t++){
// 		omp_set_nRowsum_threads(t + 1);
// 		time[t] = omp_get_wtime();

// 		for(int i =0 ; i < niter; i++){
// 			axpby(X, Y, 1.2, 0.5);
// 		}

// 		time[t] = omp_get_wtime() - time[t];

// 		printf("axpby threads = %d time = %6.3fs GFLOPS = %6.3f  speedup = %fX\n", t + 1, time[t],  totalFLOP / time[t], time[0]/time[t]);
// 	}

// }

// void test_dot( const int Nx,const int My,const int Kz, const int niter ){
	
	

// 	using namespace std;
	
// 	int N = Nx * My * Kz;

// 	vector<double> X(N), Y(N);
// 	for(int i =0 ;i < X.Size(); i++){
// 		Y[i] = X[i] = sin(static_cast<double>(i));
// 	}

// 	double totalFLOP = static_cast<double>( niter ) * static_cast<double>( N * 2 ) * 1e-9;
	
// 	printf("test dot N = %d Nx=%d My=%d Kz = %d niter = %d \n", N, Nx, My, Kz, niter);
// 	vector<double> time(4, 0);
// 	for(int t= 0; t < 4; t++){
// 		omp_set_nRowsum_threads(t + 1);
// 		time[t] = omp_get_wtime();

// 		for(int i =0 ; i < niter; i++){
// 			dot(X, Y);
// 		}

// 		time[t] = omp_get_wtime() - time[t];

// 		printf("dot threads = %d time = %6.3fs GFLOPS = %6.3f  speedup = %fX\n", t + 1, time[t],  totalFLOP / time[t], time[0]/time[t]);
// 	}

// }

void test_solver( const int Nx, const int Ny,const int Nz,
	const int Px, const int Py,const int Pz,
	double tol = 1e-8, int maxit = 1000, bool info = false){

	using namespace std;
	
	int N = Nx * Ny * Nz;

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if( rank == 0 )
		printf("test SOLVER N %d Nx %d My %d Kz %d niter = %d \n", N, Nx, Ny, Nz, 1);

	Matrix A(Nx, Ny, Nz, Px, Py, Pz);

	// A.PrintBuffers();
	// exit(0);
	double time = MPI_Wtime();
	double last_iter = Solver_BiCGSTAB(A, Nx, Ny, Nz, Px, Py, Pz, tol, maxit, info);
	double FLOP = 1e-9 * ( last_iter * ( 34.0 * N + 8.0 * A.GetNonZeroCnt()) + 9.0 * A.GetNonZeroCnt());

	time = MPI_Wtime() - time;
	if( rank == 0 )
		printf("SOLVER nodes %d time %6.3fs GFLOPS %6.3f \n", Px * Py * Pz, time, FLOP/time);

}


int main(int argc, char * argv[]) {
	
	if(argc < 10){
		cout << "use ./main Nx Ny Nz Px Py Pz tol maxit info" << endl;
		return 0;
	}

	MPI_Init(&argc, &argv);

	const int Nx = atoi(argv[1]);
	const int Ny = atoi(argv[2]);
	const int Nz = atoi(argv[3]);
	

	const int Px = atoi(argv[4]);
	const int Py = atoi(argv[5]);
	const int Pz = atoi(argv[6]);

	const double tol = atof(argv[7]);
	const int maxit = atoi(argv[8]);
	const bool info = (bool)atoi(argv[9]);

	// test_SpMV(Nx, My, Kz, 100);
	// test_axpby(Nx, My, Kz, 1000);
	// test_dot(Nx, My, Kz, 1000);

	// Matrix A(Nx, Ny, Nz, Px, Py, Pz);

	// A.PrintToTxt();
	// A.PrintOwning();
	// A.PrintBuffers();

	test_solver(Nx, Ny, Nz, Px, Py, Pz, tol, maxit, info);

	MPI_Finalize();
	return 0;
}
