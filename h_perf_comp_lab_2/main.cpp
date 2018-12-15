#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <mpi.h>

using namespace std;

void axpby(std::vector<double>& X, const std::vector<double>& Y, double a, double b){
	for(int i=0 ; i < X.size(); i++){
		X[i] = a * X[i] + b * Y[i];
	}
}

double dot(const std::vector<double>& X, const std::vector<double>& Y){

	double ans = 0;
	for(int i=0 ; i < X.size(); i++){
		ans += X[i] * Y[i];
	}
	MPI_Allreduce(MPI_IN_PLACE, &ans, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	return ans;
}

void CopyVector(const std::vector<double>& X, std::vector<double>& Y){
	for(int i=0 ; i < X.size(); i++){
		Y[i] = X[i];
	}
}

void SetVector(std::vector<double>& X, const double c){
	for(int i=0 ; i < X.size(); i++){
		X[i] = c;
	}
}

class Matrix{
public:
	Matrix(int N, int M, int K, int Px, int Py, int Pz){
		Generate(N, M, K, Px, Py, Pz);
		UpdateHalo();
	}
	// create inverse matrixspmv
	Matrix(const Matrix& other){


		int commRank, commSize;
		MPI_Comm_rank(MPI_COMM_WORLD, &commRank);
		MPI_Comm_size(MPI_COMM_WORLD, &commSize);
		
		NxLocal = other.NxLocal;
		NyLocal = other.NyLocal;
		NzLocal = other.NzLocal;
		N = other.N;
		P = other.P;
		_nRows = other._nRows;
		_nCols = other._nCols;
		_nRowsLocal = other._nRowsLocal;
		_nonZerocnt = other._nonZerocnt;

		_rowLocal = other._rowLocal;
		_col = other._col;
		_val = other._val;

		_owners = other._owners;
		_L2G = other._L2G ;
		_G2L = other._G2L;

		_nRows = other._nRows;
		_nCols = other._nCols;
		_nRowsLocal = other._nRowsLocal;

		
		for(int i=0 ;i < other._rowLocal.size() - 1; i++){
			for(int j = other._rowLocal[i]; j < other._rowLocal[i + 1]; j++){
				if( other._col[j] == _L2G[i] ){
					// cout << "FUCK " << other._val[j] << endl;
					if( other._val[j] != 0.0 )
						_val[j] = 1.0 / other._val[j];
					else
						_val[j] = 0.0;	

				}else{
					_val[j] = 0.0;
				}
			}
		}
	}

	double MakeVal(double col, double row){
		return sin(col + row + 1.0);
	}
	void Generate(int Nx, int Ny, int Nz, int Px, int Py, int Pz){
		using namespace std;

		int commRank, commSize;
		MPI_Comm_rank(MPI_COMM_WORLD, &commRank);
		MPI_Comm_size(MPI_COMM_WORLD, &commSize);

		NxLocal = Nx / Px, NyLocal = Ny/Py, NzLocal = Nz/Pz;
		N = Nx * Ny * Nz;
		P = Px * Py * Pz;

		_owners.assign( N, -1);
		for(int k = 0; k < Nz; k++){
			for(int j = 0; j < Ny; j++){
				for(int i = 0; i < Nx; i++){
					int sid = k * (Nx * Ny) + j * Nx + i;
					_owners[sid] = commRank;
				}
			}
		}

		int px, py, pz;
		{
			int slice = Px * Py;
			pz = commRank / slice;
			int pxpy = commRank % slice;
			py = pxpy / Py;
			px = pxpy % Py;
		}
		
		_nRowsLocal = NxLocal * NyLocal * NzLocal;
		
		_nonZerocnt  = 0;
		_nCols = _nRows = Nx * Ny * Nz;
		_rowLocal.resize(_nRowsLocal + 1);
		_L2G.assign( N, -1 );
		_G2L.assign( N, -1 );

		int kRange[] = { pz * NzLocal, (pz+ 1) * NzLocal };
		int yRange[] = { py * NyLocal, (py+ 1) * NyLocal };
		int xRange[] = { px * NxLocal, (px+ 1) * NxLocal };
		for(int k = kRange[0]; k < kRange[1]; k++){
			for(int j = yRange[0]; j < yRange[1]; j++){
				for(int i = xRange[0]; i < xRange[1]; i++){
					int kLocal = k - kRange[0];
					int jLocal = j - yRange[0];
					int iLocal = i - xRange[0];
					int idGlobal = k * (Nx * Ny) + j * Nx + i;

					int idLocal = kLocal * (NxLocal * NyLocal) + jLocal * NxLocal + iLocal;
					_L2G[idLocal] = idGlobal;
					_G2L[idGlobal] = idLocal;
					int cnt_cur = 
					+ int(k > 0) 
					+ int(j > 0)
					+ int(i > 0)
					+ 1
					+ int(i < Nx - 1)
					+ int(j < Ny - 1)
					+ int(k < Nz - 1);
					_rowLocal[idLocal] = _nonZerocnt;
					_nonZerocnt += cnt_cur;
				}
			}
		}

		_rowLocal.back() = _nonZerocnt;
		_col.resize(_rowLocal.back());
		_val.resize(_rowLocal.back());

		for(int k = kRange[0]; k < kRange[1]; k++){
			for(int j = yRange[0]; j < yRange[1]; j++){
				for(int i = xRange[0]; i < xRange[1]; i++){
					int rowGlobal = k * (Nx * Ny) + j * Nx + i;
					int kLocal = k - kRange[0];
					int jLocal = j - yRange[0];
					int iLocal = i - xRange[0];
					int rowLocal = kLocal * (NxLocal * NyLocal) + jLocal * NxLocal + iLocal;
					// cout << row << endl;
					int rowStartPos = _rowLocal[rowLocal];
					int cnt = 0;
					if(k > 0){
						int c = rowGlobal - Nx * Ny;
						_col[rowStartPos + cnt] = c;
						_val[rowStartPos + cnt] = MakeVal(c, rowGlobal);
						cnt++;
					}
					if(j > 0){
						int c = rowGlobal - Nx;
						_col[rowStartPos + cnt] = c;
						_val[rowStartPos + cnt] = MakeVal(c, rowGlobal);
						cnt++;
					}
					if(i > 0){
						int c = rowGlobal - 1;
						_col[rowStartPos + cnt] = c;
						_val[rowStartPos + cnt] = MakeVal(c, rowGlobal);
						cnt++;
					}
					_col[rowStartPos + cnt] = rowGlobal;
					_val[rowStartPos + cnt] = MakeVal(i, rowGlobal);
					cnt++;
					if(i < Nx - 1){
						int c = rowGlobal + 1;
						_col[rowStartPos + cnt] = c;
						_val[rowStartPos + cnt] = MakeVal(c, rowGlobal);
						cnt++;
					}
					if(j < Ny - 1){
						int c = rowGlobal + Nx;
						_col[rowStartPos + cnt] = c;
						_val[rowStartPos + cnt] = MakeVal(c, rowGlobal);
						cnt++;

					}
					if(k < Nz - 1){
						int c = rowGlobal + Nx * Ny;
						_col[rowStartPos + cnt] = c;
						_val[rowStartPos + cnt] = MakeVal(c, rowGlobal);
						cnt++;
					}
					int pos = _rowLocal[rowLocal];	
					int pp = 0;
					double sum = 0;
					for(;pos < _rowLocal[rowLocal + 1]; pos++){
						if( _col[pos] == rowGlobal ){
							pp = pos;
						}else{
							sum += fabs(_val[pos]);
						}
					}
					_val[pp] = 1.1 * sum;
				}
			}
		}

		recvFrom.resize(P);
		sendTo.resize(P);

		for( int rowLocal = 0; rowLocal < _rowLocal.size() - 1; rowLocal++){
			int rowStartPos = _rowLocal[rowLocal];
			for( int pos = rowStartPos; pos < _rowLocal[rowLocal + 1]; pos++){
				int cellGlobal = _col[pos];
				int owner = _owners[ cellGlobal ];
				if( _G2L[cellGlobal] == -1 ){
					recvFrom[ owner ].push_back( cellGlobal );
					sendTo[ owner ].push_back( _L2G[ rowLocal ] );	
				}
			}
		}

		// add  receivable (halo) cells to our matrix
		for( int owner = 0 ; owner < recvFrom.size(); owner++ ){
			if( recvFrom[owner].size() == 0 || owner == commRank){
				continue;
			}
			_G2L[owner] = _rowLocal.size() - 1;
			_L2G[_rowLocal.size() - 1] = owner;

			_rowLocal.push_back( _rowLocal.back() + recvFrom[ owner ].size() );
			for( int i =0 ; i < recvFrom[owner].size(); i++){
				_col.push_back( recvFrom[owner][i] );
				_val.push_back( 0.0 );
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

	void PrintBuffers(){

		int commRank, commSize;
		MPI_Comm_rank(MPI_COMM_WORLD, &commRank);
		MPI_Comm_size(MPI_COMM_WORLD, &commSize);
		using namespace std;
		cout << "rank " << commRank << " RowSt ";
		for(int i =0 ;i < _rowLocal.size(); i++){
			cout << _rowLocal[i] << " ";
		}
		cout << "rank " << commRank << " Row ";
		for(int i =0 ;i < _rowLocal.size(); i++){
			cout << _L2G[i] << " ";
		}
		cout << " Col ";
		for(int i =0 ;i < _rowLocal.size() - 1; i++){
			for(int j = _rowLocal[i] ; j < _rowLocal[i + 1]; j++){
				cout << _col[j]<< " ";
			}
		}
		cout << "val ";
		for(int i =0 ;i < _val.size(); i++){
			cout << _val[i] << " ";
		}
		cout << endl;
		
	}

	friend std::ostream& operator<<(std::ostream& ofs, const Matrix& m);

	virtual ~Matrix(){
	}

	void WriteToFile(const std::string& outFileName);
	
	void SpMV(const std::vector<double>& X, std::vector<double>& Y){
		SetVector(Y, 0.0);

		for(int rowLocal =0; rowLocal < _rowLocal.size() - 1; rowLocal++){
			for(int j = _rowLocal[rowLocal]; j < _rowLocal[rowLocal + 1]; j++){
				Y[ rowLocal ] += X[ _G2L[ _col[j] ] ] * _val[j];
			}
		}
		// cout << "_rowLocal ";
		// for(int rowLocal =0; rowLocal < _rowLocal.size() - 1; rowLocal++){
		// 	cout << _rowLocal[ rowLocal ] << " ";
		// }
		// cout << " _col ";
		// for(int j = 0; j < _col.size(); j++){
		// 		cout << _col[j] << " ";
		// }

		// cout << " Y: ";
		// for( int i =0 ; i < Y.size(); i++ ){
		// 	cout << Y[i] << " ";
		// }
		// cout << " _G2L ";
		// for(int j = 0; j < _G2L.size(); j++){
		// 		cout << _G2L[j] << " ";
		// }
		// cout << " _val ";
		// for(int j = 0; j < _val.size(); j++){
		// 		cout << _val[j] << " ";
		// }
		// cout << endl;
		// MPI_Finalize();
		// exit(0);
		// MPI_Allreduce( MPI_IN_PLACE, Y.data(), (int)Y.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	}

	void UpdateHalo(){
		vector< MPI_Request > req( sendTo.size() + recvFrom.size() );
		int lastREQ = 0;

		vector< double * > sendBuf( sendTo.size(), 0 );
		for( int owner = 0; owner < sendTo.size(); owner++){
			if( sendTo[owner].size() == 0)
			{
				continue;
			}
			int cnt = sendTo[owner].size();
			sendBuf[owner] = new double[cnt];
			for(int i =0 ; i < sendTo[owner].size(); i++){
				int sid = sendTo[owner][i];
				int rowGlobal = sid / N;
				int rowLocal = _G2L[rowGlobal];
				int c = sid % N;
				for( int j = _rowLocal[rowLocal]; j < _rowLocal[rowLocal + 1]; j++){
					if( _col[j] == c ){
						sendBuf[owner][i] = _val[j];
						break;
					}
				}
			}
			MPI_Isend(sendBuf[owner], cnt, MPI_DOUBLE, owner, 0, MPI_COMM_WORLD, &req[lastREQ]);
			lastREQ++;
		}

		vector< double * > recvBuf(recvFrom.size(), 0 );
		for( int owner = 0; owner < recvFrom.size(); owner++){
			if( recvFrom[owner].size() == 0)
			{
				continue;
			}
			int cnt = recvFrom[owner].size();
			recvBuf[owner] = new double[cnt];
			MPI_Irecv( recvBuf[owner], cnt, MPI_DOUBLE, owner, 0, MPI_COMM_WORLD, &req[lastREQ]);
			lastREQ++;
		}

		MPI_Waitall( lastREQ, req.data(), MPI_STATUSES_IGNORE );
		
		for( int owner = 0; owner < recvFrom.size(); owner++){
			if( recvFrom[owner].size() == 0)
			{
				continue;
			}

			for(int i =0 ; i < recvFrom[owner].size(); i++){
				int sid = recvFrom[owner][i];
				int r = sid / N;
				r = _G2L[r];
				int c = sid % N;
				for( int j = _rowLocal[r]; j < _rowLocal[r + 1]; j++){
					if( _col[j] == c ){
						_val[j] = recvBuf[owner][i];
					}
				}
			}
		}
		for( int i =0; i < sendBuf.size(); i++ ){
			if( sendBuf[i] != 0 )
				delete[] sendBuf[i];
		}
		for( int i =0; i < recvBuf.size(); i++ ){
			if( recvBuf[i] != 0 )
				delete[] recvBuf[i];
		}	
	}
	void Clean();

	int GetNonZeroCnt(){
		return _nonZerocnt;
	}

	int NxLocal, NyLocal, NzLocal;
	int N, P;
	size_t _nRows, _nCols, _nRowsLocal;
	int _nonZerocnt;
	std::vector<double> _val;
	std::vector<int> _rowLocal;
	std::vector<int> _col;
	std::vector<int> _owners;
	std::vector<int> _L2G, _G2L;
	std::vector< std::vector<int> > sendTo, recvFrom;
};

void Message(int rank, const string& s){
	if( rank == 0){
		cout << s << endl;
	}
}


int Solver_BiCGSTAB( Matrix& A,  const int Nx, const int Ny,const int Nz,
								const int Px, const int Py, const int Pz,
	double tol = 1e-8, int maxit = 1000, bool info = false){

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int N = Nx * Ny * Nz;
	int NLoc = N / Px / Py / Pz;
	// int NHalo = 
	int NLocal = NLoc;
	double mineps = 1e-15;
	double RhoMin = 1e-60;
	double Rhoi_1=1.0, alphai = 1.0, wi = 1.0, betai_1 = 1.0, 
	Rhoi_2=1.0, alphai_1 = 1.0, wi_1 = 1.0;
	Matrix DD(A);

	std::vector<double> BB(NLocal), XX(NLocal), PP(NLocal), 
		PP2(NLocal), RR(NLocal), RR2(NLocal), TT(NLocal), VV(NLocal), SS(NLocal), SS2(NLocal);

	for(int i=0 ;i < BB.size(); i++){
		BB[ i ] = sin(static_cast<double>(A._L2G[i]));
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
		if(info && rank == 0 ) 
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
		}
		else{
			betai_1=(Rhoi_1*alphai_1)/(Rhoi_2*wi_1);
			axpby(PP, RR, betai_1, 1.0); // p=r+betai_1*(p-w1*v)
			axpby(PP, VV, 1.0, -wi_1*betai_1);
		}
		
		DD.SpMV(PP,PP2);
		A.SpMV(PP2,VV);

		alphai = dot(RR2,VV);
		
		if(fabs(alphai)<RhoMin){
			// cout << "-3" << endl;
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
