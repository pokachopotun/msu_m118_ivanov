#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <mpi.h>

using namespace std;


template<typename T>
void printVector(const vector<T>& x, const string& name){
		cout << name << " ";
		for( int i =0 ;i < x.size(); i++ ){
			cout << x[i] << " ";
		}
		cout << endl;
}

void SetVector(std::vector<double>& X, const double c){
	for(int i=0 ; i < X.size(); i++){
		X[i] = c;
	}
}

int getOwner( int cellGlobal, int Nx, int Ny, int Nz, int Px, int Py, int Pz ){
	int x, y, z;
	int NxLocal = Nx / Px;
	int NyLocal = Ny / Py;
	int NzLocal = Nz / Pz;
	{
		int slice = Nx * Ny;
		z = cellGlobal / slice;
		int xy = cellGlobal % slice;
		y = xy / Nx;
		x = xy % Nx;
	}
	z /= NzLocal;
	y /= NyLocal;
	x /= NxLocal;

	return z * Px * Py + y * Px + x;
}

class Matrix{
public:
	Matrix(int N, int M, int K, int Px, int Py, int Pz){
		GenerateCSR(N, M, K, Px, Py, Pz);
		Preprocess();
	}
	// create inverse diagonal matrix
	Matrix(const Matrix& other){


		int commRank, commSize;
		MPI_Comm_rank(MPI_COMM_WORLD, &commRank);
		MPI_Comm_size(MPI_COMM_WORLD, &commSize);
		
		NxLocal = other.NxLocal;
		NyLocal = other.NyLocal;
		NzLocal = other.NzLocal;
		N = other.N;
		P = other.P;
		NLocal = other.NLocal;

		_owners = other._owners;
		GlobalRows = other.GlobalRows;
		_rowLocal = other._rowLocal;
		_col = other._col;
		_val = other._val;
		_G2L = other._G2L;

		recvFrom = other.recvFrom;
		sendTo = other.sendTo;

		
		for(int i=0 ;i < other._rowLocal.size() - 1; i++){
			for(int j = other._rowLocal[i]; j < other._rowLocal[i + 1]; j++){
				if( other._col[j] == GlobalRows[i] ){
					if( other._val[j] != 0.0 )
						_val[j] = 1.0 / other._val[j];
					else
						_val[j] = 0.0;	

				}
				else{
					_val[j] = 0.0;
				}
			}
		}
	}

	double MakeVal(double col, double row){
		return sin(col + row + 1.0);
	}
	void GenerateCSR(int Nx, int Ny, int Nz, int Px, int Py, int Pz){
		using namespace std;

		int commRank, commSize;
		MPI_Comm_rank(MPI_COMM_WORLD, &commRank);
		MPI_Comm_size(MPI_COMM_WORLD, &commSize);

		NxLocal = Nx / Px, NyLocal = Ny/Py, NzLocal = Nz/Pz;

		N = Nx * Ny * Nz;
		P = Px * Py * Pz;

		int px, py, pz;
		{
			int slice = Px * Py;
			pz = commRank / slice;
			int pxpy = commRank % slice;
			py = pxpy / Px;
			px = pxpy % Px;
		}

		int zRange[] = { pz * NzLocal, (pz+ 1) * NzLocal };
		int yRange[] = { py * NyLocal, (py+ 1) * NyLocal };
		int xRange[] = { px * NxLocal, (px+ 1) * NxLocal };

		int shifts[][3] = { {0, 0, -1}, {0, -1, 0}, {-1, 0, 0}, {0, 0, 0},
							 {1, 0, 0}, {0, 1, 0}, { 0, 0, 1}
		};

		_rowLocal.push_back( 0 );
		_owners.assign(N, -1);
		for(int k = zRange[0]; k < zRange[1]; k++){
			for(int j = yRange[0]; j < yRange[1]; j++){
				for(int i = xRange[0]; i < xRange[1]; i++){
					int kLocal = k - zRange[0];
					int jLocal = j - yRange[0];
					int iLocal = i - xRange[0];

					int rowGlobal = k * (Nx * Ny) + j * Nx + i;
					int rowLocal = kLocal * (NxLocal * NyLocal) + jLocal * NxLocal + iLocal;

					GlobalRows.push_back(rowGlobal);
					
					for( int sh = 0; sh < 7; sh++ ){
						int x = i + shifts[sh][0];
						int y = j + shifts[sh][1];
						int z = k + shifts[sh][2];

						int neighbourGlobal = z * Nx * Ny + y * Nx + x;
						if( neighbourGlobal < 0 || N <= neighbourGlobal )
							continue;

						if( !( 0 <= x && x < Nx )  || !( 0 <= y && y < Ny ) || !( 0 <= z && z < Nz ) ){
							continue;
						}

						int neighbourOwner = getOwner( neighbourGlobal, Nx, Ny, Nz, Px, Py, Pz );
						_owners[ neighbourGlobal ] = neighbourOwner;
						
						_col.push_back( neighbourGlobal );
						_val.push_back( MakeVal( neighbourGlobal, rowGlobal ) );
					}
					if( _col.size() > _rowLocal.back() ){
						_rowLocal.push_back( _col.size() );
					}
					int pp = -1;
					double sum = 0;
					for( int pos = _rowLocal[rowLocal]; pos < _rowLocal[rowLocal + 1]; pos++){
						if( _col[pos] == rowGlobal ){
							pp = pos;
						}else{
							sum += fabs(_val[pos]);
						}
					}
					if( pp != -1 )
						_val[pp] = 1.1 * sum;
				}
			}
		}
	}

	void Preprocess(){
		using namespace std;
		int commRank, commSize;
		MPI_Comm_rank(MPI_COMM_WORLD, &commRank);
		MPI_Comm_size(MPI_COMM_WORLD, &commSize);

		sendTo.resize(P);
		recvFrom.resize(P);
		_G2L.assign( N, -1 );

		for( int rowLocal = 0; rowLocal < GlobalRows.size(); rowLocal++){
			_G2L[ GlobalRows[rowLocal] ] = rowLocal;
		}
		for( int rowLocal = 0; rowLocal < _rowLocal.size() - 1; rowLocal++){
			for( int j = _rowLocal[ rowLocal ]; j < _rowLocal[ rowLocal + 1]; j++ ){
				int cell = _col[j];
				if( _G2L[ cell ] == -1 ){
					int owner = _owners[ cell ];
					// cout << owner << " ";4
					recvFrom[ owner ].push_back(cell);
					sendTo[owner].push_back( GlobalRows[ rowLocal ] );
				}
			}
		}
		for(int owner = 0; owner < recvFrom.size(); owner++ ){
			if( recvFrom[owner].size() == 0 ){
				continue;
			}
			for( int i = 0; i < recvFrom[owner].size(); i++ ){
				int cell = recvFrom[owner][i];
				_G2L[ cell ] = GlobalRows.size() + _halo.size();
				_halo.push_back( cell );
			}
		}

		for(int owner = 0; owner < sendTo.size(); owner++ ){
			if( sendTo[owner].size() == 0 ){
				continue;
			}
			for( int i = 0; i < sendTo[owner].size(); i++ ){
				int cell = sendTo[owner][i];
				_sendHalo.push_back( cell );
			}
		}
	}
	
	void SpMV( std::vector<double>& X, std::vector<double>& Y){
		SetVector(Y, 0.0);

		UpdateHalo( X );
		
		for(int rowLocal = 0; rowLocal < _rowLocal.size() - 1; rowLocal++){
			for(int j = _rowLocal[rowLocal]; j < _rowLocal[rowLocal + 1]; j++){
				Y[ rowLocal ] += X[ _G2L[ _col[j] ] ] * _val[ j ];
			}
		}

		UpdateHalo( Y );
		
	}


	double dot(const std::vector<double>& X, const std::vector<double>& Y ){

		double ans = 0;
		for( int j = 0; j < GlobalRows.size(); j++){
			ans += X[j] * Y[j];	
		}

		MPI_Allreduce(MPI_IN_PLACE, &ans, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		return ans;
	}

	void axpby(std::vector<double>& X, std::vector<double>& Y, double a, double b){
		UpdateHalo(X);
		UpdateHalo(Y);
		for(int i=0 ; i < X.size(); i++){
			X[i] = a * X[i] + b * Y[i];
		}
		UpdateHalo(X);
	}



	void CopyVector(std::vector<double>& X, std::vector<double>& Y){
		UpdateHalo(X);
		for(int i=0 ; i < X.size(); i++){
			Y[i] = X[i];
		}
	}

	void UpdateHalo(vector<double>& val){

		vector< MPI_Request > req( sendTo.size() + recvFrom.size() );
		int lastREQ = 0;

		vector< double*  > sendBuf( sendTo.size(), 0 );
		for( int owner = 0; owner < sendTo.size(); owner++){

			if( sendTo[owner].size() == 0)
			{
				continue;
			}
			int cnt = sendTo[owner].size();
			sendBuf[owner] = new double[cnt];
			for(int i =0 ; i < cnt; i++){
				int cellGlobal = sendTo[owner][i];
				int cellLocal = _G2L[cellGlobal];
				sendBuf[owner][i] = val[cellLocal];
			}
			MPI_Isend(sendBuf[owner], cnt, MPI_DOUBLE, owner, 0, MPI_COMM_WORLD, &req[lastREQ]);
			lastREQ++;
		}

		vector< double * > recvBuf(recvFrom.size(), 0 );
		// cout << ": ";
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
		// cout << ": ";
		// 		cout << "fuck" << endl;
		// MPI_Finalize();
		// exit(0);
		for( int owner = 0; owner < recvFrom.size(); owner++){
			if( recvFrom[owner].size() == 0)
			{
				continue;
			}
			// cout << owner << ": ";
			for(int i =0 ; i < recvFrom[owner].size(); i++){
				int cellGlobal = recvFrom[owner][i];
				int cellLocal = _G2L[ cellGlobal ];
				val[cellLocal] = recvBuf[owner][i];
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
	// void Clean();

	// int GetNonZeroCnt(){
	// 	return _nonZerocnt;
	// }

	int NxLocal, NyLocal, NzLocal;
	int N, P, NLocal;
	std::vector<double> _val;
	std::vector<int> _rowLocal;
	std::vector<int> _col, _halo, _sendHalo;
	std::vector<int> _owners, GlobalRows, _G2L;
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

	int commRank, commSize;
	MPI_Comm_rank(MPI_COMM_WORLD, &commRank);
	MPI_Comm_size(MPI_COMM_WORLD, &commSize);
	int N = Nx * Ny * Nz;
	int NLoc = N / Px / Py / Pz;
	int ni = A.GlobalRows.size();
	int nh = A._halo.size();
	int NLocal =  ni + nh;
	double mineps = 1e-15;
	double RhoMin = 1e-60;
	double Rhoi_1=1.0, alphai = 1.0, wi = 1.0, betai_1 = 1.0, 
	Rhoi_2=1.0, alphai_1 = 1.0, wi_1 = 1.0;
	Matrix DD(A);
	std::vector<double> BB(NLocal), XX(NLocal), PP(NLocal), 
		PP2(NLocal), RR(NLocal), RR2(NLocal), TT(NLocal), VV(NLocal), SS(NLocal), SS2(NLocal);

	for(int i=0 ;i < ni; i++){
		BB[ i ] = sin(static_cast<double>(A.GlobalRows[i]));
	}

	for(int i= 0; i < nh; i++){
		BB[ ni + i ] = sin(static_cast<double>(A._halo[i]));
	}

	A.CopyVector(BB, RR);
	A.CopyVector(BB, RR2);

	double initres = sqrt(A.dot(RR,RR));
	double eps = max(mineps, tol*initres);
	double res = initres;
	int I;
	for(I=0; I<maxit; I++){
		if(info && commRank == 0 ) 
			printf("It %d: res = %e tol=%e\n", I,res, res/initres);
		
		if(res<eps){
			break;
		}
		if(res>initres/mineps){
			return -1;
		}
		if(I==0) 
			Rhoi_1 = initres*initres;
		else 
			Rhoi_1 = A.dot(RR2,RR);
		if(fabs(Rhoi_1)<RhoMin){
			cout << "83" << endl;
			return -1;
		}
		if(I==0){
			A.CopyVector(RR, PP);
		}
		else{
			betai_1=(Rhoi_1*alphai_1)/(Rhoi_2*wi_1);
			A.axpby(PP, RR, betai_1, 1.0); // p=r+betai_1*(p-w1*v)
			A.axpby(PP, VV, 1.0, -wi_1*betai_1);
		}
				
		DD.SpMV(PP,PP2);
		A.SpMV(PP2,VV);
		alphai = A.dot(RR2,VV);
		if(fabs(alphai)<RhoMin){
			cout << "-3" << endl;
			return -3;
		}
		alphai = Rhoi_1/alphai;
		A.CopyVector(RR, SS);
		A.axpby(SS, VV, 1.0, -alphai);
		DD.SpMV(SS, SS2);
		A.SpMV(SS2, TT);
		wi = A.dot(TT, TT);

		if(fabs(wi)<RhoMin){
			cout << "112" << endl;
			return -4;
		}
		
		wi = A.dot(TT, SS) / wi;
		if(fabs(wi)<RhoMin) {
			cout << "117" << endl;
			return -5;
		}
		// x=x+alphai*p2+wi*s2
		A.axpby(XX, PP2,1.0,alphai);
		A.axpby(XX, SS2,1.0,wi);
		
		//RR=SS; // r=s-wi*t
		A.CopyVector(SS, RR);
		
		A.axpby(RR, TT, 1.0, -wi);
		alphai_1=alphai;
		Rhoi_2=Rhoi_1;
		wi_1=wi;
		res = sqrt(A.dot(RR,RR));
	}
	if(info && commRank == 0 ) 
		printf("Solver_BiCGSTAB: outres: %g\n",res);
	return I;
}

void test_SpMV( const int Nx, const int Ny,const int Nz,
	const int Px, const int Py,const int Pz, int niter){

	using namespace std;
	
	int N = Nx * Ny * Nz;
	int P = Px * Py * Pz;

	int commRank, commSize;
	MPI_Comm_rank(MPI_COMM_WORLD, &commRank);
	MPI_Comm_size(MPI_COMM_WORLD, &commSize);

	Matrix A(Nx, Ny, Nz, Px, Py, Pz);

	int ni = A.GlobalRows.size();
	int nh = A._halo.size();
	int NLocal =  ni + nh;

	std::vector<double> XX(NLocal), YY(NLocal);

	for(int i=0 ;i < ni; i++){
		XX[ i ] = sin(static_cast<double>(A.GlobalRows[i]));
	}

	for(int i= 0; i < nh; i++){
		XX[ ni + i ] = sin(static_cast<double>(A._halo[i]));
	}
	A.CopyVector(XX, YY);

	double time = MPI_Wtime();

	for( int i = 0; i < niter; i++){
		A.SpMV(XX,YY);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	time = MPI_Wtime() - time;

	if( commRank == 0 )
		printf("System MPI Operation SpM N %d Nx %d Ny %d Nz %d P %d Px %d Py %d Pz %d Time %6.3f niter %d\n", N, Nx, Ny, Nz, P, Px, Py, Pz, time, niter);
}

void test_axpby( const int Nx, const int Ny,const int Nz,
	const int Px, const int Py,const int Pz, int niter){

	using namespace std;
	
	int N = Nx * Ny * Nz;
	int P = Px * Py * Pz;

	int commRank, commSize;
	MPI_Comm_rank(MPI_COMM_WORLD, &commRank);
	MPI_Comm_size(MPI_COMM_WORLD, &commSize);

	Matrix A(Nx, Ny, Nz, Px, Py, Pz);

	int ni = A.GlobalRows.size();
	int nh = A._halo.size();
	int NLocal =  ni + nh;

	std::vector<double> XX(NLocal), YY(NLocal);

	for(int i=0 ;i < ni; i++){
		XX[ i ] = sin(static_cast<double>(A.GlobalRows[i]));
	}

	for(int i= 0; i < nh; i++){
		XX[ ni + i ] = sin(static_cast<double>(A._halo[i]));
	}
	A.CopyVector(XX, YY);

	double time = MPI_Wtime();

	for( int i = 0; i < niter; i++){
		A.axpby(XX,YY, 1.1, 1.1 );
	}

	MPI_Barrier(MPI_COMM_WORLD);

	time = MPI_Wtime() - time;

	if( commRank == 0 )
		printf("System MPI Operation xpy N %d Nx %d Ny %d Nz %d P %d Px %d Py %d Pz %d Time %6.3f niter %d\n", N, Nx, Ny, Nz, P, Px, Py, Pz, time, niter);
}

void test_dot( const int Nx, const int Ny,const int Nz,
	const int Px, const int Py,const int Pz, int niter){

	using namespace std;
	
	int N = Nx * Ny * Nz;
	int P = Px * Py * Pz;

	int commRank, commSize;
	MPI_Comm_rank(MPI_COMM_WORLD, &commRank);
	MPI_Comm_size(MPI_COMM_WORLD, &commSize);

	Matrix A(Nx, Ny, Nz, Px, Py, Pz);

	int ni = A.GlobalRows.size();
	int nh = A._halo.size();
	int NLocal =  ni + nh;

	std::vector<double> XX(NLocal), YY(NLocal);

	for(int i=0 ;i < ni; i++){
		XX[ i ] = sin(static_cast<double>(A.GlobalRows[i]));
	}

	for(int i= 0; i < nh; i++){
		XX[ ni + i ] = sin(static_cast<double>(A._halo[i]));
	}
	A.CopyVector(XX, YY);

	double time = MPI_Wtime();

	double sum = 0;
	for( int i = 0; i < niter; i++){
		A.dot(XX,YY);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	time = MPI_Wtime() - time;

	if( commRank == 0 )
		printf("System MPI Operation dot N %d Nx %d Ny %d Nz %d P %d Px %d Py %d Pz %d Time %6.3f niter %d\n", N, Nx, Ny, Nz, P, Px, Py, Pz, time, niter);
}

void test_solver( const int Nx, const int Ny,const int Nz,
	const int Px, const int Py,const int Pz,
	double tol = 1e-8, int maxit = 1000, bool info = false){

	using namespace std;
	
	int N = Nx * Ny * Nz;
	int P = Px * Py * Pz;

	int commRank, commSize;
	MPI_Comm_rank(MPI_COMM_WORLD, &commRank);
	MPI_Comm_size(MPI_COMM_WORLD, &commSize);

	Matrix A(Nx, Ny, Nz, Px, Py, Pz);

	double time = MPI_Wtime();

	Solver_BiCGSTAB(A, Nx, Ny, Nz, Px, Py, Pz, tol, maxit, false);

	MPI_Barrier(MPI_COMM_WORLD);

	time = MPI_Wtime() - time;

	if( commRank == 0 )
		printf("System MPI Operation sol N %d Nx %d Ny %d Nz %d P %d Px %d Py %d Pz %d Time %6.3f niter %d \n", N, Nx, Ny, Nz, P, Px, Py, Pz, time, 1);

}


int main(int argc, char * argv[]) {
	
	if(argc < 11){
		cout << "use ./main Nx Ny Nz Px Py Pz tol maxit info niter" << endl;
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
	const int niter = atoi(argv[10]);

	test_SpMV(Nx, Ny, Nz, Px, Py, Pz, niter);
	test_axpby(Nx, Ny, Nz, Px, Py, Pz, niter);
	test_dot(Nx, Ny, Nz, Px, Py, Pz, niter);
	test_solver(Nx, Ny, Nz, Px, Py, Pz, tol, maxit, info);

	MPI_Finalize();
	return 0;
}
