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

void axpby(std::vector<double>& X, const std::vector<double>& Y, double a, double b){
	for(int i=0 ; i < X.size(); i++){
		X[i] = a * X[i] + b * Y[i];
	}
}

double dot(const std::vector<double>& X, const std::vector<double>& Y ){

	double ans = 0;
	// double sum = 0;
	for(int i=0 ; i < X.size(); i++){
		// cout << X[i] << " " << Y[i] << endl;
		ans += X[i] * Y[i];
	}

	// MPI_Allreduce(MPI_IN_PLACE, &ans, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	// cout << "fuck" << endl;
	// MPI_Finalize();
	// exit(0);
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

	// cout << cellGlobal << " " << z << " " << y << " " << x << " " << z * Px * Py + y * Px + x << endl;
	return z * Px * Py + y * Px + x;
}

// void GetCoords( int cell, int Nx, int Ny, int Nz, int& x, int& y, int& z){
// 	int slice = Nx * Ny;
// 	z = cell / slice;
// 	int xy = cell % slice;
// 	y = xy / Nx;
// 	x = xy / Ny;
// }

class Matrix{
public:
	Matrix(int N, int M, int K, int Px, int Py, int Pz){
		Generate(N, M, K, Px, Py, Pz);
		// UpdateHalo( _val );
		// printVector( _val, "_val");
		// MPI_Finalize();
		// exit( 0 );
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
		NLocal = other.NLocal;
		Halo = other.Halo;
		_nRows = other._nRows;
		_nCols = other._nCols;
		_nRowsLocal = other._nRowsLocal;
		_nonZerocnt = other._nonZerocnt;

		_rowLocal = other._rowLocal;
		_haloStart = other._haloStart;
		_col = other._col;
		_val = other._val;

		// _owners = other._owners;
		_L2G = other._L2G ;
		_G2L = other._G2L;
 
		_nRows = other._nRows;
		_nCols = other._nCols;
		_nRowsLocal = other._nRowsLocal;

		recvFrom = other.recvFrom;
		sendTo = other.sendTo;

		
		for(int i=0 ;i < other._rowLocal.size() - 1; i++){
			for(int j = other._rowLocal[i]; j < other._rowLocal[i + 1]; j++){
				if( other._col[j] == _L2G[i] ){
					// cout << "FUCK " << other._val[j] << endl;
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
	void Generate(int Nx, int Ny, int Nz, int Px, int Py, int Pz){
		using namespace std;

		int commRank, commSize;
		MPI_Comm_rank(MPI_COMM_WORLD, &commRank);
		MPI_Comm_size(MPI_COMM_WORLD, &commSize);

		NxLocal = Nx / Px, NyLocal = Ny/Py, NzLocal = Nz/Pz;
		NLocal = NxLocal * NyLocal * NzLocal;
		Halo = NLocal;
		N = Nx * Ny * Nz;
		P = Px * Py * Pz;

		// _owners.assign( N, -1);
		// for(int pk = 0; pk < Pz; pk++){
		// 	for(int pj = 0; pj < Py; pj++){
		// 		for(int pi = 0; pi < Px; pi++){
					
		// 			int pid = pk * Px * Py + pj * Px + pi;

		// 			for(int k = pk * NzLocal; k < ( pk + 1)* NzLocal; k++){
		// 				for(int j = pj * NyLocal; j < ( pj + 1)* NyLocal; j++){
		// 					for(int i = pi * NxLocal; i < ( pi + 1) * NxLocal; i++){
		// 						int sid = k * (Nx * Ny) + j * Nx + i;
		// 						_owners[sid] = pid;
		// 					}
		// 				}
		// 			}
		// 		}
		// 	}
		// }
		// cout << " commRank: ";
		// for( int i = 0; i < _owners.size(); i++ ){
		// 	cout << _owners[i] << " ";
		// }
		// cout << endl;
		// MPI_Finalize();
		// exit(0);

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
		// _rowLocal.resize(_nRowsLocal + 1);
		_L2G.assign( N, -1 );
		_G2L.assign( N, -1 );

		int zRange[] = { pz * NzLocal, (pz+ 1) * NzLocal };
		int yRange[] = { py * NyLocal, (py+ 1) * NyLocal };
		int xRange[] = { px * NxLocal, (px+ 1) * NxLocal };
		// for(int k = kRange[0]; k < kRange[1]; k++){
		// 	for(int j = yRange[0]; j < yRange[1]; j++){
		// 		for(int i = xRange[0]; i < xRange[1]; i++){
		// 			int kLocal = k - kRange[0];
		// 			int jLocal = j - yRange[0];
		// 			int iLocal = i - xRange[0];
		// 			int idGlobal = k * (Nx * Ny) + j * Nx + i;

		// 			int idLocal = kLocal * (NxLocal * NyLocal) + jLocal * NxLocal + iLocal;

		// 			int cnt_cur = 
		// 			+ int(k > 0) 
		// 			+ int(j > 0)
		// 			+ int(i > 0)
		// 			+ 1
		// 			+ int(i < Nx - 1)
		// 			+ int(j < Ny - 1)
		// 			+ int(k < Nz - 1);
		// 			_rowLocal[idLocal] = _nonZerocnt;
		// 			_nonZerocnt += cnt_cur;
		// 		}
		// 	}
		// }

		// _rowLocal.back() = _nonZerocnt;
		// _col.resize(_rowLocal.back());
		// _val.resize(_rowLocal.back());
		// cout << _val.size() << endl;
		// printVector( _rowLocal, "rowLocal");
		// MPI_Finalize();
		// exit( 0 );
		int shifts[][3] = { {0, 0, -1}, {0, -1, 0}, {-1, 0, 0}, {0, 0, 0},
							 {1, 0, 0}, {0, 1, 0}, { 0, 0, 1}
		};
		// vector< int > shifts[3];
		// shifts.push_back(  );
		// shifts.push_back( -Nx );
		// shifts.push_back( -1 );
		// shifts.push_back( 0 );
		// shifts.push_back( 1 );
		// shifts.push_back( Nx );
		// shifts.push_back( Nx * Ny );

		recvFrom.resize(P);
		sendTo.resize(P);

		_rowLocal.push_back( 0 );
		_haloStart.push_back( 0 );
		for(int k = zRange[0]; k < zRange[1]; k++){
			for(int j = yRange[0]; j < yRange[1]; j++){
				for(int i = xRange[0]; i < xRange[1]; i++){
					int kLocal = k - zRange[0];
					int jLocal = j - yRange[0];
					int iLocal = i - xRange[0];


					int rowGlobal = k * (Nx * Ny) + j * Nx + i;
					int rowLocal = kLocal * (NxLocal * NyLocal) + jLocal * NxLocal + iLocal;

					_L2G[rowLocal] = rowGlobal;
					_G2L[rowGlobal] = rowLocal;
					
					// _rowLocal[ rowLocal ] = _rowLocal[ rowLocal - 1 ];
					// cout << row << endl;
					// int rowStartPos = _rowLocal[rowLocal];
					// int cnt = 0;
					for( int sh = 0; sh < 7; sh++ ){
						int x = i + shifts[sh][0];
						int y = j + shifts[sh][1];
						int z = k + shifts[sh][2];
						if( !( 0 <= x && x < Nx ) ){
							continue;
						}
						if( !( 0 <= y && y < Ny ) ){
							continue;
						}
						if( !( 0 <= z && z < Nz ) ){
							continue;
						}
						int neighbourGlobal = z * Nx * Ny + y * Nx + x;
						if( neighbourGlobal < 0 || N <= neighbourGlobal )
							continue;

						if( commRank == getOwner( neighbourGlobal, Nx, Ny, Nz, Px, Py, Pz ) ){
							_col.push_back( neighbourGlobal );
							_val.push_back( MakeVal( neighbourGlobal, rowGlobal ) );
						}
					}
					int cnt = 0;

					for( int sh = 0; sh < 7; sh++ ){
						int x = i + shifts[sh][0];
						int y = j + shifts[sh][1];
						int z = k + shifts[sh][2];
						if( !( 0 <= x && x < Nx ) ){
							continue;
						}
						if( !( 0 <= y && y < Ny ) ){
							continue;
						}
						if( !( 0 <= z && z < Nz ) ){
							continue;
						}
						int neighbourGlobal = z * Nx * Ny + y * Nx + x;
						if( neighbourGlobal < 0 || N <= neighbourGlobal )
							continue;
						
						int owner = getOwner( neighbourGlobal, Nx, Ny, Nz, Px, Py, Pz );
						if( commRank != owner ){
							Halo++;
							cnt++;
							_col.push_back( neighbourGlobal );
							_val.push_back( MakeVal( neighbourGlobal, rowGlobal ) );

							recvFrom[ owner ].push_back( neighbourGlobal );
							sendTo[ owner ].push_back( rowGlobal );	

							_G2L[neighbourGlobal] = _col.size() - 1;
							_L2G[_col.size() - 1] = neighbourGlobal;
						}
					}


					// if(k > 0){
					// 	int c = rowGlobal - Nx * Ny;
					// 	_col[rowStartPos + cnt] = c;
					// 	_val[rowStartPos + cnt] = MakeVal(c, rowGlobal);
					// 	cnt++;
					// }
					// if(j > 0){
					// 	int c = rowGlobal - Nx;
					// 	_col[rowStartPos + cnt] = c;
					// 	_val[rowStartPos + cnt] = MakeVal(c, rowGlobal);
					// 	cnt++;
					// }
					// if(i > 0){
					// 	int c = rowGlobal - 1;
					// 	_col[rowStartPos + cnt] = c;
					// 	_val[rowStartPos + cnt] = MakeVal(c, rowGlobal);
					// 	cnt++;
					// }
					// _col[rowStartPos + cnt] = rowGlobal;
					// _val[rowStartPos + cnt] = MakeVal(rowGlobal, rowGlobal);
					// cnt++;
					// if(i < Nx - 1){
					// 	int c = rowGlobal + 1;
					// 	_col[rowStartPos + cnt] = c;
					// 	_val[rowStartPos + cnt] = MakeVal(c, rowGlobal);
					// 	cnt++;
					// }
					// if(j < Ny - 1){
					// 	int c = rowGlobal + Nx;
					// 	_col[rowStartPos + cnt] = c;
					// 	_val[rowStartPos + cnt] = MakeVal(c, rowGlobal);
					// 	cnt++;

					// }
					// if(k < Nz - 1){
					// 	int c = rowGlobal + Nx * Ny;
					// 	_col[rowStartPos + cnt] = c;
					// 	_val[rowStartPos + cnt] = MakeVal(c, rowGlobal);
					// 	cnt++;
					// }

					if( _col.size() > _rowLocal.back() ){
						_rowLocal.push_back( _col.size() );
						_haloStart.push_back( _col.size() - cnt );
						// printVector( _rowLocal, "rowLocal");
					}
					int pp = -1;
					double sum = 0;
					for( int pos = _rowLocal[rowLocal] ;pos < _rowLocal[rowLocal + 1]; pos++){
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
		// cout << commRank << ": ";
		// // for( int i = 0 ; i < sendTo.size(); i++ ){
		// // 	printVector( sendTo[i], "sendTO");
		// // }
		// cout << endl;
		// // printVector( _L2G, "L2G");
		// printVector( _G2L, "G2L");
		// printVector( _rowLocal, "_rowLocal");
		
		// printVector( _col, "_col");
		// printVector( _val, "_val");
		// MPI_Finalize();
		// exit(0);



		// for( int rowLocal = 0; rowLocal < _rowLocal.size() - 1; rowLocal++){
		// 	int rowStartPos = _rowLocal[rowLocal];
		// 	for( int pos = rowStartPos; pos < _rowLocal[rowLocal + 1]; pos++){
		// 		int cellGlobal = _col[pos];

		// 		int owner = getOwner( cellGlobal, Nx, Ny, Nz, Px, Py, Pz );

		// 		if( _G2L[cellGlobal] == -1 ){

		// 			Halo++;
		// 		}
		// 	}
		// }

		// // add  receivable (halo) cells to our matrix
		// for( int owner = 0 ; owner < recvFrom.size(); owner++ ){
		// 	if( recvFrom[owner].size() == 0 || owner == commRank){
		// 		continue;
		// 	}
		// 	// cout << commRank << " " << owner << ": ";

		// 	_rowLocal.push_back( _rowLocal.back() + recvFrom[ owner ].size() );
		// 	for( int i =0 ; i < recvFrom[owner].size(); i++){
		// 		int cellGlobal = recvFrom[owner][i];
		// 		_G2L[cellGlobal] = _col.size();
		// 		_L2G[_col.size()] = cellGlobal;
		// 		// _col.push_back( cellGlobal );
		// 		// _val.push_back( 0.0 );
		// 	}
		// 	// cout << endl;

		// }

		// printVector( _L2G, "L2G");
		// printVector( _G2L, "G2L");
		// printVector( _rowLocal, "_rowLocal");
		
		// printVector( _col, "_col");
		// printVector( _val, "_val");
		// MPI_Finalize();
		// exit(0);
		// MPI_Finalize();
		// 	exit(0);
		
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
			// cout << _L2G[ rowLocal ] << " += ";
			for(int j = _rowLocal[rowLocal]; j < _rowLocal[rowLocal + 1]; j++){
				Y[ _L2G[rowLocal] ] += X[ _col[j] ] * _val[ j ];
				// cout << " + " << _col[j]<< " * " << _col[j] ;
			}
			// cout << endl;
		}
		MPI_Allreduce( MPI_IN_PLACE, Y.data(), (int)Y.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
		// UpdateHalo(Y);
		// printVector(X, "X");
		// printVector(_val, "val");
		// printVector(Y, "Y");
		// // printVector( Y, "Y");
		// MPI_Finalize();
		// exit(0);
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
		
	}

	void UpdateHalo(vector<double>& val){
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
				int cellGlobal = sendTo[owner][i];
				int cellLocal = _G2L[cellGlobal];
				sendBuf[owner][i] = val[cellLocal];
				// for( int j = _rowLocal[cellLocal]; j < _rowLocal[cellLocal + 1]; j++){
				// 	if( _col[j] == cellGlobal ){
				// 		sendBuf[owner][i] = _val[j];
				// 	}
				// }
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
				// cout << recvBuf[owner][i] << " ";
				// // int r = sid / N;
				// // if( rowLocal >= _rowLocal.size()){
				// // 	cout << "FUCK" << endl;
				// // 		MPI_Finalize();
				// // 		exit(0);
				// // }
				// // r = _G2L[r];
				// // int c = sid % N;
				// // cout << rowLocal << " ";
				// for( int j = _rowLocal[rowLocal]; j < _rowLocal[rowLocal + 1]; j++){
				// 	if( _col[j] == sid ){
						
				// 		// cout << recvBuf[owner][i] << endl;
				// 	}
				// }
			}

		}

		// cout << endl;
		// printVector( val, "VAL");
		// MPI_Finalize();
		// exit(0);	// cout << endl;
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
	int N, P, NLocal, Halo;
	size_t _nRows, _nCols, _nRowsLocal;
	int _nonZerocnt;
	std::vector<double> _val;
	std::vector<int> _rowLocal, _haloStart;
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
	// printVector( A._rowLocal, "_rowLocal");
	// printVector( A._val, "_val");
	// printVector( A._col, "_col");
	// MPI_Finalize();
	// exit(0);
	// cout << "Nloc " << NLoc << endl;
	// cout << " row local size " <<  A._rowLocal.size() << endl;
	// int NHalo = 
	// int NLocal = A._rowLocal.size() - 1;
	int NLocal = N;
	double mineps = 1e-15;
	double RhoMin = 1e-60;
	double Rhoi_1=1.0, alphai = 1.0, wi = 1.0, betai_1 = 1.0, 
	Rhoi_2=1.0, alphai_1 = 1.0, wi_1 = 1.0;
	Matrix DD(A);
	// DD.PrintBuffers();
	// MPI_Finalize();
	// exit(0);

	std::vector<double> BB(NLocal), XX(NLocal), PP(NLocal), 
		PP2(NLocal), RR(NLocal), RR2(NLocal), TT(NLocal), VV(NLocal), SS(NLocal), SS2(NLocal);


	// printVector( A._L2G, "L2G");
	// MPI_Finalize();
	// exit(0);
	for(int i=0 ;i < BB.size(); i++){
		BB[ i ] = sin( static_cast<double>(i + 1));// sin(static_cast<double>(A._L2G[i]));
		// BB[ i ] = A._L2G[i];
	}
	// printVector(BB, "BB");
	// printVector(A._L2G, "L2G");
	// MPI_Finalize();
	// exit(0);
	int nit;

	CopyVector(BB, RR);
	CopyVector(BB, RR2);
	// printVector(RR, "RR");
	// printVector(RR2, "RR2");

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

		// MPI_Finalize();
		// exit(0);

		alphai = Rhoi_1/alphai;

		// printVector(RR, "RR");
		// SS=RR; // s=r-alphai*v
		CopyVector(RR, SS);
		// printVector(SS, "SS");
		// cout << alphai << endl;
		axpby(SS, VV, 1.0, -alphai);

		// printVector(SS, "SS");

		DD.SpMV(SS, SS2);


		
		A.SpMV(SS2, TT);

		// cout << "TT ";
		// for( int i =0 ;i < TT.size(); i++ ){
		// 	cout << TT[i] << " ";
		// }
		// cout << endl;

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
