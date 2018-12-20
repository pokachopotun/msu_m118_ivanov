#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <omp.h>
#include <cmath>

using namespace std;


void axpby(std::vector<double>& X, const std::vector<double>& Y, double a, double b){

	#pragma omp parallel for
	for(int i=0 ; i < X.size(); i++){
		X[i] = a * X[i] + b * Y[i];
	}
}

double dot(const std::vector<double>& X, const std::vector<double>& Y){

	double ans = 0;
	#pragma omp parallel for reduction(+ : ans)
	for(int i=0 ; i < X.size(); i++){
		ans += X[i] * Y[i];
	}
	return ans;
}

void CopyVector(const std::vector<double>& X, std::vector<double>& Y){
	#pragma omp parallel for
	for(int i=0 ; i < X.size(); i++){
		Y[i] = X[i];
	}
}

void SetVector(std::vector<double>& X, const double c){
	#pragma omp parallel for
	for(int i=0 ; i < X.size(); i++){
		X[i] = c;
	}
}


class Matrix{
public:
	Matrix(int N, int M, int K){
		Generate(N, M, K);
	}

	Matrix(const Matrix& rhs){
		_n = rhs._n;
		_m = rhs._m;
		
		_row.resize(_n + 1);
		_col.resize(_n);
		_val.resize(_n);

		#pragma omp parallel for
		for(int i =0 ; i <= _n; i++){
			_row[i] = i;
		}

		#pragma omp parallel for
		for(int i =0 ; i < _n; i++){
			_col[i] = i;
		}

		#pragma omp parallel for
		for(int i=0 ;i < rhs._n; i++){
			for(int j = rhs._row[i]; j < rhs._row[i + 1]; j++){
				if( rhs._col[j] == i){
					_val[i] = 1.0 / rhs._val[j];
					break;
				}
			}
		}
		_cnt = (int)_row.size() - 1;
	}

	double MakeVal(double col, double row){
		return sin(col + row + 1.0);
	}
	void Generate(int N, int M, int K){
		using namespace std;
		// cout << "Generation: " << endl;
		_cnt  = 0;
		_m = _n = N * M * K;
		_row.resize(_n + 1);
		for(int k = 0; k < K; k++){
			for(int j = 0; j < M; j++){
				for(int i = 0; i < N; i++){
					int id = k * (N * M) + j * N + i;
					// cout << id << endl;
					int cnt_cur = 
					+ int(k > 0) 
					+ int(j > 0)
					+ int(i > 0)
					+ 1
					+ int(i < N - 1)
					+ int(j < M - 1)
					+ int(k < K - 1);
					_row[id] = _cnt;
					_cnt += cnt_cur;
				}
			}
		}
		_row[_n] = _cnt;
		_col.resize(_cnt);
		_val.resize(_cnt);

		for(int k = 0; k < K; k++){
			for(int j = 0; j < M; j++){ 
				for(int i = 0; i < N; i++){
					int row = k * (N * M) + j * N + i;
					// cout << row << endl;
					int pos = _row[row];
					int cnt = 0;
					if(k > 0){
						int c = row - N * M;
						_col[pos + cnt] = c;
						_val[pos + cnt] = MakeVal(c, row);
						cnt++;
					}
					if(j > 0){
						int c = row - N;
						_col[pos + cnt] = c;
						_val[pos + cnt] = MakeVal(c, row);
						cnt++;
					}
					if(i > 0){
						int c = row - 1;
						_col[pos + cnt] = c;
						_val[pos + cnt] = MakeVal(c, row);
						cnt++;
					}
					_col[pos + cnt] = row;
					_val[pos + cnt] = MakeVal(i, row);
					cnt++;
					if(i < N - 1){
						int c = row + 1;
						_col[pos + cnt] = c;
						_val[pos + cnt] = MakeVal(c, row);
						cnt++;
					}
					if(j < M - 1){
						int c = row + N;
						_col[pos + cnt] = c;
						_val[pos + cnt] = MakeVal(c, row);
						cnt++;

					}
					if(k < K - 1){
						int c = row + N * M;
						_col[pos + cnt] = c;
						_val[pos + cnt] = MakeVal(c, row);
						cnt++;
					}
					pos = _row[row];	
					int pp = 0;
					double sum = 0;
					for(;pos < _row[row + 1]; pos++){
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
		for(int i =0 ;i < _n; i++){
			int pos = _row[i];
			for(int j =0 ; j < _m; j++){
				if(pos < _row[i + 1] && _col[pos] == j){
					cout << _val[pos];
					pos++;
				}else{
					cout << " ";
				}
				cout << " ";
			}
			cout << "," << endl;
		}
	}

	void PrintBuffers(){
		using namespace std;
		cout << "Row " << endl;
		for(int i =0 ;i < _n + 1; i++){
			cout << _row[i] << " ";
		}
		cout << endl;
		cout << "val " << endl;
		for(int i =0 ;i < _cnt; i++){
			cout << _val[i] << " ";
		}
		cout << endl;
		for(int i =0 ;i < _n; i++){
			for(int j = _row[i] ; j < _row[i + 1]; j++){
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
	
	void SpMV(const std::vector<double>& X, std::vector<double>& Y){
		SetVector(Y, 0.0);
		#pragma omp parallel for
		for(int i =0; i < _n; i++){
			for(int j = _row[i]; j < _row[i + 1]; j++){
				Y[i] += X[_col[j]] * _val[j];
			}
		}
	}  
	void Clean();

	int GetNonZeroCnt(){
		return _cnt;
	}
private:
	size_t _n, _m;
	int _cnt;
	std::vector<double> _val;
	std::vector<int> _row;
	std::vector<int> _col;
};


int Solver_BiCGSTAB( Matrix& A,  const int Nx, const int My,const int Kz, int nt,
	double tol = 1e-8, int maxit = 1000, bool info = false){
	omp_set_num_threads(nt);
	int N = Nx * My * Kz;
	double mineps = 1e-15;
	double RhoMin = 1e-60;
	double Rhoi_1=1.0, alphai = 1.0, wi = 1.0, betai_1 = 1.0, 
	Rhoi_2=1.0, alphai_1 = 1.0, wi_1 = 1.0;
	Matrix DD(A);
	// DD.PrintToTxt();
	std::vector<double> BB(N), XX(N), PP(N), 
		PP2(N), RR(N), RR2(N), TT(N), VV(N), SS(N), SS2(N);

	#pragma omp parallel for
	for(int i=0 ;i < BB.size(); i++){
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
			// cout << "RR" << hperf::dot(RR, RR) << endl;
			// cout << "PP" << hperf::dot(PP, PP) << endl;
		}
		else{
			betai_1=(Rhoi_1*alphai_1)/(Rhoi_2*wi_1);
			axpby(PP, RR, betai_1, 1.0); // p=r+betai_1*(p-w1*v)
			axpby(PP, VV, 1.0, -wi_1*betai_1);
		}

		DD.SpMV(PP,PP2);
		A.SpMV(PP2,VV);

		// DD.PrintToTxt();
		// cout << "PP2 " << hperf::dot(PP2, PP2) << endl;
		
		alphai = dot(RR2,VV);
		
		if(fabs(alphai)<RhoMin){
			// cout <<  fabs(alphai)  << " " << RhoMin << endl;
			// cout << "99" << endl;
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
	if(info) 
		printf("Solver_BiCGSTAB: outres: %g\n",res);
	return I;
}

void test_SpMV( const int Nx,const int My,const int Kz, const int nt, const int niter ){

	using namespace std;
	int N = Nx * My * Kz;
	Matrix A(Nx, My, Kz);
	vector<double> X(N), Y(N);
	omp_set_num_threads( nt );
	for(int i =0 ;i < X.size(); i++){
		Y[i] = X[i] = sin(static_cast<double>(i));
	}
	
	double time = omp_get_wtime();

	for(int i =0 ; i < niter; i++){
		A.SpMV(X, Y);
	}

	time = omp_get_wtime() - time;

	printf("System OpenMP Operation SpM N %d Nx %d Ny %d Nz %d P %d Px %d Py %d Pz %d Time %6.3f niter %d \n", N, Nx, My, Kz, nt, 1, 1, 1, time, niter);

}

void test_axpby( const int Nx,const int My,const int Kz, int nt, const int niter ){

	using namespace std;
	int N = Nx * My * Kz;

	vector<double> X(N), Y(N);
	for(int i =0 ;i < X.size(); i++){
		Y[i] = X[i] = sin(static_cast<double>(i));
	}

	omp_set_num_threads(nt);

	double time = omp_get_wtime();

	for(int i =0 ; i < niter; i++){
		axpby(X, Y, 1.2, 0.5);
	}

	time = omp_get_wtime() - time;

	printf("System OpenMP Operation xpy N %d Nx %d Ny %d Nz %d P %d Px %d Py %d Pz %d Time %6.3f niter %d \n", N, Nx, My, Kz, nt, 1, 1, 1, time, niter);

}

void test_dot( const int Nx,const int My,const int Kz, int nt, const int niter ){

	using namespace std;
	
	int N = Nx * My * Kz;

	omp_set_num_threads(nt);

	vector<double> X(N), Y(N);
	for(int i =0 ;i < X.size(); i++){
		Y[i] = X[i] = sin(static_cast<double>(i));
	}
		
	double time = omp_get_wtime();

	for(int i =0 ; i < niter; i++){
		dot(X, Y);
	}

	time = omp_get_wtime() - time;

	printf("System OpenMP Operation dot N %d Nx %d Ny %d Nz %d P %d Px %d Py %d Pz %d Time %6.3f niter %d \n", N, Nx, My, Kz, nt, 1, 1, 1, time, niter);

}

void test_solver( const int Nx, const int My,const int Kz, int nt,
	double tol = 1e-8, int maxit = 1000, bool info = false){

	using namespace std;
	
	int N = Nx * My * Kz;
	
	// printf("test SOLVER N = %d Nx=%d My=%d Kz = %d niter = %d \n", N, Nx, My, Kz, 1);
	omp_set_num_threads(nt);

	Matrix A(Nx, My, Kz);
	
	double time = omp_get_wtime();

	double last_iter = Solver_BiCGSTAB(A, Nx, My, Kz, nt, tol, maxit, 0);

	time = omp_get_wtime() - time;

	printf("System OpenMP Operation sol N %d Nx %d Ny %d Nz %d P %d Px %d Py %d Pz %d Time %6.3f niter %d \n", N, Nx, My, Kz, nt, 1, 1, 1, time, 1);

}


int main(int argc, char * argv[]) {
	
	if(argc < 9){
		cout << "use ./main N M K num_threads tol maxit info niter" << endl;
		return 0;
	}
	const int Nx = atoi(argv[1]);
	const int My = atoi(argv[2]);
	const int Kz = atoi(argv[3]);
	const int nt = atoi(argv[4]);
	const double tol = atof(argv[5]);
	const int maxit = atoi(argv[6]);
	const bool info = (bool)atoi(argv[7]);
	const int niter = atoi(argv[8]);

	test_SpMV(Nx, My, Kz, nt, niter);
	test_axpby(Nx, My, Kz, nt, niter);
	test_dot(Nx, My, Kz, nt, niter);

	test_solver(Nx, My, Kz, nt, tol, maxit, info);

	return 0;
}
