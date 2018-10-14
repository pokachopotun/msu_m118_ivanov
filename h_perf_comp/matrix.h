#ifndef MATRIX_H_20092018
#define MATRIX_H_20092018
#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <chrono>
#include <limits>
#include <cmath>

class Matrix{
public:
	Matrix(int N, int M, int K){
		Generate(N, M, K);
	}

	size_t BufferSize() const;

	// Matrix(size_t n, size_t m, char type, char* buf = nullptr)
	// 	: _n(n), _m(m), _type(type)
	// {

	// }

	double MakeVal(double col, double row){
		return 1.0;//sin(col + row + 1.0);
	}
	void Generate(int N, int M, int K){
		using namespace std;
		// cout << "Generation: " << endl;
		_cnt  = 0;
		_m = _n = N * M * K;
		_row = (int*) calloc(_n + 1, sizeof(int));
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
		_col = (int*) calloc(_cnt, sizeof(int));
		_val = (double*) calloc(_cnt, sizeof(double));

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
			cout << endl;
		}
	}

	void PrintBuffers(){
		using namespace std;
		for(int i =0 ;i < _n + 1; i++){
			cout << _row[i] << " ";
		}
		cout << endl;
		for(int i =0 ;i < _n; i++){
			for(int j = _row[i] ; j < _row[i + 1]; j++){
				cout << _col[j]<< " ";
			}
			cout << endl;
		}
	}
	// void ReadInTxt(const std::string& fileName);

	// void MakeUno();

	size_t MatrixSize() const;

	size_t Cols() const;
	
	size_t Rows() const;

	friend std::ostream& operator<<(std::ostream& ofs, const Matrix& m);

	virtual ~Matrix(){
	}

	void WriteToFile(const std::string& outFileName);

	friend Matrix Multiply(const Matrix& lhs, const Matrix& rhs, const std::string& mode);
	friend bool Compare(const Matrix& lhs, const Matrix& rhs);
	
	//TODO
	size_t GetPos(size_t row, size_t col) const;

	void Clean();
private:
	size_t _n, _m;
	int _cnt;
	double* _val;
	int* _row;
	int* _col;
};

// std::ostream& operator<<(std::ostream& ofs, const Matrix& m);

// Matrix Multiply(const Matrix& lhs, const Matrix& rhs, const std::string& mode);

// bool Compare(const Matrix& lhs, const Matrix& rhs);

#endif
