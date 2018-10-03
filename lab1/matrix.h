#ifndef MATRIX_H_20092018
#define MATRIX_H_20092018
#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <chrono>
#include <limits>

class Matrix{
public:
	Matrix(const std::string& fileName){
		std::ifstream ifs(fileName, std::ios::binary);

		ifs.read(reinterpret_cast<char *>( &_type ), sizeof(char));
		ifs.read(reinterpret_cast<char *>( &_n ), sizeof(size_t)/sizeof(char));
		ifs.read(reinterpret_cast<char *>( &_m ), sizeof(size_t)/sizeof(char));
		size_t bufSize = BufferSize();
		// std::cout << bufSize << std::endl;
		_charBuf = new char[bufSize];
		ifs.read(_charBuf, bufSize);
	}

	size_t BufferSize() const;

	Matrix(size_t n, size_t m, char type, char* buf = nullptr)
		: _n(n), _m(m), _type(type)
	{
		size_t bufSize = BufferSize();
		_charBuf = new char[bufSize];
		if(buf != nullptr){
			for(size_t i = 0 ;i < bufSize; i++){
				_charBuf[i] = buf[i];
			}
		}else{
			for(size_t i = 0 ;i < bufSize; i++){
				_charBuf[i] = 0;
			}
		}
	}

	void MakeRandom();

	void ReadInTxt(const std::string& fileName);

	void MakeUno();

	size_t MatrixSize() const;

	size_t Cols() const;
	
	size_t Rows() const;

	friend std::ostream& operator<<(std::ostream& ofs, const Matrix& m);

	virtual ~Matrix(){
	}

	void WriteToFile(const std::string& outFileName);

	friend Matrix Multiply(const Matrix& lhs, const Matrix& rhs, const std::string& mode);
	friend bool Compare(const Matrix& lhs, const Matrix& rhs);
	
	size_t getPos(size_t row, size_t col) const;

	void clean();
protected:
	size_t _n, _m;
	char _type;
	char* _charBuf;
};

class FloatMatrix : public Matrix {
public:
	FloatMatrix(const Matrix& other) :  Matrix(other) {
		floatBuf = reinterpret_cast<float* >(_charBuf);
	}
	float& operator[](size_t pos){
		return floatBuf[pos];
	}
	const float& operator[](size_t pos) const{
		return floatBuf[pos];
	}
	float& operator()(size_t row, size_t col){
		return floatBuf[getPos(row, col)];
	}
	const float& operator()(size_t row, size_t col) const{
		return floatBuf[getPos(row, col)];
	}

private:
	float* floatBuf;
};

class DoubleMatrix :  public Matrix{
public:
	DoubleMatrix(const Matrix& other) : Matrix(other){
		doubleBuf = reinterpret_cast<double* >(_charBuf);
	}
	double& operator[](size_t pos){
		return doubleBuf[pos];
	}
	double& operator()(size_t row, size_t col){
		return doubleBuf[getPos(row, col)];
	}
	const double& operator[](size_t pos) const{
		return doubleBuf[pos];
	}
	const double& operator()(size_t row, size_t col) const{
		return doubleBuf[getPos(row, col)];
	}

private:
	double* doubleBuf;
};



class MyTimer{
public:
	 static std::chrono::high_resolution_clock::time_point t1;
	 static std::chrono::high_resolution_clock::time_point t2;
	 static double GetSeconds(){
	 	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
	 	return time_span.count();
	 }
	 static void Tick(){
	 	MyTimer::t1 = std::chrono::high_resolution_clock::now();
	 }
	 static void Tock(){
	 	MyTimer::t2 = std::chrono::high_resolution_clock::now();
	 }
};


std::ostream& operator<<(std::ostream& ofs, const Matrix& m);

FloatMatrix MultiplyFloat(const FloatMatrix& lhs, const FloatMatrix& rhs, const std::string& mode);

DoubleMatrix MultiplyDouble(const DoubleMatrix& lhs, const DoubleMatrix& rhs, const std::string& mode);

Matrix Multiply(const Matrix& lhs, const Matrix& rhs, const std::string& mode);

bool CompareDouble(const DoubleMatrix& lhs, const DoubleMatrix& rhs);

bool CompareFloat(const FloatMatrix& lhs, const FloatMatrix& rhs);

bool Compare(const Matrix& lhs, const Matrix& rhs);

#endif