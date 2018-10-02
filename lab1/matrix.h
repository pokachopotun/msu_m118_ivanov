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

	size_t BufferSize() const{
		if( _type == 'f' ){
			return MatrixSize() * sizeof(float)/sizeof(char);
		}else{
			return MatrixSize() * sizeof(double)/sizeof(char);
		}
	}

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

	void MakeRandom(){
		using namespace std;
		random_device rd;
		mt19937 gen(rd());
		uniform_real_distribution<float> dis(0, 1);
		if(_type == 'f'){
			float * floatBuf = reinterpret_cast< float * >(_charBuf);
			for(size_t i = 0 ; i < Rows(); i++){
				for(size_t j = 0; j < Cols(); j++){
					size_t pos = i * Cols() + j;
					floatBuf[pos] = static_cast<float>( dis(gen) );
				}	
			}
		}else{
			double * doubleBuf = reinterpret_cast< double * >(_charBuf);
			for(size_t i=0 ; i < Rows(); i++){
				for(size_t j = 0; j < Cols(); j++){
					size_t pos = i * Cols() + j;
					doubleBuf[pos] = static_cast<double>( dis(gen) );
				}	
			}
		}
	}

	void ReadInTxt(const std::string& fileName){
		using namespace std;
		ifstream txt(fileName);
		if(_type == 'f'){
			float * floatBuf = reinterpret_cast< float * >(_charBuf);
			for(size_t i = 0 ; i < Rows(); i++){
				for(size_t j = 0; j < Cols(); j++){
					size_t pos = i * Cols() + j;
					txt >> floatBuf[pos];
				}	
			}
		}else{
			double * doubleBuf = reinterpret_cast< double * >(_charBuf);
			for(size_t i=0 ; i < Rows(); i++){
				for(size_t j = 0; j < Cols(); j++){
					size_t pos = i * Cols() + j;
					txt >> doubleBuf[pos];
				}	
			}
		}
	}

	void MakeUno(){
		using namespace std;
		if(_type == 'f'){
			float * floatBuf = reinterpret_cast< float * >(_charBuf);
			for(size_t i = 0 ; i < Rows(); i++){
				for(size_t j = 0; j < Cols(); j++){
					size_t pos = i * Cols() + j;
					if(i == j){
						floatBuf[pos] = 1.0f;
					}else{
						floatBuf[pos] = 0.0f;
					}
				}	
			}
		}else{
			double * doubleBuf = reinterpret_cast< double * >(_charBuf);
			for(size_t i=0 ; i < Rows(); i++){
				for(size_t j = 0; j < Cols(); j++){
					size_t pos = i * Cols() + j;
					if(i == j){
						doubleBuf[pos] = 1.0;
					}else{
						doubleBuf[pos] = 0.0;
					}
				}	
			}
		}
	}

	size_t MatrixSize() const{
		return Cols() * Rows();
	}

	size_t Cols() const{
		return _m;
	}
	
	size_t Rows() const{
		return _n;
	}

	friend std::ostream& operator<<(std::ostream& ofs, const Matrix& m);

	virtual ~Matrix(){
		//delete[] _charBuf;
	}

	void WriteToFile(const std::string& outFileName){
		std::ofstream ofs(outFileName, std::ios::binary);
		ofs.write(&_type, sizeof(char));
		ofs.write(reinterpret_cast< char* >( &_n ), sizeof(size_t)/sizeof(char));
		ofs.write(reinterpret_cast< char* >( &_m ), sizeof(size_t)/sizeof(char));
		ofs.write(_charBuf, BufferSize() );
	}

	friend Matrix Multiply(const Matrix& lhs, const Matrix& rhs, const std::string& mode);
	friend bool Compare(const Matrix& lhs, const Matrix& rhs);
	
	size_t getPos(size_t row, size_t col) const {
		return row * Cols() + col;
	}

	void clean(){
		delete _charBuf;
	}

protected:
	size_t _n, _m;
	char _type;
	char* _charBuf;
};

class FloatMatrix : public Matrix {
public:
	FloatMatrix(const Matrix& other) :  Matrix(other) {
		// _n = other._n;
		// _m = other._m;
		// _type = other._type;
		// _charBuf = other._charBuf;
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
		// _n = other._n;
		// _m = other._m;
		// _type = other._type;
		// _charBuf = other._charBuf;
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


std::ostream& operator<<(std::ostream& ofs, const Matrix& m){
	if(m._type == 'f'){
		float * floatBuf = reinterpret_cast<float *> (m._charBuf);
		for(size_t i = 0 ; i < m.Rows(); i++){
			for(size_t j = 0; j < m.Cols(); j++){
				size_t pos = i * m.Cols() + j;
				ofs << floatBuf[pos] << " ";
			}	
			ofs << std::endl;
		}
	}else{
		double * doubleBuf = reinterpret_cast<double *> (m._charBuf);
		for(size_t i=0 ; i < m.Rows(); i++){
			for(size_t j = 0; j < m.Cols(); j++){
				size_t pos = i * m.Cols() + j;
				ofs << doubleBuf[pos] << " ";
			}	
			ofs << std::endl;
		}
	}
	return ofs;
}


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

std::chrono::high_resolution_clock::time_point MyTimer::t1 = std::chrono::high_resolution_clock::now();
std::chrono::high_resolution_clock::time_point MyTimer::t2 = std::chrono::high_resolution_clock::now();

FloatMatrix MultiplyFloat(const FloatMatrix& lhs, const FloatMatrix& rhs, const std::string& mode){
	
	Matrix tmp(lhs.Rows(), rhs.Cols(), 'f', nullptr);
	FloatMatrix res(tmp);
	// std::cout << "MultiplyFloat init done" << std::endl;
	// std::cout << res.Rows() << " " << res.Cols() << std::endl;
	if(mode == "ijk"){
		
		MyTimer::Tick();
		
		for(size_t i =0 ; i < lhs.Rows(); i++){
			for(size_t j = 0; j < rhs.Cols(); j++){
				for(size_t k = 0; k < lhs.Cols(); k++){
					res(i, j) += lhs(i, k) * rhs(k, j);
				}
			}
		}

		MyTimer::Tock();
		
  		return res;
	}
	if(mode == "jik"){

		MyTimer::Tick();
		for(size_t j = 0; j < rhs.Cols(); j++){
			for(size_t i =0 ; i < lhs.Rows(); i++){
				for(size_t k = 0; k < lhs.Cols(); k++){
					res(i, j) += lhs(i, k) * rhs(k, j);
				}
			}
		}
		MyTimer::Tock();
		return res;
	}
	if(mode == "ikj"){
		MyTimer::Tick();
		for(size_t i =0 ; i < lhs.Rows(); i++){
			for(size_t k = 0; k < lhs.Cols(); k++){
				float r = lhs(i, k);
				for(size_t j = 0; j < rhs.Cols(); j++){	
					res(i, j) += r * rhs(k, j);
				}
			}
		}
		MyTimer::Tock();
		return res;
	}
	if(mode == "kij"){
		MyTimer::Tick();
		for(size_t k = 0; k < lhs.Cols(); k++){
			for(size_t i =0 ; i < lhs.Rows(); i++){
				float r = lhs(i, k);
				for(size_t j = 0; j < rhs.Cols(); j++){	
					res(i, j) += r * rhs(k, j);
				}
			}
		}
		MyTimer::Tock();
		return res;
	}
	if(mode == "jki"){
		MyTimer::Tick();
		for(size_t j = 0; j < rhs.Cols(); j++){	
			for(size_t k = 0; k < lhs.Cols(); k++){
				float r = rhs(k, j);
				for(size_t i =0 ; i < lhs.Rows(); i++){
					res(i, j) += lhs(i, k) * r;
				}
			}
		}
		MyTimer::Tock();
		return res;
	}
	if(mode == "kji"){
		MyTimer::Tick();
		for(size_t k = 0; k < lhs.Cols(); k++){
			for(size_t j = 0; j < rhs.Cols(); j++){	
				float r = rhs(k, j);
				for(size_t i =0 ; i < lhs.Rows(); i++){
					res(i, j) += lhs(i, k) * r;
				}
			}
		}
		MyTimer::Tock();
		return res;
	}
	return res;
}

DoubleMatrix MultiplyDouble(const DoubleMatrix& lhs, const DoubleMatrix& rhs, const std::string& mode){
	Matrix tmp(lhs.Rows(), rhs.Cols(), 'd', nullptr);
	DoubleMatrix res(tmp);
	// std::cout << "MultiplyFloat init done" << std::endl;
	// std::cout << res.Rows() << " " << res.Cols() << std::endl;
	if(mode == "ijk"){
		MyTimer::Tick();
		for(size_t i =0 ; i < lhs.Rows(); i++){
			for(size_t j = 0; j < rhs.Cols(); j++){
				for(size_t k = 0; k < lhs.Cols(); k++){
					res(i, j) += lhs(i, k) * rhs(k, j);
				}
			}
		}
		MyTimer::Tock();
		return res;
	}
	if(mode == "jik"){
		MyTimer::Tick();
		for(size_t j = 0; j < rhs.Cols(); j++){
			for(size_t i =0 ; i < lhs.Rows(); i++){
				for(size_t k = 0; k < lhs.Cols(); k++){
					res(i, j) += lhs(i, k) * rhs(k, j);
				}
			}
		}
		MyTimer::Tock();
		return res;
	}
	if(mode == "ikj"){
		MyTimer::Tick();
		for(size_t i =0 ; i < lhs.Rows(); i++){
			for(size_t k = 0; k < lhs.Cols(); k++){
				double r = lhs(i, k);
				for(size_t j = 0; j < rhs.Cols(); j++){	
					res(i, j) += r * rhs(k, j);
				}
			}
		}
		MyTimer::Tock();
		return res;
	}
	if(mode == "kij"){
		MyTimer::Tick();
		for(size_t k = 0; k < lhs.Cols(); k++){
			for(size_t i =0 ; i < lhs.Rows(); i++){
				double r = lhs(i, k);
				for(size_t j = 0; j < rhs.Cols(); j++){	
					res(i, j) += r * rhs(k, j);
				}
			}
		}
		MyTimer::Tock();
		return res;
	}
	if(mode == "jki"){
		MyTimer::Tick();
		for(size_t j = 0; j < rhs.Cols(); j++){	
			for(size_t k = 0; k < lhs.Cols(); k++){
				double r = rhs(k, j);
				for(size_t i =0 ; i < lhs.Rows(); i++){
					res(i, j) += lhs(i, k) * r;
				}
			}
		}
		MyTimer::Tock();
		return res;
	}
	if(mode == "kji"){
		MyTimer::Tick();
		for(size_t k = 0; k < lhs.Cols(); k++){
			for(size_t j = 0; j < rhs.Cols(); j++){	
				double r = rhs(k, j);
				for(size_t i =0 ; i < lhs.Rows(); i++){
					res(i, j) += lhs(i, k) * r;
				}
			}
		}
		MyTimer::Tock();
		return res;
	}
	return res;
}


Matrix Multiply(const Matrix& lhs, const Matrix& rhs, const std::string& mode){
	char type = std::min(lhs._type, rhs._type);
	// std::cout << mode << std::endl;
	if(type == 'f'){
		return MultiplyFloat( FloatMatrix(lhs), FloatMatrix(rhs), mode);
	}else{
		return MultiplyDouble( DoubleMatrix(lhs), DoubleMatrix(rhs), mode);
	}
}


bool CompareDouble(const DoubleMatrix& lhs, const DoubleMatrix& rhs)
{
		for(size_t i =0 ; i < lhs.Rows(); i++){
			for(size_t j = 0; j < lhs.Cols(); j++){
				//std::cout << "looping" << std::endl;
				//std::cout << fabs( lhs(i, j) - rhs(i, j) ) << std::endl;
				if( fabs( lhs(i, j) - rhs(i, j) >= 1e-6 ) )
						return false;
			}
		}
		return true;
}

bool CompareFloat(const FloatMatrix& lhs, const FloatMatrix& rhs)
{
		for(size_t i =0 ; i < lhs.Rows(); i++){
			for(size_t j = 0; j < lhs.Cols(); j++){
				//std::cout << "looping" << std::endl;
				//std::cout << fabs( lhs(i, j) - rhs(i, j) ) << std::endl;
				if( fabs( lhs(i, j) - rhs(i, j) ) >= 1e-6 )
						return false;
			}
		}
		return true;
}

bool Compare(const Matrix& lhs, const Matrix& rhs){
	if(lhs._type != rhs._type || lhs.Cols() != rhs.Cols() || lhs.Rows() != rhs.Rows() ){
		return false;
	}
	char type = lhs._type;
	if(type == 'f'){
		return CompareFloat( FloatMatrix(lhs), FloatMatrix(rhs) );
	}else{
		return CompareDouble( DoubleMatrix(lhs), DoubleMatrix(rhs) );
	}
}


#endif