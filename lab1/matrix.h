#ifndef MATRIX_H_20092018
#define MATRIX_H_20092018
#include <iostream>
#include <fstream>
#include <string>
#include <random>

class Matrix{
public:
	Matrix(const std::string& fileName){
		std::ifstream ifs(fileName, std::ios::binary);

		ifs.read(reinterpret_cast<char *>( &_type ), sizeof(char));
		ifs.read(reinterpret_cast<char *>( &_n ), sizeof(size_t)/sizeof(char));
		ifs.read(reinterpret_cast<char *>( &_m ), sizeof(size_t)/sizeof(char));
		size_t bufSize = BufferSize();
		std::cout << bufSize << std::endl;
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

class DoubleMatrix : Matrix{
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

FloatMatrix MultiplyFloat(const FloatMatrix& lhs, const FloatMatrix& rhs, const std::string& mode){
	Matrix tmp(lhs.Rows(), rhs.Cols(), 'f', nullptr);
	FloatMatrix res(tmp);
	std::cout << "MultiplyFloat init done" << std::endl;
	std::cout << res.Rows() << " " << res.Cols() << std::endl;
	if(mode == "ijk"){
		for(size_t i =0 ; i < lhs.Rows(); i++){
			for(size_t j = 0; j < rhs.Cols(); j++){
				for(size_t k = 0; k < lhs.Cols(); k++){
					res(i, j) += lhs(i, k) * rhs(k, j);
				}
			}
		}
		return res;
	}
	if(mode == "jik"){
		for(size_t i =0 ; i < lhs.Rows(); i++){
			for(size_t j = 0; j < rhs.Cols(); j++){
				for(size_t k = 0; k < lhs.Cols(); k++){
					res(i, j) += lhs(i, k) * rhs(k, j);
				}
			}
		}
		return res;
	}
	if(mode == "ikj"){
		for(size_t i =0 ; i < lhs.Rows(); i++){
			for(size_t k = 0; k < lhs.Cols(); k++){
				double r = lhs(i, k);
				for(size_t j = 0; j < rhs.Cols(); j++){	
					res(i, j) += r * rhs(k, j);
				}
			}
		}
		return res;
	}
	if(mode == "kij"){
		for(size_t i =0 ; i < lhs.Rows(); i++){
			for(size_t k = 0; k < lhs.Cols(); k++){
				double r = lhs(i, k);
				for(size_t j = 0; j < rhs.Cols(); j++){	
					res(i, j) += r * rhs(k, j);
				}
			}
		}
		return res;
	}
	if(mode == "jki"){
		for(size_t j = 0; j < rhs.Cols(); j++){	
			for(size_t k = 0; k < lhs.Cols(); k++){
				double r = rhs(k, j);
				for(size_t i =0 ; i < lhs.Rows(); i++){
					res(i, j) += lhs(i, k) * r;
				}
			}
		}
		return res;
	}
	if(mode == "kji"){
		for(size_t j = 0; j < rhs.Cols(); j++){	
			for(size_t k = 0; k < lhs.Cols(); k++){
				double r = rhs(k, j);
				for(size_t i =0 ; i < lhs.Rows(); i++){
					res(i, j) += lhs(i, k) * r;
				}
			}
		}
		return res;
	}
	return res;
}

FloatMatrix MultiplyDouble(const DoubleMatrix& lhs, const DoubleMatrix& rhs, const std::string& mode){
	Matrix tmp(lhs.Rows(), rhs.Cols(), 'd', nullptr);
	DoubleMatrix res(tmp);
	std::cout << "MultiplyFloat init done" << std::endl;
	std::cout << res.Rows() << " " << res.Cols() << std::endl;
	iif(mode == "ijk"){
		for(size_t i =0 ; i < lhs.Rows(); i++){
			for(size_t j = 0; j < rhs.Cols(); j++){
				for(size_t k = 0; k < lhs.Cols(); k++){
					res(i, j) += lhs(i, k) * rhs(k, j);
				}
			}
		}
		return res;
	}
	if(mode == "jik"){
		for(size_t i =0 ; i < lhs.Rows(); i++){
			for(size_t j = 0; j < rhs.Cols(); j++){
				for(size_t k = 0; k < lhs.Cols(); k++){
					res(i, j) += lhs(i, k) * rhs(k, j);
				}
			}
		}
		return res;
	}
	if(mode == "ikj"){
		for(size_t i =0 ; i < lhs.Rows(); i++){
			for(size_t k = 0; k < lhs.Cols(); k++){
				double r = lhs(i, k);
				for(size_t j = 0; j < rhs.Cols(); j++){	
					res(i, j) += r * rhs(k, j);
				}
			}
		}
		return res;
	}
	if(mode == "kij"){
		for(size_t i =0 ; i < lhs.Rows(); i++){
			for(size_t k = 0; k < lhs.Cols(); k++){
				double r = lhs(i, k);
				for(size_t j = 0; j < rhs.Cols(); j++){	
					res(i, j) += r * rhs(k, j);
				}
			}
		}
		return res;
	}
	if(mode == "jki"){
		for(size_t j = 0; j < rhs.Cols(); j++){	
			for(size_t k = 0; k < lhs.Cols(); k++){
				double r = rhs(k, j);
				for(size_t i =0 ; i < lhs.Rows(); i++){
					res(i, j) += lhs(i, k) * r;
				}
			}
		}
		return res;
	}
	if(mode == "kji"){
		for(size_t j = 0; j < rhs.Cols(); j++){	
			for(size_t k = 0; k < lhs.Cols(); k++){
				double r = rhs(k, j);
				for(size_t i =0 ; i < lhs.Rows(); i++){
					res(i, j) += lhs(i, k) * r;
				}
			}
		}
		return res;
	}
	return res;
}


Matrix Multiply(const Matrix& lhs, const Matrix& rhs, const std::string& mode){
	char type = std::min(lhs._type, rhs._type);
	std::cout << mode << std::endl;
	if(type == 'f'){
		return MultiplyFloat( FloatMatrix(lhs), FloatMatrix(rhs), mode);
	}else{
		return MultiplyDouble( DoubleMatrix(lhs), DoubleMatrix(rhs), mode);
	}
}


#endif