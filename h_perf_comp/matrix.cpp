#include <matrix.h>
#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <chrono>
#include <limits>

size_t Matrix::BufferSize() const{
	return MatrixSize();
}

// void Matrix::ReadInTxt(const std::string& fileName){
// 	using namespace std;
// 	ifstream txt(fileName);
// 	for(size_t i = 0 ; i < Rows(); i++){
// 		for(size_t j = 0; j < Cols(); j++){
// 			size_t pos = i * Cols() + j;
// 			txt >> _buf[pos];
// 		}	
// 	}
// }

// void Matrix::MakeUno(){
// 	for(size_t i = 0 ; i < Rows(); i++){
// 		for(size_t j = 0; j < Cols(); j++){
// 			//TODO correct calculation
// 			size_t pos = i * Cols() + j;
// 			if(i == j){
// 				_buf[pos] = 1.0f;
// 			}else{
// 				_buf[pos] = 0.0f;
// 			}
// 		}	
// 	}
// }

size_t Matrix::MatrixSize() const{
	return Cols() * Rows();
}

size_t Matrix::Cols() const{
	return _m;
}
	
size_t Matrix::Rows() const{
	return _n;
}

//TODO
// void Matrix::WriteToFile(const std::string& outFileName){
// 	std::ofstream ofs(outFileName, std::ios::binary);
// 	ofs.write(&_type, sizeof(char));
// 	ofs.write(reinterpret_cast< char* >( &_n ), sizeof(size_t)/sizeof(char));
// 	ofs.write(reinterpret_cast< char* >( &_m ), sizeof(size_t)/sizeof(char));
// 	ofs.write(_charBuf, BufferSize() );
// }
	
// size_t Matrix::GetPos(size_t row, size_t col) const {
// 	return _row_i[row];
// }

// void Matrix::clean(){
// 	delete[] _buf;
// }

// std::ostream& operator<<(std::ostream& ofs, const Matrix& m){
// 	for(size_t i=0 ; i < m.Rows(); i++){
// 		for(size_t j = 0; j < m.Cols(); j++){
// 			size_t pos = i * m.Cols() + j;
// 			ofs << _buf[pos] << " ";
// 		}	
// 		ofs << std::endl;
// 	}
// 	return ofs;
// }

// Matrix Multiply(const Matrix& lhs, const Matrix& rhs){
// 	//TODO
	
// }


// Matrix Multiply(const Matrix& lhs, const Matrix& rhs, const std::string& mode){
// 	//TODO
// }

// bool Compare(const Matrix& lhs, const Matrix& rhs){
// 	//TODO
// }