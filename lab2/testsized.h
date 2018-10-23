#include <string>



class TestSized{
public:
	TestSized(const int n):N(n){
		_a = new float[N * N];
		_b = new float[N * N];
		_c = new float[N * N];
	}
	int Multiply(int bSize, const std::string& mode){
		if(mode == "ijk"){
			for(int i = 0;i < N; i+=bSize){
				for(int j=0 ; j < N; j+= bSize){
					for(int k = 0 ; k < N; k+=bSize){
						//block
						for(int x = i; x < std::min(i + bSize, N); x++){
							for(int y = j; y < std::min(j + bSize, N); y++){
								for(int z = k; z < std::min(k + bSize, N); z++){
									c(x, y) += a(x, z) * b(z, y);
								}
							}
						}
					}
				}
			}
			return 0;
		}
		if(mode == "ikj"){
			for(int i = 0;i < N; i+=bSize){
				for(int k = 0 ; k < N; k+=bSize){
					for(int j=0 ; j < N; j+= bSize){
						//block
						for(int x = i; x < std::min(i + bSize, N); x++){
							for(int z = k; z < std::min(k + bSize, N); z++){
								float r = a(x, z); 
								for(int y = j; y < std::min(j + bSize, N); y++){
									c(x, y) += r * b( z, y );
								}
							}
						}
					}
				}
			}
			return 0;
		}
		return 1;
	}

	float& a(const int x, const int y){
		return _a[x * N + y];
	}
	float& b(const int x, const int y){
		return _b[x * N + y];
	}

	float& c(const int x, const int y){
		return _c[x * N + y];
	}

	~TestSized(){
		delete[] _a;
		delete[] _b;
		delete[] _c;
	}

private:
	const int N;
	float *_a;
	float *_b;
	float *_c;
};