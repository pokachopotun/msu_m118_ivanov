#include <iostream>
#include <fstream>
#include <random>
#include <string>
#include <algorithm>

#include <cstdlib>

using namespace std;

int main(int argc, char* argv[]){
	if(argc < 2){
		cout << "use ./genmap mapping.txt" << endl;
		return 0;
	}

	string outputFilename(argv[1]);

	ofstream ofs(outputFilename);

	vector< vector< int > > a;
	for(int i = 0;i < 8; i++){
		for(int j =0 ; j < 8; j++){
			for(int k = 0; k < 8; k++){
				std::vector<int> v;
				v.push_back(i);
				v.push_back(j);
				v.push_back(k);
				a.push_back(v);
			}
		}
	}
 
    std::random_device rd;
    std::mt19937 g(rd());
 
    std::shuffle(a.begin(), a.end(), g);
 
    for(int i = 0; i < a.size(); i++){
		for(int j = 0; j < 3; j++){
			ofs << a[i][j] << " ";
		}
		ofs << 0;
		ofs << endl;
	}
	return 0;
}
