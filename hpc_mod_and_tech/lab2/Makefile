polus: main.cpp
	mpixlC -O2 main.cpp -o polus -fopenmp -lpthread

bluegene: main.cpp
	mpixlC -O2 main.cpp -o bluegene -fopenmp -lpthread

local: main.cpp
	mpic++ -O3 main.cpp -o local -fopenmp -lpthread

test: local
	mpirun -n 4 ./local 2 2 1000 1000 2 0.001 0.001 0.0000001 
