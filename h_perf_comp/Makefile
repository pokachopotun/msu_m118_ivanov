all: main

main: main.cpp
	g++ -I./ -O3 -std=c++11 -fopenmp -o main main.cpp

clean:
	rm -f main print compare gen fromtxt timeit *.o *.gch *.time *.dat *.out *.svg
