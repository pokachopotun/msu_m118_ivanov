all: main gen print compare mul

main: main.cpp
	mpic++ -O3 main.cpp -o main

gen: gen.cpp
	g++ -std=c++11 -O3 gen.cpp -o gen

print: print.cpp
	g++ -std=c++11 -O3 print.cpp -o print

compare: compare.cpp
	g++ -std=c++11 -O3 compare.cpp -o compare

mul: mul.cpp
	g++ -std=c++11 -O3 mul.cpp -o mul

clean:
	rm -f main print compare gen mul *.dat *.out *.in