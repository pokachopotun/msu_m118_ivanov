all: main print gen compare

runall:
	./gen
	./main
	./print
	./compare
	cat output.txt
main:
	g++ main.cpp -o main
print:
	g++ print.cpp -o print
gen:
	g++ -std=c++11 gen.cpp -o gen
compare:
	g++ compare.cpp -o compare
clean:
	rm main print compare gen
