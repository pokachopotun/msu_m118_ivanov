all: main print gen compare


define timeit
 #    ./gen f 3 3 $(1)
	# ./gen f 3 3 $(2)
	./timeit
	# ./print $(1)
	# ./print $(2)
	# ./print $(3)
endef

runall:
	./gen d 11 5 generated.dat
	./main 
	./print generated.dat
	./compare

timeit:
	g++ -I./ -std=c++11 matrix.h timeit.cpp -o timeit
print:
	g++ -I./ -std=c++11 matrix.h print.cpp -o print
gen:
	g++ -I./ -std=c++11 matrix.h gen.cpp -o gen
compare:
	g++ -I./ -std=c++11 matrix.h compare.cpp -o compare
test:
	$(call test_once,E.in,F.in,G.out)
	# ./gen 4 3 3 A.in
	# ./gen f 3 3 B.in
	# ./main A.in B.in C.out
	# ./print A.in
	# ./print B.in
	# ./print C.out
clean:
	rm main print compare gen
