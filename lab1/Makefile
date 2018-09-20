all: main print gen compare


define test_once
    ./gen f 3 3 $(1)
	./gen f 3 3 $(2)
	./main $(1) $(2) $(3)
	./print $(1)
	./print $(2)
	./print $(3)
endef

runall:
	./gen d 11 5 generated.dat
	./main 
	./print generated.dat
	./compare

main:
	g++ -I./ -std=c++11 matrix.h main.cpp -o main
print:
	g++ -I./ -std=c++11 matrix.h print.cpp -o print
gen:
	g++ -I./ -std=c++11 matrix.h gen.cpp -o gen
compare:
	g++ compare.cpp -o compare
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
