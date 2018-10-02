all: main print gen compare fromtxt timeit

define report_size
    ./timeit f $(1) $(1) ijk > $(1).ijk.time
	./timeit f $(1) $(1) jik > $(1).jik.time
	./timeit f $(1) $(1) kij > $(1).kij.time
	./timeit f $(1) $(1) ikj > $(1).ikj.time
	./timeit f $(1) $(1) jki > $(1).jki.time
	./timeit f $(1) $(1) kji > $(1).kji.time
	cat $(1).*.time > $(1)x$(1).time.dat
	gnuplot -e "set terminal svg size 600,500; set output '$(1).plot.svg'; plot '$(1)x$(1).time.dat' w lp t '$(1)x$(1)'"
endef

define test_size
	#test $(1)x$(1)
	./main A_$(1)x$(1).test B_$(1)x$(1).test $(1)x$(1).test.out
	./compare C0_$(1)x$(1).test $(1)x$(1).test.out
endef

define gen_test_size
	#gen test $(1)x$(1)
	./gen f $(1) $(1) A_$(1)x$(1).test
	./gen f $(1) $(1) B_$(1)x$(1).test
	#A
	./print A_$(1)x$(1).test
	#B
	./print B_$(1)x$(1).test
endef


timeit:
	g++ -I./ -O3 -std=c++11 matrix.h timeit.cpp -o timeit
print:
	g++ -I./ -std=c++11 matrix.h print.cpp -o print
gen:
	g++ -I./ -std=c++11 matrix.h gen.cpp -o gen
fromtxt:
	g++ -I./ -std=c++11 matrix.h fromtxt.cpp -o fromtxt
compare:
	g++ -I./ -std=c++11 matrix.h compare.cpp -o compare
main:
	g++ -I./ -std=c++11 matrix.h main.cpp -o main
test:
	$(call test_size,3)
	$(call test_size,4)
	$(call test_size,5)

gen_test:
	$(call gen_test_size,4)	
	$(call gen_test_size,5)

report:
	$(call report_size,50)
	$(call report_size,200)
	$(call report_size,300)
	$(call report_size,500)
	$(call report_size,1000)
clean:
	rm main print compare gen timeit *.time *.dat *.out