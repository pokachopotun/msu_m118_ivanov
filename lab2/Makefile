all: timeit flops totins L1_cache L2_cache


timeit: timeit.cpp testsized.h
	g++ -I./ -O3 -std=c++11 testsized.h timeit.cpp -lpapi -o timeit

flops: flops.cpp testsized.h
	g++ -I./ -O3 -std=c++11 testsized.h flops.cpp -lpapi -o flops

totins: totins.cpp testsized.h
	g++ -I./ -O3 -std=c++11 testsized.h totins.cpp -lpapi -o totins

L1_cache: L1_cache.cpp testsized.h
	g++ -I./ -O3 -std=c++11 testsized.h L1_cache.cpp -lpapi -o L1_cache

L2_cache: L2_cache.cpp testsized.h
	g++ -I./ -O3 -std=c++11 testsized.h L2_cache.cpp -lpapi -o L2_cache





define submit_one
	mpisubmit.pl --stdout=$(1).$(2).$(3).out --stderr=$(1).$(2).$(3).err L2_cache $(1) $(2) $(3)
endef

define submit_bSize
	$(call submit_one,1000,$(1),ijk)
	$(call submit_one,1000,$(1),ikj)

	$(call submit_one,2000,$(1),ijk)
	$(call submit_one,2000,$(1),ikj)
	
	$(call submit_one,3000,$(1),ijk)
	$(call submit_one,3000,$(1),ikj)
	
	$(call submit_one,4000,$(1),ijk)
	$(call submit_one,4000,$(1),ikj)

	$(call submit_one,5000,$(1),ijk)
	$(call submit_one,5000,$(1),ikj)
endef

submit: flops
	$(call submit_bSize,836)
	$(call submit_bSize,32)

clean:
	rm -f timeit totins flops L1_cache L2_cache *.o *.gch *.time *.dat *.out *.svg

# define report_bsize_mode
# 	./timeit $(1) $(2) ijk > $(1).$(2).ijk.time
# 	cat $(1).$(2).*.time > $.time.dat
# 	gnuplot -e "set terminal svg size 600,500; set output '$(1).plot.svg'; plot '$(1)x$(1).time.dat' w lp t '$(1)x$(1)'"
# endef

# define report_
# 	./timeit $(1) $(2) ijk > $(1).$(2).ijk.time
# 	./timeit $(1) $(2) ikj > $(1).$(2).jkj.time
# 	cat $(1).$(2).*.time > $(1)x$(1).time.dat
# 	gnuplot -e "set terminal svg size 600,500; set output '$(1).plot.svg'; plot '$(1)x$(1).time.dat' w lp t '$(1)x$(1)'"
# endef
