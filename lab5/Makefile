all: main gen print compare mul

define test_size
		$(call submit_one,32,$(1),$(2))
		$(call submit_one,64,$(1),$(2))
		$(call submit_one,128,$(1),$(2))
		$(call submit_one,256,$(1),$(2))
		$(call submit_one,512,$(1),$(2))
endef

define submit_one
		mpisubmit.bg -n $(1) --stdout=$(1).$(2)x$(3).out ./main $(2)x$(3).in $(3)x1.in $(1)_$(2)x$(3).matrix.out
endef

define gen_size
		mpisubmit.bg -n 1 ./gen $(1) $(2) $(1)x$(2).in
		mpisubmit.bg -n 1 ./gen $(2) 1 $(2)x1.in
endef

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

genMap: genMap.cpp
	g++ -std=c++11 -O3 genMap.cpp -o genMap

clean:
	rm -f main print compare gen mul *.dat *.out *.in

generate: gen
	$(call gen_size,512,512)
	$(call gen_size,1024,1024)
	$(call gen_size,2048,2048)
	$(call gen_size,4096,4096)
	$(call gen_size,4096,1024)
	$(call gen_size,1024,4096)

submit: main
	$(call test_size,512,512)
	$(call test_size,1024,1024)
	$(call test_size,2048,2048)
	$(call test_size,4096,4096)
	$(call test_size,4096,1024)
	$(call test_size,1024,4096)
