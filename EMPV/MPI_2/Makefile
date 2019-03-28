all: main

define test_size
		$(call submit_one,32,$(1),$(2))
		$(call submit_one,64,$(1),$(2))
		$(call submit_one,128,$(1),$(2))
		$(call submit_one,256,$(1),$(2))
		$(call submit_one,512,$(1),$(2))
endef

define submit
		mpisubmit.bg -n $(1) -w 00:03:00 --stdout=$(3).n.$(2).p.$(1).out ./main 0 10 5 0.5 $(2) 
endef

define gen_size
		mpisubmit.bg -n 1 ./gen $(1) $(2) $(1)x$(2).in
		mpisubmit.bg -n 1 ./gen $(2) 1 $(2)x1.in
endef

main: main.cpp
	mpic++ -O3 -std=c++11 main.cpp -o main -lmpi

clean:
	rm -f main *.dat *.out

series_2: main
	$(call submit,1,10000000,2)
	$(call submit,64,10000000,2)
	$(call submit,256,10000000,2)
	$(call submit,1024,10000000,2)

series_1: main
	$(call submit,1,1,1)
	$(call submit,1,10,1)
	$(call submit,1,100,1)
	$(call submit,1,1000,1)
	$(call submit,1,10000,1)
	$(call submit,1,100000,1)
	$(call submit,1,1000000,1)
	$(call submit,1,10000000,1)

	$(call submit,64,1,1)
	$(call submit,64,10,1)
	$(call submit,64,100,1)
	$(call submit,64,1000,1)
	$(call submit,64,10000,1)
	$(call submit,64,100000,1)
	$(call submit,64,1000000,1)
	$(call submit,64,10000000,1)

	$(call submit,256,1,1)
	$(call submit,256,10,1)
	$(call submit,256,100,1)
	$(call submit,256,1000,1)
	$(call submit,256,10000,1)
	$(call submit,256,100000,1)
	$(call submit,256,1000000,1)
	$(call submit,256,10000000,1)

	$(call submit,1024,1,1)
	$(call submit,1024,10,1)
	$(call submit,1024,100,1)
	$(call submit,1024,1000,1)
	$(call submit,1024,10000,1)
	$(call submit,1024,100000,1)
	$(call submit,1024,1000000,1)
	$(call submit,1024,10000000,1)

series_3: main
	$(call submit,1,1000,3)
	$(call submit,64,64000,3)
	$(call submit,256,256000,3)
	$(call submit,1024,1024000,3)
	
