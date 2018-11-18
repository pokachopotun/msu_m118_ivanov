all: main

define submit_one
	mpisubmit.pl -n $(1) -w 00:01 --stdout=$(2).$(3).$(1).out --stderr=$(1).$(2).$(3).err main $(2) $(3) $(4)
endef

define submit_n
	$(call submit_one,$(1),1000,2000,$(1)_1e3.txt)
	$(call submit_one,$(1),10000,20000,$(1)_1e4.txt)
	$(call submit_one,$(1),100000,200000,$(1)_1e5.txt)
	$(call submit_one,$(1),1000000,2000000,$(1)_1e6.txt)
	$(call submit_one,$(1),10000000,20000000,$(1)_1e7.txt)
	$(call submit_one,$(1),100000000,200000000,$(1)_1e8.txt)
endef

run: main
	$(call submit_n,16)
	$(call submit_n,32)
	$(call submit_n,64)
	$(call submit_n,128)
	$(call submit_n,256)
	$(call submit_n,512)
	$(call submit_n,1024)

main: main.cpp
	g++ -std=c++11 main.cpp -o main -lpthread
