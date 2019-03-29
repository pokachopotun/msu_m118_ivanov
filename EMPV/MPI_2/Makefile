all: main

define submit
	mpisubmit.pl -p ${1} --stdout=${1}.rozenbrooke.out ./main -- 1000 20 100000 onepoint rozenbrooke 10 0.2
	mpisubmit.pl -p ${1} --stdout=${1}.rastrigin.out ./main -- 1000 20 100000 onepoint rastrigin 10 0.2
	mpisubmit.pl -p ${1} --stdout=${1}.square.out ./main -- 1000 20 100000 onepoint square 10 0.2	
endef

submit: main
	$(call submit,1)
	$(call submit,2)
	$(call submit,4)
	$(call submit,8)

main: main.cpp
	mpixlC -O3 -std=c++11 main.cpp -o main

clean:
	rm -f main *.dat *.out *.err
