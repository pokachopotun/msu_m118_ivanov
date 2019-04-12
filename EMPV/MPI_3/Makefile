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

define one_glyder
	mpirun -n 10 ./local 10 ${1} ${1}.out glyder
	python visualize.py ${1}.out 10
endef

glyder: local
	$(call one_glyder,0)
	$(call one_glyder,1)
	$(call one_glyder,2)
	$(call one_glyder,3)
	$(call one_glyder,4)
	$(call one_glyder,5)
	$(call one_glyder,6)
	$(call one_glyder,7)
	$(call one_glyder,8)
	$(call one_glyder,9)
	$(call one_glyder,10)
	$(call one_glyder,11)	
	
	
	
	
	
	

main: main.cpp
	mpixlC -O3 main.cpp -o main

local: main.cpp
	mpic++ -O3 main.cpp -o local -lmpi

clean:
	rm -f main *.dat *.out *.err *.tar
