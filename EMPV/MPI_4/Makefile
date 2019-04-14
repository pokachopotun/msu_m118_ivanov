all: main

define submit
	mpisubmit.bg -n ${1} --stdout=${1}.1e1.out ./main -- 10 40 null glyder
	mpisubmit.bg -n ${1} --stdout=${1}.1e2.out ./main -- 100 400 null glyder
	mpisubmit.bg -n ${1} --stdout=${1}.1e3.out ./main -- 1000 4000 null glyder
	mpisubmit.bg -n ${1} --stdout=${1}.0.5e3.out ./main -- 500 2000 null glyder
endef

submit: main
	$(call submit,16)
	$(call submit,64)
	$(call submit,256)
	$(call submit,1024)

define one_glyder
	mpirun -n 4 ./local 5 ${1} ${1}.out glyder
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
