all: main

define submit_one
	mpisubmit.bg -n ${1} --stdout=${1}.${2}.out ./main 20000 20 ${2}
endef

define submit_size
	$(call submit_one,{1},1)
	$(call submit_one,{1},2)
	$(call submit_one,{1},3)
endef

submit: main
	$(call submit_size,2)
	$(call submit_size,4)
	$(call submit_size,8)
	$(call submit_size,16)

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
	mpixlcxx_r -O2 main.cpp -o main

local: main.cpp
	mpic++ -O3 main.cpp -o local -lmpi

clean:
	rm -f main *.dat *.out *.err *.tar
