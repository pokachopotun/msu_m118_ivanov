polus: main.cpp
	mpixlC -O2 main.cpp -o polus 

local: main.cpp
	mpic++ -O2 main.cpp -o local 

clean:
	rm -f main print compare gen fromtxt timeit *.o *.gch *.time *.dat *.out *.svg

define submit_mpi_job
	./local $(1) $(2) $(3) $(4) 0.00000001 1000 0 $(5)
endef

define submit_mpi_job_polus
	mpisubmit.pl -p $(4) --stdout mpi.$(6).$(1).$(1).$(1).$(4).out --stderr mpi.$(6).$(1).$(1).$(1).$(4).err ./local -- $(1) $(2) $(3) $(7) $(8) $(9) 0.00000001 1000 0 $(5)
endef

mode_a: polus
	$(call submit_mpi_job_polus, 128, 128, 128, 1, 10, mode_a,1,1,1)
	$(call submit_mpi_job_polus, 128, 128, 128, 2, 10, mode_a,1,1,2)
	$(call submit_mpi_job_polus, 128, 128, 128, 4, 10, mode_a,1,2,2)
	$(call submit_mpi_job_polus, 128, 128, 128, 8, 10, mode_a,2,2,2)


mode_b: polus
	$(call submit_mpi_job_polus, 64, 64, 64, 1, 10, mode_b,1,1,1)
	$(call submit_mpi_job_polus, 64, 64, 64, 2, 10, mode_b,1,1,2)
	$(call submit_mpi_job_polus, 64, 64, 64, 4, 10, mode_b,1,2,2)
	$(call submit_mpi_job_polus, 64, 64, 64, 8, 10, mode_b,2,2,2)

	$(call submit_mpi_job_polus, 128, 128, 128, 1, 10, mode_b,1,1,1)
	$(call submit_mpi_job_polus, 128, 128, 128, 2, 10, mode_b,1,1,2)
	$(call submit_mpi_job_polus, 128, 128, 128, 4, 10, mode_b,1,2,2)
	$(call submit_mpi_job_polus, 128, 128, 128, 8, 10, mode_b,2,2,2)


	$(call submit_mpi_job_polus, 256, 256, 256, 1, 10, mode_b,1,1,1)
	$(call submit_mpi_job_polus, 256, 256, 256, 2, 10, mode_b,1,1,2)
	$(call submit_mpi_job_polus, 256, 256, 256, 4, 10, mode_b,1,2,2)
	$(call submit_mpi_job_polus, 256, 256, 256, 8, 10, mode_b,2,2,2)

mode_c: polus
	#10^4
	$(call submit_mpi_job_polus, 32, 32, 32, 1, 10, mode_a,1,1,1)
	$(call submit_mpi_job_polus, 32, 32, 64, 2, 10, mode_a,1,1,2)
	$(call submit_mpi_job_polus, 32, 64, 64, 4, 10, mode_a,1,2,2)
	$(call submit_mpi_job_polus, 64, 64, 64, 8, 10, mode_a,2,2,2)	

	#10^5
	$(call submit_mpi_job_polus, 64, 64, 64, 1, 10, mode_a,1,1,1)
	$(call submit_mpi_job_polus, 64, 64, 128, 2, 10, mode_a,1,1,2)
	$(call submit_mpi_job_polus, 64, 128, 128, 4, 10, mode_a,1,2,2)
	$(call submit_mpi_job_polus, 128, 128, 128, 8, 10, mode_a,2,2,2)

	#10^6
	$(call submit_mpi_job_polus, 128, 128, 128, 1, 10, mode_c,1,1,1)
	$(call submit_mpi_job_polus, 128, 128, 256, 2, 10, mode_c,1,1,2)
	$(call submit_mpi_job_polus, 128, 256, 256, 4, 10, mode_c,1,2,2)
	$(call submit_mpi_job_polus, 256, 256, 256, 8, 10, mode_c,2,2,2)


