#include <mpi.h>

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <ctime>

using namespace std;

double frand(double a, double b)
{
	return a+(b-a)*(rand()/double(RAND_MAX));
}

long do_walk(long a, long b, long x, double p, double& t)
{
	static const long maxSteps = 1e6;
	long step = 0;
	while( x>a && x<b && step < maxSteps )
	{
		if( frand(0,1)<p )
			x += 1;
		else
			x -= 1;
		t += 1.0;
		step += 1;
	}
	return x;
}

int main(int argc, char** argv)
{

	MPI_Init(&argc, &argv);

	long a = atoi(argv[1]);
	long b = atoi(argv[2]);
	long x = atoi(argv[3]);
	double p = atof(argv[4]);
	long nGlobal = atoi(argv[5]);


	long size, rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);



	MPI_Barrier(MPI_COMM_WORLD);
	double time = MPI_Wtime();
	// presume that division is done evenly
	long nLocal = nGlobal / size;

	double totalTime = 0.0;
	long bCount = 0;

	for( long i=0; i< nLocal; i++ )
	{
		long out = do_walk(a, b, x, p, totalTime);
		if( out == b ) {
			bCount += 1;
		}
	}

	if( rank == 0){
		MPI_Reduce( MPI_IN_PLACE, &totalTime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
		MPI_Reduce( MPI_IN_PLACE, &bCount, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD );
	}
	else{
		MPI_Reduce( &totalTime, 0, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
		MPI_Reduce( &bCount, 0, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD );
	}

	time = MPI_Wtime() - time;
	
	if( rank == 0 ){
		printf("a %d b %d x %d p %6.2f nExperiments %d nProc %d bCount %d bProbability %6.2f avgTime %6.2f runTime %f\n", 
			a, b, x, p, nGlobal, size, bCount, 1.0*bCount/nGlobal, totalTime / nGlobal, time);
	}

	MPI_Finalize();

	return 0;
}