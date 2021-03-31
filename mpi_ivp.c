#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 3

typedef double (ivp_function) (double, double*); // function(t, state)
typedef double (* ivp_function_ptr) (double, double*); // *function(t, state)
ivp_function f01, f02, f03;

int main(int argc, char *argv[]){

	double t = nan("1");

	if (argc > 1) {
		t = atoi(argv[1]);
	}

	if (isnan(t)) {
		return EXIT_FAILURE;
	}

	int ntasks, rank;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);



	MPI_Finalize();

	return EXIT_SUCCESS;
}

double f01 (double t, double * state) {

}
double f02 (double t, double * state) {
}
double f03 (double t, double * state) {
}
