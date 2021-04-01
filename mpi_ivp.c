#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ivp.h"

#define N 2


ivp_function f01, f02, f03;

const ivp_function_ptr inst[4] = {f01, f02, f03, f02};

int main(int argc, char *argv[]){

	double t = nan("1");

	if (argc > 1) {
		t = atoi(argv[1]);
	}

	if (isnan(t)) {
		return EXIT_FAILURE;
	}

	unsigned nump = N;
	int ntasks, rank;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == 0) {
		double c[] = { -2., 1., 1.};
		unsigned count[] = {2, 2};
		unsigned nconst[] = {0, 1, 1, 1};
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}
// P1
double f01 (double t, double * state, size_t n, double c[]) {
	return sin(t) * state[0];
}
double f02 (double t, double * state, size_t n, double c[]) {
	return c[0] * t * state[1];
}

//P2
double f03 (double t, double * state, size_t n, double c[]) {
	return c[0] * exp(t);
}
/* double f04 (double t, double * state, size_t n, double param[]) { */
/*     return t * state[1]; */
/* } */
