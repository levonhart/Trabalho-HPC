#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <omp.h>
#include <math.h>

#include "ivp.h"

void runge_kutta(const double t,const ivp_function_ptr y[], unsigned count[],
		unsigned size, double c[], int nconst[], int nump) {
	int ntasks, rank;
	double * constants;
	int * nconstants;
	double * buffer;
	double * dv;
	double * x0, * x1;
	size_t bufsiz, constants_size;

	MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	printf("Process: %d\tt=%lf\n", rank, t);

	if (rank == 0) {
		int dconst[nump], displs[nump], csize[nump];
		MPI_Request request[nump];

		nconstants = malloc(count[0]*sizeof(int));
		memcpy(nconstants, nconst, count[0]*sizeof(int));
		csize[0] = 0;
		for (size_t i = 0; i < count[0]; i++) {
			csize[0] += nconstants[i];
		}
		constants_size = csize[0];

		constants = malloc(constants_size*sizeof(double));
		memcpy(constants, c, constants_size*sizeof(double));

		dconst[0] = 0;
		for (int i = 1; i < nump; i++) {
			dconst[i] += count[i-1];
			MPI_Send(nconst+dconst[i], count[i], MPI_INT, i, 1, MPI_COMM_WORLD);
		}
		displs[0] = 0;
		for (int i = 1; i < nump; i++) {
			csize[i] = 0;
			for (size_t j = 0; j < count[i]; j++) {
				csize[i] += nconst[ dconst[i] + j ];
			}
			displs[i] = displs[i-1] + csize[i-1];
			MPI_Isend(c+displs[i], csize[i], MPI_DOUBLE, i, 2,
					MPI_COMM_WORLD, &request[i]);
		}
	} else if (rank < nump) {
		MPI_Request request;
		MPI_Status status;
		nconstants = malloc(count[rank]*sizeof(int));
		MPI_Recv(nconstants, count[rank], MPI_INT, 0, 1,
				MPI_COMM_WORLD, &status);
		constants_size = 0;
		for (size_t i = 0; i < count[rank]; i++) {
			constants_size += nconstants[i];
		}
		constants = malloc(constants_size*sizeof(double));
		MPI_Irecv(constants, constants_size, MPI_DOUBLE, 0, 2,
				MPI_COMM_WORLD, &request);
		MPI_Wait(&request, &status);
	}
	else {
		fprintf(stderr, "Processo %d ocioso.\n", rank);
		return;
	}
	printf("Process: %d\t"
			"c[] =", rank);
	for (size_t i = 0; i < constants_size; i++) {
		printf("    %0.3lf%s", constants[i],
				i+1 == constants_size ? "\n" : "");
	}

	dv = malloc(count[rank]*sizeof(double));
	x0 = malloc(count[rank]*sizeof(double));
	x1 = malloc(count[rank]*sizeof(double));

#pragma omp parallel for
	for (int i = 0; i < count[rank]; i++) {
		
		dv[i] = +INFINITY;
		while (dv[i] > MIN_DV) {
			dv[i] = 0.;
		}
	}

	free(constants);
	free(nconstants);
	free(dv);
	free(x0);
	free(x1);
}
