#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include <stdlib.h>
#include <unistd.h>

#define Prob 0.5
int main( int argc, char *argv[])
{
	int D[24][3] = {
		14,12,0,
		7,11,0,
		1,12,0,
		9,11,0,
		14,5,0,
		2,4,0,
		6,10,0,
		8,3,0,
		4,13,0,
		15,16,0,
		5,6,0,
		13,16,0,
		12,3,0,
		11,15,0,
		8,10,0,
		2,7,0,
		14,1,0,
		2,15,0,
		6,8,0,
		4,9,0,
		5,3,0,
		1,3,0,
		9,15,0,
		9,16,0
	};
	int myrank, size;
	MPI_Status status;
	MPI_Init( &argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);



	
	MPI_Finalize();
	return 0;
}
