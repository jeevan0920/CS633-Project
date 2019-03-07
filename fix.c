#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include <stdlib.h>
#include <unistd.h>

#define Prob 0.5
#define ne 24
#define ts 3
int main( int argc, char *argv[])
{
	int D[ne][ts] = {
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

	int myrank, size, i;
	MPI_Status status;
	MPI_Init( &argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	int last_node_size = ne/size + ne % size, node_size = ne/size;
	int *Data;
	if ( myrank == size-1 ){
		Data = (int *) malloc( sizeof(int) * last_node_size *ts );
	}
	else{
		Data = (int *) malloc( sizeof(int) * node_size * ts );
	}
	if(myrank == 0){
		
		int *sendcnts = (int *) malloc( sizeof(int) * size );
		for( i=0; i<size; i++)
			sendcnts[i] = node_size * ts ;
			sendcnts[size-1] = last_node_size * ts ;
		for ( i=0 ; i<size; i++)
			printf("%d ", sendcnts[i]);
		int *displs = (int *) malloc ( sizeof(int) * size);
		for( i=1,displs[0] = 0; i<size; i++)
			displs[i] = displs[i-1] + node_size * ts;
		for ( i=0 ; i<size; i++)
                        printf("%d ", displs[i]);

		MPI_Scatterv( D, sendcnts, displs, MPI_INT, Data, node_size * ts, MPI_INT, myrank , MPI_COMM_WORLD ); 
	}
	else{
		if( myrank == size-1 )
			MPI_Scatterv( NULL, NULL, NULL, MPI_INT, Data, last_node_size * ts , MPI_INT, 0 , MPI_COMM_WORLD );
		else
			MPI_Scatterv( NULL, NULL, NULL, MPI_INT, Data, node_size * ts , MPI_INT, 0 , MPI_COMM_WORLD );
	}
	if( myrank == 0 ){
		for( i=0; i<node_size * ts ; i += ts) 
			printf("%d %d %d\n", Data[i], Data[i+1], Data[i+2] );

	}
	  if( myrank == 1 ){
                for( i=0; i<node_size  * ts; i += ts)
                        printf("%d %d %d\n", Data[i], Data[i+1], Data[i+2] );

        }
  	if( myrank == 2 ){
                for( i=0; i<node_size * ts; i += ts)
                        printf("%d %d %d\n", Data[i], Data[i+1], Data[i+2] );

        }
  	if( myrank == 3 ){
                for( i=0; i<node_size * ts; i += ts)
                        printf("%d %d %d\n", Data[i], Data[i+1], Data[i+2] );

        }

	MPI_Finalize();
	return 0;
}
