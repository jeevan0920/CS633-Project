#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include <stdlib.h>
#include <unistd.h>

#define Prob 0.5
#define V 17
#define ne 24
#define ts 3
int find_leaders( int *leaders, int *Bcast_leaders, int node_vert , int Base); 
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
	
	int last_node_size = ne/size + ne % size, node_size = ne/size, last_node_vert = V/size + V % size, node_vert = V/size;
	/*
	 * last_node_size Number of edges at (n-1) th node
	 * Data is a int pointer that will hold all the edges list for
	 * the particular node.
	 * MPI_Scatterv will done at 0th node process that means.. 0th
	 * node will read all the data from the file and scatterv to all
	 * th nodes.
	 */
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
/*	
	// print all the edges at every node  
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
*/
	/*
	 * Nodes  : This will tells that parent of the node with respect to 
	 * 		this node
	 * Leaders : This will contain leaders list till now
	 * Bcast_leaders : This will be broadcasted to all other nodes.
	 * Num_of_leaders : Num_of_leaders is selected num of leaders;
	 */
	int *Nodes, Base = myrank * node_vert, *Leaders, *Bcast_leaders;
        if ( myrank == size-1 ){
                Nodes = (int *) malloc( sizeof(int) * last_node_vert );
                Leaders = (int *) malloc( sizeof(int) * last_node_vert );
                Bcast_leaders = (int *) malloc( sizeof(int) * last_node_vert );
        	for( i=0; i<last_node_vert; i++){
			Leaders[i] = -1;
			Nodes[i] = -1;	
		}
	}
        else{
                Nodes = (int *) malloc( sizeof(int) * node_vert );	
                Leaders = (int *) malloc( sizeof(int) * node_vert );
              	Bcast_leaders = (int *) malloc( sizeof(int) * node_vert );
		for( i=0; i<node_vert; i++){
			Nodes[i] = -1;
			Leaders[i] = -1;
		}
        }
	
	int Num_of_leaders;
	if( myrank == size-1 ){
		Num_of_leaders = find_leaders( Leaders, Bcast_leaders, last_node_vert, Base);
		for(i=0; i<last_node_vert ; i++)
			printf("%d ", Leaders[i]);
		for( i=0; i<Num_of_leaders ; i++)
			printf(" %d ", Bcast_leaders[i]);
	}
	else{
		Num_of_leaders = find_leaders( Leaders, Bcast_leaders, node_vert, Base);
	}
	if(myrank == 1){
		 for(i=0; i<node_vert ; i++)
                        printf("%d ", Leaders[i]);
                for( i=0; i<Num_of_leaders ; i++)
                        printf(" %d ", Bcast_leaders[i]);
	}
	MPI_Finalize();
	return 0;
}
/*
 * find_leaders function will find the leaders the the vertices that that node hold
 * I am taking every node and genearating random probability if that probablity is
 * greater than our desired probablity then it will consider that node has LEADER
 */
int find_leaders( int *Leaders, int *Bcast_leaders, int node_vert, int Base){
	int i, count = 0, rand_value;
	srand( getpid() );
	for( i = 0; i< node_vert; i++){
		rand_value = (int)rand()  % 100;
		float rand_prob = (float) rand_value/100;
		printf("\n %f ", (float)rand_value/100 );
		if( rand_prob >= Prob && Leaders[i] ){
			Leaders[i] = 0;
			Bcast_leaders[count++] = Base + i;
		}
	}
	return count;
}
