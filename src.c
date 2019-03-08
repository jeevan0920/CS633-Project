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
int find_contract_edges( int size, int *Data, int node_size, int *Alltocounts, int *Total_Leaders );
void find_alltoallv_data( int size, int *Data, int node_size, int *Alltocounts, int *Total_Leaders, int *AlltoData);
void find_cummulative_displs( int size, int *Alltocounts, int *sdispls);
void setup_leader_contraction(int *Alltoall_count, int *Alltoall_data, int *Nodes );
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

		free(sendcnts);
		free(displs);
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
	 * Total_Leadres : This list contain all the leaders present round;
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
	int *Total_Leaders;
	Total_Leaders = (int *) malloc( sizeof(int) * V );

	int *recvcnts = (int *) malloc( sizeof(int) * size );
	for( i=0; i < size; i++)
        	recvcnts[i] = node_vert ;
        	recvcnts[size-1] = last_node_vert;
        int *displs = (int *) malloc ( sizeof(int) * size);
	for( i=1,displs[0] = 0; i<size; i++)
        	displs[i] = displs[i-1] + node_vert;

	if( myrank == size-1 ){
		Num_of_leaders = find_leaders( Leaders, Bcast_leaders, last_node_vert, Base);
		for(i=0; i<last_node_vert ; i++)
			printf("%d %d\n", Leaders[i], Bcast_leaders[i]);
		
		MPI_Allgatherv( Bcast_leaders, last_node_vert, MPI_INT, Total_Leaders, recvcnts, displs, MPI_INT, MPI_COMM_WORLD);  
	
		free(recvcnts);
		free(displs);

		for( i=0 ; i<V; i++)
                printf("%d ", Total_Leaders[i]);
	}
	else{
		Num_of_leaders = find_leaders( Leaders, Bcast_leaders, node_vert, Base);
		MPI_Allgatherv( Bcast_leaders, node_vert, MPI_INT, Total_Leaders, recvcnts, displs, MPI_INT, MPI_COMM_WORLD);

	}
	MPI_Barrier ( MPI_COMM_WORLD );
	printf("All should reach Barrier here\n");

	int *Alltocounts, *AlltoData, Contracted_edges;
	Alltocounts = (int *) malloc( sizeof(int) * size );
	for( i=0; i<size; i++)
	{
		Alltocounts[i] = 0;
	}
	if( myrank == size-1)
	{
		Contracted_edges = find_contract_edges(size, Data, last_node_size, Alltocounts, Total_Leaders);
		AlltoData = (int *) malloc( sizeof(int) * Contracted_edges * 2 );
		find_alltoallv_data(size, Data, last_node_size, Alltocounts, Total_Leaders, AlltoData);
		printf("After edges selection ");
	}
	else
	{
		Contracted_edges = find_contract_edges(size, Data, node_size, Alltocounts, Total_Leaders);
		AlltoData = (int *) malloc( sizeof(int) * Contracted_edges * 2 );
		find_alltoallv_data(size, Data, node_size, Alltocounts, Total_Leaders, AlltoData);

	}
	
	int j;
	for( j=0; j< size; j++ ){
		MPI_Barrier(MPI_COMM_WORLD);
		if( myrank == j )
		{
			for( i=0; i<size; i++)
                        	printf("%d ",Alltocounts[i]);

		printf("\n");	
		}
	}
	int *Alltoall_count;
	Alltoall_count = (int *) malloc( sizeof(int) * size );
	MPI_Alltoall( Alltocounts, 1 , MPI_INT, Alltoall_count, 1 , MPI_INT, MPI_COMM_WORLD);

        for( j=0; j< size; j++ ){
                MPI_Barrier(MPI_COMM_WORLD);
                if( myrank == j )
                {
                        for( i=0; i<size; i++)
                                printf("%d ",Alltoall_count[i]);

                printf("\n");
                }
        }
	
	int *sdispls, *rdispls, *Alltoall_data, count = 0;
	sdispls = (int *) malloc( sizeof(int) * size );
	rdispls = (int *) malloc( sizeof(int) * size );
	for( i=0; i<size; i++)
	{
		count += Alltoall_count[i];
	}
	Alltoall_data = (int *) malloc( sizeof(int) * count * 2 );
	find_cummulative_displs( size, Alltocounts, sdispls);
	find_cummulative_displs( size, Alltoall_count, rdispls);
	for( j=0; j< size; j++ ){
                MPI_Barrier(MPI_COMM_WORLD);
                if( myrank == j )
                {
                        for( i=0; i<size; i++)
                                printf("%d %d  ", sdispls[i], rdispls[i] );

                printf("\n");
                }
        }
	MPI_Barrier(MPI_COMM_WORLD);
	for(i=0; i<size; i++)
	{
		Alltocounts[i] *= 2;
		Alltoall_count[i] *= 2;
	}
	for( j=0; j< size; j++ ){
                MPI_Barrier(MPI_COMM_WORLD);
                if( myrank == j )
                {
		 printf(" Contracted edges : %d ", Contracted_edges );

			 for(i=0; i < Contracted_edges * 2; i += 2)
                         printf("%d %d  ", AlltoData[i], AlltoData[i+1] );

                printf("\n");
                }
        }

	MPI_Alltoallv( AlltoData, Alltocounts, sdispls, MPI_INT, Alltoall_data, Alltoall_count, rdispls, MPI_INT, MPI_COMM_WORLD);

	for( j=0; j< size; j++ ){
                MPI_Barrier(MPI_COMM_WORLD);
                if( myrank == j )
                {
			int k = 0;
			for( i=0; i<size; i++)
				k += Alltoall_count[i];
                 printf(" Recieved Edges  : %d ", k/2 );

                         for(i=0; i < k; i += 2)
                         printf("%d %d  ", Alltoall_data[i], Alltoall_data[i+1] );

                printf("\n");
                }
        }

	setup_leader_contraction( Alltoall_count, Alltoall_data, Nodes );	
	MPI_Finalize();
	return 0;
}
void setup_leader_contraction(int *Alltoall_count, int *Alltoall_data, int *Nodes ){

}
/*
 * This function will find the send displacements , Recevie displacements
 * to MPI_Alltoallv() function 
 */
void find_cummulative_displs( int size, int *Alltocounts, int *sdispls)
{
	int i, count = 0;
	for( i=0; i < size; i++ )
	{
		 count += Alltocounts[i];
                 sdispls[i]= (count - Alltocounts[i]) * 2;
                 //printf("%d  ", cum_count[i] );
	}
	if( Alltocounts[ size-1 ] == 0 )
	{
		sdispls[ size-1 ] = sdispls [ size -2 ] ;
	}
}
/*
 * find_leaders function will find the leaders the the vertices that that node hold
 * I am taking every node and genearating random probability if that probablity is
 * greater than our desired probablity then it will consider that node has LEADER
 */
int find_leaders( int *Leaders, int *Bcast_leaders, int node_vert, int Base)
{
	int i, count = 0, rand_value;
	srand( getpid() );
	for( i = 0; i< node_vert; i++){
		rand_value = (int)rand()  % 100;
		float rand_prob = (float) rand_value/100;
		if( rand_prob >= Prob && Leaders[i] ){
			Bcast_leaders[i] = Base + i;
			count += 1;
		}
		else
			Bcast_leaders[i] = 0;
	}
	return count;
}
/*
 * find_contract_edges function will find the how many edges are contracted  for each
 * node.
 * By this we can allocate memory for sending this data to all other nodes so that
 * they will update nodes's parent pointer( Means that vertex has been contracted )
 */
int find_contract_edges(int size, int *Data, int node_size, int *Alltocounts , int *Total_Leaders)
{
	int i, node, count = 0;
	for( i=0 ; i<node_size * ts; i += ts)
	{
		if(!Data[i+2])
		{
			if( Total_Leaders[ Data[i] ] && Total_Leaders[ Data[i+1] ] )
			{
				// if both vertices are leaders nothing to do
			}
			else if ( !Total_Leaders[ Data[i] ] && !Total_Leaders[ Data[i+1] ] )
			{
				// if both are not leaders nothing to do
			}
			else
			{
				count += 1;
				if( Total_Leaders[ Data[i] ] )
				{
					node = Data[i + 1] / size;
					if (node > size-1 )
						node = size - 1; 
					Alltocounts[node] += 1;
				}
				else
				{
					node = Data[i] / size;
					if (node > size-1)
						node = size - 1;
                                        Alltocounts[node] += 1;
				}
			}
		}
	}
	return count;
}
/*
 * This function will set up data format to send all other node about 
 * edge contracted information 
 * AlltoData : It is the data to be sent. here only that memory been allocated.
 */
void find_alltoallv_data( int size, int *Data, int node_size, int *Alltocounts, int *Total_Leaders, int *AlltoData)
{
	int i, node, count=0, local_count[size],cum_count[size], index ;
	for( i=0,cum_count[0] = 0; i<size;  i++)
	{
		 count += Alltocounts[i];
		 local_count[i] = 0;
		 cum_count[i]= (count - Alltocounts[i]) * 2;
		 //printf("%d  ", cum_count[i] );
	}

	for( i=0 ; i<node_size * ts; i += ts)
        {
                if(!Data[i+2])
                {
                        if( Total_Leaders[ Data[i] ] && Total_Leaders[ Data[i+1] ] )
                        {
                                // if both vertices are leaders nothing to do
                        }
                        else if ( !Total_Leaders[ Data[i] ] && !Total_Leaders[ Data[i+1] ] )
                        {
                                // if both are not leaders nothing to do
                        }
                        else
                        {
                                if( Total_Leaders[ Data[i] ] )
                                {
                                        node = Data[i + 1] / size;
                                        if (node > size-1 )
                                                node = size - 1;
					index = cum_count[node] + local_count[node] * 2 ;
                                        AlltoData[index] = Data[i+1];
                                        AlltoData[index +1] = Data[i];

                                        local_count[node] += 1;
					Data[i+2] = 1;

                                }
                                else
                                {
                                        node = Data[i] / size;
					if (node > size-1)
						node = size - 1;
					index = cum_count[node] + local_count[node] * 2 ;
                                        AlltoData[index] = Data[i];
                                        AlltoData[index +1] = Data[i + 1];
                                        
					local_count[node] += 1;
					Data[i+2] = 1;
                                }
                        }
                }
        }
	/*	
 	for( i=0; i<count * 2; i += 2)
	{
		printf("%d %d\n", AlltoData[i], AlltoData[i+1] );
	}
	*/
}
