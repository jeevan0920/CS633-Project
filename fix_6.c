#include <stdio.h>
#include <string.h> 
#include "mpi.h"
#include <stdlib.h>
#include <unistd.h>
#include <time.h>

#define Prob 0.5
#define ts 3

//global variables 
int V; //to hold number of vertices
int ne; //to hold number of edges



int CURRENT_ITERATE = 1;
int find_leaders( int *leaders, int *Bcast_leaders, int node_vert , int Base);
int find_contract_edges( int size, int node_vert, int *Data, int node_size, int *Alltocounts, int *Total_Leaders );
void find_alltoallv_data( int size, int node_vert, int *Data, int node_size, int *Alltocounts, int *Total_Leaders, int *AlltoData);
void find_cummulative_displs( int size, int *Alltocounts, int *sdispls);
void setup_leader_contraction(int *Alltoall_count, int *Alltoall_data, int *Nodes, int Base , int size);
unsigned int relink_edges( int node_size, int *Total_Clist, int *Data);
void Broadcast_Leaders(int myrank, int size, int node_vert, int last_node_vert, int *Leaders, int *Bcast_Leaders, int *Total_Leaders, int Base);

int main( int argc, char *argv[])
{

	if(argc != 4){
		printf("Invalid arguments\n");
		printf("Usage : executable path_to_dataset  num_edges num_vertices\n");
		return -1;
	}	

	int num_verts = atoi(argv[3]); //number of vertices
	int num_edges = atoi(argv[2]); //number of edges
	char *file_path = argv[1]; //path to input data set
	int filesize = num_edges*3*sizeof(int); //file size in bytes
	int num_ints ;//number of integers to be read by process


	V = num_verts;
	ne = num_edges;

	/*	
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
		9,16,0,
		7,9,0,
	};
	*/



	int myrank, size, i;
	MPI_Status status;
	MPI_Init( &argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	//============================================= parallel io reading ============================================

	//for handling parallel IO read
	MPI_File fh;
        MPI_Offset offset;
	
	MPI_File_open(MPI_COMM_WORLD,file_path,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);
	
	//calculating number of integers to be read by each process
	num_ints = 3*( num_edges / size);
	//assign remaining integers to last process
	if(myrank == size-1){
        	num_ints  = num_ints + 3*(num_edges % size);
	}
	
	

	int last_node_size = ne/size + ne % size, node_size = ne/size, last_node_vert = V/size + V % size, node_vert = V/size;
	/*
	 * last_node_size Number of edges at (n-1) th node
	 * Data is a int pointer that will hold all the edges list for
	 * the particular node.
	 * MPI_Scatterv will be done at 0th node process that means.. 0th
	 * node will read all the data from the file and scatterv to all
	 * the nodes.
	 */
	int *Data;
	if ( myrank == size-1 ){
		Data = (int *) malloc( sizeof(int) * last_node_size *ts );
	}
	else{
		Data = (int *) malloc( sizeof(int) * node_size * ts );
	}

	/*
	if(myrank == 0){
		
		int *sendcnts = (int *) malloc( sizeof(int) * size );
		for( i=0; i<size; i++)
			sendcnts[i] = node_size * ts ;
			sendcnts[size-1] = last_node_size * ts ;
		int *displs = (int *) malloc ( sizeof(int) * size);
		for( i=1,displs[0] = 0; i<size; i++)
			displs[i] = displs[i-1] + node_size * ts;
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
	*/

	//calculating offset to read from the file for current process
	offset = (num_edges/size)*3*sizeof(int)*myrank;

        //printf("## process %d offset is %lld \n",myrank,offset);
	//printf("## process %d number of integers assigned is %d\n",myrank,num_ints);
	
	//parallel reading from  file from differnt offsets from each process
        MPI_File_read_at(fh,offset,Data,num_ints,MPI_INT, &status);
	
	//printing the data read from each process
	//for(int i=0;i<num_ints;i += 3)
        //        printf("process %d read %d %d %d \n",myrank,Data[i],Data[i+1],Data[i+2]);


	//closing the file 
	MPI_File_close(&fh);

	//=============================   end of parallel io reading   ===========================================================	
	
	

		

	// print all the edges at every node
	
	if( myrank == 0 ){
		for( i=0; i<node_size * ts ; i += ts) 
			printf("%d %d\n", Data[i], Data[i+1]);	
			//printf("%d %d %d\n", Data[i], Data[i+1], Data[i+2] );
	}
	  if( myrank == 1 ){
                for( i=0; i<node_size  * ts; i += ts)
			printf("%d %d\n", Data[i], Data[i+1]);                        
			//printf("%d %d %d\n", Data[i], Data[i+1], Data[i+2] );

        }
  	if( myrank == 2 ){
                for( i=0; i<last_node_size * ts; i += ts)
			printf("%d %d\n", Data[i], Data[i+1]);
                        //printf("%d %d %d\n", Data[i], Data[i+1], Data[i+2] );

        }
	/*
	 * Nodes  : Contains the leader for each vertex
	 * Leaders : This will contain leaders list till now
	 * Bcast_leaders : This will be broadcasted to all other nodes.
	 * Num_of_leaders : Num_of_leaders is selected num of leaders;
	 * Total_Leadres : This list contain all the leaders present round;
	 */
	int *Nodes, Base = myrank * node_vert, *Leaders, *Bcast_leaders;
        unsigned int edges_Left, edges_Next_iter = ne;
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
	int sample = 4;
while ( edges_Next_iter )	
{
	last_node_size = ne/size + ne % size, node_size = ne/size, last_node_vert = V/size + V % size, node_vert = V/size;
	Broadcast_Leaders(myrank, size, node_vert, last_node_vert, Nodes, Bcast_leaders, Total_Leaders, Base);
	

	MPI_Barrier ( MPI_COMM_WORLD );
	//printf("All should reach Barrier here\n");

	/*
	 * Alltocounts : This will count the edges to be sent to each node
	 * AlltoData : This will store the edges to be sent to each node 
	 * Contracted_edges : This will tells us to how many edges has been contracted
	 */
	int *Alltocounts, *AlltoData, Contracted_edges;
	Alltocounts = (int *) malloc( sizeof(int) * size );
	for( i=0; i<size; i++)
	{
		Alltocounts[i] = 0;
	}
	if( myrank == size-1)
	{	
		node_vert = last_node_vert;
		node_size = last_node_size;

	}
	Contracted_edges = find_contract_edges(size, node_vert, Data, node_size, Alltocounts, Total_Leaders);
	AlltoData = (int *) malloc( sizeof(int) * Contracted_edges * 2 );
	find_alltoallv_data(size, node_vert, Data, node_size, Alltocounts, Total_Leaders, AlltoData);
	
	/* 
	 * This will print Alltocounts of every node means How many edges has been contracted and whom to sent.
	 */
	int j;
	
	for( j=0; j< size; j++ ){
		MPI_Barrier(MPI_COMM_WORLD);
		if( myrank == j )
		{
			if( myrank == 0)
	                {
        	                printf("\nNumber of possible edges contracting to every  node\n at process %d :", myrank);
                	}
			for( i=0; i<size; i++)
                        	printf("%d ",Alltocounts[i]);
		if( myrank < size-1 )
		printf("\n at process %d :", myrank + 1);	
		}
	}
	
	/*
	 * Alltoall_count : This will collect number of edges from every node ( Means the vertex belongs to This node
	 * has been contracted some where else ) 
	 * MPI_Alltoall : This will collect all the counts they are sending in MPI_Alltoallv
	 */
	int *Alltoall_count;
	Alltoall_count = (int *) malloc( sizeof(int) * size );
	MPI_Alltoall( Alltocounts, 1 , MPI_INT, Alltoall_count, 1 , MPI_INT, MPI_COMM_WORLD);
	
	 /*
         * This will print Alltoall_count  of every node (Means how many edges each node is sent to this node)
         */
	
        for( j=0; j< size; j++ ){
                MPI_Barrier(MPI_COMM_WORLD);
                if( myrank == j )
                {
			 if( myrank == 0)
	                {
        	                printf("\n Number of possible contracting edges recieved \n at process %d :", myrank);
                	}
                        for( i=0; i<size; i++)
                                printf("%d ",Alltoall_count[i]);
		if( myrank < size-1 )
                printf("\n at process %d :", myrank + 1);
                }
        }
	
	/*
	 * sdispls : sending offsets to other nodes
	 * rdispls : Recveing offset form other nodes
	 */
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

	/*
	 * Printing sdispls, rdispls here to check
	 */
	for( j=0; j< size; j++ ){
                MPI_Barrier(MPI_COMM_WORLD);
                if( myrank == j )
                {
                        for( i=0; i<size; i++);
                              //  printf("%d %d  ", sdispls[i], rdispls[i] );
               // printf("\n");
                }
        }
	MPI_Barrier(MPI_COMM_WORLD);

	/*
	 * Alltocounts , Alltoall_count should be interms of how many MPI_INTEGERS
	 * each edge will contain two Vertexes so. multiplied by 2. 
	 */
	for(i=0; i<size; i++)
	{
		Alltocounts[i] *= 2;
		Alltoall_count[i] *= 2;
	}

	/*
	 * printing All the contracted edges 
	 */
	
	for( j=0; j< size; j++ ){
                MPI_Barrier(MPI_COMM_WORLD);
                if( myrank == j )
                {
			 printf(" Contracted edges : %d\t ", Contracted_edges );
			 for(i=0; i < Contracted_edges * 2; i += 2)
                         printf("%d %d  ", AlltoData[i], AlltoData[i+1] );

                printf("\n");
                }
        }
	
	/*
	 * MPI_Alltoallv will get all the data that belongs to respective nodes. so that it will mark as a contracted vertex.
	 */
	MPI_Alltoallv( AlltoData, Alltocounts, sdispls, MPI_INT, Alltoall_data, Alltoall_count, rdispls, MPI_INT, MPI_COMM_WORLD);

	setup_leader_contraction( Alltoall_count, Alltoall_data, Nodes , Base, size);	
	
	/*
	 * After setup_leader_contraction checking what are edges are being contracted.
	 */
	
  	for( j=0; j< size; j++ ){
                MPI_Barrier(MPI_COMM_WORLD);
                if( myrank == j )
                {
                        int k = 0;
                        for( i=0; i<size; i++)
                                k += Alltoall_count[i];
                 printf("\n Recieved Edges  : %d \t", k/2 );

                         for(i=0; i < k; i += 2)
                         printf("%d %d  ", Alltoall_data[i], Alltoall_data[i+1] );
		if( myrank == size -1 )
			node_vert = last_node_vert;
                printf("\n Parent Vertex: %d\t", node_vert );
			for(i=0; i < node_vert; i++)
			{
				printf("%d %d  ", Base + i, Nodes[i] );
			}
		printf("\n");
                }
        }
	
	MPI_Barrier(MPI_COMM_WORLD);
	int *recvcnts = (int *) malloc( sizeof(int) * size );
        last_node_size = ne/size + ne % size, node_size = ne/size, last_node_vert = V/size + V % size, node_vert = V/size;
	for( i=0; i < size; i++)
                recvcnts[i] = node_vert ;
                recvcnts[size-1] = last_node_vert;
        int *displs = (int *) malloc ( sizeof(int) * size);
        for( i=1,displs[0] = 0; i<size; i++)
                displs[i] = displs[i-1] + node_vert;
	int *Total_Clist = (int *) malloc( sizeof(int) * V);
	if( myrank == size-1 )
	{
		node_vert = last_node_vert;
		node_size = last_node_size;
	}
	
	for( j=0; j< size; j++)
        {
                MPI_Barrier( MPI_COMM_WORLD );
                if( myrank == j)
                {	
			printf(" Vertex contracted information \t ");
                        for( i=0; i< node_vert ; i++ )
                        {
                                printf("%d ", Nodes[i] );
                        }
                        printf("\n");
                }

        }
	
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Allgatherv( Nodes, node_vert, MPI_INT, Total_Clist, recvcnts, displs, MPI_INT, MPI_COMM_WORLD);
	
	for( j=size-1; j>=0; j--)
	{
		MPI_Barrier( MPI_COMM_WORLD );
		if( myrank == j)
                {
			if( myrank == size-1)
	                {
                	        printf("Gatherv Contracting vertex information :\n");
        	        }
			for( i=0; i< V; i++ )
			{
				printf("%d-%d: ",i, Total_Clist[i] );
			}
			printf("\n");
		}

	}
	
	//relink_edges( node_size, Total_Clist, Data );
 	for( j=0; j< size; j++ ){
                MPI_Barrier(MPI_COMM_WORLD);
                if(myrank == j)
                {
                        for( i=0; i<node_size * ts; i += ts)
                        {
                             // printf("%d %d %d \n", Data[i], Data[i+1], Data[i+2] );
                        }
                }
        }
	MPI_Barrier( MPI_COMM_WORLD);
	for(j=0; j< size; j++)
	{
		MPI_Barrier(MPI_COMM_WORLD);
		if(myrank == j)
		{
			edges_Left = relink_edges( node_size, Total_Clist, Data );
		}
	}
	for( j=0; j< size; j++ ){
                MPI_Barrier(MPI_COMM_WORLD);
		if(myrank == j)
		{
			for( i=0; i<node_size * ts; i += ts)
			{
	//			printf("%d %d %d \n", Data[i], Data[i+1], Data[i+2] );
			}
		}
        }
	for(j=0; j< size; j++)
        {
                MPI_Barrier(MPI_COMM_WORLD);
                if(myrank == j)
                {
                	printf("%d \n", edges_Left);
		}
        }

	MPI_Allreduce( &edges_Left, &edges_Next_iter, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	free(Alltocounts);
	free(AlltoData);
	free(Alltoall_count);
	free(sdispls);
	free(rdispls);
	free(Alltoall_data);
	free(recvcnts);
	free(displs);
	free(Total_Clist);
	sample --;
}
	MPI_Finalize();
	return 0;
}
void Broadcast_Leaders(int myrank, int size, int node_vert, int last_node_vert, int *Nodes, int *Bcast_leaders, int *Total_Leaders, int Base)
{
	  /*
         * recvcnts : It is recv_counts while doing Allgatherv
         * displs : where we put the data
         * MPI_Allgatherv : It is for broadcastint all the leaders it have to 
         * all other nodes.
         */

	int Num_of_leaders, i=0;
	int *recvcnts = (int *) malloc( sizeof(int) * size );
        for( i=0; i < size; i++)
                recvcnts[i] = node_vert ;
                recvcnts[size-1] = last_node_vert;
        int *displs = (int *) malloc ( sizeof(int) * size);
        for( i=1,displs[0] = 0; i<size; i++)
                displs[i] = displs[i-1] + node_vert;
	
	 if( myrank == size-1 ){
                Num_of_leaders = find_leaders( Nodes, Bcast_leaders, last_node_vert, Base);

                /* for(i=0; i<last_node_vert ; i++)
                {
                        printf("%d %d\n", Nodes[i], Bcast_leaders[i]);
                } */
                MPI_Allgatherv( Bcast_leaders, last_node_vert, MPI_INT, Total_Leaders, recvcnts, displs, MPI_INT, MPI_COMM_WORLD);


                printf( "\n==================================================================\n");
                printf("All Leaders \n");
                for( i=0 ; i<V; i++)
                {
                        printf("%d ", Total_Leaders[i]);
                }
                printf( "\n==================================================================\n After contraction Edges list : \n");
        }

        else{
                Num_of_leaders = find_leaders( Nodes, Bcast_leaders, node_vert, Base);
                MPI_Allgatherv( Bcast_leaders, node_vert, MPI_INT, Total_Leaders, recvcnts, displs, MPI_INT, MPI_COMM_WORLD);

        }
	free(recvcnts);
	free(displs);
}
/*
 * re allocate edges
 */
unsigned int relink_edges( int node_size, int *Total_Clist, int *Data)
{
	int i, node;
       	unsigned int count=0;
        for( i=0 ; i<node_size * ts; i += ts)
        {
                if(!Data[ i+2 ])
                {
			//printf("Before \t %d-->%d  %d--> %d  %d \n",Data[i], Total_Clist[ Data[i] ], Data[i+1], Total_Clist[ Data[i+1] ], Data[i+2] );
                        if( Total_Clist[ Data[i] ] == -1  &&  Total_Clist[ Data[i+1] ] == -1 )
                        {
                                // if both vertices are not yet been contracted.
                        }
                        else if ( Total_Clist[ Data[i] ] != -1 && Total_Clist[ Data[i+1] ] != -1 )
                        {
                                // if both vertice are non-leaders
				Data[ i ] = Total_Clist[ Data[i] ];
				Data[ i+1 ] = Total_Clist[ Data[i+1] ] ;
				 if( Data[ i ] == Data [ i+ 1 ] )
                                {
                                        Data[ i+2 ] = 1;
                                }
                        }
                        else
                        {	
				// one leader one non-leader
				
				if( Total_Clist[ Data[i] ] == -1 )
				{
					Data[ i+1 ] = Total_Clist[ Data[i+1] ] ;
				}
				else
				{
					Data [ i ] = Total_Clist[ Data[i] ];
				}
				if( Data[ i ] == Data [ i+ 1 ] )
				{
					Data[ i+2 ] = 1;
				}
                        }
                }
		if( !Data[ i + 2 ])
		{	count += 1;
			printf("%d %d\n", Data[i], Data[i+1] );
		}
        }
	return count;
}

/*
 * This function will contract the vertex based on Alltoall_data( Recived by all
 * other nodes )
 * It will mark a parent_vertex in Node list;
 */
void setup_leader_contraction(int *Alltoall_count, int *Alltoall_data, int *Nodes, int Base, int size )
{
	 int k = 0, index, i ;
         for( i=0; i<size; i++)
                 k += Alltoall_count[i];
         for(i=0; i < k; i += 2){
		index = Alltoall_data[i] - Base;
		Nodes[ index ] = Alltoall_data[ i + 1 ] ;
 	 }
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
int find_leaders( int *Nodes, int *Bcast_leaders, int node_vert, int Base)
{
	int i, count = 0, rand_value;
	srand( getpid() + CURRENT_ITERATE );
	for( i = 0; i< node_vert; i++){
		rand_value = (int)rand()  % 100;
		float rand_prob = (float) rand_value/100;
		if( rand_prob >= Prob && Nodes[i] == -1 ){
			Bcast_leaders[i] = Base + i;
			count += 1;
		}
		else
			Bcast_leaders[i] = 0;
	}
	CURRENT_ITERATE += 1;
	return count;
}
/*
 * find_contract_edges function will find the how many edges are contracted  for each
 * node.
 * By this we can allocate memory for sending this data to all other nodes so that
 * they will update nodes's parent pointer( Means that vertex has been contracted )
 */
int find_contract_edges(int size,int node_vert, int *Data, int node_size, int *Alltocounts , int *Total_Leaders)
{
	int i, node, count = 0;
	node_vert = V/size;
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
					node = Data[i + 1] / node_vert;
					if (node > size-1 )
						node = size - 1; 
					Alltocounts[node] += 1;
				}
				else
				{
					node = Data[i] / node_vert;
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
void find_alltoallv_data( int size, int node_vert, int *Data, int node_size, int *Alltocounts, int *Total_Leaders, int *AlltoData)
{
	int i, node, count=0, local_count[size],cum_count[size], index ;
	node_vert = V/size;
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
                                        node = Data[i + 1] / node_vert;
                                        if (node > size-1 )
                                                node = size - 1;
					index = cum_count[node] + local_count[node] * 2 ;
                                        AlltoData[index] = Data[i+1];
                                        AlltoData[index +1] = Data[i];

                                        local_count[node] += 1;
			//		Data[i+2] = 1;

                                }
                                else
                                {
                                        node = Data[i] / node_vert;
					if (node > size-1)
						node = size - 1;
					index = cum_count[node] + local_count[node] * 2 ;
                                        AlltoData[index] = Data[i];
                                        AlltoData[index +1] = Data[i + 1];
                                        
					local_count[node] += 1;
			//		Data[i+2] = 1;
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
