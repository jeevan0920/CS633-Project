#include <stdio.h>
#include <string.h> 
#include "mpi.h"
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <stdio.h>

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
	FILE *op_file = fopen("output.txt","w"); //output file which contains the unque color assigned to each vertex



	V = num_verts;
	ne = num_edges;


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
	int *Vertex_Color;
	if( myrank == 0)
	{
		Vertex_Color = (int *) malloc(sizeof(int) * V );
	}
	int Num_of_leaders,j;
	int *Total_Leaders;
	Total_Leaders = (int *) malloc( sizeof(int) * V );
	int sample = 0;
	// Optimize Allocation & Deallocation of buffers
	int *Alltocounts, *Alltoall_count, *sdispls, *rdispls;
	Alltocounts = (int *) malloc( sizeof(int) * size );
	Alltoall_count = (int *) malloc( sizeof(int) * size );
	sdispls = (int *) malloc( sizeof(int) * size );
        rdispls = (int *) malloc( sizeof(int) * size );
	int *recvcnts = (int *) malloc( sizeof(int) * size );
	int *displs = (int *) malloc ( sizeof(int) * size);
	int *Total_Clist = (int *) malloc( sizeof(int) * V);

while ( edges_Next_iter )	
{
	last_node_size = ne/size + ne % size, node_size = ne/size, last_node_vert = V/size + V % size, node_vert = V/size;
	Broadcast_Leaders(myrank, size, node_vert, last_node_vert, Nodes, Bcast_leaders, Total_Leaders, Base);
	

	/*
	 * Alltocounts : This will count the edges to be sent to each node
	 * AlltoData : This will store the edges to be sent to each node 
	 * Contracted_edges : This will tells us to how many edges has been contracted
	 */
	int  *AlltoData, Contracted_edges;
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
	 * Alltoall_count : This will collect number of edges from every node ( Means the vertex belongs to This node
	 * has been contracted some where else ) 
	 * MPI_Alltoall : This will collect all the counts they are sending in MPI_Alltoallv
	 */
	MPI_Alltoall( Alltocounts, 1 , MPI_INT, Alltoall_count, 1 , MPI_INT, MPI_COMM_WORLD);
	
	 
	/*
	 * sdispls : sending offsets to other nodes
	 * rdispls : Recveing offset form other nodes
	 */
	int *Alltoall_data, count = 0;
	for( i=0; i<size; i++)
	{
		count += Alltoall_count[i];
	}
	Alltoall_data = (int *) malloc( sizeof(int) * count * 2 );
	find_cummulative_displs( size, Alltocounts, sdispls);
	find_cummulative_displs( size, Alltoall_count, rdispls);

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
	 * MPI_Alltoallv will get all the data that belongs to respective nodes. so that it will mark as a contracted vertex.
	 */
	MPI_Request Recv_request[i], Send_request[i];
        MPI_Status status;
        for(i=0; i< size; i++)
        {
                MPI_Isend(AlltoData + sdispls[i] , Alltocounts[i], MPI_INT, i, myrank, MPI_COMM_WORLD, Send_request+i);
        }
        for(i=0; i<size; i++)
        {
                MPI_Irecv(Alltoall_data + rdispls[i] , Alltoall_count[i], MPI_INT, i, i, MPI_COMM_WORLD, Recv_request + i);
        }
        for(i=0; i<size;i++)
        {
                MPI_Wait(Send_request + i, &status);
        }
        for(i=0; i<size; i++)
        {
                MPI_Wait(Recv_request + i, &status);
        }

	//MPI_Alltoallv( AlltoData, Alltocounts, sdispls, MPI_INT, Alltoall_data, Alltoall_count, rdispls, MPI_INT, MPI_COMM_WORLD);

	setup_leader_contraction( Alltoall_count, Alltoall_data, Nodes , Base, size);	
	
	/*
	 * After setup_leader_contraction checking what are edges are being contracted.
	 */
	
        last_node_size = ne/size + ne % size, node_size = ne/size, last_node_vert = V/size + V % size, node_vert = V/size;
	for( i=0; i < size; i++)
                recvcnts[i] = node_vert ;
                recvcnts[size-1] = last_node_vert;
        for( i=1,displs[0] = 0; i<size; i++)
                displs[i] = displs[i-1] + node_vert;
	if( myrank == size-1 )
	{
		node_vert = last_node_vert;
		node_size = last_node_size;
	}
	
	for(i=0; i< node_vert; i++)
	{
		Total_Clist[i+Base] = Nodes[i];
	}
	MPI_Allgatherv( MPI_IN_PLACE, node_vert, MPI_INT, Total_Clist, recvcnts, displs, MPI_INT, MPI_COMM_WORLD);
	
	
	for(j=0; j< size; j++)
	{
		if(myrank == j)
		{
			edges_Left = relink_edges( node_size, Total_Clist, Data );
		}
	}
	
	MPI_Allreduce( &edges_Left, &edges_Next_iter, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	
	 if(edges_Left == 0)
        {
                if(myrank == 0)
                {
                                for( i=0; i< V; i++ )
                                {
                                        //printf("%d-%d: ",i, Total_Clist[i] );
                                        j = i;
                                        while( Total_Clist[j] != -1 )
                                        {
                                                j= Total_Clist[j];
                                        }
                                        Vertex_Color[i] = j;
                                }
                }
        }


	free(AlltoData);
	free(Alltoall_data);
	sample ++;
	if(myrank == 0)
	{
		printf(" Round %d is Completed \n",sample);
	}
}
		
	if(myrank == 0)
	{
		for(i=0; i<V; i++)
		{
			fprintf(op_file,"%d\t%d\n",i, Vertex_Color[i]);
		}
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
		
		for(i=0 ; i < last_node_vert ; i++)
		{	
			Total_Leaders[Base + i] = Bcast_leaders[i];
		}	
                MPI_Allgatherv( MPI_IN_PLACE, last_node_vert, MPI_INT, Total_Leaders, recvcnts, displs, MPI_INT, MPI_COMM_WORLD);

		/*
                printf( "\n==================================================================\n");
                printf("All Leaders \n");
                for( i=0 ; i<V; i++)
                {
                        printf("%d ", Total_Leaders[i]);
                }
                printf( "\n==================================================================\n After contraction Edges list : \n");
		*/
        }

        else{
		  for(i=0 ; i < last_node_vert ; i++)
                 {
                         Total_Leaders[Base + i] = Bcast_leaders[i];
                 }

                Num_of_leaders = find_leaders( Nodes, Bcast_leaders, node_vert, Base);
                MPI_Allgatherv( MPI_IN_PLACE, node_vert, MPI_INT, Total_Leaders, recvcnts, displs, MPI_INT, MPI_COMM_WORLD);

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
			//printf("%d %d\n", Data[i], Data[i+1] );
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
