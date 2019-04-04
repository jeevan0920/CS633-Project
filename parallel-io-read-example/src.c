/*
* Sample code for performing paralle io reading
*/


#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>


int main(int argc , char *argv[]){
	
	int myrank,size,nints,filesize = atoi(argv[1]);
	int nprocs;
	int* buf;

	
	MPI_Status status;
	MPI_Init(&argc,&argv);

	MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);

	nprocs = size;

	MPI_File fh;
	MPI_Offset offset;

	MPI_File_open(MPI_COMM_WORLD,"./graph_dataset",MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);
	nints = filesize / (nprocs*sizeof(int));

	buf = (int*)malloc(sizeof(int) * nints);


	offset = myrank*nints*sizeof(int);
	printf("## process %d offset is %d \n",myrank,offset);
	MPI_File_read_at(fh,offset,buf,nints,MPI_INT, &status);

	for(int i=0;i<nints;i++)
		printf("process %d read %d \n",myrank,buf[i]);

	
	MPI_File_close(&fh);
	MPI_Finalize();
	return 0;
		
	
}
