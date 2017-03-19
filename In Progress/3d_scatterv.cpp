#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "mpi.h"

#define MASTER 0

int allocate3D(double ****arr3D, int l,int m,int n);
void deallocate3D(double*** arr3D,int l,int m);

int main(int argc, char **argv){
	double*** global;
	 double*** local;
	//Size and processor count values.
	int size_x, size_y, size_z;
	int proc_x, proc_y, proc_z;
	int procsize_x, procsize_y, procsize_z, procsize_total, size_total;

	int rank, size;

	if(argc<=6){
		printf("Usage: ./o size-x size-y size-z proc-x proc-y proc-z\n");
		return 1;
	}

	size_x = atoi(argv[1]);
	size_y = atoi(argv[2]);
	size_z = atoi(argv[3]);

	proc_x = atoi(argv[4]);
	proc_y = atoi(argv[5]);
	proc_z = atoi(argv[6]);

	procsize_x = size_x / proc_x;
	procsize_y = size_y / proc_y;
	procsize_z = size_z / proc_z;
	procsize_total = procsize_x*procsize_y*procsize_z;
	size_total = size_x*size_y*size_z;


	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);


	if (size != (proc_x*proc_y*proc_z)) {
		printf("Enter correct dimensions.\n");
		MPI_Abort(MPI_COMM_WORLD,1);
	}

	if (rank == MASTER){
		allocate3D(&global, size_x,size_y,size_z);

		for (int i=0; i<size_x; i++) {
			for (int j=0; j<size_y; j++)
				for (int k=0; k<size_z; k++)
					global[i][j][k] = k*size_x*size_y+j*size_x+i;
		}		

		for (int i=0; i<size_x; i++) {
			for (int j=0; j<size_y; j++){
				for (int k=0; k<size_z; k++)
					printf("%3.1f  ", global[i][j][k] );
				printf("\n");
			}
			printf("\n");
		}
	}

	// create the local array which we'll process
	allocate3D(&global, procsize_x, procsize_y, procsize_z);

	// create a datatype to describe the subarrays of the global array
	int sizes[3] = {size_x, size_y, size_z};
	int subsizes[3] = {procsize_x, procsize_y, procsize_z};
	int starts[3] = {0, 0, 0};

	MPI_Datatype type, subarrtype;
	MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &type);
	MPI_Type_create_resized(type, 0, procsize_x*sizeof(double), &subarrtype);
	MPI_Type_commit(&subarrtype);

	double *globalptr = NULL;
	if (rank == MASTER) 
		globalptr = &global[0][0][0];

	//scatter the array

	int sendcounts[size_total];
	int displs[size_total];

	if (rank == MASTER){
		for (int i=0; i<size_total; i++)
			sendcounts[i] = 1;
		int disp = 0;
		for (int i=0; i<size_x; i++) {
			for (int j=0; j<size_y; j++){
				for (int k=0; k<size_z; k++){
					displs[k*size_x*size_y+j*size_x+i] = disp;
					disp += 1;
				}
			}
		}		
	}

	MPI_Scatterv(globalptr, sendcounts, displs, subarrtype, &(local[0][0][0]),
				 size_total/procsize_total, MPI_DOUBLE,
				 0, MPI_COMM_WORLD);

	MPI_Finalize();
	return 0;
}

int allocate3D(double ****arr3D, int l,int m,int n){
	
	int i,j,k;

	*arr3D = (double***)malloc(l * sizeof(double **));

	for(i=0;i<l;i++){
		(*arr3D)[i] = (double**)malloc(m * sizeof(double*));

		for(j=0;j<m;j++){
			(*arr3D)[i][j] = (double*)malloc(n*sizeof(double));
		}
	}

	return 0;
}

void deallocate3D(double*** arr3D,int l,int m)
{
	int i,j;

	for(i=0;i<l;i++)
	{
		for(int j=0;j<m;j++)
		{
			free(arr3D[i][j]);
		}
		free(arr3D[i]);
	}
	free(arr3D);
}