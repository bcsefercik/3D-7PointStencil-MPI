#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define SIZE_X  12
#define SIZE_Y  8

int main() {

    MPI_Init(NULL, NULL);
    int comm_sz, my_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    char i;

    int array[SIZE_X*SIZE_Y];
    const int yDim=2;
    const int xDim=3;
    const int sub_XSize = SIZE_Y/yDim;  /* number of rows in _block_ */
    const int sub_YSize = SIZE_X/xDim; /* number of cols in _block_ */

    if (my_rank == 0) {
        for (int i=0; i<SIZE_X*SIZE_Y; i++) {
            array[i] = i;
        }
    }


    int sub_array[sub_XSize*sub_YSize];

    int receive_up[sub_XSize];
    int receive_down[sub_XSize];
    int receive_left[sub_YSize];
    int receive_right[sub_YSize];
    int send_left[sub_YSize];
    int send_right[sub_YSize];


    MPI_Datatype blocktype;
    MPI_Datatype blocktype2;

    MPI_Type_vector(sub_XSize, sub_YSize, SIZE_X, MPI_INT, &blocktype2);
    MPI_Type_create_resized(blocktype2, 0, sizeof(int), &blocktype);
    MPI_Type_commit(&blocktype);

    int offset[yDim*xDim];
    int counts[yDim*xDim];
    for (int i=0; i<yDim; i++) {
        for (int j=0; j<xDim; j++) {
            offset[i*xDim+j] = i*SIZE_X*sub_XSize+j*sub_YSize;
            counts[i*xDim+j] = 1;
        }
    }

    MPI_Scatterv(array, counts, offset, blocktype, sub_array, sub_XSize*sub_YSize, MPI_INT, 0, MPI_COMM_WORLD);

    for (int rank=0; rank<comm_sz; rank++) {
        if (rank == my_rank) {
            printf("Rank = %d\n", my_rank);
            if (my_rank == 0) {
                printf("Global matrix: \n");
                for (int i=0; i<SIZE_Y; i++) {
                    for (int j=0; j<SIZE_X; j++) {
                        printf("%3d ",array[i*SIZE_X+j]);
                    }
                    printf("\n");
                }
            }
            printf("Local Matrix:\n");
            for (int i=0; i<sub_XSize; i++) {
                for (int j=0; j<sub_YSize; j++) {
                    printf("%3d ",sub_array[i*sub_YSize+j]);
                }
                printf("\n");
            }
            printf("\n");
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Request reqs[8];
    //communication - further will be combined with execution to improve performance

    for (int rank=0; rank<comm_sz; rank++) {
        if (rank == my_rank) {
            printf("Rank = %d\n", my_rank);
            if (my_rank <xDim) {
                MPI_Isend(&sub_array[(sub_XSize-1)*sub_YSize], sub_XSize, MPI_INT, my_rank+xDim, my_rank, MPI_COMM_WORLD, &reqs[0]);
                MPI_Irecv(&receive_up[0], sub_XSize, MPI_INT, my_rank+xDim, my_rank+xDim, MPI_COMM_WORLD, &reqs[1]);
            }

            if (my_rank >=xDim) {
                MPI_Isend(&sub_array[0], sub_XSize, MPI_INT, my_rank-xDim, my_rank, MPI_COMM_WORLD, &reqs[2]);
                MPI_Irecv(&receive_down[0], sub_XSize, MPI_INT, my_rank-xDim, my_rank-xDim, MPI_COMM_WORLD, &reqs[3]);
            }
            if ((my_rank+1)%xDim!=0) {
                for (int i=0; i<sub_XSize; i++) {
                        send_right[i] = sub_array[i*sub_YSize+sub_XSize];
                    }
                MPI_Isend(&send_right[0], sub_XSize, MPI_INT, my_rank+1, 99+my_rank, MPI_COMM_WORLD, &reqs[4]);
                MPI_Irecv(&receive_right[0], sub_XSize, MPI_INT, my_rank+1, 99+my_rank+1, MPI_COMM_WORLD, &reqs[5]);
            }
            if ((my_rank)%xDim!=0) {
                for (int i=0; i<sub_XSize; i++) {
                        send_left[i] = sub_array[i*sub_YSize];
                    }
                MPI_Isend(&send_left[0], sub_XSize, MPI_INT, my_rank-1, 999+my_rank, MPI_COMM_WORLD, &reqs[6]);
                MPI_Irecv(&receive_left[0], sub_XSize, MPI_INT, my_rank-1, 999+my_rank-1, MPI_COMM_WORLD, &reqs[7]);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    // sub_array[i*sub_YSize+j] += sub_array[j+(i+1)*sub_YSize] + sub_array[j+(i-1)*sub_YSize] + sub_array[j+1+i*sub_YSize] + sub_array[j-1+i*sub_YSize];
                        sub_array[i*sub_YSize+j] += receive_up[i] + sub_array[j+(i-1)*sub_YSize] + sub_array[j+1+i*sub_YSize] + sub_array[j-1+i*sub_YSize];


    //execution - will be performed considering sensitive situations

    for (int rank=0; rank<comm_sz; rank++) {
        if (rank == my_rank) {
            printf("Rank = %d\n", my_rank);
            for (int i=0; i<sub_XSize; i++) {
                for (int j=0; j<sub_YSize; j++) {
                    if (my_rank <xDim && (i==sub_XSize-1)){
                        if (my_rank%xDim==0)
                        {
                            sub_array[i*sub_YSize+j] += receive_up[i] + sub_array[j+(i-1)*sub_YSize] + receive_right[i] + sub_array[j-1+i*sub_YSize];
                        }
                        else if (my_rank+1)%xDim==0)
                        {
                            sub_array[i*sub_YSize+j] += receive_up[i] + sub_array[j+(i-1)*sub_YSize] + sub_array[j+1+i*sub_YSize] + receive_left[i];
                        }
                        else
                            sub_array[i*sub_YSize+j] += receive_up[i] + sub_array[j+(i-1)*sub_YSize] + receive_right[i] + receive_left[i];
                    }
                }
                    else
                        if (my_rank%xDim==0)
                        {
                            sub_array[i*sub_YSize+j] += sub_array[j+(i+1)*sub_YSize] + receive_down[i] + receive_right[i] + sub_array[j-1+i*sub_YSize];
                        }
                        else if (my_rank+1)%xDim==0)
                        {
                            sub_array[i*sub_YSize+j] += sub_array[j+(i+1)*sub_YSize] + receive_down[i] + sub_array[j+1+i*sub_YSize] + receive_left[i];
                        }
                        else
                            sub_array[i*sub_YSize+j] += sub_array[j+(i+1)*sub_YSize] + receive_down[i] + receive_right[i] + receive_left[i];
                }
            }
        MPI_Barrier(MPI_COMM_WORLD);
        }
    }


    MPI_Finalize();

    return 0;
}