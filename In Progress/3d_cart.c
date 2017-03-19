#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

void setAmbient(double*** T0, int x_dim, int y_dim, int z_dim, double ambient_T, double interior_T);
void communicationBetweenThreads(double*** x, int size_tot_x, int size_tot_y, int size_tot_z, int* neighbor, MPI_Comm comm3d, MPI_Datatype matrix_type_oxz, MPI_Datatype matrix_type_oxy, MPI_Datatype matrix_type_oyz , int me, int* xs, int* ys, int* zs, int* xe, int* ye, int* ze);
void perform7PointOperation(double*** x0, double*** x, int size_x, int size_y, int size_z, int me, int* xs, int* ys, int* zs, int* xe, int* ye, int* ze);
void processToMap(int me, int *xs, int *xe, int *ys, int *ye, int *zs, int *ze, int xcell, int ycell, int zcell, int x_domains, int y_domains, int z_domains);


int main()
{

    double*** T;
    double*** T0; // a pointer that points to an array of pointers apparently ??
    double* T_alloc;
    double* T0_alloc;

    int size_x, size_y, size_z; // initializing larger array (target data)
    size_x = 12; size_y = 12; size_z = 12; // sizes of the larger array
    int x_domains, y_domains, z_domains;
    int dims[3];// initializing subdomain number in each dimension (# of dimensions being 3)
    x_domains = 3; y_domains = 3; z_domains = 3; // keep in mind # of processes must be equal to their multiplication (here being 9)
    dims[0]=x_domains;
    dims[1]=y_domains;
    dims[2]=z_domains;  // storing # of domains in each direction into another array
    int xcell, ycell, zcell, size_tot_x, size_tot_y, size_tot_z;
    xcell=(size_x/x_domains); ycell=(size_y/y_domains); zcell=(size_z/z_domains); // sizes of each cell

    size_tot_x = size_x + 2*x_domains + 2;
    size_tot_y = size_y + 2*y_domains + 2;
    size_tot_z = size_z + 2*z_domains + 2;

    T_alloc = malloc(size_tot_x * size_tot_y * size_tot_z * sizeof(*T_alloc));
    T0_alloc = malloc(size_tot_x * size_tot_y * size_tot_z * sizeof(*T0_alloc));
    double* Tfinal;
    double* Ttemp;
    Tfinal = malloc(size_x*size_y*size_z*sizeof(double));
    Ttemp = malloc(xcell*ycell*zcell*sizeof(double));


    T = malloc(size_tot_x*sizeof(double**));
    T0 = malloc(size_tot_x*sizeof(double**));

    int i, j;
    for (i = 0; i < size_tot_x; i++) {
        T[i] = malloc(size_tot_y * sizeof(**T));
        T0[i] = malloc(size_tot_y * sizeof(**T0));
        for (j = 0; j < size_tot_y; j++) {
            T[i][j] = T_alloc;
            T0[i][j] = T0_alloc;
            T_alloc += size_tot_z;
            T0_alloc += size_tot_z;
        }
    }
    int sizes[3], subsizes1[3], subsizes2[3], subsizes3[3];
    sizes[0]=size_tot_x; sizes[1]=size_tot_y; sizes[2]=size_tot_z;
    subsizes1[0]=xcell; subsizes1[1]=1; subsizes1[2]=zcell; // Create matrix data type to communicate on vertical Oxz plan //
    subsizes2[0]=1; subsizes2[1]=ycell; subsizes2[2]=zcell; // Create matrix data type to communicate on vertical Oyz plan //
    subsizes3[0]=xcell; subsizes3[1]=ycell; subsizes3[2]=1; // Create matrix data type to communicate on vertical Oxy plan //


    // time parameters
    double total_time = 15.0; //(seconds)
    double dt = 1.0; //(second) update steps

    // temperature paremeters
    double ambient_T = 80; //(degrees celsius)
    double interior_T = 25; //(degrees celsius) initial state of the computation space


    // defining communication variables - explain further what they all do
    int comm_sz, rank;
    int starts[3]; // starting points of the sub arrays
    for (int i = 0; i < 3; i++) {starts[i] = 0;}

    MPI_Datatype matrix_type_oxz, matrix_type_oxy, matrix_type_oyz;

    MPI_Comm MPI_COMM_WORLD_3D; //Cartesian communication world

    // setting boundary looping for each processor dimension to false, setting ranking reordering to false
    int loops[3];
    int reorder;
    loops[0] = 0;
    loops[1] = 0;
    loops[2] = 0;
    reorder = 0;

    // defining neighboring relations (may have to put after init !!!!!)
    int neighbor[6];
    for (i = 0; i < 6; i++){neighbor[i] = MPI_PROC_NULL;} // consider east-west  0,1 - north-south being 2,3 - etc.

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD,&comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    int *xs, *ys, *zs, *xe, *ye, *ze;
    xs = malloc(comm_sz*sizeof(int));
    xe = malloc(comm_sz*sizeof(int));
    ys = malloc(comm_sz*sizeof(int));
    ye = malloc(comm_sz*sizeof(int));
    zs = malloc(comm_sz*sizeof(int));
    ze = malloc(comm_sz*sizeof(int));

    processToMap(rank, xs, xe, ys, ye, zs, ze, xcell, ycell, zcell, x_domains, y_domains, z_domains);

    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, loops, reorder, &MPI_COMM_WORLD_3D);
    // defining exchange between processes in new cartesian communication world (shift by one)
    MPI_Cart_shift(MPI_COMM_WORLD_3D,0,1,&neighbor[0],&neighbor[1]); // x dimension: West - East neigbors //
    MPI_Cart_shift(MPI_COMM_WORLD_3D,1,1,&neighbor[2],&neighbor[3]); // y dimension: North -  South neigbors //
    MPI_Cart_shift(MPI_COMM_WORLD_3D,2,1,&neighbor[4],&neighbor[5]); // z dimension: In - Out neigbors //

    MPI_Type_create_subarray(3, sizes, subsizes1, starts, MPI_ORDER_C, MPI_DOUBLE, &matrix_type_oxz);
    MPI_Type_commit(&matrix_type_oxz);

    MPI_Type_create_subarray(3, sizes, subsizes2, starts, MPI_ORDER_C, MPI_DOUBLE, &matrix_type_oyz);
    MPI_Type_commit(&matrix_type_oyz);

    MPI_Type_create_subarray(3, sizes, subsizes3, starts, MPI_ORDER_C, MPI_DOUBLE, &matrix_type_oxy);
    MPI_Type_commit(&matrix_type_oxy);

    setAmbient(T0, size_tot_x, size_tot_y, size_tot_z, ambient_T, interior_T);

    MPI_Barrier( MPI_COMM_WORLD_3D ) ;

    communicationBetweenThreads(T0, size_tot_x, size_tot_y, size_tot_z, neighbor, MPI_COMM_WORLD_3D, matrix_type_oxz, matrix_type_oxy, matrix_type_oyz, rank, xs, ys, zs, xe, ye, ze);

    int t;
    for (t = 0; t < total_time; t=t+dt)
    {
        perform7PointOperation(T0, T, size_tot_x, size_tot_y, size_tot_z, rank, xs, ys, zs, xe, ye, ze);
        communicationBetweenThreads(T0, size_tot_x, size_tot_y, size_tot_z, neighbor, MPI_COMM_WORLD_3D, matrix_type_oxz, matrix_type_oxy, matrix_type_oyz, rank, xs, ys, zs, xe, ye, ze);
        MPI_Barrier( MPI_COMM_WORLD_3D ) ;
    }

    MPI_Gather(Ttemp,xcell*ycell*zcell, MPI_DOUBLE, Tfinal, xcell*ycell*zcell, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank==0){} // write data to file for master thread //

    for (i=0;i<=size_tot_x-1;i++){
        free(T[i]);
        free(T0[i]);}
    free(T);
    free(T0);
    free(T_alloc);
    free(T0_alloc);
    free(Tfinal);
    free(Ttemp);
    free(xs);
    free(xe);
    free(ys);
    free(ye);
    free(zs);
    free(ze);

    MPI_Type_free(&matrix_type_oxz);
    MPI_Type_free(&matrix_type_oxy);
    MPI_Type_free(&matrix_type_oyz);

    MPI_Finalize();

    return 0;


}


void setAmbient(double*** T0, int x_dim, int y_dim, int z_dim, double ambient_T, double interior_T){
  int i, j, k;
  /* Setup temperature on edges */
    for(j=0; j<=y_dim-1; j++)
    { for(k=0; k<=z_dim-1; k++)
      {       T0[0][j][k] = ambient_T;
        T0[x_dim-1][j][k] = ambient_T; }
    }

    for(i=0; i<=x_dim-1; i++)
    { for(k=0; k<=z_dim-1; k++)
      {       T0[i][0][k] = ambient_T;
        T0[i][y_dim-1][k] = ambient_T; }
    }

    for(i=0; i<=x_dim-1; i++)
    { for(j=0; j<=y_dim-1; j++)
      {       T0[i][j][0] = ambient_T;
        T0[i][j][z_dim-1] = ambient_T; }
    }


  /* Setup temp2_init inside */
  for(i=1; i<=x_dim-2; i++)
  { for(j=1; j<=y_dim-2; j++)
    {  for(k=1; k<=z_dim-2; k++)
           T0[i][j][k] = interior_T;
    }
  }}



void communicationBetweenThreads(double*** x, int size_tot_x, int size_tot_y, int size_tot_z, int* neighbor, MPI_Comm comm3d, MPI_Datatype matrix_type_oxz, MPI_Datatype matrix_type_oxy, MPI_Datatype matrix_type_oyz , int me, int* xs, int* ys, int* zs, int* xe, int* ye, int* ze) {

    int S=3, E=0, N=2, W=1, Zd=4, Zu=5;
    int flag;
    MPI_Status status;

    /********* North/South communication ************************************/
    flag = 1;
    /* Send my boundary to North and receive from South */
    MPI_Sendrecv(&x[xs[me]][ys[me]][zs[me]], 1, matrix_type_oxz ,neighbor[N],flag, &x[xs[me]][ye[me]+1][zs[me]], 1, matrix_type_oxz, neighbor[S], flag, comm3d, &status);

    /* Send my boundary to South and receive from North */
    MPI_Sendrecv(&x[xs[me]][ye[me]][zs[me]], 1, matrix_type_oxz,neighbor[S],flag, &x[xs[me]][ys[me]-1][zs[me]], 1, matrix_type_oxz,neighbor[N], flag, comm3d, &status);

    /********* Est/West communication ***************************************/
    flag = 2;
    /* Send my boundary to Est and receive from West */
    MPI_Sendrecv(&x[xe[me]][ys[me]][zs[me]], 1, matrix_type_oyz, neighbor[E], flag, &x[xs[me]-1][ys[me]][zs[me]], 1, matrix_type_oyz, neighbor[W], flag, comm3d, &status);

    /* Send my boundary to West and receive from Est */
    MPI_Sendrecv(&x[xs[me]][ys[me]][zs[me]], 1, matrix_type_oyz, neighbor[W], flag, &x[xe[me]+1][ys[me]][zs[me]], 1, matrix_type_oyz, neighbor[E], flag, comm3d, &status);

    /********* Zdown/Zup communication **************************************/
    flag = 3;
    /* Send my boundary to Zup and receive from Zdown */
    MPI_Sendrecv(&x[xs[me]][ys[me]][ze[me]], 1, matrix_type_oxy, neighbor[Zu], flag, &x[xs[me]][ys[me]][zs[me]-1], 1, matrix_type_oxy, neighbor[Zd], flag, comm3d, &status);

    /* Send my boundary to Zdown and receive from Zup */
    MPI_Sendrecv(&x[xs[me]][ys[me]][zs[me]], 1, matrix_type_oxy, neighbor[Zd], flag, &x[xs[me]][ys[me]][ze[me]+1], 1, matrix_type_oxy, neighbor[Zu], flag, comm3d, &status);}

void perform7PointOperation(double*** x0, double*** x, int size_x, int size_y, int size_z, int me, int* xs, int* ys, int* zs, int* xe, int* ye, int* ze){
  int i, j, k;

  /* Perform an explicit update on the points within the domain */
  for(k=zs[me];k<=ze[me];k++)
  { for(j=ys[me];j<=ye[me];j++)
    { for(i=xs[me];i<=xe[me];i++)
      { x[i][j][k] = (x0[i-1][j][k] + x0[i+1][j][k] +
                     x0[i][j-1][k] + x0[i][j+1][k] +
                     x0[i][j][k-1] + x0[i][j][k+1] +
                     x0[i][j][k])/7;
      }
    }
  }

  /* copy back the computed value : x  <-- x^(n+1)               */
  /*                                x0 <-- x^n                   */
  /*   and compute at the same time the 2_norm of the "residual" */
  for(k=zs[me];k<=ze[me];k++)
  { for(j=ys[me];j<=ye[me];j++)
    { for(i=xs[me];i<=xe[me];i++)
      {
        x0[i][j][k] = x[i][j][k];
      }
    }
  }}