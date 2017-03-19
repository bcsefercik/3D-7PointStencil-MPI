#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>

#define MASTER 0
#define OUTER 1.0
#define INNER -1.0

void usage(char **argv);
void initValues(double* array, int width, int height, int depth, double inner_temp, double outer_temp);
double getElement(double* array, int x, int y, int z, int width, int height, int depth);
double getElement(double* array, int offset, int x, int y, int z, int width, int height, int depth);
void setElement(double* array, double value, int x, int y, int z, int width, int height, int depth);
void computation(double *result, double *temp, int offset, int width, int height, int depth);
void print_3d(double *array, int width, int height, int depth);
double average_array(double *array, int width, int height, int depth);
double getTime();

int main(int argc, char **argv) {

	int grid_width, grid_height, grid_depth, num_iterations;
	int ghost_width, ghost_height, ghost_depth, ghost_size;

	double *ghost_temp, *ghost_result;


	double start_time, end_time, total_time;


	if (argc < 5)
		usage(argv);

	if((grid_width = atoi(argv[1])) <= 0 ||
	   (grid_height = atoi(argv[2])) <= 0 ||
	   (grid_depth = atoi(argv[3])) <= 0 ||
	   (num_iterations = atoi(argv[4])) <= 0)
		usage(argv);

	ghost_width = grid_width + 2;
	ghost_height = grid_height +2;
	ghost_depth = grid_depth + 2;

	ghost_size = ghost_width * ghost_height * ghost_depth;


	ghost_temp = (double *) malloc(ghost_size * sizeof(double));
	ghost_result = (double *) malloc(ghost_size * sizeof(double));

	if(!ghost_temp || !ghost_result){
		fprintf(stderr, "error: unable to allocate memory\n");
		exit(1);
	}

	initValues(ghost_result, ghost_width, ghost_height, ghost_depth, INNER, OUTER);

	int  num = 0;
	start_time = getTime();
	while(num<num_iterations){
		num++;

		computation(ghost_result, ghost_temp, 0, ghost_width, ghost_height, ghost_depth);

	}
	
	end_time = getTime();
	total_time = end_time - start_time;
	double average = average_array(ghost_result,ghost_width,ghost_height,ghost_depth);


	//print_3d(ghost_result,ghost_width,ghost_height,ghost_depth);

	printf("Average        : %.2f \n", average);
	printf("Time           : %f \n", total_time);

	free(ghost_temp);
	free(ghost_result);
	return 0;
}

void usage(char **argv)
{
    fprintf(stderr, "Usage: %s <width> <height> <depth> <num_iterations>\n", argv[0]);
    exit(1);
}

void initValues(double* array, int width, int height, int depth, double inner_temp, double outer_temp){
	for(int i=1; i<(width-1); i++){
		for(int j=1; j<(height-1); j++){
			for(int k=1; k<(depth-1); k++){
				setElement(array, inner_temp, i, j, k, width, height, depth);
			}
		}
	}

	for(int j=0; j<height; j++){
		for(int k=0; k<depth; k++){
			setElement(array, outer_temp, 0, j, k, width, height, depth);
			setElement(array, outer_temp, width-1, j, k, width, height, depth);
		}
	}

	for(int i=0; i<width; i++){
		for(int k=0; k<depth; k++){
			setElement(array, outer_temp, i, 0, k, width, height, depth);
			setElement(array, outer_temp, i, height-1, k, width, height, depth);
		}
	}

	for(int i=0; i<width; i++){
		for(int j=0; j<height; j++){
			setElement(array, outer_temp, i, j, 0, width, height, depth);
			setElement(array, outer_temp, i, j, depth-1, width, height, depth);
		}
	}
}

double getElement(double* array, int x, int y, int z, int width, int height, int depth){
	int index = x*height*depth + y*depth + z;
	return (array)[index];
}

double getElement(double* array, int offset, int x, int y, int z, int width, int height, int depth){
	int index = x*height*depth + y*depth + z + offset;
	return (array)[index];
}

void setElement(double* array, double value, int x, int y, int z, int width, int height, int depth){
	int index = x*height*depth + y*depth + z;
	array[index] = value;
}

void computation(double *result, double *temp, int offset, int width, int height, int depth){
	int index; 
	for(int i = 1; i<(width-1); i++){
		for(int j = 1; j<(height-1); j++){
			for(int k = 1; k<(depth-1); k++){
				index = offset + i*height*depth + j*depth + k;
				double a = result[index];
				for(int b = 0; b<130; b++){
					a = (a+getElement(result, offset, i, j, k, width, height, depth)+
								getElement(result, offset, i-1, j, k, width, height, depth)+
								getElement(result, offset, i+1, j, k, width, height, depth)+
								getElement(result, offset, i, j-1, k, width, height, depth)+
								getElement(result, offset, i, j+1, k, width, height, depth)+
								getElement(result, offset, i, j, k-1, width, height, depth)+
								getElement(result, offset, i, j, k+1, width, height, depth))/8.0;
				}
				result[index] = a;
			}
		}
	}
}

void print_3d(double *array, int width, int height, int depth){
	for(int k=0; k<depth; k++){
		for(int j=0; j<height; j++){
			for(int i=0; i<width; i++){
				printf("%3.3f  ", getElement(array, i, j, k, width, height, depth));
			}
			printf("\n");
		}
		printf("\n\n");
	}
}

double average_array(double *array, int width, int height, int depth){
	double average = 0;
	for(int i=0; i<width; i++){
		for(int j=0; j<height; j++){
			for(int k=0; k<depth; k++){
				average += getElement(array, i, j, k, width, height, depth);
			}
		}
	}
	return average/(width*height*depth);
}

double getTime(){
	const double kMicro = 1.0e-6;
	struct timeval TV;
	
	const int RC = gettimeofday(&TV, NULL);
	if(RC == -1)
	{
		printf("ERROR: Bad call to gettimeofday\n");
		return(-1);
	}
	return( ((double)TV.tv_sec) + kMicro * ((double)TV.tv_usec) );
}
