#include <mpi.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>

// upper bound for random number
#define UPPER_BOUND 10

// declear the functions
void allocate_memory(int w, int l, float **mx_flat, float ***mx);
void initialize_mx(int w, int l, float ***mx);
void print_mx(int w, int l, float ***mx);



int main(int argc, char **argv) {
	float *init_mx_flat; // one dimension version of the initial matrix
	float **init_mx; // two dimension version of the initial matrix
	float *init_max_flat_copy; // the copy of the one dimension version of the initial matrix
	float **init_mx_copy;	// the copy of the two dimension version of the initial matrix
	float *result_unsorted_flat; // one dimension version of the unsorted result
	float **result_unsorted; // two dimension version of the unsorted result
	float *result_flat; // one dimension version of the sorted result
	float **result; // two dimension version of the sorted result

	// one dimension version of the task matrix (one pair, a and b)
	float *mx_flat_a;	
	float *mx_flat_b;
	// two dimension version of the task matrix (one pair, a and b)
	float **mx_a;	
	float **mx_b;

	int n, m; // the size of the matrix, which is given by users 

	int myid;
	int numprocs;
	MPI_Status status;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	// check the command line arguments
	if (argc != 3) {
		if (myid == 0) {
			printf("Wrong number of arguments.\n");
			printf("Please enter the command in the following format:\n");
			printf("mpirun -np [proc num] pairwise_inner_products [the number of row] [the number of column]\n");
			printf("Note: [proc num] should be odd number; [the number of rwo] / [proc num] = 0");
			printf("\n");
		}
		goto EXIT;
	}

	// parse the command line arguments
	n = atoi(argv[1]);
	m = atoi(argv[2]);

	// check the validity of the arguments
	if ((numprocs % 2 != 1) || (n % numprocs != 0)) {
		if (myid == 0) {
			printf("Illegal arguments.\n");
			printf("Please enter the command in the following format:\n");
			printf("mpirun -np [proc num] pairwise_inner_products [the number of row] [the number of column]\n");
			printf("Note: [proc num] should be odd number; [the number of rwo] / [proc num] = 0");
			printf("\n");
		}
		goto EXIT;
	}

	if (myid == 0) {
		allocate_memory(n, m, &init_mx_flat, &init_mx);
		allocate_memory(n, m, &init_max_flat_copy, &init_mx_copy);

		initialize_mx(n, m, &init_mx);

		printf("The initial grid: \n");
		print_mx(n, m, &init_mx);

		
		// if there is only one process, the sequential computation is performed
		if (numprocs == 1) {

		}
	}

EXIT:
	MPI_Finalize();
	return 0;
}

// memory allocation for grid
void allocate_memory(int w, int l, float **mx_flat, float ***mx) {
	int count = w * l;
	*mx_flat = (float *)malloc(sizeof(float) * count);
	*mx = (float **)malloc(sizeof(float *) * w);
	int i;

	for (i = 0; i < w; i++) {
		(*mx)[i] = &((*mx_flat)[i * l]);
	}
}

// initialize the matrix
void initialize_mx(int w, int l, float ***mx) {
	time_t s;
	srand((unsigned)time(&s));
	int i, j;

	for (i = 0; i < w; i++) {
		for (j = 0; j < l; j++) {
			(*mx)[i][j] = ((float)rand() / (float)(RAND_MAX)) * UPPER_BOUND;
			(*mx)[i][j] = (int)(100.0 * (*mx)[i][j] + 0.5) / 100.0;
		}
	}
}

// print the matrix
void print_mx(int w, int l, float ***mx) {
	int i, j;

	for (i = 0; i < w; i++) {
		for (j = 0; j < l; j++) {
			printf("%.2f ", (*mx)[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}