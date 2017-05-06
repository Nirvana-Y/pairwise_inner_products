#include <mpi.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>

// upper bound for random number
#define UPPER_BOUND 10

// declear the functions
void allocate_memory(int w, int l, float **mx_flat, float ***mx);
void allocate_memory_result(int w, int l, float **result_flat, float ***result);
void initialize_mx(int w, int l, float ***mx);
void print_mx(int w, int l, float ***mx);
void sequential_computation(int w, int l, float ***mx);
void sort_result(int w, int l, int numprocs, float ***unsorted_result, float ***result);



int main(int argc, char **argv) {
	float *init_mx_flat; // one dimension version of the initial matrix
	float **init_mx; // two dimension version of the initial matrix
	float *init_mx_flat_copy; // the copy of the one dimension version of the initial matrix
	float **init_mx_copy;	// the copy of the two dimension version of the initial matrix
	float *unsorted_result_flat; // one dimension version of the unsorted result
	float **unsorted_result; // two dimension version of the unsorted result
	float *result_flat; // one dimension version of the sorted result
	float **result; // two dimension version of the sorted result

	// one dimension version of the task matrixs (one pair, a and b)
	float *mx_a_flat;	
	float *mx_b_flat;
	// two dimension version of the task matrixs (one pair, a and b)
	float **mx_a;	
	float **mx_b;
	
	float inner_product;
	int row; // the number of row distributed to each process
	int n, m; // the size of the matrix, which is given by users 

	int myid;
	int numprocs;
	MPI_Status status;
	MPI_Request request;

	int i, j, k, u;

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

	row = n / numprocs;

	// allocate memory for the task matrixs for every processes
	allocate_memory(row, m, &mx_a_flat, &mx_a);
	allocate_memory(row, m, &mx_b_flat, &mx_b);

	if (myid == 0) {
		allocate_memory(n, m, &init_mx_flat, &init_mx);
		allocate_memory_result(n, m, &result_flat, &result);
		allocate_memory(row * numprocs, row * ((numprocs - 1) / 2 + 1) - 1, &unsorted_result_flat, &unsorted_result);

		initialize_mx(n, m, &init_mx);
		printf("The initial grid: \n");
		print_mx(n, m, &init_mx);
	
		// if there is only one process, the sequential computation is performed
		if (numprocs == 1) {
			sequential_computation(n, m, &init_mx);
			goto EXIT;
		}

		// distribute the tasks to all the processes, process 0 and other processes
		memcpy(mx_a_flat, init_mx_flat, sizeof(float) * row * m);

		for (i = 1; i < numprocs; i++) {
			MPI_Send(&init_mx_flat[row * m * i], row * m, MPI_FLOAT, i, 1, MPI_COMM_WORLD);
		}

	}
	else {
		allocate_memory(row, row * ((numprocs - 1) / 2 + 1) - 1, &unsorted_result_flat, &unsorted_result);

		// receive the tasks from process 0
		MPI_Recv(mx_a_flat, row * m, MPI_FLOAT, 0, 1, MPI_COMM_WORLD, &status);
	}

	for (i = 0; i < (row - 1); i++) {
		for (j = (i + 1); j < row; j++) {
			inner_product = 0;
			for (k = 0; k < m; k++) {
				inner_product = inner_product + mx_a[i][k] * mx_a[j][k];
			}
			unsorted_result[i][j - i - 1] = inner_product;
		}
	}

	MPI_Isend(mx_a_flat, row * m, MPI_FLOAT, (myid + numprocs - 1) % numprocs, 1, MPI_COMM_WORLD, &request);
	for (u = 0; u < (numprocs - 1) / 2; u++) {
		MPI_Recv(mx_b_flat, row * m, MPI_FLOAT, (myid + 1) % numprocs, 1, MPI_COMM_WORLD, &status);
		for (i = 0; i < row; i++) {
			for (j = 0; j < row; j++) {
				inner_product = 0;
				for (k = 0; k < m; k++) {
					inner_product = inner_product + mx_a[i][k] * mx_b[j][k];
				}
				unsorted_result[i][(row - 1 - i) + j + u * row] = inner_product;
			}
		}
		if (u < (numprocs - 1) / 2 - 1) {
			MPI_Isend(mx_b_flat, row * m, MPI_FLOAT, (myid + numprocs - 1) % numprocs, 1, MPI_COMM_WORLD, &request);
		}
	}

	if (myid == 0) {
		for (i = 1; i < numprocs; i++) {
			MPI_Recv(&unsorted_result_flat[row * (row * ((numprocs - 1) / 2 + 1) - 1) * i], row * (row * ((numprocs - 1) / 2 + 1) - 1), MPI_FLOAT, i, 1, MPI_COMM_WORLD, &status);
		}

		sort_result(n, m, numprocs, &unsorted_result, &result);

		printf("The parallel computation result: \n");
		for (i = 0; i < (n - 1); i++) {
			for (j = 0; j < (n - i - 1); j++) {
				printf("%.2f ", result[i][j]);
			}
			printf("\n");
		}
	}
	else {
		MPI_Send(unsorted_result_flat, row * (row * ((numprocs - 1) / 2 + 1) - 1), MPI_FLOAT, 0, 1, MPI_COMM_WORLD, &status);
	}


EXIT:
	MPI_Finalize();
	return 0;
}

// memory allocation for matrix
void allocate_memory(int w, int l, float **mx_flat, float ***mx) {
	int count = w * l;
	*mx_flat = (float *)malloc(sizeof(float) * count);
	*mx = (float **)malloc(sizeof(float *) * w);
	int i;

	for (i = 0; i < w; i++) {
		(*mx)[i] = &((*mx_flat)[i * l]);
	}
}

// memory allocation for result
void allocate_memory_result(int w, int l, float **result_flat, float ***result) {
	int count = w * (w - 1) / 2;
	int i, j;
	*result_flat = (float *)malloc(sizeof(float) * count);
	*result = (float **)malloc(sizeof(float *) * (w - 1));
	
	int *index = (int *)malloc(sizeof(int) * (w - 1));
	index[0] = 0;

	for (i = 1; i < (w - 1); i++) {
		index[i] = 0;
		for (j = (w - 1); j > (w - 1 -i); j--) {
				index[i] = index[i] + j;
		}
	}

	for (i = 0; i < (w - 1); i++) {
		(*result)[i] = &((*result_flat)[index[i]]);
	}
}

// initialize the matrix with the random float number between 0 and UPPER_BOUND round to 2 decimal places
void initialize_mx(int w, int l, float ***mx) {
	time_t s;
	srand((unsigned)time(&s));
	int i, j;

	for (i = 0; i < w; i++) {
		for (j = 0; j < l; j++) {
			// (*mx)[i][j] = ((float)rand() / (float)(RAND_MAX)) * UPPER_BOUND;
			(*mx)[i][j] = rand() % 10;
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

// sequential program for pairwise inner products
void sequential_computation(int w, int l, float ***mx) {
	float inner_product;
	float *result_flat;
	float **result;
	int i, j, k;

	allocate_memory_result(w, l, &result_flat, &result);

	for (i = 0; i < (w - 1); i++) {
		for (j = (i + 1); j < w; j++) {
			inner_product = 0;
			for (k = 0; k < l; k++) {
				inner_product = inner_product + (*mx)[i][k] * (*mx)[j][k];
			}
			result[i][j - i - 1] = inner_product;
		}
	}

	printf("The sequential computation result: \n");
	for (i = 0; i < (w - 1); i++) {
		for (j = 0; j < (w - i - 1); j++) {
			printf("%.2f ", result[i][j]);
		}
		printf("\n");
	}
}

// sort the unsorted result
void sort_result(int w, int l, int numprocs, float ***unsorted_result, float ***result) {
	int row = w / numprocs;
	float *mediate_mx_flat;
	float **mediatae_mx;
	int i, j, k, u;

	allocate_memory(row * (numprocs - 1) / 2, row * ((numprocs - 1) / 2 + 1) - 2, &mediate_mx_flat, &mediatae_mx);

	for (i = 0; i < row * (numprocs + 1) / 2; i++) {
		for (j = 0; j < (row * ((numprocs - 1) / 2 + 1) - 1); j++) {
			(*result)[i][j] = (*unsorted_result)[i][j];
		}
	}

	for (i = row * (numprocs + 1) / 2; i < w - 1; i++) {
		k = 0;
		for (j = 0; j < ((numprocs - 1) / 2 * row - 1 - k); j++) {
			(*result)[i][j] = (*unsorted_result)[i][j];
		}
		k++;
	}

	k = 0;
	for (i = row * (numprocs + 1) / 2; i < w; i++) {
		u = 0;
		for (j = ((numprocs - 1) / 2 * row - 1 - k); j < (row * ((numprocs - 1) / 2 + 1) - 1); j++) {
			mediatae_mx[k][u] = (*unsorted_result)[i][j];
			u++;
		}
		k++;
	}

	for (i = 0; i < row * (numprocs - 1) / 2; i++) {
		k = i % row;
		u = i / row * row;
		for (j = (row * ((numprocs - 1) / 2 + 1) - 1) - k; j < w - i - 1; j++) {
			(*result)[i][j] = mediatae_mx[u][i];
			u++;
		}
	}
}