#ifndef _MATRIX_H
#define _MATRIX_H

#define MOD 10007

// Function for freeing the memory used by a matrix
void free_matrix(int rows, int **matrix);

// Function for allocating the necessary memory of a matrix
int **alloc_matrix(int rows, int columns);

// Function for reading the elements of a matrix and calculating the
// sum of them.
int **read_matrix(int rows, int columns, int *sum);

// Struct for saving the number of elements and sum of elements in a matrix
typedef struct {
	int rows, columns, sum_elements;
} dimensions;

// Function for loading the matrix in the memory and storing it inside
// the array of all matrices.
void load_matrix(int ****all_matrices, dimensions **matrices_dim,
				 int *nr_matrices);

// Function for printing the elements of a matrix
void print_matrix(int ***all_matrices, dimensions *matrices_dim,
				  int nr_matrices);

// Function for resizing the matrix
void resize_matrix(int ***all_matrices,
				   dimensions *matrices_dim, int nr_matrices);

// Function for sorting all matrices
void sort_matrices(int ***all_matrices, dimensions *matrices_dim,
				   int nr_matrices);

// Function for transposing a matrix
void transpose_matrix(int ***all_matrices, dimensions *matrices_dim,
					  int nr_matrices);

// Function for deleting a matrix
void delete_matrix(int ****all_matrices,
				   dimensions **matrices_dim, int *nr_matrices);

// Function for multiplying 2 matrices
void multiply_matrices(int ****all_matrices, dimensions **matrices_dim,
					   int *nr_matrices);

// Function for printing elements of a matrix
void print_dimensions(dimensions *matrices_dim, int nr_matrices);

// Struct for partial products of the 8 sub-matrices
typedef struct {
	int **P1, **P2, **P3, **P4, **P5, **P6, **P7;
} strassen_product;

// Function for calculating the sum of 2 matrices
int **add_matrices(int **first_matrix, int **second_matrix, int rows);

// Function for subtracting one matrix from another
int **subtract_matrices(int **first_matrix, int **second_matrix, int rows);

// Function for calculating the partial products of the 8 sub-matrices
// recursively
strassen_product find_products(int **first_matrix11, int **first_matrix12,
							   int **first_matrix21, int **first_matrix22,
							   int **second_matrix11, int **second_matrix12,
							   int **second_matrix21, int **second_matrix22,
							   int k);

// Function for calculating the final matrix
int **find_quarters(strassen_product all_p, int k);

// Function for multiplying 2 matrices using Strassen's algorithm
int **strassen(int **first_matrix, int **second_matrix, int n);

// Function for adding the result of Strassen's multiplication algorithm
// in the array of matrices
void add_strassen(int ****all_matrices, dimensions **matrices_dim,
				  int *nr_matrices);

#endif
