#include<stdio.h>
#include<stdlib.h>
#include "my_octave.h"

int main(void)
{
	// Declaration of data needed in the program:
	// array of matrices and number of matrices
	int ***all_matrices, nr_matrices = 0;

	// Declaration of array which stores dimensions of all matrices
	dimensions *matrices_dim;

	// Reading the input and calling functions to solve the problem
	char c = getc(stdin);
	while (c != 'Q') {
		switch (c) {
		// Loading a matrix in the array
		case 'L':
			load_matrix(&all_matrices, &matrices_dim, &nr_matrices);
			break;
		// Printing dimensions of a matrix
		case 'D':
			print_dimensions(matrices_dim, nr_matrices);
			break;
		// Printing elements of a matrix
		case 'P':
			print_matrix(all_matrices, matrices_dim, nr_matrices);
			break;
		// Resizing a matrix
		case 'C':
			resize_matrix(all_matrices, matrices_dim, nr_matrices);
			break;
		// Multiplying two matrices
		case 'M':
			multiply_matrices(&all_matrices, &matrices_dim, &nr_matrices);
			break;
		// Sorting all matrices
		case 'O':
			sort_matrices(all_matrices, matrices_dim, nr_matrices);
			break;
		// Transposing a matrix
		case 'T':
			transpose_matrix(all_matrices, matrices_dim, nr_matrices);
			break;
		// Deleting a matrix from the array
		case 'F':
			delete_matrix(&all_matrices, &matrices_dim, &nr_matrices);
			break;
		// Using Strassen's algorithm to multiply 2 matrices
		case 'S':
			add_strassen(&all_matrices, &matrices_dim, &nr_matrices);
			break;
		// The input is wrong and doesn't match any defined command
		default:
			printf("Unrecognized command\n");
			break;
		}
		getc(stdin);
		c = getc(stdin);
	}

	// Freeing the memory of initialized matrices and arrays
	for (int i = 0; i < nr_matrices; i++)
		free_matrix(matrices_dim[i].rows, all_matrices[i]);
	free(matrices_dim);
	free(all_matrices);

	return 0;
}

// Function for freeing the memory used by a matrix
void free_matrix(int rows, int **matrix)
{
	for (int i = 0; i < rows; i++)
		free(matrix[i]);
	free(matrix);
}

// Function for allocating the necessary memory of a matrix
int **alloc_matrix(int rows, int columns)
{
	int **matrix;
	matrix = malloc(rows * sizeof(int *));
	for (int i = 0; i < rows; i++)
		matrix[i] = malloc(columns * sizeof(int));
	return matrix;
}

// Function for reading the elements of a matrix and calculating the
// sum of them.
int **read_matrix(int rows, int columns, int *sum)
{
	int **matrix;
	matrix = alloc_matrix(rows, columns);
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < columns; j++) {
			scanf("%d", (*(matrix + i) + j));
			*sum += *(*(matrix + i) + j) % MOD;
		}
	}
	return matrix;
}

// Function for loading the matrix in the memory and storing it inside
// the array of all matrices.
void load_matrix(int ****all_matrices, dimensions **matrices_dim,
				 int *nr_matrices)
{
	// Allocating memory in arrays
	if (*nr_matrices == 0) {
		*all_matrices = malloc(1 * sizeof(int **));
		*matrices_dim = malloc(1 * sizeof(dimensions));
	} else {
		*all_matrices = realloc(*all_matrices, ((*nr_matrices) + 1) *
							   sizeof(int **));
		*matrices_dim = realloc(*matrices_dim, ((*nr_matrices) + 1) *
							   sizeof(dimensions));
	}

	// Initializing a new matrix with number of rows, columns
	// and sum of elements
	int **matrix, rows, columns, sum = 0;
	scanf("%d%d", &rows, &columns);
	matrix = read_matrix(rows, columns, &sum);
	if (!matrix) {
		exit(-1);
	} else {
		// Storing the elements in arrays
		(*all_matrices)[*nr_matrices] = matrix;
		(*matrices_dim)[*nr_matrices].rows = rows;
		(*matrices_dim)[*nr_matrices].columns = columns;
		sum = sum % MOD;
		if (sum < 0)
			sum += MOD;
		// Storing the sum in the array as to shorten the execution
		// time when sorting the matrices
		(*matrices_dim)[*nr_matrices].sum_elements = sum;
	}
	(*nr_matrices)++;
}

// Function for printing the elements of a matrix
void print_matrix(int ***all_matrices, dimensions *matrices_dim,
				  int nr_matrices)
{
	int position;
	scanf("%d", &position);

	// Checking if the index is correct
	if (position >= nr_matrices || position < 0) {
		printf("No matrix with the given index\n");
	} else {
		for (int i = 0; i < matrices_dim[position].rows; i++) {
			for (int j = 0; j < matrices_dim[position].columns; j++)
				printf("%d ", all_matrices[position][i][j]);
			printf("\n");
		}
	}
}

// Function for resizing the matrix
void resize_matrix(int ***all_matrices,
				   dimensions *matrices_dim, int nr_matrices)
{
	// Initializing variables for details of the new resized matrix
	int number_resized_rows, number_resized_columns, position;
	int *resized_rows, *resized_columns;

	// Reading data (number of resized rows and columns, their indices)
	scanf("%d", &position);
	scanf("%d", &number_resized_rows);
	resized_rows = malloc(number_resized_rows * sizeof(int));
	for (int i = 0; i < number_resized_rows; i++)
		scanf("%d", resized_rows + i);
	scanf("%d", &number_resized_columns);
	resized_columns = malloc(number_resized_columns * sizeof(int));
	for (int i = 0; i < number_resized_columns; i++)
		scanf("%d", resized_columns + i);

	// Checking if the index is correct
	if (position >= nr_matrices || position < 0) {
		printf("No matrix with the given index\n");
	} else {
		// Initializing a new matrix with number of rows, columns
		// and sum of elements
		int **new_matrix, **curr = all_matrices[position];
		new_matrix = alloc_matrix(number_resized_rows, number_resized_columns);
		int sum = 0;

		// Copying the elements in the new resized matrix and calculating
		// their sum
		for (int i = 0; i < number_resized_rows; i++) {
			for (int j = 0; j < number_resized_columns; j++) {
				new_matrix[i][j] = curr[resized_rows[i]][resized_columns[j]];
				sum += new_matrix[i][j] % MOD;
				sum %= MOD;
			}
		}
		sum %= MOD;
		if (sum < 0)
			sum += MOD;

		// Freeing the memory used for the old matrix
		free_matrix(matrices_dim[position].rows, all_matrices[position]);

		// Storing the new resized matrix
		all_matrices[position] = new_matrix;
		matrices_dim[position].rows = number_resized_rows;
		matrices_dim[position].columns = number_resized_columns;
		matrices_dim[position].sum_elements = sum;
	}

	// Freeing arrays used for storing indeces of rows and columns
	free(resized_rows);
	free(resized_columns);
}

// Function for sorting all matrices
void sort_matrices(int ***all_matrices, dimensions *matrices_dim,
				   int nr_matrices)
{
	int **aux;
	dimensions aux_struct;
	for (int i = 0; i < nr_matrices - 1; i++) {
		for (int j = i + 1; j < nr_matrices; j++) {
			if (matrices_dim[i].sum_elements > matrices_dim[j].sum_elements) {
				aux = all_matrices[i];
				all_matrices[i] = all_matrices[j];
				all_matrices[j] = aux;

				aux_struct = matrices_dim[i];
				matrices_dim[i] = matrices_dim[j];
				matrices_dim[j] = aux_struct;
			}
		}
	}
}

// Function for transposing a matrix
void transpose_matrix(int ***all_matrices, dimensions *matrices_dim,
					  int nr_matrices)
{
	int position;
	scanf("%d", &position);

	// Checking if the index is correct
	if (position >= nr_matrices || position < 0) {
		printf("No matrix with the given index\n");
	} else {
		// Initializing a new matrix with number of rows and columns
		int **new_matrix, number_rows, number_columns;
		number_rows = matrices_dim[position].columns;
		number_columns = matrices_dim[position].rows;
		new_matrix = alloc_matrix(number_rows, number_columns);

		// Copying the elements in the new transposed matrix
		for (int i = 0; i < number_columns; i++) {
			for (int j = 0; j < number_rows; j++)
				new_matrix[j][i] = all_matrices[position][i][j];
		}

		// Freeing the memory used for the old matrix
		free_matrix(number_columns, all_matrices[position]);

		// Storing the new transposed matrix
		matrices_dim[position].columns = number_columns;
		matrices_dim[position].rows = number_rows;
		all_matrices[position] = new_matrix;
	}
}

// Function for deleting a matrix
void delete_matrix(int ****all_matrices,
				   dimensions **matrices_dim, int *nr_matrices)
{
	int position;
	scanf("%d", &position);

	// Checking if the index is correct
	if (position >= *nr_matrices || position < 0) {
		printf("No matrix with the given index\n");
		return;
	}

	// Freeing the matrix and removing the pointer from the array
	int number_rows = (*matrices_dim)[position].rows;
	free_matrix(number_rows, (*all_matrices)[position]);

	// Moving the remaining matrices one position to the left
	for (int i = position + 1; i < *nr_matrices; i++) {
		(*all_matrices)[i - 1] = (*all_matrices)[i];
		(*matrices_dim)[i - 1] = (*matrices_dim)[i];
	}

	// Last pointer is now NULL, as the matrices were moved to the left
	(*all_matrices)[(*nr_matrices) - 1] = NULL;

	// Updating dimensions and allocated memory for arrays
	(*all_matrices) = realloc((*all_matrices), ((*nr_matrices) - 1)
							  * sizeof(int **));
	(*matrices_dim) = realloc((*matrices_dim), ((*nr_matrices) - 1)
							  * sizeof(dimensions));
	(*nr_matrices)--;
}

// Function for multiplying 2 matrices
void multiply_matrices(int ****all_matrices, dimensions **matrices_dim,
					   int *nr_matrices)
{
	int first_position, second_position;
	scanf("%d%d", &first_position, &second_position);

	// Checking if the indices are correct
	if (first_position >= *nr_matrices || second_position >=
		*nr_matrices || first_position < 0 || second_position < 0) {
		printf("No matrix with the given index\n");
		return;
	}

	// Checking if the dimensions are compatible
	if ((*matrices_dim)[first_position].columns !=
		(*matrices_dim)[second_position].rows) {
		printf("Cannot perform matrix multiplication\n");
		return;
	}

	// Initializing a new matrix with number of rows, columns
	// and sum of elements
	int **result, number_rows, number_columns, number_aux;
	int sum = 0;
	number_rows = (*matrices_dim)[first_position].rows;
	number_columns = (*matrices_dim)[second_position].columns;
	number_aux = (*matrices_dim)[second_position].rows;
	result = alloc_matrix(number_rows, number_columns);

	*all_matrices = realloc(*all_matrices, ((*nr_matrices) + 1) *
							sizeof(int **));
	*matrices_dim = realloc(*matrices_dim, ((*nr_matrices) + 1) *
							sizeof(dimensions));

	// Using the algorithm for multiplying 2 matrices
	for (int i = 0; i < number_rows; i++) {
		for (int j = 0; j < number_columns; j++) {
			result[i][j] = 0;
			for (int k = 0; k < number_aux; k++)
				result[i][j] += (((*all_matrices)[first_position][i][k] % MOD) *
								((*all_matrices)[second_position][k][j] % MOD))
								% MOD;
			result[i][j] %= MOD;

			// Calculating sum for elements of the new matrix
			if (result[i][j] < 0)
				result[i][j] = result[i][j] + MOD;
			sum += result[i][j];
			sum %= MOD;
		}
	}

	// Storing the result of multiplying 2 matrices as another matrix
	(*all_matrices)[*nr_matrices] = result;

	// Storing the dimensions and sum of elements for the new matrix
	(*matrices_dim)[*nr_matrices].rows = number_rows;
	(*matrices_dim)[*nr_matrices].columns = number_columns;
	sum = sum % MOD;
	if (sum < 0)
		sum += MOD;
	(*matrices_dim)[*nr_matrices].sum_elements = sum;
	(*nr_matrices)++;
}

// Function for printing elements of a matrix
void print_dimensions(dimensions *matrices_dim, int nr_matrices)
{
	int position;
	scanf("%d", &position);

	// Checking if the index is correct
	if (position >= nr_matrices || position < 0)
		printf("No matrix with the given index\n");
	else
		printf("%d %d\n", matrices_dim[position].rows,
			   matrices_dim[position].columns);
}

// Function for calculating the sum of 2 matrices
int **add_matrices(int **first_matrix, int **second_matrix, int rows)
{
	// Initializing a temporary matrix
	int **temp_matrix = alloc_matrix(rows, rows);

	// Using the algorithm for sum of 2 matrices
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < rows; j++) {
			temp_matrix[i][j] = (first_matrix[i][j] % MOD) +
								(second_matrix[i][j] % MOD);
			temp_matrix[i][j] = temp_matrix[i][j] % MOD;
			if (temp_matrix[i][j] < 0)
				temp_matrix[i][j] += MOD;
		}
	}

	// Returning the resulted matrix
	return temp_matrix;
}

// Function for subtracting one matrix from another
int **subtract_matrices(int **first_matrix, int **second_matrix, int rows)
{
	// Initializing a temporary matrix
	int **temp_matrix = alloc_matrix(rows, rows);

	// Using the algorithm for subtracting one matrix from another
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < rows; j++) {
			temp_matrix[i][j] = first_matrix[i][j] - second_matrix[i][j];
			temp_matrix[i][j] = temp_matrix[i][j] % MOD;
			if (temp_matrix[i][j] < 0)
				temp_matrix[i][j] += MOD;
		}
	}

	// Returning the resulted matrix
	return temp_matrix;
}

// Function for calculating the partial products of the 8 sub-matrices
// recursively
strassen_product find_products(int **first_matrix11, int **first_matrix12,
							   int **first_matrix21, int **first_matrix22,
							   int **second_matrix11, int **second_matrix12,
							   int **second_matrix21, int **second_matrix22,
							   int k)
{
	// Using Strassen's algorithm
	strassen_product all_p;
	int **aux1, **aux2;

	aux1 = subtract_matrices(second_matrix12, second_matrix22, k);
	all_p.P1 = strassen(first_matrix11, aux1, k);
	free_matrix(k, aux1);

	aux1 = add_matrices(first_matrix11, first_matrix12, k);
	all_p.P2 = strassen(aux1, second_matrix22, k);
	free_matrix(k, aux1);

	aux1 = add_matrices(first_matrix21, first_matrix22, k);
	all_p.P3 = strassen(aux1, second_matrix11, k);
	free_matrix(k, aux1);

	aux1 = subtract_matrices(second_matrix21, second_matrix11, k);
	all_p.P4 = strassen(first_matrix22, aux1, k);
	free_matrix(k, aux1);

	aux1 = add_matrices(first_matrix11, first_matrix22, k);
	aux2 = add_matrices(second_matrix11, second_matrix22, k);
	all_p.P5 = strassen(aux1, aux2, k);
	free_matrix(k, aux1);
	free_matrix(k, aux2);

	aux1 = subtract_matrices(first_matrix12, first_matrix22, k);
	aux2 = add_matrices(second_matrix21, second_matrix22, k);
	all_p.P6 = strassen(aux1, aux2, k);
	free_matrix(k, aux1);
	free_matrix(k, aux2);

	aux1 = subtract_matrices(first_matrix11, first_matrix21, k);
	aux2 = add_matrices(second_matrix11, second_matrix12, k);
	all_p.P7 = strassen(aux1, aux2, k);
	free_matrix(k, aux1);
	free_matrix(k, aux2);

	// Returning the resulted 8 sub-matrices as a struct
	return all_p;
}

// Function for calculating the final matrix
int **find_quarters(strassen_product all_p, int k)
{
	// Using Strassen's algorithm
	int **aux1, **aux2;

	aux1 = add_matrices(all_p.P5, all_p.P4, k);
	aux2 = add_matrices(aux1, all_p.P6, k);

	int **final_matrix = alloc_matrix(k * 2, k * 2);
	int **final_matrix11 = subtract_matrices(aux2, all_p.P2, k);
	free_matrix(k, aux2);
	free_matrix(k, aux1);

	int **final_matrix12 = add_matrices(all_p.P1, all_p.P2, k);
	int **final_matrix21 = add_matrices(all_p.P3, all_p.P4, k);
	aux1 = add_matrices(all_p.P5, all_p.P1, k);
	aux2 = subtract_matrices(aux1, all_p.P3, k);
	int **final_matrix22 = subtract_matrices(aux2, all_p.P7, k);
	free_matrix(k, aux2);
	free_matrix(k, aux1);

	for (int i = 0; i < k; i++) {
		for (int j = 0; j < k; j++) {
			final_matrix[i][j] = final_matrix11[i][j];
			final_matrix[i][j + k] = final_matrix12[i][j];
			final_matrix[k + i][j] = final_matrix21[i][j];
			final_matrix[k + i][k + j] = final_matrix22[i][j];
		}
	}

	// Freeing the allocated memory
	free_matrix(k, final_matrix11);
	free_matrix(k, final_matrix12);
	free_matrix(k, final_matrix21);
	free_matrix(k, final_matrix22);

	// Returning the resulted final matrix
	return final_matrix;
}

// Function for multiplying 2 matrices using Strassen's algorithm
int **strassen(int **first_matrix, int **second_matrix, int n)
{
	// Elementary case of the algorithm
	if (n == 1) {
		int **final_matrix = alloc_matrix(1, 1);
		final_matrix[0][0] = first_matrix[0][0] * second_matrix[0][0];
		return final_matrix;
	}

	// Initializing the final (result) matrix
	int **final_matrix;

	// Splitting the matrix in 8 sub-matrices
	int k = n / 2;
	int **first_matrix11 = alloc_matrix(k, k);
	int **first_matrix12 = alloc_matrix(k, k);
	int **first_matrix21 = alloc_matrix(k, k);
	int **first_matrix22 = alloc_matrix(k, k);
	int **second_matrix11 = alloc_matrix(k, k);
	int **second_matrix12 = alloc_matrix(k, k);
	int **second_matrix21 = alloc_matrix(k, k);
	int **second_matrix22 = alloc_matrix(k, k);

	// Filling the sub-matrices
	for (int i = 0; i < k; i++) {
		for (int j = 0; j < k; j++) {
			first_matrix11[i][j] = first_matrix[i][j];
			first_matrix12[i][j] = first_matrix[i][k + j];
			first_matrix21[i][j] = first_matrix[k + i][j];
			first_matrix22[i][j] = first_matrix[k + i][k + j];
			second_matrix11[i][j] = second_matrix[i][j];
			second_matrix12[i][j] = second_matrix[i][k + j];
			second_matrix21[i][j] = second_matrix[k + i][j];
			second_matrix22[i][j] = second_matrix[k + i][k + j];
		}
	}

	// Calculating the partial products of the 8 sub-matrices recursively
	strassen_product all_p = find_products(first_matrix11, first_matrix12,
										   first_matrix21, first_matrix22,
										   second_matrix11, second_matrix12,
										   second_matrix21, second_matrix22, k);

	// Calculating the final matrix
	final_matrix = find_quarters(all_p, k);

	// Freeing the allocated memory
	free_matrix(k, first_matrix11);
	free_matrix(k, first_matrix12);
	free_matrix(k, first_matrix21);
	free_matrix(k, first_matrix22);
	free_matrix(k, second_matrix11);
	free_matrix(k, second_matrix12);
	free_matrix(k, second_matrix21);
	free_matrix(k, second_matrix22);
	free_matrix(k, all_p.P1);
	free_matrix(k, all_p.P2);
	free_matrix(k, all_p.P3);
	free_matrix(k, all_p.P4);
	free_matrix(k, all_p.P5);
	free_matrix(k, all_p.P6);
	free_matrix(k, all_p.P7);

	// Returning the resulted matrix
	return final_matrix;
}

// Function for adding the result of Strassen's multiplication algorithm
// in the array of matrices
void add_strassen(int ****all_matrices, dimensions **matrices_dim,
				  int *nr_matrices)
{
	int first_position, second_position;
	scanf("%d%d", &first_position, &second_position);

	// Checking if the index is correct
	if (first_position >= *nr_matrices || second_position >=
		*nr_matrices || first_position < 0 || second_position < 0) {
		printf("No matrix with the given index\n");
		return;
	}

	// Checking if the dimensions are compatible
	if ((*matrices_dim)[first_position].rows !=
		(*matrices_dim)[second_position].rows) {
		printf("Cannot perform matrix multiplication\n");
		return;
	}

	// Initializing a new matrix with number of rows, columns
	// and sum of elements
	int n = (*matrices_dim)[first_position].rows, sum = 0;
	int **final_matrix;
	final_matrix = strassen((*all_matrices)[first_position],
							(*all_matrices)[second_position], n);
	*all_matrices = realloc(*all_matrices, ((*nr_matrices) + 1) *
							sizeof(int **));
	*matrices_dim = realloc(*matrices_dim, ((*nr_matrices) + 1) *
							sizeof(dimensions));

	// Storing the result of multiplying 2 matrices as another matrix
	(*all_matrices)[*nr_matrices] = final_matrix;

	// Storing the dimensions of elements for the new matrix
	(*matrices_dim)[*nr_matrices].rows = n;
	(*matrices_dim)[*nr_matrices].columns = n;

	// Calculating the sum of elements for the new matrix
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			sum += *(*(final_matrix + i) + j) % MOD;
	sum = sum % MOD;
	if (sum < 0)
		sum += MOD;

	// Storing the sum of elements for the new matrix
	(*matrices_dim)[*nr_matrices].sum_elements = sum;
	(*nr_matrices)++;
}
