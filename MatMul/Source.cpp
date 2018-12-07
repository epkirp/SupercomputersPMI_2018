#include "stdafx.h"
#include "Header.h"

void upper_triangle_write_file(const char* fname)
{
	int** A = new int* [N];
	for (size_t i(0); i < N; i++)
		A[i] = new int[N];

	for (size_t i(0); i < N; ++i)
		for (size_t j(0); j < N; ++j)
			if (i <= j) A[i][j] = rand() % 10;
			else A[i][j] = 0;

	FILE* file = fopen(fname, "w");

	for (size_t i(0); i < N; ++i)
	{
		for (size_t j(0); j < N; ++j)
			fprintf(file, "%4d", A[i][j]);
		fprintf(file, "\n");
	}
	fclose(file);
	matrix_free(A);
}

void lower_triangle_write_file(const char* fname)
{
	int** A = new int*[N];
	for (size_t i(0); i < N; i++)
		A[i] = new int[N];

	for (size_t i(0); i < N; ++i)
		for (size_t j(0); j < N; ++j)
		{
			if (i >= j) A[i][j] = rand() % 10;
			else A[i][j] = 0;
		}

	FILE* file = fopen(fname, "w");

	for (size_t i(0); i < N; ++i)
	{
		for (size_t j(0); j < N; ++j)
			fprintf(file, "%4d", A[i][j]);
		fprintf(file, "\n");
	}
	fclose(file);
	matrix_free(A);

}

void matrix_free(int **A)
{
	for (size_t i(0); i < N; i++)
		delete[] A[i];

	delete[] A;
}

int* lower_triangle_convert_to_vec(const size_t sz_b, int** A)
{
	size_t S(N / sz_b);
	int* vec = new int[(S + 1)*S / 2 * sz_b*sz_b];
	size_t i = 0, j = 0, t0 = 0, t1 = 0, k = 0;

	for (size_t t0(0); t0 < S; ++t0)
		for (size_t t1(t0); t1 < S; ++t1)
			for (size_t i(0); i < sz_b; ++i)
				for (size_t j(0); j < sz_b; ++j)
						vec[k++] = A[i + t1 * sz_b][j + t0 * sz_b];
	return vec;
}

int* upper_triangle_convert_to_vec(const size_t sz_b, int** A)
{
	size_t S(N / sz_b);
	int* vec = new int[(S + 1)*S / 2 * sz_b*sz_b];
	size_t i = 0, j = 0, t0 = 0, t1 = 0, k = 0;

	for (size_t t0(0); t0 < S; ++t0)
		for (size_t t1(t0); t1 < S; ++t1)
			for (size_t i(0); i < sz_b; ++i)
				for (size_t j(0); j < sz_b; ++j)
					vec[k++] = A[i + t0 * sz_b][j + t1 * sz_b];

	return vec;
}

int** read_file(const char* fname)
{
	int** A = new int*[N];
	for (size_t i(0); i < N; i++)
		A[i] = new int[N];

	FILE *file;
	file = fopen(fname, "r");

	for (size_t i (0); i<N; ++i)
		for (size_t j(0); j<N; ++j)
			if (!feof(file)) fscanf(file, "%d", &A[i][j]);
	fclose(file);
	return A;
}

int* mulpiplication(int* A, int* B, size_t sz_b)
{
	int* C = new int[N*N];
	for (size_t i(0); i < N*N; ++i)
		C[i] = 0;

	double t1, t2;
	t1 = omp_get_wtime();
	size_t S = N / sz_b;

	for (size_t i(0); i<S; ++i)
		for (size_t j(0); j<S; ++j)
			for (size_t k(j); k<S; ++k)
			{
				int *a = A + (i * (S + 1) - (i + 1)*i / 2 + (k - i)) * sz_b*sz_b,
					*b = B + (j * (S + 1) - (j + 1)*j / 2 + (k - j)) * sz_b*sz_b,
					*c = C + i * N*sz_b + j * sz_b*sz_b;

					for (size_t ii(0); ii < sz_b; ++ii)
						for (size_t jj(0); jj < sz_b; ++jj)
							for (size_t kk(0); kk < sz_b; ++kk)
								c[ii*sz_b + jj] += a[ii*sz_b + kk] * b[kk* sz_b + jj ];
			}

	t2 = omp_get_wtime();

	FILE* file = fopen("time.txt", "a");
	fprintf(file, "%f\n", (t2 - t1));
	fclose(file);

	return C;
}

int* mulpiplication_parallel(int* A, int* B, size_t sz_b)
{
	int* C = new int[N*N];
	for (size_t i(0); i < N*N; ++i)
		C[i] = 0;

	double t1, t2;
	t1 = omp_get_wtime();
	
#pragma omp parallel num_threads(4) 
	{
		size_t S = N / sz_b;
#pragma omp for  schedule(static)
		for (int i(0); i < S; ++i)
			for (int j(0); j < S; ++j)
				for (int k(j); k < S; ++k)
				{
					int *a = A + (i * (S + 1) - (i + 1)*i / 2 + (k - i)) * sz_b*sz_b,
						*b = B + (j * (S + 1) - (j + 1)*j / 2 + (k - j)) * sz_b*sz_b,
						*c = C + i * N*sz_b + j * sz_b*sz_b;

					for (int ii(0); ii < sz_b; ++ii)
						for (int jj(0); jj < sz_b; ++jj)
							for (int kk(0); kk < sz_b; ++kk)
								c[ii*sz_b + jj] += a[ii*sz_b + kk] * b[kk* sz_b + jj];
				}
	}
	t2 = omp_get_wtime();;

	FILE* file = fopen("time_parallel.txt", "a");
	fprintf(file, "%f\n", (t2 - t1));
	fclose(file);

	return C;
}

int* mulpiplication_parallel_divide(int* A, int* B, size_t sz_b)
{
	int* C = new int[N*N];
	for (size_t i(0); i < N*N; ++i)
		C[i] = 0;

	double t1, t2;
	t1 = omp_get_wtime();

#pragma omp parallel num_threads(4)
	{
		size_t S = N / sz_b;
#pragma omp for  schedule(static)
		for (int i(0); i < S; ++i)
			for (int j(0); j < S; ++j)
				for (int k(j); k < S; ++k)
				{
					int *a = get_block(A + (i * (S + 1) - (i + 1)*i / 2 + (k - i)) * sz_b*sz_b, sz_b),
						*b = get_block(B + (j * (S + 1) - (j + 1)*j / 2 + (k - j)) * sz_b*sz_b, sz_b),
						*c = C + i * N*sz_b + j * sz_b*sz_b;

					for (int ii(0); ii < sz_b; ++ii)
						for (int jj(0); jj < sz_b; ++jj)
							for (int kk(0); kk < sz_b; ++kk)
								c[ii*sz_b + jj] += a[ii*sz_b + kk] * b[kk* sz_b + jj];
					delete[] a;
					delete[] b;
				}
	}
	t2 = omp_get_wtime();;

	FILE* file = fopen("time_parallel_divide.txt", "a");
	fprintf(file, "%f\n", (t2 - t1));
	fclose(file);

	return C;
}

int* get_block(int* A, int sz_b) 
{
	int* res = new int[sz_b*sz_b];

	for (int i = 0; i < sz_b*sz_b; ++i)
		res[i] = A[i];

	return res;
}