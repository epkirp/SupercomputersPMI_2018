// ExtraTask.cpp: определяет точку входа для консольного приложения.
#include "stdafx.h"
#include "Header.h"

int main(int argc, char**argv)
{
	upper_triangle_write_file("upper.txt");
	int** u_mtrx = read_file("upper.txt");
	lower_triangle_write_file("lower.txt");
	int** l_mtrx = read_file("lower.txt");

	for (size_t sz_b(1) ; sz_b <= N; ++sz_b)
	{
		if (N%sz_b == 0)
		{
			FILE* dfile = fopen("dividers.txt", "a");
			fprintf(dfile, "%d\n", sz_b);
			printf("\nblock size: %d\n", sz_b);

			int* A = upper_triangle_convert_to_vec(sz_b, u_mtrx);
			int* B = lower_triangle_convert_to_vec(sz_b, l_mtrx);
			int* C1 = mulpiplication(A, B, sz_b);
			delete[] C1;
			int* C2 = mulpiplication_parallel(A, B, sz_b);
			delete[] C2;
			int* C3 = mulpiplication_parallel_divide(A, B, sz_b);
			delete[] C3;

			delete[] A;
			delete[] B;

			fclose(dfile);
		}
	}
	matrix_free(u_mtrx);
	matrix_free(l_mtrx);
	printf("\n");
	return 0;
}

