#pragma once
#pragma warning(disable : 4996)
#include "stdafx.h"

#include <stdio.h>
#include <stdlib.h>
#include <ctime>
#include "omp.h"

const int N = 360;

void upper_triangle_write_file(const char* fname);
void lower_triangle_write_file(const char* fname);
int** read_file(const char* fname);
int* lower_triangle_convert_to_vec(const size_t sz_b, int** A);
int* upper_triangle_convert_to_vec(const size_t sz_b, int** A);
void matrix_free(int **A);

int* mulpiplication(int* A, int* B, const size_t sz_b);
int* mulpiplication_parallel(int* A, int* B, const size_t sz_b);
int* mulpiplication_parallel_divide(int* A, int* B, const size_t sz_b);
int* get_block(int* A, int sz_b);
