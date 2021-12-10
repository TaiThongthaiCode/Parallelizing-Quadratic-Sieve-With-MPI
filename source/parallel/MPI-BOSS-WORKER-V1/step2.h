/*

Tai Thongthai and Tarang Saluja

"The measure of man is what he does with power" - Cardi B
*/

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <gmp.h>
#include <stdarg.h>
#include <obstack.h>
#include <mpi.h>

#include "step1.h"

using namespace std;

#pragma once

struct polynomial_element {
    mpz_t poly;
};


prime_element * load(mpz_t N, int fbs);
polynomial_element * generate_sieving_interval(mpz_t N, int pes, mpz_t T);
void sieve(polynomial_element * SI, prime_element * FB);
void sieving_step(prime_element *FB, mpz_t N, int fbs, int rank, MPI_Status status, int block_size, int num_proc);
int **alloc_2d_int(int rows, int cols);

void prime_divide(polynomial_element* SI, int** power_storage, int block_size, int size_FB, int smallest, int prime, int* counter, int i);
unsigned long prime_find_min(int block_size, mpz_t a, mpz_t p, mpz_t min, mpz_t T, mpz_t r, mpz_t idx, int rank);

void reduce_and_transpose( int** relations, int** power_storage, int block_size, int size_FB);

//void solve_matrix();
