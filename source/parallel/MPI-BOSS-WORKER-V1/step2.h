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
polynomial_element * generate_sieving_interval(mpz_t N, int pes);
void sieve(polynomial_element * SI, prime_element * FB);
void sieving_step(polynomial_element *SI, prime_element *FB, mpz_t N, polynomial_element *SI_SAVE, int fbs, int pes, int rank, MPI_Status status, int block_size, int num_proc);
int **alloc_2d_int(int rows, int cols);

void prime_divide(polynomial_element* SI, int** power_storage, int size_SI, int size_FB, int smallest, int prime, int* counter, int i, int block_size);
int prime_find_min(int size_SI, mpz_t a, mpz_t p, mpz_t min, mpz_t T, mpz_t r, mpz_t idx, int base_index, int block_size, int rank);

void solve_matrix();
