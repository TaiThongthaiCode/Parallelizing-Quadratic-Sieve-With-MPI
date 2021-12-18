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

/*
struct holding mpz_t of polynomial
*/
struct polynomial_element {
    mpz_t poly;
};

/*
Description: load in the factor base elements and the corresponding values
            of a and b which correspond to the smallest values a and b such
            that a^2 - T is 0 mod p and b^2 - T is 0 mod p.
Params: (1) mpz_t number N (product of two primes);
        (2) integer factor base size (fbs)
Return: A pointer to prime_element array. The struct prime_element holds ints,
        prime p, and their solutions a and b.
*/
prime_element * load(mpz_t N, int fbs);


/*
Description: step where we repeatedly divide our sieving interval using our
              prime factor base until we have the required number of relations
              and writing the relevant data to file. Rank 1 is the boss node,
               and all other ranks are the worker nodes. Only the workers sieve,
                the boss compiles results
Params: (1) Pointer to prime_element array (the factorbase);
        (2) mpz_t N;
        (3) int factorbasesize (fbs);
        (4) int rank of the current node
        (5) MPI_Status object, status
        (6) int block_size: the chunk size that we generate our temporary
         sieving sub intervals on
        (7) int num_proc: the number of processes allocated
Return: Nothing
Output [All Text Files]: (1) Power_Matrix.txt (matrix representation of
                              prime powers of found relations);
                          (2) Bit_Matrix.txt (Power_Matrix mod 2);
                          (3) Smooth_Numbers.txt (Numbers we have found to be
                          smooth)
*/
void sieving_step(prime_element *FB, mpz_t N, int fbs, int rank, 
                MPI_Status status, int block_size, int num_proc);

/*
Description: Receives expo matrix and smooth nums (creates bit matrix
from expo matrix)from workers and saves them to file. Unpacks smooth nums
before saving to file.
Params: (1) int pointer to counter
        (2) int factor base size
        (3) int pointer to need_more (does the boss need more relations)
        (4) ofstream reference to expo_matrix_file
        (5) ofstream reference to bit_matrix_file
        (6) ofstream reference to smooth_num_file
        (7) int pointer to dead_processes
Return: none
Output [All Text Files]: (1) Power_Matrix.txt (matrix representation of
                              prime powers of found relations);
                          (2) Bit_Matrix.txt (Power_Matrix mod 2);
                          (3) Smooth_Numbers.txt (Numbers we have found to
                          be smooth)
*/
void master_unpack_save(int* total_counter, int size_FB, int* need_more, 
                    ofstream& expo_matrix_file, ofstream& bit_matrix_file, 
                    ofstream& smooth_num_file, int* dead_processes);

/*
Description: Helper function to sieving step: the bulk of the work
              of the worker's sieving task is done here.
Params: (1) int double pointer to power_storage 2D array
        (2) int pointer to counter
        (3) int block size of our temporary sieving subinterval
        (4) mpz_t N our key
        (5) mpz_t T from polynomial Q(T) = T^2 - N
        (6) int size of factor base
        (7) prime_element pointer to our factor base array
        (8) int rank of our current node
        (9) polynomial_element pointer to our sieving subinterval array
Return: Nothing
*/
void worker_sieves(int** power_storage, int* counter, 
                int block_size, mpz_t N, mpz_t T, int size_FB,
                prime_element* FB, int* relations_amt, int rank, 
                polynomial_element* SI);

/*
Description: Generates our sieving interval given a size and starting bound.
Params: (1) mpz_t number N (product of two primes);
        (2) int size_SSI (interval size);
        (3) mpz_t number T (interval starting bound starting bound).
Return: A pointer to the polynomial_element array. The polynomial_element
        struct holds an mpz_t type (we named "poly"). This is done
        due to weird memory constraints of mpz_t types
*/
polynomial_element * generate_sieving_interval(mpz_t N, int pes, mpz_t T);

/*
Description: Function for worker to pack smooth_nums as strings and then
            send it to the boss node, using MPI_send
Params: (1) int pointer to relations amount
        (2) int size of factor base
        (3) double int pointer to relations
        (4) string pointer (array of strings) to smooth nums
Return: Nothing
*/
void worker_pack_send(int* relations_amt, int size_FB, 
                        int** relations, string* smooth_nums);

/*
Description: Allocates a contiguous 2D matrix given a row and column size
Params: (1) int row size
        (2) int column size
Returns: A int double pointer to a two dimensional contigous array of ints
Taken verbatim from
https://coderedirect.com/questions/388175/mpi-matrix-
multiplication-with-dynamic-allocation-seg-fault
*/
int **alloc_2d_int(int rows, int cols);


/*
Description: Helper function used in sieving_step. Here we divide
            elements in the sieving interval by primes when
            applicable
Params: (1) polynomial_element array pointer SI;
        (2) int double pointer to 2D-power_storage array;
        (3) int the block size of our temporary sieving subinterval
        (4) int factorbase size (size_FB)
        (5) int smallest index for which this prime divides
            this sieving element
        (6) int prime (prime of interest)
        (7) int pointer counter (keeps track of number of relations)
        (8) int i (index of the prime of interest)

Return [Explicit]: Nothing
Return [Implicit]: [Note: We write in place using passed pointers]:
                    (1) polynomial_element array pointer SI;
                    (2) int double pointer to 2D-power_storage array;
                    (3) int pointer counter (keeping track of number
                    of relations)
*/
void prime_divide(polynomial_element* SI, int** power_storage, 
            int block_size, int size_FB, int smallest, int prime, 
            int* counter, int i);


/*
Description: Helper function: for each prime, figure out the
              smallest polynomial expressed as (a+pk)^2 - N
Params: (1) int block size
        (2) mpz_t a (solution to tonelli shanks (Q(x) is divisble
            by p when evaluated at numbers of form a+kp and b+kp)
        (3) mpz_t prime p
        (4) mpz_t min
        (5) mpz_t T from Q(T)
        (6) mpz_t r
        (7) mpz_t idx
        (8) int rank of current node

Return: unsigned long: smallest polynomial expressed as (a+pk)^2 - N
*/
unsigned long prime_find_min(int block_size, mpz_t a, mpz_t p, 
                        mpz_t min, mpz_t T, mpz_t r, mpz_t idx, int rank);

/*
Description: Helper function: filter out polynomials that we have
            not found to be relations. Also changes the shape
            representation of our relations 2D array
Params: (1) string array of smooth numbers
        (2) int 2D array of relations (represented as prime powers)
        (3) int 2D array of power_storage
        (4) int block size of our temporary sieving sub_interval
        (5) int size of factor base
        (6) polynomial_element array (our sieving sub_interval)
Return [Explicit]: Nothing
Return [Implicit]: (1) smooth_nums array
                   (2) relations 2D matrix
*/
void reduce_and_transpose(string* smooth_nums, int** relations, 
                        int** power_storage, int block_size, 
                        int size_FB, polynomial_element* SI);


/*
Description: Helper function for worker nodes sieve - packs smooth nums
              into a string to prepare for MPI send
Params: (1) int pointer to string_length
        (2) string array of smooth_nums
        (3) amount of relations found
Return: string - packed smooth nums string
*/
string pack(int* string_length, string* smooth_nums, int relations_amt);


//void solve_matrix();
