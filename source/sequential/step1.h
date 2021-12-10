/*
Tai Thongthai and Tarang Saluja

"To live is to risk it all, otherwise you might as well be an inert chunk of molecules" - Richard Thaler

*/
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <gmp.h>
#include <stdarg.h>
#include <obstack.h>


using namespace std;

#pragma once

struct prime_element {
    unsigned long int p, a, b;
};

int getprimes(int l, mpz_t N, prime_element * primes, int fbs);
void shanktonellis(mpz_t N, prime_element *prime);
//int convertx2e(int x, int& e);
// int powmod(int base, int exponent, int mod);

int mpz_sqrtm(mpz_ptr rop, mpz_t a, mpz_t q);

vector<mpz_t> removeNotQR(vector<mpz_t> primes, mpz_t n);
