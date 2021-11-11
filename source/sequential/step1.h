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

vector<mpz_t> getprimes(mpz_t l);
vector<mpz_t> getprimesold(mpz_t l, mpz_t N);
vector<mpz_t> removeNotQR(vector<mpz_t> primes, mpz_t n);
