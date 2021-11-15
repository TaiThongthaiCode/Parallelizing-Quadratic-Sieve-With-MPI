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

#include "step1.h"

using namespace std;

#pragma once

struct polynomial_element {
    mpz_t poly;
};

prime_element * load(mpz_t N);
polynomial_element * generate_sieving_interval(mpz_t N);
void sieve(polynomial_element * SI, prime_element * FB);