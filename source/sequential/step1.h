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
    int p, a, b;
};

void getprimes(int l, mpz_t N, prime_element * primes);
void shanktonellis(mpz_t N, mpz_t p);
vector<mpz_t> removeNotQR(vector<mpz_t> primes, mpz_t n);


  long s = 0;
  long q = p - 1;
  while ((q & 1) == 0) { q /= 2; ++s; }
  if (s == 1) {
    long r = pow_mod(n, (p+1)/4, p);
    if ((r * r) % p == n) return r;
    return 0;
  }
  // Find the first quadratic non-residue z by brute-force search
  long z = 1;
  while (pow_mod(++z, (p-1)/2, p) != p - 1);
  long c = pow_mod(z, q, p);
  long r = pow_mod(n, (q+1)/2, p);
  long t = pow_mod(n, q, p);
  long m = s;
  while (t != 1) {
    long tt = t;
    long i = 0;
    while (tt != 1) {
      tt = (tt * tt) % p;
      ++i;
      if (i == m) return 0;
    }
    long b = pow_mod(c, pow_mod(2, m-i-1, p-1), p);
    long b2 = (b * b) % p;
    r = (r * b) % p;
    t = (t * b2) % p;
    c = b2;
    m = i;
  }
  if ((r * r) % p == n) return r;
  return 0;

