/*
Tai Thongthai and Tarang Saluja

"If I had 8 hours to chop a tree, I would spend the first 6 hours sharpening my axe" - Harlene Quinzell

*/

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <gmp.h>
#include <stdarg.h>
#include <obstack.h>

#include "step1.h"

using namespace std;

int main(int argc, char *argv[]){

    //setting N
    mpz_t N;
    mpz_init(N);
    mpz_set_str(N, "69", 10);



    size_t j = mpz_sizeinbase (N, 10);

    int size = static_cast<int>(j);

    int fbs, l;
    if (size < 25) {
      fbs = 150;
    } else if (size < 50) {
      fbs = 1560;
    } else if (size < 75) {
      fbs = 6000;
    } else if (size < 100) {
      fbs = 60000;
    }
    l = 2*fbs*log(2*fbs);

    prime_element * primes = (prime_element *)calloc(fbs, sizeof(prime_element));
    getprimes(l, N, primes);


    //writing a factorbase file
    ofstream fb;
    fb.open ("factorbase.txt");

    for (int i = 0; i < fbs; i++){
      fb << primes[i].p <<" " << primes[i].a << " " << primes[i].b << endl;
    }

    fb.close();


    return 0;
}


/*
Sieve of erasthones
*/
void getprimes(int l, mpz_t N, prime_element * primes){

    int idx = 0;

    mpz_t pp;
    mpz_init(pp);

    bool *truth = new bool[l+1];
    for (int i = 0; i < l+1; i++){
        truth[i] = true;
    }
    truth[0] = false;
    truth[1] = false;
    for (int i = 0; i*i < l; i++){
        if (truth[i]){
            for (int j=i*i; j<=l; j=j+i){
                truth[j] = false;
            }
        }
    }
    for (int i = 1; i<l+1; i++){
        if (truth[i]){
            mpz_set_ui(pp, i);

            //if qr then shank tonelli 
            if (mpz_legendre(N, pp) != -1){
              primes[idx].p = i;
              shanktonellis(N, &primes[idx]);
              idx++;
            }
        }
    }
}


void shanktonellis(mpz_t N, prime_element *prime){

 int res_int; 
 mpz_t b, q, res;
 mpz_init(b);
 mpz_init(q);
 mpz_init(res);
 mpz_set(b, N);
 mpz_set_ui(q, prime->p);
 mpz_sqrtm(res, b, q);


 res_int = mpz_get_ui(res);
 prime->a = res_int;
 prime->b = prime->p - prime->a;

 return;
  
}


int mpz_sqrtm(mpz_ptr rop, mpz_t a, mpz_t q)
{
  mpz_t g, temp, t, gInv, qDiv, h, b;
  int i, s, e, y;

  if (mpz_legendre(a,q) == -1) {
    return 0;
  }

  mpz_init(g); mpz_init(temp);

  mpz_sub_ui(temp, q, 1);

  while (mpz_legendre( g, q ) != -1)
    {
      mpz_random(g, 2);
      mpz_mod(g, g, temp);
      mpz_add_ui(g, g, 1);
    }

  mpz_init_set(t, q);
  mpz_sub_ui(t, t, 1);
  s = mpz_scan1(t, 0);
  mpz_tdiv_q_2exp(t, t, s);

  e = 0;

  mpz_init(gInv);
  if (!mpz_invert(gInv, g, q))
    return 0;

  mpz_init(qDiv);
  mpz_init(h);
  for (i = 2; i <= s; i++)
    {
      mpz_powm_ui(temp, gInv, e, q);
      mpz_mul(h, a, temp);
      mpz_sub_ui(temp, q, 1);
      mpz_tdiv_q_2exp(qDiv, temp, i);
      mpz_powm(temp, h, qDiv, q);
      if (mpz_cmp_ui(temp, 1 ) != 0)
{
 y = 1 << (i - 1);
 e += y;
}
    }

  mpz_powm_ui(temp, gInv, e, q);
  mpz_mul(h, a, temp);

  mpz_init(b);
  mpz_add_ui(t, t, 1);
  mpz_tdiv_q_2exp(t, t, 1);
  mpz_powm(h, h, t, q);
  mpz_powm_ui(g, g, (int)e/2, q);
  mpz_mul(b, g, h);
  mpz_mod(b, b, q);
  mpz_set(rop, b);

  mpz_clear(g); mpz_clear(temp); mpz_clear(t); mpz_clear(gInv); mpz_clear(qDiv); mpz_clear(h); mpz_clear(b);
  return 1;
}

