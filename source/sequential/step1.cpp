/*
Tai Thongthai and Tarang Saluja

"If I had 8 hours to chop a tree, I would spend the first 6 hours sharpening my axe" - Harlene Quinzell

*/

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
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

    for (int i = 0; i < fbs; i++){
      cout << primes[i].p << endl;
    }


    return 0;
}


/*
Sieve of erectus
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
            if (mpz_legendre(N, pp)){
              primes[idx].p = i;
              idx++;
            }
        }
    }
}


void shanktonellis(mpz_t N, mpz_t p){
  int s = 0;
  int q = p - 1;
  while ((q & 1) == 0) {
    q /= 2;
    s++
  }

  if (s==1) {
    int r = 
  }

}


vector<mpz_t> removeNotQR(vector <mpz_t>primes, mpz_t n){
  vector <mpz_t> ok_primes;
//   ok_primes.push_back(primes[0]);
//   cout << "size " << primes.size() << endl;


//   mpz_t p;
//   mpz_init(p);



//   for (int i = 1; i < primes.size(); i++){

//     if (1){
//       // ok_primes.push_back(p);
//     }
//   }

  return ok_primes;
}
