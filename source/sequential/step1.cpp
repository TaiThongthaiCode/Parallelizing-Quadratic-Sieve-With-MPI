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
    mpz_init_set_ui(N, 69);

    // setting e
    mpf_t e;
    mpf_init(e);

    double param;
    float te;
    param = 1;
    te = exp (param);

    mpf_set_d(e, te);

    //setting k

    //setting l (Don't know how to get optimal value due to log constraint)
    mpz_t l;
    mpz_init(l);
    mpz_mul(l, N, N); //setting l = N*N

    vector<mpz_t> primes = getprimes(l);


    return 0;
}

vector<mpz_t> getprimes(mpz_t l){
    vector<mpz_t> primes;

    mpz_t currprime;
    mpz_init(currprime);
    mpz_set_ui(currprime, 2);

    //base case, if the same
    while (mpz_cmp_ui(l, 0) != 0) {
        primes.push_back(currprime); //TODO COMPILER DOESN"T LIKE THIS LINE
        mpz_sub_ui(l, l, 1); //l = l -1
        mpz_nextprime(currprime, currprime); // get next prime and set to currprime

    }

    return primes;

}


vector<mpz_t> getprimesold(mpz_t l, mpz_t N){
    vector<mpz_t> primes;
    // bool * truth = new bool[l+1];
    // for (int i = 0; i < l+1; i++){
    //     truth[i] = true;
    // }
    // truth[0] = false;
    // truth[1] = false;
    // for (int i = 0; i*i < l; i++){
    //     if (truth[i]){
    //         for (int j=i*i; j<=l; j=j+i){
    //             truth[j] = false;
    //         }
    //     }
    // }
    // for (int i = 1; i<l+1; i++){
    //     if (truth[i]){
    //         primes.push_back(i);
    //     }
    // }
    return primes;
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
