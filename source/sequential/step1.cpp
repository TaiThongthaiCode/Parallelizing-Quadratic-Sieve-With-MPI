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
    mpz_t N;
    mpz_init(N);

    double param, e;
    param = 1;
    e = exp (param);
    double k = pow (e, sqrt(2*log(N) * log(log(N)) ) );
    int l = 2*k *log(2*k) + 1;
    cout << l << endl;
    vector<mpz_t> primes = getprimes(l, N);
    // for (int i=0; i<primes.size(); i++){
    //     cout << primes[i] << endl;
    // }

    vector<mpz_t> ok_primes = removeNotQR(primes, N);

    cout << "results" << endl;
    for (int i=0; i<ok_primes.size(); i++){
        cout << primes[i] << endl;
    }




    return 0;
}


vector<mpz_t> getprimes(int l, mpz_t N){
    vector<mpz_t> primes;
    bool * truth = new bool[l+1];
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
            primes.push_back(i);
        }
    }
    return primes;
}


vector<mpz_t> removeNotQR(vector <mpz_t>primes, mpz_t n){
  vector <mpz_t> ok_primes;
  ok_primes.push_back(primes[0]);
  cout << "size " << primes.size() << endl;


  mpz_t p;
  mpz_init(p);



  for (int i = 1; i < primes.size(); i++){

    if (1){
      // ok_primes.push_back(p);
    }
  }

  return ok_primes;
}
