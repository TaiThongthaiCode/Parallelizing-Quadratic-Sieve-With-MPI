/*

Tai Thongthai and Tarang Saluja

"Life is not a problem to be solved, but a reality to be experienced" - Benjamin Tennyson 10;

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



#include "step2.h"



int main(int argc, char *argv[]){

    mpz_t N;
    mpz_init_set_ui(N, 69);

    polynomial_element * SI = generate_sieving_interval(N);
    prime_element * FB = load(N);


    int size = 80000;
    for (int i = 0; i < 150; i++) {
        cout << FB[i].p << endl;
        cout << FB[i].a << endl;
        cout << FB[i].b << endl;
    }

    return 0;
}

prime_element * load(mpz_t N){

    size_t j = mpz_sizeinbase (N, 10);
    int size = static_cast<int>(j);
    int fbs;
    if (size < 25) {
      fbs = 150;
    } else if (size < 50) {
      fbs = 1560;
    } else if (size < 75) {
      fbs = 6000;
    } else if (size < 100) {
      fbs = 60000;
    }

    prime_element * FB = (prime_element *)calloc(fbs, sizeof(prime_element));

    string line;
    ifstream myfile ("factorbase.txt");
    if (myfile.is_open()){
        int idx = 0;
        while ( getline (myfile,line) ){
          char str_array[line.length()];
          strcpy(str_array, line.c_str());
          char* number = strtok(str_array, " ");
          int item = 0;
          while (number != NULL){
            if (item == 0){
                FB[idx].p = atoi(number);
                item++;
            } else if (item == 1) {
                FB[idx].a = atoi(number);
                item++;
            } else if (item == 2) {
                FB[idx].b = atoi(number);
                item = 0;
            }
            number = strtok(NULL, " ");
        }
          idx++;
        }
        myfile.close();
    }

    else cout << "Unable to open file";

    return FB;

}

polynomial_element * generate_sieving_interval(mpz_t N){

    int size = 80000;
    polynomial_element * SI = new polynomial_element[size];

    mpz_t T, Tsq, res;
    mpz_init(Tsq);
    mpz_init(res);
    mpz_init_set_ui(T, 1);
    mpz_root(T, N, 2); // T = sqrt(N)
    mpz_add_ui(T, T, 1); //Buffer T by one to ensure non negativity

    for (int i = 0; i < size; i++){
        //res = T^2 - N
        mpz_pow_ui(Tsq, T, 2);
        mpz_sub(res, Tsq, N);

        mpz_init(SI[i].poly);
        mpz_set(SI[i].poly, res);

        mpz_add_ui(T, T, 1);

    }

    return SI;
}

void* sieving_step(polynomial_element *SI, prime_element *FB, mpz_t N){
  mpz_t T, p, a, b, diff, idx1, r1, min1, min2;
  mpz_init_set_ui(T, 1);
  mpz_root(T, N, 2); // T = sqrt(N)
  mpz_add_ui(T, T, 1); //Buffer T by one to ensure non negativity

  int counter = 0;

  for (int i = 1; i < prime_element.length(); i++){
    //convert p, a, b to primes
    mpz_init_set_ui(p, prime_element[i].p);
    mpz_init_set_ui(a, prime_element[i].a);
    mpz_init_set_ui(b, prime_element[i].b);


    for (int j; j < polynomial_element.length(); j++){
      //for each prime, figure out the smallest polynomial expressed as (a+pk)^2 - N
      mpz_set_ui(idx1, j);
      mpz_add(idx1, idx1, T);
      mpz_add_ui(idx1, idx1, 1);
      mpz_sub(idx1, idx1, a);
      mpz_mod(r1, idx1, p);
      if (r1 == 0){
        min1 = idx1;
      }
      break;
    }

    for (int j; j < polynomial_element.length(); j++){
      mpz_set_ui(idx1, j);
      mpz_add(idx1, idx1, T);
      mpz_add_ui(idx1, idx1, 1);
      mpz_sub(idx1, idx1, b);
      mpz_mod(r1, idx1, p);

      if (r1 == 0){
        min2 = idx1;
      }
      break;
    }

    int step = unsigned long int mpz_get_ui (const mpz_t p);
    int init1 = unsigned long int mpz_get_ui (const mpz_t a);
    int init2 = unsigned long int mpz_get_ui (const mpz_t b);

    for (int j = init1; j < polynomial_element.length(); j = j + step){
      if (prime_element[i].p != 1){
          prime_element[i].p = prime_element[i].p/p;
          if (prime_element[i].p == 1){
            counter += 1;
          }
      }
    }

    for (int j = init2; j < polynomial_element.length(); j = j + step){
      if (prime_element[i].p != 1){
          prime_element[i].p = prime_element[i].p/p;
          if (prime_element[i].p == 1){
            counter += 1;
          }
      }
    }

    if (counter >= prime_element.length() + 10){
      break;
    }
  }


}
