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

    int pes = 80000;
    polynomial_element * SI = generate_sieving_interval(N);

    size_t j = mpz_sizeinbase (N, 10);
    int size = static_cast<int>(j);
    int fbs;

    string line;
    ifstream myfile ("fb_size.txt");
    if (myfile.is_open()){
      getline(myfile, line);
      char str_array[line.length()];
      strcpy(str_array, line.c_str());
      fbs = atoi(str_array);
    }

    prime_element * FB = load(N, fbs);


    // for (int i = 0; i < 150; i++) {
    //     cout << FB[i].p << "-"<< FB[i].a << "-"<< FB[i].b << endl;
    // }

    // for (int i = 0; i < 80000; i++){
    //   cout << i << ":" << SI[i].poly << endl;
    // }

    int** relations = sieving_step(SI, FB, N, fbs, pes);
    for (int i = 0; i < pes; i++){
        if (relations[i][fbs] == 1){
          cout << relations[i][fbs] << endl;
        }
    }

    return 0;
}

prime_element * load(mpz_t N, int fbs){


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

int** sieving_step(polynomial_element *SI, prime_element *FB, mpz_t N, int fbs, int pes){
  mpz_t T, p, a, b, idx, r, min1, min2;
  int size_FB = fbs;
  int size_SI = pes;
  int power;

  mpz_init(idx);
  mpz_init(r);
  mpz_init(min1);
  mpz_init(min2);
  mpz_init_set_ui(T, 1);
  mpz_root(T, N, 2); // T = sqrt(N)
  mpz_add_ui(T, T, 1); //Buffer T by one to ensure non


  int** power_storage = new int*[size_SI];
  for (int i = 0; i < size_SI; i++){
    power_storage[i] = new int[size_FB + 1];
    for (int j = 0; j< size_FB + 1; j++){
      power_storage[i][j] = 0;
    }
  }


  int counter = 0;

  for (int i = 0; i < size_FB; i++){
    //convert p, a, b to primes

    if (FB[i].p == 0){
      break;
    }

    mpz_init_set_ui(p, FB[i].p);
    mpz_init_set_ui(a, FB[i].a);
    mpz_init_set_ui(b, FB[i].b);

    cout << i << endl;

    // cout << p << "-"<< a << "-"<< b<< endl;

    for (int j = 0; j < size_SI; j++){
      //for each prime, figure out the smallest polynomial expressed as (a+pk)^2 - N
      mpz_set_ui(idx, j);
      mpz_add(idx, idx, T);
      mpz_sub(idx, idx, a);
      //int test = mpz_cmp_ui(p, 0);

      //cout << "idx:" << idx << endl;
      mpz_mod(r, idx, p);

      int q = mpz_cmp_ui(r, 0);

      if (q == 0){
        mpz_set_ui(min1, j);
        // cout <<  << endl;
        //cout << min1 << endl;
        break;
      }
    }


    //for each prime, figure ou the smallest polynomial expressed as (b+pk)^2 - N
    for (int j = 0; j < size_SI; j++){
      //for each prime, figure out the smallest polynomial expressed as (a+pk)^2 - N
      mpz_set_ui(idx, j);
      mpz_add(idx, idx, T);
      mpz_sub(idx, idx, b);
      //int test = mpz_cmp_ui(p, 0);

      //cout << "idx:" << idx << endl;
      mpz_mod(r, idx, p);

      int q = mpz_cmp_ui(r, 0);

      if (q == 0){
        mpz_set_ui(min2, j);
        // cout <<  << endl;
        //cout << min2 << endl;
        break;
      }
    }
    //prepare for for loop
    int step = mpz_get_ui (p);
    int init1 =  mpz_get_ui (min1);
    int init2 = mpz_get_ui (min2);

    // cout << "step:" << p << "-" << "init1:" << init1 << "-" << "init2:" << init2 << endl;

    //starting with smallest polynomial exprsesed as (a+pk)^2 - N, iterate through (a+ pk')^2 - N
    mpz_t res;
    mpz_init(res);

    for (int j = init1; j < size_SI; j = j + step){
      // cout << "idx:" << j << endl;
      int q = mpz_cmp_ui(SI[j].poly, 1);
      if (q != 0){
        power = 0;
        q = mpz_divisible_ui_p(SI[j].poly, step);
        // cout << "Sieving Polynomial" << SI[j].poly << "step: " << step << endl;
        while (q != 0){
          mpz_divexact_ui(SI[j].poly, SI[j].poly, step);
          power += 1;
          q = mpz_divisible_ui_p(SI[j].poly, step);
        }
        // cout << power << endl;
        power_storage[i][j] = power;
        q = mpz_cmp_ui(SI[j].poly, 1);
        if (q == 0){
          counter += 1;
          power_storage[size_FB][j] = 1;
          // cout << "We are here!" << endl;
        }
      }
    }

    for (int j = init2; j < size_SI; j = j + step){
      // cout << "idx:" << j << endl;
      int q = mpz_cmp_ui(SI[j].poly, 1);
      if (q != 0){
        power = 0;
        q = mpz_divisible_ui_p(SI[j].poly, step);
        // cout << "Sieving Polynomial" << SI[j].poly << "step: " << step << endl;
        while (q != 0){
          mpz_divexact_ui(SI[j].poly, SI[j].poly, step);
          power += 1;
          q = mpz_divisible_ui_p(SI[j].poly, step);
        }
        // cout << power << endl;
        power_storage[i][j] = power;
        q = mpz_cmp_ui(SI[j].poly, 1);
        if (q == 0){
          counter += 1;
          power_storage[size_FB][j] = 1;
          // cout << "We are here!" << endl;
        }
      }
    }

    cout << "count:" << counter << endl;


    // //if there are enough relations, it is time to return
    if (counter >= size_FB + 10){
      return power_storage;
    }
  }
  return power_storage;


}
