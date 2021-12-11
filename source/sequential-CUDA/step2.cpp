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

#include </usr/include/eigen3/Eigen/Dense>

#include "step2.h"

using namespace Eigen;

int main(int argc, char *argv[]){

    //Set value of N
    mpz_t N;
    mpz_init(N);
    mpz_set_str(N, "10141516181932272496105647768491411342341", 10);

    size_t j = mpz_sizeinbase (N, 10);
    int size = static_cast<int>(j);

    //current size of sieving interval
    //int pes = 386*size*size -23209.3*size + 2352768;
    int pes = 800000;
    //polynomial_element * SI = generate_sieving_interval(N, pes);
    //polynomial_element * SISAVE = generate_sieving_interval(N, pes);

    int relation_count = 0;

    //grab size of factor base from file
    int fbs;
    string line;
    ifstream myfile ("fb_size.txt");
    if (myfile.is_open()){
      getline(myfile, line);
      char str_array[line.length()];
      strcpy(str_array, line.c_str());
      fbs = atoi(str_array);
    }

    //fill in factor base
    prime_element * FB = load(N, fbs);


    for (int i = 0; i < fbs; i++) {
        cout << FB[i].p << "-"<< FB[i].a << "-"<< FB[i].b << endl;
    }


    //Write complete columns as rows into a text file
    sieving_step(FB, N, fbs, pes);

    return 0;
}

//load in the factor base elements and the corresponding values of a and b which correspond to the smallest values a and b such that a^2 - T is 0 mod p and b^2 - T is 0 mod p.
prime_element * load(mpz_t N, int fbs){

    prime_element * FB = (prime_element *)calloc(fbs, sizeof(prime_element));

    string line;
    ifstream myfile ("factorbase.txt");
    if (myfile.is_open()){
        int idx = 0;
        //go through each line and get the required data
        while ( getline (myfile,line) ){
          char str_array[line.length()];
          strcpy(str_array, line.c_str());
          //Use strtok to get one nuber at a time
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
          //iterate index, so that it is stored in the correct place
          idx++;
        }
        myfile.close();
    }

    else cout << "Unable to open file";

    return FB;
}

//create the sieving iterval of size 80, 000
polynomial_element * generate_sieving_interval(mpz_t N, int size_SSI, mpz_t T){

    int size = size_SSI;
    polynomial_element * SI = new polynomial_element[size];

    mpz_t Tsq, res;
    mpz_init(Tsq);
    mpz_init(res);


    //Evaluate for 80,000 values and add to array.
    for (int i = 0; i < size_SSI; i++){
        //res = T^2 - N
        mpz_pow_ui(Tsq, T, 2);
        mpz_sub(res, Tsq, N);

        mpz_init(SI[i].poly);
        mpz_set(SI[i].poly, res);

        mpz_add_ui(T, T, 1);
    }
    return SI;
}


//step where we repeatedly divide until we have the required number of relations
void sieving_step(prime_element *FB, mpz_t N, int fbs, int pes){

  //intialize values and recompute value for T
  mpz_t T, T_hold, p, a, b, idx, r, min1, min2;
  int size_FB = fbs;
  int size_SI = pes;
  int size_SSI = 80000;
  int power;

  mpz_init(idx);
  mpz_init(r);
  mpz_init(min1);
  mpz_init(min2);
  mpz_init(T_hold);
  mpz_init_set_ui(T, 1);
  mpz_root(T, N, 2); // T = sqrt(N)
  mpz_add_ui(T, T, 1); //Buffer T by one to ensure non
  mpz_set(T_hold, T);
  string temp;
  temp = mpz_get_str(NULL, 10, T);
  cout << "The value of T is " <<  temp << endl;

  int counter = 0;


  ofstream smooth_num_store;
  smooth_num_store.open ("Smooth_Num.txt");

  ofstream power_matrix_store;
  power_matrix_store.open("Power_Matrix.txt");

  ofstream bit_matrix_store;
  bit_matrix_store.open("Bit_Matrix.txt");


  while (counter < size_FB + 10){
    cout << "Next iter" << endl;
    temp = mpz_get_str(NULL, 10, T);
    cout << "The value of T is " <<  temp << endl;
    //initialize matrix which stores max power of a prime which divides any given polynomial evaluation
    int** power_storage = new int*[size_FB+1];
    for (int i = 0; i < size_FB+1; i++){
      power_storage[i] = new int[size_SSI];
      for (int j = 0; j< size_SSI; j++){
        power_storage[i][j] = 0;
      }
    }

    // counting number of relation

    polynomial_element * SI = generate_sieving_interval(N, size_SSI, T);
    mpz_set(T, T_hold);
    polynomial_element * SI_SAVE = generate_sieving_interval(N, size_SSI, T);
    mpz_set(T, T_hold);


    for (int i = 0; i < size_FB; i++){
      //convert p, a, b to mpz types
      mpz_init_set_ui(p, FB[i].p);
      mpz_init_set_ui(a, FB[i].a);
      mpz_init_set_ui(b, FB[i].b);

      //find smallest indices such that the polynomial evaluation at that index is divisble by p


      int init1 = prime_find_min(size_SI, a, p, min1, T_hold, r, idx);
      int init2 = prime_find_min(size_SI, b, p, min2, T_hold, r, idx);
      string temp;
      temp = mpz_get_str(NULL, 10, T_hold);

      // cout << "We find init1: "  << init1 << " and init2: " << init2 << "for the prime " << p << "and value " << temp << endl;

      //prepare for for loop
      int step = mpz_get_ui (p);
      mpz_t res;
      mpz_init(res);

      //go ahead and do all of the divisions
      prime_divide(SI, power_storage, size_SSI, size_FB, init1, step, &counter, i);
      prime_divide(SI, power_storage, size_SSI, size_FB, init2, step, &counter, i);

      // cout << "count:" << counter << endl;


      // //if there are enough relations, it is time to return
      if (counter >= size_FB + 10){
        cout << "We are here!" << endl;
        break;
        //return power_storage;
      }
    }

    for (int j = 0; j < size_SSI; j++){
      if (power_storage[size_FB][j] == 1){
        temp = mpz_get_str(NULL, 10, SI_SAVE[j].poly);
        smooth_num_store << temp << endl;
        for (int i = 0; i < size_FB; i++){
          power_matrix_store << power_storage[i][j];
          int val = power_storage[i][j];
          val = power_storage[i][j] % 2;
          bit_matrix_store << val;
        }
        power_matrix_store << endl;
        bit_matrix_store << endl;
      }
    }



    delete SI;
    delete SI_SAVE;

    mpz_add_ui(T, T, size_SSI);
    mpz_set(T_hold, T);

    for (int i = 0; i < size_FB + 1; i++){
      delete[] power_storage[i];
    }
    delete[] power_storage;

  }

  if (counter < size_FB + 1){
    cout << "Not enough relations found" << endl;
    exit(0);
  }
  //return power_storage;


}

//divide all the numbers by the prime when applicable
void prime_divide(polynomial_element* SI, int** power_storage, int size_SI, int size_FB, int smallest, int prime, int* counter, int i){
    int power = 0;
    int step = prime;

    //iterate through all numbers which will be divisble by prime
    for (int j = smallest; j < size_SI; j = j + step){

      //keep dividing until 1
      int q = mpz_cmp_ui(SI[j].poly, 1);
      if (q != 0){
        power = 0;
        q = mpz_divisible_ui_p(SI[j].poly, step);
        while (q != 0){
          mpz_divexact_ui(SI[j].poly, SI[j].poly, step);
          power += 1;
          q = mpz_divisible_ui_p(SI[j].poly, step);
        }
        //store the power
        power_storage[i][j] += power;

        q = mpz_cmp_ui(SI[j].poly, 1);
        if (q == 0){
          *counter += 1; //iterate counter if it has now been reduced to 1
          power_storage[size_FB][j] = 1;
          if (*counter >= size_FB + 10){
            return;
          }
        }
      }
    }
}


int prime_find_min(int size_SI, mpz_t a, mpz_t p, mpz_t min, mpz_t T, mpz_t r, mpz_t idx){
  for (int j = 0; j < size_SI; j++){
      //for each prime, figure out the smallest polynomial expressed as (a+pk)^2 - N
      mpz_set_ui(idx, j);
      mpz_add(idx, idx, T);
      mpz_sub(idx, idx, a);
      //check if i congruent to a - T mod p
      mpz_mod(r, idx, p);
      int q = mpz_cmp_ui(r, 0);
      if (q == 0){
        mpz_set_ui(min, j);
        break;
      }
    }
    return mpz_get_ui (min);
}



void solve_matrix(){
  MatrixXf A = MatrixXf::Random(3, 2);
  VectorXf b = VectorXf::Random(3);


}