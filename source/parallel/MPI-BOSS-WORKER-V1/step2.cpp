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
#include <mpi.h>

#include </usr/include/eigen3/Eigen/Dense>

#include "step2.h"


int main(int argc, char *argv[]){
    int block_size = 0;
    unsigned int rank = 0;
    unsigned long num_proc= 0;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, (int*) &num_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, (int*) &rank);
    //Set value of N
    mpz_t N;
    mpz_init(N);
    mpz_set_str(N, "10141516181932272496105647768491411342341", 10);

    size_t j = mpz_sizeinbase (N, 10);
    int size = static_cast<int>(j);

    //current size of sieving interval
    //int pes = 386*size*size -23209.3*size + 2352768;
    int pes = 80000;
    block_size = 80000/(10*num_proc);
    polynomial_element * SI = generate_sieving_interval(N, pes);
    polynomial_element * SISAVE = generate_sieving_interval(N, pes);

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
    int** relations = sieving_step(SI, FB, N, fbs, pes, rank, status, block_size, num_proc);


    //count number of relations to create array of polynomials which have been factorized
    for (int i = 0; i < pes; i++){
      if (relations[fbs][i] == 1){
        relation_count += 1;
      }
    }
    cout << relation_count << endl;
    cout << "PES: " << pes << endl;
    polynomial_element * SISUB = new polynomial_element[relation_count];

    //initialize exponent matrix
    int count = 0;


    //Fill out the  polynomial array
    int sieve_number = 0;
    for (int i = 0; i < pes; i++){
      if (relations[fbs][i] == 1){
        mpz_set(SISUB[count].poly, SISAVE[i].poly);
        count += 1;
      }
    }

    ofstream fb;
    fb.open ("Smooth_Num.txt");
    for (int i = 0; i < count; i++){
      fb << mpz_get_str(NULL, 10, SISUB[i].poly) << endl;
    }
    fb.close();

    fb.open ("Power_Matrix.txt");
    for (int i = 0; i < pes; i++){
        if (relations[fbs][i] == 1){
          for (int j = 0; j < fbs; j++){
            fb << relations[j][i];
          }
          fb << endl;
        }
    }
    fb.close();

    fb.open ("Expo_Matrix.txt");
    for (int i = 0; i < pes; i++){
        if (relations[fbs][i] == 1){
          for (int j = 0; j < fbs; j++){
            relations[j][i] = relations[j][i] % 2;
            fb << relations[j][i];
          }
          fb << endl;
        }
    }
    fb.close();

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
polynomial_element * generate_sieving_interval(mpz_t N, int pes){

    int size = pes;
    polynomial_element * SI = new polynomial_element[size];

    //Find smallest value of T such that T^2 - N >= 0
    mpz_t T, Tsq, res;
    mpz_init(Tsq);
    mpz_init(res);
    mpz_init_set_ui(T, 1);
    mpz_root(T, N, 2); // T = sqrt(N)
    mpz_add_ui(T, T, 1); //Buffer T by one to ensure non negativity

    //Evaluate for 80,000 values and add to array.
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


//step where we repeatedly divide until we have the required number of relations
int** sieving_step(polynomial_element *SI, prime_element *FB, mpz_t N, int fbs, int rank, int pes, MPI_Status status, int block_size, int num_proc){

  //intialize values and recompute value for T
  mpz_t T, p, a, b, idx, r, min1, min2;
  int size_FB = fbs;
  int size_SI = pes;
  int power;

  int continue_sieving = 1;
  int block_base;
  int block_relations = (size_SI * block_size)/(num_proc * (size_FB + 10)) + 1;

  if (rank == 0){
    cout << "We are in the rank 0" << endl;
    int total_counter = 0;
    int new_relations = 0;
    while (total_counter < size_FB + 10){
      MPI_Recv(&new_relations, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
      total_counter += new_relations;

      int location = status.MPI_SOURCE;
      if (total_counter > size_FB + 10){
        continue_sieving = 0;
        MPI_Send(&continue_sieving, 1, MPI_INT, location, 0, MPI_COMM_WORLD);
      } else{
        continue_sieving = 1;
        MPI_Send(&continue_sieving, 1, MPI_INT, location, 0, MPI_COMM_WORLD);
      }
    }
  } else{

    block_base = (rank-1)*block_size;

    while (continue_sieving == 1){

      //initialize matrix which stores max power of a prime which divides any given polynomial evaluation
      int** power_storage = new int*[size_FB+1];
      for (int i = 0; i < size_FB+1; i++){
        power_storage[i] = new int[block_size];
        for (int j = 0; j< block_size; j++){
          power_storage[i][j] = 0;
        }
      }

      // counting number of relations
      int counter = 0;

      for (int i = 0; i < size_FB; i++){
        //convert p, a, b to mpz types
        mpz_init_set_ui(p, FB[i].p);
        mpz_init_set_ui(a, FB[i].a);
        mpz_init_set_ui(b, FB[i].b);

        //find smallest indices such that the polynomial evaluation at that index is divisble by p
        int init1 = prime_find_min(size_SI, a, p, min1, T, r, idx, block_base, block_size);
        int init2 = prime_find_min(size_SI, b, p, min2, T, r, idx, block_base, block_size);

        //prepare for for loop
        int step = mpz_get_ui (p);
        mpz_t res;
        mpz_init(res);

        //go ahead and do all of the divisions
        prime_divide(SI, power_storage, size_SI, size_FB, init1, step, &counter, i, block_size);
        prime_divide(SI, power_storage, size_SI, size_FB, init2, step, &counter, i, block_size);

        cout << "count:" << counter << endl;


        // //if there are enough relations, it is time to return
        if (counter >= block_relations){
          cout << "We are here!" << endl;
          int size = (size_FB + 1) * block_size;
          MPI_Send(&power_storage, size. MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        }
      }
    }
  }
  return NULL;








  // if (counter < size_FB + 1){
  //   cout << "Not enough relations found" << endl;
  //   exit(0);
  // }
  // //return power_storage;


}

//divide all the numbers by the prime when applicable
void prime_divide(polynomial_element* SI, int** power_storage, int size_SI, int size_FB, int smallest, int prime, int* counter, int i, int block_size){
    int power = 0;
    int step = prime;

    //iterate through all numbers which will be divisble by prime
    for (int j = smallest; j < smallest + block_size; j = j + step){

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
        }
      }
    }
}


int prime_find_min(int size_SI, mpz_t a, mpz_t p, mpz_t min, mpz_t T, mpz_t r, mpz_t idx, int base_index, int block_size){
  for (int j = base_index; j < base_index + block_size; j++){
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
