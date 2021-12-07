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
    unsigned long num_proc = 0;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, (int*) &num_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, (int*) &rank);
    //Set value of N
    mpz_t N;
    mpz_init(N);
    mpz_set_str(N, "103603", 10);

    size_t j = mpz_sizeinbase (N, 10);
    int size = static_cast<int>(j);

    //current size of sieving interval
    //int pes = 386*size*size -23209.3*size + 2352768;
    int pes = 80000;
    block_size = 80000/(10*(num_proc-1));
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


    // for (int i = 0; i < fbs; i++) {
    //     // cout << FB[i].p << "-"<< FB[i].a << "-"<< FB[i].b << endl;
    // }

    cout << "Before sieve" << endl;
    //Write complete columns as rows into a text file
    sieving_step(SI, FB, N, SISAVE, fbs, pes, rank, status, block_size, num_proc);

    cout << "AFter sieve" << endl;

    // for (int i = 0; i < pes; i++){
    //   if (relations[fbs][i] == 1){
    //     relation_count += 1;
    //   }
    // }
    // cout << relation_count << endl;
    // cout << "PES: " << pes << endl;
    // polynomial_element * SISUB = new polynomial_element[relation_count];

    //initialize exponent matrix
    // int count = 0;


    //Fill out the  polynomial array
    // int sieve_number = 0;
    // for (int i = 0; i < pes; i++){
    //   if (relations[fbs][i] == 1){
    //     mpz_set(SISUB[count].poly, SISAVE[i].poly);
    //     count += 1;
    //   }
    // }

    // ofstream fb;
    // fb.open ("Smooth_Num.txt");
    // for (int i = 0; i < count; i++){
    //   fb << mpz_get_str(NULL, 10, SISUB[i].poly) << endl;
    // }
    // fb.close();
    //
    // fb.open ("Power_Matrix.txt");
    // for (int i = 0; i < pes; i++){
    //     if (relations[fbs][i] == 1){
    //       for (int j = 0; j < fbs; j++){
    //         fb << relations[j][i];
    //       }
    //       fb << endl;
    //     }
    // }
    // fb.close();
    //
    // fb.open ("Expo_Matrix.txt");
    // for (int i = 0; i < pes; i++){
    //     if (relations[fbs][i] == 1){
    //       for (int j = 0; j < fbs; j++){
    //         relations[j][i] = relations[j][i] % 2;
    //         fb << relations[j][i];
    //       }
    //       fb << endl;
    //     }
    // }
    // fb.close();

    MPI_Finalize();
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
void sieving_step(polynomial_element *SI, prime_element *FB, mpz_t N, polynomial_element *SI_SAVE, int fbs, int pes, int rank, MPI_Status status, int block_size, int num_proc){

  //intialize values and recompute value for T
  mpz_t T, p, a, b, idx, r, min1, min2, poly;
  int size_FB = fbs;
  int size_SI = pes;
  int power;
  int size;
  int continue_sieving = 1;
  int block_base;

  //Find smallest value of T such that T^2 - N >= 0
  mpz_t Tsq, res;
  mpz_init(Tsq);
  mpz_init(res);
  mpz_init(a);
  mpz_init(b);
  mpz_init(p);
  mpz_init(idx);
  mpz_init(r);
  mpz_init(min1);
  mpz_init(min2);
  mpz_init(poly);
  mpz_init_set_ui(T, 1);
  mpz_root(T, N, 2); // T = sqrt(N)
  mpz_add_ui(T, T, 1); //Buffer T by one to ensure non negativity

  int** power_storage;
  power_storage = alloc_2d_int(size_FB + 1, block_size);

  if (rank == 0){

    ofstream smooth_num_file;
    smooth_num_file.open ("Smooth_Num.txt");

    ofstream expo_matrix_file;
    expo_matrix_file.open ("Expo_Matrix.txt");

    ofstream bit_matrix_file;
    bit_matrix_file.open("Bit_Matrix.txt");


    cout << "We are in the rank 0" << endl;
    int total_counter = 0;
    int new_relations = 0;
    MPI_Recv(&new_relations, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
    int location = status.MPI_SOURCE;

    int** relations_storage;
    relations_storage = alloc_2d_int(new_relations, size_FB);
    size = new_relations * size_FB;
    MPI_Recv(&relations_storage[0][0], size, MPI_INT, location, 0, MPI_COMM_WORLD, &status);

    cout << "The goods has been received" << endl;



   //  while (total_counter < size_FB + 10){
   //    MPI_Recv(&new_relations, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
   //    int location = status.MPI_SOURCE;
   //    cout << "We find location : " << location << endl;
   //    total_counter += new_relations;
   //    cout << new_relations << endl;
   //
   //
   //    size = block_size * (size_FB + 1);
   //    MPI_Recv(&power_storage[0][0], size, MPI_INT, location, 0, MPI_COMM_WORLD, &status);
   //
   //
   //    //change name to relations later
   //    for (int i = 0; i < block_size; i++){
   //      if (power_storage[fbs][i] == 1){
   //        smooth_num_file << mpz_get_str(NULL, 10, SI_SAVE[i+block_base].poly) << endl;
   //        for (int j = 0; j < fbs; j++){
   //          expo_matrix_file << power_storage[j][i];
   //          cout << power_storage[j][i];
   //          power_storage[j][i] = power_storage[j][i] % 2;
   //          bit_matrix_file << power_storage[j][i];
   //        }
   //        cout << endl;
   //        expo_matrix_file << endl;
   //        bit_matrix_file << endl;
   //      }
   //    }
   //
   //    // if (total_counter > size_FB + 10){
   //    //   continue_sieving = 0;
   //    //   MPI_Send(&continue_sieving, 1, MPI_INT, location, 0, MPI_COMM_WORLD);
   //    // } else{
   //    //   continue_sieving = 1;
   //    //   MPI_Send(&continue_sieving, 1, MPI_INT, location, 0, MPI_COMM_WORLD);
   //    // }
   //
   // }
   // smooth_num_file.close();
   // expo_matrix_file.close();
   // bit_matrix_file.close();
   // cout << "We here" << endl;
 } else{

   // int val = rank;
   // MPI_Send(&val, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
   //
   // int** temp_array = alloc_2d_int(size_FB + 1, block_size);
   // for (int i = 0; i < size_FB + 1; i++){
   //   for (int j = 0; j < block_size; j++){
   //     temp_array[i][j] = rank;
   //     temp_array[size_FB][j] = 1;
   //   }
   // }
   //
   //
   // if (rank == 1 || rank == 3){
   //   for (int j = 0; j < block_size; j++){
   //       temp_array[size_FB][j] = 0;
   //   }
   // }
   //
   //
   //
   // size = (size_FB + 1)*block_size;
   //
   // MPI_Send(&temp_array[0][0], size, MPI_INT, 0, 0, MPI_COMM_WORLD);



    block_base = (rank-1)*block_size;
    cout << block_size << endl;

    while (continue_sieving == 1){


      //initialize matrix which stores max power of a prime which divides any given polynomial evaluation
      int** power_storage;
      power_storage = alloc_2d_int(size_FB + 1, block_size);

      // counting number of relations
      int counter = 0;


      ofstream prime_divide_result;
      prime_divide_result.open("Prime_Divide.txt");


      for (int i = 0; i < size_FB; i++){

        //convert p, a, b to mpz types
        mpz_init_set_ui(p, FB[i].p);
        mpz_init_set_ui(a, FB[i].a);
        mpz_init_set_ui(b, FB[i].b);

        //cout << "Before finding minimum primes" << endl;
        //find smallest indices such that the polynomial evaluation at that index is divisble by p
        unsigned long init1 = prime_find_min(size_SI, a, p, min1, T, r, idx, block_base, block_size, rank);

        unsigned long init2 = prime_find_min(size_SI, b, p, min2, T, r, idx, block_base, block_size, rank);

        //  ofstream smooth_num_file;
        //  smooth_num_file.open ("Smooth_Num.txt");

        prime_divide_result << "We find init1: " << init1 << " and init2: " << init2 << "for prime " << p << endl;

        //prepare for for loop
        int step = mpz_get_ui (p);
        mpz_t res;
        mpz_init(res);

        //go ahead and do all of the divisions
          if (init1 < size_SI + 1){
            prime_divide(SI, power_storage, size_SI, size_FB, init1, step, &counter, i, block_size, block_base);
          }

          if (init2 < size_SI + 1){
            prime_divide(SI, power_storage, size_SI, size_FB, init2, step, &counter, i, block_size, block_base);
          }




        //cout << "count:" << counter << endl;
        }
        prime_divide_result.close();

        int relations_amt = 0;
        for (int j = 0; j < block_size; j++){
          if (power_storage[size_FB][j] == 1){
            relations_amt += 1;
          }
        }

        string* smooth_nums = new string[relations_amt];
        int** relations = alloc_2d_int(relations_amt, size_FB);

        // ofstream reduced_relations;
        // reduced_relations.open("reduce_relations.txt");

        // ofstream smooth_nums_test;
        // smooth_nums_test.open("smooth_nums_test.txt");

        reduce_and_transpose(smooth_nums, relations, power_storage, block_size, size_FB, SI_SAVE, block_base);
        //
        // for (int i = 0; i < relations_amt; i++){
        //   for (int j = 0; j < size_FB; j++){
        //     cout << relations[i][j];
        //   }
        //   cout << endl;
        // }

        // for (int f= 0; f < relations_amt; f++){
        //   cout << smooth_nums[f] << endl;
        // }
        //
        // smooth_nums_test.close();
        // // reduced_relations.close();


        int size = (size_FB) * relations_amt;
        int ret = MPI_Send(&relations_amt, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        ret = MPI_Send(&relations[0][0], size, MPI_INT, 0, 0, MPI_COMM_WORLD);
        //cout << "Return val: " << ret << endl;
      //  cout << "Sent something" << endl;
        //MPI_Recv(&continue_sieving, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        block_base = block_base + num_proc * block_size;
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    cout << "YEH KYA HAIN" << endl;
    //return power_storage;
  }

int **alloc_2d_int(int rows, int cols) {
    int *data = (int *)malloc(rows*cols*sizeof(int));
    int **array= (int **)malloc(rows*sizeof(int*));
    for (int i=0; i<rows; i++)
        array[i] = &(data[cols*i]);
    return array;
}

//divide all the numbers by the prime when applicable
void prime_divide(polynomial_element* SI, int** power_storage, int size_SI, int size_FB, int smallest, int prime, int* counter, int i, int block_size, int block_base){
    int power = 0;
    int step = prime;

    //iterate through all numbers which will be divisble by prime
    for (int j = smallest; j < block_base + block_size; j = j + step){

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


unsigned long prime_find_min(int size_SI, mpz_t a, mpz_t p, mpz_t min, mpz_t T, mpz_t r,
   mpz_t idx, int base_index, int block_size, int rank){
     unsigned long temp = size_SI+1;

  for (unsigned long j = base_index; j < base_index + block_size; j++){
      //for each prime, figure out the smallest polynomial expressed as (a+pk)^2 - N
      mpz_set_ui(idx, j);
      mpz_add(idx, idx, T);
      mpz_sub(idx, idx, a);

      //check if i congruent to a - T mod p
      mpz_mod(r, idx, p);
      int q = mpz_cmp_ui(r, (unsigned long) 0);
      if (q == 0){
        mpz_set_ui(min, j);
        temp = mpz_get_ui (min);
        break;
      }
    }
    return temp;
}

void reduce_and_transpose(string* smooth_nums, int** relations, int** power_storage, int block_size, int size_FB, polynomial_element *SI, int block_base){
  cout << "Reducing and transposing" << endl;

  int s_idx = 0;
  for (int j = 0; j < block_size; j++){
    if (power_storage[size_FB][j] == 1){
      //cout << "In this loop" << endl;
      char tmp[1024];
      mpz_get_str(tmp, 10, SI[j+block_base].poly);
      smooth_nums[s_idx] = tmp;
      for (int idx = 0; idx < size_FB; idx++){
        relations[s_idx][idx] = power_storage[idx][j];
      }
      s_idx++;
    }
  }
}
