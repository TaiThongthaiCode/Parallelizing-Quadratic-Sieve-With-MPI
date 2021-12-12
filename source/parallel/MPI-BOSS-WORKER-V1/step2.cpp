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

    if (argc != 2) {
      cout << "Wrong amount of arguments. Please enter desired N" << endl;
      exit(1);
    }

    char* pq = argv[1];
    // cout << "Step2 now running for N=" << pq << endl;

    unsigned int rank = 0;
    unsigned long num_proc = 0;
    MPI_Status status;
    int block_size = 3;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, (int*) &num_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, (int*) &rank);
    //Set value of N
    mpz_t N;
    mpz_init(N);
    mpz_set_str(N, pq, 10);
    prime_element * FB;

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
    FB = load(N, fbs);

    //Bosses and workers interact to obtain smooth numbers
    //and factorizations
    sieving_step(FB, N, fbs, rank, status, block_size, num_proc);

    MPI_Finalize();
    // cout << "Step2 now running for N=" << pq << endl;
    return 0;
}


/*
Description: load in the factor base elements and the corresponding values
            of a and b which correspond to the smallest values a and b such
            that a^2 - T is 0 mod p and b^2 - T is 0 mod p.
Params: (1) mpz_t number N (product of two primes);
        (2) integer factor base size (fbs)
Return: A pointer to prime_element array. The struct prime_element holds ints,
        prime p, and their solutions a and b.
*/
prime_element * load(mpz_t N, int fbs){

    int idx, item;
    char * number;
    prime_element * FB = (prime_element *)calloc(fbs, sizeof(prime_element));

    string line;
    ifstream myfile ("factorbase.txt");
    if (myfile.is_open()){
        idx = 0;
        //go through each line and get the required data
        while ( getline (myfile,line) ){
          char str_array[line.length()];
          strcpy(str_array, line.c_str());
          //Use strtok to get one nuber at a time
          number = strtok(str_array, " ");
          item = 0;
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
        delete [] number;
    }

    else cout << "Unable to open file";

    return FB;
}


/*
Description: Generates our sieving interval given a size and starting bound.
Params: (1) mpz_t number N (product of two primes);
        (2) int size_SSI (interval size);
        (3) mpz_t number T (interval starting bound starting bound).
Return: A pointer to the polynomial_element array. The polynomial_element
        struct holds an mpz_t type (we named "poly"). This is done
        due to weird memory constraints of mpz_t types
*/
polynomial_element * generate_sieving_interval(mpz_t N, int pes, mpz_t T){

    int size = pes;
    polynomial_element * SI = new polynomial_element[size];

    //Find smallest value of T such that T^2 - N >= 0
    mpz_t Tsq, res;
    mpz_init(Tsq);
    mpz_init(res);

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

void master_unpack_save(int* total_counter, int size_FB, int* need_more, ofstream& expo_matrix_file, ofstream& bit_matrix_file, ofstream& smooth_num_file, int* dead_processes){
  int new_relations, location, size, bit_val, packed_str_length;
  MPI_Status status;
  int** relations_storage;
  char* packed_smooth_nums_m;


  // STEP 1 REC
  MPI_Recv(&new_relations, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
  *total_counter += new_relations;
  location = status.MPI_SOURCE;
  size = new_relations * size_FB;
  relations_storage = alloc_2d_int(new_relations, size_FB);

  //STEP 2 REC
  MPI_Recv(&relations_storage[0][0], size, MPI_INT, location, 0, MPI_COMM_WORLD, &status);

  if (*need_more == 1){
    for (int i = 0; i < new_relations; i++){
      for (int j = 0; j < size_FB; j++){
        expo_matrix_file << relations_storage[i][j];
        bit_val = relations_storage[i][j] % 2;
        bit_matrix_file << bit_val;
      }
      expo_matrix_file << endl;
      bit_matrix_file << endl;
    }
  }

  // STEP 3
  MPI_Recv(&packed_str_length, 1, MPI_INT, location, 0, MPI_COMM_WORLD, &status);
  // cout << "got length" << endl;
  packed_smooth_nums_m = new char[packed_str_length];
  MPI_Recv(&packed_smooth_nums_m[0], packed_str_length, MPI_CHAR, location, 0, MPI_COMM_WORLD, &status);
  // cout << "got strs" << endl;

  if (*need_more == 1){
    for (int i = 0; i < packed_str_length; i++){
      if (packed_smooth_nums_m[i] != '|' && packed_smooth_nums_m[i] != '\0') {
        smooth_num_file << packed_smooth_nums_m[i];
      } else if (packed_smooth_nums_m[i] == '|') {
        smooth_num_file << endl;
      }
    }
  }

  //set need more to zero, once enough relations exist
  if (*total_counter >= size_FB + 10){
    *need_more = 0;
  }


  MPI_Send(need_more, 1, MPI_INT, location, 0, MPI_COMM_WORLD);

  // cout << location << endl;

  if (*need_more == 0){
    *dead_processes += 1;
    cout << "We are here" << endl;
  }
  cout << "We are not here" << endl;


  delete [] packed_smooth_nums_m;
  free(relations_storage[0]);
  free(relations_storage);

  // cout << location << endl;
}

void worker_sieves(int** power_storage, int* counter, int block_size, mpz_t N, mpz_t T, int size_FB, prime_element* FB, int* relations_amt, int rank, polynomial_element* SI){
  mpz_t p, a, b, r, idx, min1, min2;
  mpz_init(p);
  mpz_init(a);
  mpz_init(b);
  mpz_init(r);
  mpz_init(idx);
  mpz_init(min1);
  mpz_init(min2);
  int step;


  for (int i = 0; i < size_FB; i++){
      //convert p, a, b to mpz types
      mpz_set_ui(p, FB[i].p);
      mpz_set_ui(a, FB[i].a);
      mpz_set_ui(b, FB[i].b);

      unsigned long init1 = prime_find_min(block_size, a, p, min1, T, r, idx, rank);
      unsigned long init2 = prime_find_min(block_size, b, p, min2, T, r, idx, rank);
      step = mpz_get_ui (p);

      //go ahead and do all of the divisions
      if (init1 < block_size + 1){
          prime_divide(SI, power_storage, block_size, size_FB, init1, step, counter, i);
      }
      if (init2 < block_size + 1){
          prime_divide(SI, power_storage, block_size, size_FB, init2, step, counter, i);
      }

      if(*counter >= size_FB + 10){
        break;
      }
    }


    for (int j = 0; j < block_size; j++){
      if (power_storage[size_FB][j] == 1){
        *relations_amt += 1;
      }
    }
    int rels = *relations_amt;

}

//step where we repeatedly divide until we have the required number of relations
void sieving_step(prime_element *FB, mpz_t N, int fbs, int rank, MPI_Status status, int block_size, int num_proc){


  //intialize values and recompute value for T
  mpz_t T, T_hold, poly, Tsq, res;
  int size_FB = fbs;
  int continue_sieving = 1;
  int power, size;

  if (rank == 0){
    int need_more = 1;
    int dead_processes = 0;
    int total_counter = 0;

    ofstream smooth_num_file;
    smooth_num_file.open ("Smooth_Num.txt");
    ofstream expo_matrix_file;
    expo_matrix_file.open ("Expo_Matrix.txt");
    ofstream bit_matrix_file;
    bit_matrix_file.open("Bit_Matrix.txt");

    while (dead_processes < num_proc - 1){
      master_unpack_save(&total_counter, size_FB, &need_more, expo_matrix_file, bit_matrix_file, smooth_num_file, &dead_processes);

    }
    smooth_num_file.close();
    expo_matrix_file.close();
    bit_matrix_file.close();

 } else{
  mpz_init(Tsq);
  mpz_init(res);
  mpz_init(T_hold);
  mpz_init(poly);

   int block_offset = (rank - 1)* block_size;
   //Find smallest value of T such that T^2 - N >= 0
    mpz_init_set_ui(T, 1);
    mpz_root(T, N, 2); // T = sqrt(N)
    mpz_add_ui(T, T, 1); //Buffer T by one to ensure non negativity
    mpz_set(T_hold, T);
    mpz_add_ui(T, T, block_offset);
    mpz_set(T_hold, T);

    while (continue_sieving == 1){

        cout << "Sieve start" << endl;

        int** relations;
        int counter = 0;
        int relations_amt = 0;
        int** power_storage = new int*[size_FB+1];
        string* smooth_nums;

        for (int i = 0; i < size_FB+1; i++){
          power_storage[i] = new int[block_size];
          for (int j = 0; j< block_size; j++){
            power_storage[i][j] = 0;
          }
        }

        polynomial_element * SI = generate_sieving_interval(N, block_size, T);
        mpz_set(T, T_hold);
        polynomial_element * SI_SAVE = generate_sieving_interval(N, block_size, T);
        mpz_set(T, T_hold);

        // cout << "Time to sieve!" << endl;
        worker_sieves(power_storage, &counter, block_size, N, T, size_FB, FB, &relations_amt, rank, SI);


        relations = alloc_2d_int(relations_amt, size_FB);
        smooth_nums = new string[relations_amt];



        reduce_and_transpose(smooth_nums, relations, power_storage, block_size, size_FB, SI_SAVE);

        worker_pack_send(&relations_amt, size_FB, relations, smooth_nums);


        for (int i = 0; i < size_FB + 1; i++){
          delete[] power_storage[i];
        }
        delete[] power_storage;
        free(relations[0]);  // remove the pool
        free(relations);     // remove the pointers
        delete [] smooth_nums;

        int offset = block_size * (num_proc - 1);

        mpz_add_ui(T, T, offset);
        mpz_set(T_hold, T);
        MPI_Recv(&continue_sieving, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        cout << continue_sieving << endl;
      }
    }
  }

void worker_pack_send(int* relations_amt, int size_FB, int** relations, string* smooth_nums){
  string packed_smooth_nums;
  int rels = *relations_amt;
  int* packed_length = new int;
  char *packed;

  //STEP 1 SEND
  int ret = MPI_Send(&rels, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);


  //STEP2
  int size = (size_FB) * rels;
  ret = MPI_Send(&relations[0][0], size, MPI_INT, 0, 0, MPI_COMM_WORLD);


  //STEP3 SENDING SMOOTH NUMS
  packed_smooth_nums = pack(packed_length, smooth_nums, rels);
  packed = new char[*packed_length];
  strcpy(packed, packed_smooth_nums.c_str());
  ret = MPI_Send(packed_length, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
  // cout << "sent length" << endl;
  // cout << *packed_length << endl;
  ret = MPI_Send(&packed[0], *packed_length, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
  // cout << "sent strs" << endl;


  delete[] packed;
  delete packed_length;


}

int **alloc_2d_int(int rows, int cols) {
    int *data = (int *)malloc(rows*cols*sizeof(int));
    int **array= (int **)malloc(rows*sizeof(int*));
    for (int i=0; i<rows; i++)
        array[i] = &(data[cols*i]);
    return array;
}

//divide all the numbers by the prime when applicable
void prime_divide(polynomial_element* SI, int** power_storage, int block_size, int size_FB, int smallest, int prime, int* counter, int i){
    int power = 0;
    int step = prime;

    //iterate through all numbers which will be divisble by prime
    for (int j = smallest; j < block_size; j = j + step){

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


unsigned long prime_find_min(int block_size, mpz_t a, mpz_t p, mpz_t min, mpz_t T, mpz_t r, mpz_t idx, int rank){
     unsigned long temp = block_size+1;


  for (unsigned long j = 0; j < block_size; j++){
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

void reduce_and_transpose(string* smooth_nums, int** relations, int** power_storage, int block_size, int size_FB, polynomial_element* SI){
  string temp;

  int s_idx = 0;
  for (int j = 0; j < block_size; j++){
    if (power_storage[size_FB][j] == 1){
      temp = mpz_get_str(NULL, 10, SI[j].poly);
      smooth_nums[s_idx] = temp;
      for (int idx = 0; idx < size_FB; idx++){
        relations[s_idx][idx] = power_storage[idx][j];
      }
      s_idx++;
    }
  }
}

string pack(int* string_length, string* smooth_nums, int relations_amt){

  *string_length = 0;
  string packed_smooth_nums = "";
  for (int i = 0; i < relations_amt; i++){
    packed_smooth_nums = packed_smooth_nums + smooth_nums[i] + "|";
    *string_length = *string_length + smooth_nums[i].length() + 1;
  }
  *string_length = *string_length + 1;

  return packed_smooth_nums;

}
