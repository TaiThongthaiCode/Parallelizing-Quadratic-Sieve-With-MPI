/*

Tai Thongthai and Tarang Saluja

"The philosophers have only interpreted the world in various ways. The point, however, is to change it."  - Pikachu;

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


#include "step3.h"


int main(int argc, char* arv[]){

  //Set value of N
  mpz_t N;
  mpz_init_set_ui(N, 10062757);

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
  //prime_element * FB = load(N, fbs);



  return 0;
}


//TODO

//(1)  Insert the matrix row by row from the output of step2.cpp

//(2) Use algorithm to get solutions over GF(2)

//(3) For each solution (1) Find the number T_1, T_2 ... T_k and the square root of the product of primes (2) Subtract take two numbers and take gcd with N

//Once solution is found, quit!
