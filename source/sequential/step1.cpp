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

#include "step1.h"

using namespace std;

int main(int argc, char *argv[]){

    int N = 69;

    double param, e;
    param = 1;
    e = exp (param);
    cout << e <<endl;

    double k = pow (e, sqrt(2*log(N) * log(log(N)) ) );

    int l = 2*k *log(2*k) + 1;

    cout << l << endl;

    vector<int> primes = getprimes(l, N);


    for (int i=0; i<primes.size(); i++){
        cout << primes[i] << endl;
    }

    return 0;

}


vector<int> getprimes(int l, int N){

    vector<int> primes;

    bool * truth = new bool[l+1];
    for (int i = 0; i < l+1; i++){
        truth[i] = true;
    }

    truth[0] = false;
    truth[1] = false;

    for (int i = 0; i*i < l; i++){
        if (truth[i]){
            for (int j=i*i; j<=N; j=j+i){
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


