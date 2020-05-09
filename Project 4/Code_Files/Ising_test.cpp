#include<iostream>
#include <armadillo>
#include <vector>
#include "time.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#define EPS 3.0e-14
#define MAXIT 10

#include <random>
#include <chrono>

using namespace std;
using namespace arma;
using namespace chrono;


void Ising_Func(int N, int L, int test);
vector<int> readvalues(string file);

/*
int Ising_Test(){

    int L = 2;
    int N = 100;

    Ising_Func(N,L,1);

    return 0;

} // end of test
*/
void test_dE(){
    unsigned seed = system_clock::now().time_since_epoch().count();
    mt19937_64 generator (seed);
    mat S = mat(5,5,fill::zeros);
    for(uword i=0;i<5;i++){
        for (uword j=0;j<5;j++){
            double r = generate_canonical< double, 128 > (generator);
            if (r < 0.5){S(i,j) = -1;}
            else{S(i,j) = 1;}
        }
    }
    cout << S << endl;
    double dE = 0;
    for(uword i=0;i<5;i++){
        for (uword j=0;j<5;j++){
            dE = 2*S(i,j)*(S(i,j-1) + S(j,j+1) + S(i-1,j) + S(i+1,j));
            cout << dE << endl;
        }
    }

    return;
}
