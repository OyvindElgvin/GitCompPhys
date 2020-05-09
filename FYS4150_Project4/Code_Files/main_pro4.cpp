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
#include <mpi.h>

#include "pro4_functions.h"

using namespace std;
using namespace arma;


//void Ising_Func(vec T, int L, int N, int test);
//int Ising_Test();
void test_dE();
vector<int> readvalues(string file);


int main(){

    //void test_dE();

    // Running test for L = 2 and N = 100, T = 1


    //vec T = vec("1.0");


    //Ising_Func_Para_e(T,2,100000000,"Results_4b","order",1,10000000,"no probability");

    MPI_Init(NULL, NULL);
    // Starts clock for general time
    clock_t start, finish;
    start = clock();

    Task_e();


    finish = clock();
    double timeused = double (finish - start)/(CLOCKS_PER_SEC );
    timeused = timeused / 3600;
    cout << setiosflags(ios::showpoint | ios::uppercase);
    cout << setprecision(10) << setw(20) << "Hours used for the total calculation = " << timeused << endl;
    MPI_Finalize ();

} // end of main
