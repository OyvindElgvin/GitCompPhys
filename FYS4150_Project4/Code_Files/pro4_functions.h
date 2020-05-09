#ifndef PRO4_FUNCTIONS_H
#define PRO4_FUNCTIONS_H

#include<iostream>
#include "armadillo"
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


//void changing_state(mt19937_64 generator, uword i, uword j, vec dE, vec P,int L, mat& spin_matrix, double& Energy, double& Mmoment,int& AC);
void changing_state(double inverse_period, uword i, uword j, vec dE, vec P,int L, mat& spin_matrix, double& Energy, double& Mmoment,int& AC);
vector<int> readvalues(string file);
void Ising_Func_Para(vec T,int L,int N,string file,string order,int test,int stepsize,string probability);
void Ising_Func_Para_e(vec T,int L,int N,string file,string order,int test,int stepsize,string probability);
void initialize(int n_spins,mat& spin_matrix, double& E, double& M,string order);
void Ising_Func_Para_mpi_e(vec T,int L,int N,string file,string order,int test,int stepsize,string probability);


void Task_b();
void Task_c();
void Task_d();
void Task_e();


#endif // PRO4_FUNCTIONS_H
