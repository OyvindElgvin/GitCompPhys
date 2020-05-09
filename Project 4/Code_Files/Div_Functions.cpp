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

#include "pro4_functions.h"

using namespace std;
using namespace arma;
using namespace chrono;

// Read a Coloumn Vector from a file
vector<int> readvalues(string file){
    string line;
    vector<int> N_values;
    ifstream values;
    values.open(file,ios::in);
    while (getline(values,line)){
        if (atoi(line.c_str()) != 0){
        N_values.push_back(atoi(line.c_str()));
        }
       else{}
   }
    values.close();
    return N_values;
    }


// inline function for periodic boundary conditions
inline int periodic(int i, int limit, int add) {
    return (i+limit+add) % (limit);
}



// function to initialise energy, spin matrix and magnetization
void initialize(int n_spins,mat& spin_matrix, double& E, double& M,string order){
    // Set up spin matrix
    if (order == "order"){
        spin_matrix.ones();
    }
    else{
        spin_matrix.randu();
        for (uword i = 0; i<n_spins; i++){
            for (uword j = 0; j<n_spins; j++){
                if(spin_matrix(i,j) < 0.5){spin_matrix(i,j) = -1.;}
                else{spin_matrix(i,j) = 1.;}
            }
        }
    }

    // Set up inital magnetic moment and energy
    for (uword i = 0; i<n_spins; i++){
        for (uword j = 0; j<n_spins; j++){
            M += spin_matrix(i,j);
            E -= spin_matrix(i,j)*
                    (spin_matrix(periodic(i,n_spins,-1),j) +
                     spin_matrix(i,periodic(j,n_spins,-1)));
        }
    }
} // End initialize

// function for changing energy state
//void changing_state(mt19937_64 generator, uword i, uword j, vec dE, vec P,int L,mat& spin_matrix, double& Energy, double& Mmoment,int& AC){
void changing_state(double inverse_period, uword i, uword j, vec dE, vec P,int L,mat& spin_matrix, double& Energy, double& Mmoment,int& AC){


    // Calculating energy difference for flipped spin
    double delta_E = 2*spin_matrix(i,j) * (spin_matrix(i,periodic(j,L,-1)) + spin_matrix(i,periodic(j,L,1)) + spin_matrix(periodic(i,L,-1),j) + spin_matrix(periodic(i,L,1),j));

    // Find probability exp(-b dE) from precalculated array
    uvec index = find(dE == delta_E);
    double Prob = P(index(0));

    // Sample rule: Check if new state is accepted
    //double r = generate_canonical< double, 128 > (generator);
    double r = double(rand())*inverse_period;
    //cout << dE(index(0)) << "     " << Prob << "      " << r << endl;
    if (r < Prob){
        spin_matrix(i,j) = -spin_matrix(i,j);
        Energy = Energy + delta_E;
        Mmoment = Mmoment + 2*spin_matrix(i,j);
        AC += 1;
    }

} // end function changing_state


vec T_vector(double min, double max,double step_length){
    int num_steps = static_cast<int>((max-min)/step_length);
    vec T = vec(num_steps,fill::zeros);
    T(0) = min;
    for (uword i = 1;i<num_steps;i++){
        T(i) = T(i-1) + step_length;
    }
    return T;
} // end function T_vector


void Task_b(){
    vector<int> Nvalues = readvalues("Pro4b_Nvalues.txt");
    int L = 2;
    vec T = vec("1.0");
    int stepsize = 100000;

    for (int i = 0;i<Nvalues.size();i++){
    Ising_Func_Para(T,L,Nvalues[i],"Results_4b","order",1,stepsize,"no probability");
    }
    return;
} // end Task_b

void Task_c(){
    vector<int> Nvalues = readvalues("Pro4c_Nvalues.txt");
    int L = 20;
    int stepsize = 10000;
    vec T = vec("1.0 2.4");

    for (int i = 0;i<Nvalues.size();i++){
        Ising_Func_Para(T,L,Nvalues[i],"Results_4c_order","order",0,stepsize,"no probability");
        Ising_Func_Para(T,L,Nvalues[i],"Results_4c_random","random",0,stepsize,"no probability");
    }
    return;
} // end Task_c


void Task_d(){
    vector<int> Nvalues = readvalues("Pro4d_Nvalues.txt");
    int L = 20;
    vec T = vec("1.0 2.4");

    for (int i = 0;i<Nvalues.size();i++){
        int stepsize = pow(10,log10(Nvalues[i])-1);
        Ising_Func_Para(T,L,Nvalues[i],"Results_4d","random",2,stepsize,"probability");
    }
    return;
} // end Task_d

void Task_e(){
    vector<int> Nvalues = readvalues("Pro4e_Nvalues.txt");
    vec T = T_vector(2.2,2.3,0.0005);
    vec L = vec("40 60 80 100");
    for (int i = 0;i<Nvalues.size();i++){
        int stepsize = pow(10,log10(Nvalues[i])-1);
        for (uword j= 0;j<L.n_elem;j++){
            cout << "L = " << L(j) << endl;
            Ising_Func_Para_mpi_e(T,L(j),Nvalues[i],"Results_4e_2","order",0,stepsize,"no probability");
        }//end of L loop
    }//end of N loop
    return;
}// end Task_e

