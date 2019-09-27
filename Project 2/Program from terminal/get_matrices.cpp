
#include <iostream>
#include <cmath>
#include <armadillo>
#include <fstream>
#include <string>
#include <iomanip>
#include <chrono>
#include <ctime>



using namespace std;
using namespace arma;


// 2b

mat get_mat_A_2b(int N){

    // setting up the tridiagonal matrix A
    mat A;
    A = zeros(N,N);
    double h = 1.0 / N;


    // filling in A
    for ( int i = 0; i < N; i++){
    	A(i,i) = 2.0/(h*h);
    	if (i > 0)
    	{
    		A(i,i-1) = -1.0/(h*h);
    	}
    	if (i <= N-2)
    	{
    		A(i,i+1) = -1.0/(h*h);
    	}
    }
	return A;
}

// 2d

mat get_ham_2d(int Dim, double Rmax){


    int i, j, Orbital;
    double Rmin, Step, DiagConst, NondiagConst, OrbitalFactor;
    Rmin = 0.0;
    //Rmax = 8.0;
    //Orbital = 0;
    //Dim = 40;

    mat Hamiltonian = zeros<mat>(Dim,Dim);
    Step = Rmax / Dim;
    DiagConst = 2.0 / (Step*Step);
    NondiagConst = -1.0 / (Step*Step);
    //OrbitalFactor = Orbital * (Orbital + 1.0);

    vec rho(Dim);
    vec potential(Dim);
    for (int i = 0; i< Dim; i++){
        rho(i) = Rmin + (i+1)*Step;     // take boundaries into considerations
        potential(i) = rho(i) * rho(i);     // V_i
        //w(i) = potential(r(i)); //+ OrbitalFactor/(r(i) * r(i));
    }


    Hamiltonian(0,0) = DiagConst + potential(0);
    Hamiltonian(0,1) = NondiagConst;
    for (int i = 1; i < Dim-1; i++){
        Hamiltonian(i,i-1) = NondiagConst;
        Hamiltonian(i,i) = DiagConst + potential(i);
        Hamiltonian(i,i+1) = NondiagConst;
    }
    Hamiltonian(Dim-1,Dim-2) = NondiagConst;
    Hamiltonian(Dim-1,Dim-1) = DiagConst + potential(Dim-1);

	return Hamiltonian;
}

// 2e

mat get_ham_2e(int Dim, double omega, double Rmax){

    int i, j, Orbital;
    double Rmin, Step, DiagConst, NondiagConst, OrbitalFactor;
    Rmin = 0.0;
    //Rmax = 8.0;
    //omega = 0.01;
    //Dim = 40;
    mat R2(Dim,Dim, fill::zeros);
    mat Hamiltonian = zeros<mat>(Dim,Dim);
    Step = Rmax / Dim;
    DiagConst = 2.0 / (Step*Step);
    NondiagConst = -1.0 / (Step*Step);

    vec rho(Dim);
    vec potential(Dim);
    for (int i = 0; i< Dim; i++){
        rho(i) = Rmin + (i+1)*Step;     // take boundaries into considerations
        potential(i) = omega*omega*rho(i)*rho(i) + 1/rho(i);
        //w(i) = potential(r(i)); //+ OrbitalFactor/(r(i) * r(i));
    }

    Hamiltonian(0,0) = DiagConst + potential(0);
    Hamiltonian(0,1) = NondiagConst;
    for (int i = 1; i < Dim-1; i++){
        Hamiltonian(i,i-1) = NondiagConst;
        Hamiltonian(i,i) = DiagConst + potential(i);
        Hamiltonian(i,i+1) = NondiagConst;
    }
    Hamiltonian(Dim-1,Dim-2) = NondiagConst;
    Hamiltonian(Dim-1,Dim-1) = DiagConst + potential(Dim-1);

	return Hamiltonian;
}



mat test_maxoff_mat(int n){
    mat test_maxoff_mat(n,n,fill::eye);
    test_maxoff_mat(1,3) = 3;
	return test_maxoff_mat;

}
