#include <iostream>
#include <cmath>
#include <armadillo>
#include <fstream>
#include <string>
#include <iomanip>
#include <chrono>
#include <ctime>
//#include "jacobi.h"


using namespace std;
using namespace arma;

// decleration of functions

double maxoffdiag(mat A, int & k, int & l, int n);
void rotate (mat &A, mat &R, int k, int l, int n);
void jacobi (mat A, mat R, int n);
double potential(double);
void output(double, double, int, vec& );



int main(int argc, char* argv[])
{
    // setting ut the dimension, N
    int N = (atoi(argv[1]));

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

    // setting up an empty R matrix
    mat R(N,N, fill::zeros);



    // checking analytic eigenpairs with armadillo
    chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();

    vec eigenvalues;
    mat eigenvectors;
    eig_sym(eigenvalues, eigenvectors, A);
    //cout << "Eigenvalues for A = " << endl << eigenvalues << endl
    //    << "eigenvector for A = "<< endl
    //    << eigenvectors << endl;

    chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
    chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double> >(t2 - t1);
    cout << "Armadillo used " << time_span.count() << " seconds.";
    cout << endl;




    // checking the time for the jacobi method
    chrono::high_resolution_clock::time_point t3 = chrono::high_resolution_clock::now();

    jacobi(A,R,N);

    chrono::high_resolution_clock::time_point t4 = chrono::high_resolution_clock::now();
    chrono::duration<double> time_span2 = chrono::duration_cast<chrono::duration<double> >(t4 - t3);
    cout << "Jacobi used " << time_span2.count() << " seconds.";
    cout << endl;

    double diff_time = time_span2.count() - time_span.count();
    cout << "the differens is = " << diff_time << endl;




    // 2d, Quantum dots in three dimensions, one electron

    int i, j, Dim, Orbital;
    double Rmin, Rmax, Step, DiagConst, NondiagConst, OrbitalFactor;
    Rmin = 0.0;
    Rmax = 10.0;
    Orbital = 0;
    Dim = 400;
    mat Hamiltonian = zeros<mat>(Dim,Dim);
    Step = Rmax / Dim;
    DiagConst = 2.0 / (Step*Step);
    NondiagConst = -1.0 / (Step*Step);
    OrbitalFactor = Orbital * (Orbital + 1.0);

    vec r(Dim);
    vec w(Dim);
    for (int i = 0; i< Dim; i++){
        r(i) = Rmin + (i+1)*Step;
        w(i) = potential(r(i)) + OrbitalFactor/(r(i) * r(i));
    }





    Hamiltonian(0,0) = DiagConst + w(0);
    Hamiltonian(0,1) = NondiagConst;
    for (int i = 1; i < Dim-1; i++){
        Hamiltonian(i,i-1) = NondiagConst;
        Hamiltonian(i,i) = DiagConst + w(i);
        Hamiltonian(i,i+1) = NondiagConst;
    }

    Hamiltonian(Dim-1,Dim-2) = NondiagConst;
    Hamiltonian(Dim-1,Dim-1) = DiagConst + w(Dim-1);

    vec Eigval(Dim);
    eig_sym(Eigval, Hamiltonian);

    output(Rmin , Rmax, Dim, Eigval);

    return 0;
}





// functions



// function for finding the max value of all elements
double maxoffdiag(mat A, int & k, int & l, int N){
    double max = 0.0;
    for (int i = 0; i < N; ++i){
        for (int j = i+1; j < N; ++j){
            if (fabs(A(i,j)) > max){
                max = fabs(A(i,j));
                l = i;
                k = j;
            }
        }
    }
    return max;
}



// Function for the rotation

void rotate (mat & A, mat & R, int k, int l, int n){
    double s, c;
    if (A(k,l) != 0.0){
        double t, tau;
        tau = (A(l,l) - A(k,k)) / (2*A(k,l));
        if (tau >= 0){
            t = 1.0/(tau + sqrt(1.0 + tau*tau));
        }
        else {
            t = -1.0/(-tau + sqrt(1.0 + tau*tau));
        }
        c = 1/sqrt(1 + t*t);
        s = c*t;
    }
    else {
        c = 1.0;
        s = 0.0;
    }

    // replacing the k and l elements in the matrix
    double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
    a_kk = A(k,k);
    a_ll = A(l,l);
    A(k,k) = a_kk*c*c - 2.0*A(k,l)*c*s + a_ll*s*s;
    A(l,l) = a_kk*s*s + 2.0*A(k,l)*c*s + a_ll*c*c;
    A(k,l) = 0.0;
    A(l,k) = 0.0;
    for (int i = 0; i < n; i++){
        if (i != k && i != l){
            a_ik = A(i,k);
            a_il = A(i,l);
            A(i,k) = a_ik*c - a_il*s;
            A(k,i) = A(i,k);
            A(i,l) = a_il*c + a_ik*s;
            A(l,i) = A(i,l);
        }
        r_ik = R(i,k);
        r_il = R(i,l);
        R(i,k) = r_ik*c - r_il*s;
        R(i,l) = r_il*c + r_ik*s;
    }
    return;
}



// making the eigenvector matrix R

void jacobi (mat A, mat R, int n){

    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            if (i == j){
                R(i,j) = 1.0;
            }
            else {
                R(i,j) = 0.0;
            }
        }
    }

    // Loop for rotating the matrix
    int k, l;
    double tol = 1e-8;
    int iterations = 0.0;
    double max_nr_itera = double (n) * double (n) * double (n);
    double max_offdiag = maxoffdiag (A, k, l, n);


    while ( max_offdiag > tol && double (iterations) < max_nr_itera){
        max_offdiag = maxoffdiag (A, k, l, n);
        rotate (A, R, k, l, n);
        iterations++;
    }
    //cout << "# iterations = " << iterations << endl;
    //cout << "A = " << endl << A << endl;
    //cout << "R = " << endl << R << endl;

    vec A_diag = A.diag();
    A_diag.save("2b_A", arma_ascii);
    R.save("2b_R", arma_ascii);


    return;
}

// funtion for the potential
double potential(double x){
    return x*x;
}


void output(double Rmin , double Rmax, int Dim, vec& d)
{
    int i;
    cout << "RESULTS:" << endl;
    cout << setiosflags(ios::showpoint | ios::uppercase);
    cout <<"Rmin = " << setw(15) << setprecision(8) << Rmin << endl;
    cout <<"Rmax = " << setw(15) << setprecision(8) << Rmax << endl;
    cout <<"Number of steps = " << setw(4) << Dim << endl;
    cout << "Five lowest eigenvalues:" << endl;
    for(i = 0; i < 5; i++) {
        cout << setw(15) << setprecision(8) << d[i] << endl;
    }
}
