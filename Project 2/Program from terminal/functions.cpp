
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

void jacobi (string filename, mat A, mat R, int n){

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
    A_diag.save(filename + "A_eigvalues", arma_ascii);
    R.save(filename + "R", arma_ascii);


    return;
}


void output(string filename, double Rmin , double Rmax, int Dim, vec& d){
    int i;
    cout << "RESULTS:" << filename << endl;
    cout << setiosflags(ios::showpoint | ios::uppercase);
    cout <<"Rmin = " << setw(15) << setprecision(8) << Rmin << endl;
    cout <<"Rmax = " << setw(15) << setprecision(8) << Rmax << endl;
    cout <<"Number of steps = " << setw(4) << Dim << endl;
    cout << "Five lowest eigenvalues with armadillo:" << endl;
    for(i = 0; i < 4; i++) {
        cout << setw(15) << setprecision(8) << d[i] << endl;
    }
}
