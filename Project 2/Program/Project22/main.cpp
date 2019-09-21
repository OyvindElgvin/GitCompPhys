#include <iostream>
#include <cmath>
#include <armadillo>
#include <fstream>
#include <string>
#include <iomanip>
//#include "jacobi.h"


using namespace std;
using namespace arma;

// decleration of functions

double maxoffdiag(double** A, int* k, int* l, int n);
void rotate (double** A, double** R, int k, int l, int n);
void Jacobi (double** A, double** R, int n);



int main(int argc, char* argv[] )
{
    // setting ut the dimension, N
    unsigned int N = unsigned (atoi(argv[2]));


    // setting up an A nxn random matrix to diagonalize
    mat B = randn<mat>(N,N);
    cout << "matrix A: " << endl << B << endl;
    mat A_sym = B.t()*B;
    cout << "matrix A_sym: " << endl << A_sym << endl;



    // analytically diagonalizing with armadillo
    mat A_diagonalized;
    A_diagonalized = diagmat(A_sym,-1) + diagmat(A_sym,0) + diagmat(A_sym,1);
    cout << "arma diag A_sym: " << endl << A_diagonalized << endl;



    // setting up an empty R matrix
    mat R(N,N, fill::zeros);
    cout << "matrix R: " << endl << R << endl;



    // generating the eigenvalues in armadillo
    vec eigenvalues;
    mat eigenvectors;
    eig_sym(eigenvalues, eigenvectors, A_sym);
    cout << eigenvalues << endl;



    mat C;
    int k,l;
    B = maxoffdiag(  noe her   ,k,l,N);

    return 0;
}








// Function for finding the max value of all elements
double maxoffdiag(double** A, int* k, int* l, int n){
    double max = 0.0;
    for (int i = 0; i < n; i++){
        for (int j = i+1; j < n; j++){
            if (fabs(A[i][j]) > max){
                max = fabs(A[i][j]);
                *l = i;
                *k = j;
            }
        }
    }
    return max;
}

// Function for the rotation
void rotate (double** A, double** R, int k, int l, int n){
    double s, c;
    if (A[k][l] != 0.0){
        double t, tau;
        tau = (A[l][l] - A[k][k]) / (2*A[k][l]);
        if (tau > 0){
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
    a_kk = A[k][k];
    a_ll = A[l][l];
    A[k][k] = a_kk*c*c - 2.0*A[k][l]*c*s + a_ll*s*s;
    A[l][l] = a_ll*c*c + 2.0*A[k][l]*c*s + a_kk*s*s;
    A[k][l] = 0.0;
    A[l][k] = 0.0;
    for (int i = 0; i < n; i++){
        if (i != k && i != l){
            a_ik = A[i][k];
            a_il = A[i][l];
            A[i][k] = a_ik*c - a_il*s;
            A[k][i] = A[i][k];
            A[i][l] = a_il*c + a_ik*s;
            A[l][i] = A[i][l];
        }
        r_ik = R[i][k];
        r_il = R[i][l];
        R[i][k] = r_ik*c - r_il*s;
        R[i][l] = r_il*c - r_ik*s;
    }
    return;
}

// making the eigenvector matrix R
void Jacobi (double** A, double** R, int n){

    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            if (i == j){
                R[i][j] = 1.0;
            }
            else {
                R[i][j] = 0.0;
            }
        }
    }

    // Loop for rotating the matrix
    int k, l;
    double tol = 1e-8;
    int iterations = 0.0;
    double max_nr_itera = double (n) * double (n) * double (n);
    double max_offdiag = maxoffdiag (A, &k, &l, n);

    while ( max_offdiag > tol && double (iterations) < max_nr_itera){
        max_offdiag = maxoffdiag (A, &k, &l, n);
        rotate (A, R, k, l, n);
        iterations++;
    }
    cout << " # iterations = " << iterations << endl;
    return;
}


