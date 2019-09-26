#include <iostream>
#include <cmath>
#include <armadillo>
#include <fstream>
#include <string>
#include <iomanip>
#include <chrono>
#include <ctime>
#include "functions.cpp"
#include "get_matrices.cpp"

using namespace std;
using namespace arma;

// decleration of functions

double maxoffdiag(mat A, int & k, int & l, int n);
void rotate (mat &A, mat &R, int k, int l, int n);
void jacobi (mat A, mat R, int n);
double potential(double);
void output(string, double, double, int, vec& );




int main(int argc, char* argv[])
{



    // 2b

    // setting up the dimension, N
    int N = (atoi(argv[1]));

    // generating matrix A for 2b
    mat A = get_mat_A_2b(N);

    // setting up an empty R matrix
    mat R(N,N, fill::zeros);




    // checking analytic eigenpairs with armadillo
    // measuring the time it takes
    chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();

    vec eigenvalues;
    mat eigenvectors;
    eig_sym(eigenvalues, eigenvectors, A);

    chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
    chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double> >(t2 - t1);
    cout << "Armadillo used " << time_span.count() << " seconds.";
    cout << endl;




    // checking the time for the jacobi method for A in 2b
    chrono::high_resolution_clock::time_point t3 = chrono::high_resolution_clock::now();

    jacobi("2b",A,R,N);

    chrono::high_resolution_clock::time_point t4 = chrono::high_resolution_clock::now();
    chrono::duration<double> time_span2 = chrono::duration_cast<chrono::duration<double> >(t4 - t3);
    cout << "Jacobi used " << time_span2.count() << " seconds.";
    cout << endl;

    double diff_time = time_span2.count() - time_span.count();
    cout << "the differens in 2b is = " << diff_time << endl;






    // 2c tests
/*

    test1(1,1);
    cout << A_test << endl;

*/





    // 2d, Quantum dots in three dimensions, one electron


    int Dim = 100;      // dimension of the matrix
    int Rmin = 0.0;     // rho_min
    int Rmax = 8;       // rho_max
    mat R2(Dim,Dim, fill::zeros);
    mat ham = get_ham_2d(Dim, Rmax);    // Setting up the matrix for one electron

    // finding the analytical solution with armadillo
    vec Eigval(Dim);
    eig_sym(Eigval, ham);
    output("2d", Rmin, Rmax, Dim, Eigval);


    // finding the appoximate solution with jacobi for one electron
    jacobi("2d",ham,R2,Dim);

    // 2d end





    // 2e, Quantum dots in three dimensions, two electrons


    double omega = 5;
    mat ham_2e = get_ham_2e(Dim, omega, Rmax);

    // finding the analytical solution with armadillo
    vec Eigval_2e(Dim);
    eig_sym(Eigval_2e, ham_2e);
    output("2e", Rmin, Rmax, Dim, Eigval_2e);


    // finding the appoximate solution with jacobi
    jacobi("2e",ham_2e,R2,Dim);



    // 2e end


    return 0;
}








// functions


/*
// 2c function

void test1(){
    int t = 4;
    mat A_test;
    // filling in A
    for ( int i = 0; i < t; i++){
        A_test(i,i) = 2.0;
        if (i > 0)
        {
            A_test(i,i-1) = -1.0;
        }
        if (i <= t-2)
        {
            A_test(i,i+1) = -1.0;
        }
    }
    vec test_eigval(n);
    mat test_eigvec(n);
    eig_gen(test_eigval, test_eigvec, A_test);

}
*/
