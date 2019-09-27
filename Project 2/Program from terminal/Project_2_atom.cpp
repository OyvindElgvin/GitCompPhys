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
//void rotate (mat &A, mat &R, int k, int l, int n);
void jacobi (mat A, mat R, int n);

void output(string, double, double, int, vec& );
mat test_maxoff_mat(int n);




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

    // makes a identity matrix with one element with value 3
    // and tests if the maxoffdiag function picks it out.
    int k, l;
    mat test1 = test_maxoff_mat(5);    // matrix with 3 as highest element value
    //cout << test1 << endl;
    double testen = maxoffdiag(test1, k, l, 5);
    cout << "Test functions" << endl
         << "Test 1" << endl
         << "if the test prints out 3, it works " << endl
         << "if not, it doesn't" << endl
         << testen << endl;


    // checks if the orthogonal principal is preserved
    mat vv = randu<mat>(5,5);
    mat v = orth(vv);

    mat vT = trans(v);
    mat delta = vT * v;
    cout << "delta = " << endl
         << delta << endl;


    mat u = randu<mat>(5,5);
    mat uT = trans(u);

    mat w = u * v;
    mat wT = vT * uT;
    mat wTw = vT * uT * u * v;



    //cout << "Test 2"<< endl
    //     << "delta = " << endl
    //     << delta << endl
    //     << "wTw = " << endl
    //     << wTw << endl;





    // 2d, Quantum dots in three dimensions, one electron


    int Dim = 100;      // dimension of the matrix
    int Rmin = 0.0;     // rho_min
    int Rmax = 8;       // rho_max
    mat R2(Dim,Dim, fill::zeros);
    mat ham = get_ham_2d(Dim, Rmax);    // Setting up the matrix for one electron

    // finding the analytical solution with armadillo
    vec Eigval(Dim);
    eig_sym(Eigval, ham);
    output("2d", Rmin, Rmax, Dim, Eigval); // prints out eigenvalues


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
