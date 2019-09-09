#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <string>
#include <iomanip>
using namespace std;
using namespace arma;

ofstream ofile;

int main(int argc, char* argv[]){

    string filename;
    // checks the number of arguments
    if (argc <= 1){
        cout << "write filename, plus an integer n for the nxn matrise" << endl;
        exit(1);
    }

    else{
        filename = string (argv[2]);
        filename = "n=" + filename;

    }

    string fileout = filename;

    ofile.open(fileout);

    int n = atoi(argv[2]);

    // vector lengths

    double* x = new double[n+2];
    double* a = new double[n+2];
    double* c = new double[n+2];
    double* b = new double[n+2];
    double* f = new double[n+2];
    double* v_general = new double[n+2];
    double* exact = new double[n+2];
    double* b_tilde = new double[n+2];
    double* relativeError_general = new double[n+1];


    double h = 1 / double (n+1); // step size

    // grid points
    for (int i = 0; i < n+1; i++){
        x[i] = i*h;
    }

    // Defining the endpoints for x
    x[0] = 0;
    x[n+1] = 1;


    // Exact solution
    for (int i = 0; i < n+1; i++){
        exact[i] = 1 - (1-exp(-10))*x[i] - exp(-10*x[i]);
    }

    // function
    for (int i = 0; i < n+1; i++){
        f[i] = 100*exp(-10*x[i]);
    }

    // vector elements of b
    for (int i = 0; i < n+1; i++){
        b[i] = 2;
    }

    // vector elements of a, c
    for (int i = 0; i < n+1; i++){
        a[i] = -1;
        c[i] = -1;
    }

    // vector elements for b_tilde
    for (int i = 0; i < n+1; i++){
        b_tilde[i] = h*h*f[i];
    }



    // Starts clock for general time
    clock_t start, finish;
    start = clock();

    //forward general substitution
    for (int i = 2; i < n+1; i++){
        b[i] = b[i] - a[i-1]*c[i-1] / b[i-1];
        b_tilde[i] = b_tilde[i] - b_tilde[i-1]*c[i-1] / b[i-1];
    }

    //backwards general substitution
    v_general[n] = b_tilde[n] / b[n];

    for (int i = n-1; i > 0; i--){
        v_general[i] = (b_tilde[i] - c[i]*v_general[i+1]) / b[i];
    }

    finish = clock();
        double timeused = double (finish - start)/(CLOCKS_PER_SEC );
        cout << setiosflags(ios::showpoint | ios::uppercase);
        cout << setprecision(10) << setw(20) << "Time used for the general algorithm = " << timeused  << endl;




    // specialized algorithm



    double* b_spes = new double[n+2]; //  spes diagonalen
    double* b_tilde_spes = new double[n+2]; // spes solution
    double* v_spes = new double[n+2]; // spes toplitz



    // flops outside the clock
    for (int i = 1; i < n+1; i++){
        b_spes[i] = (i+1) / double (i);

    }

    // vector elements for b_tilde_spes
    for (int i = 1; i < n+1; i++){
        b_tilde_spes[i] = h*h*f[i];
    }

    start = clock();

    // Specialised Gauss
    // forward specialized substitution
    for (int i = 2; i < n+1; i++){
        b_tilde_spes[i] = b_tilde_spes[i] + b_tilde_spes[i-1] / (b_spes[i-1]);
    }

    // backwards specialized substitution
    v_spes[n] = b_tilde_spes[n] / b_spes[n];
    for (int i = n-1; i > 0; i--){
        v_spes[i] = (b_tilde_spes[i] + v_spes[i+1]) / (b_spes[i]);
    }


    finish = clock();
        double timeused2 = double (finish - start)/(CLOCKS_PER_SEC );
        cout << setiosflags(ios::showpoint | ios::uppercase);
        cout << setprecision(10) << setw(20) << "Time used for the specialized algorithm = " << timeused2  << endl;



    // LU decomposition

    double* y = new double[n+2];
    double* b_LU = new double[n+2];
    double* v_LU = new double[n+2];

    mat A;
    A = zeros(n+2,n+2);

    for (int i = 0; i < n+2; i++){
        A(i,i) = 2;
        if (i > 0)
        {
            A(i,i-1) = -1;
        }
        if (i+1 <= n+1)
        {
            A(i,i+1) = -1;
        }
    }

    mat L, U;

    start = clock();

    lu(L,U,A);

    //y = solve( L * b_tilde );
    //v_LU = solve( U * y );


    finish = clock();
        double timeused3 = double (finish - start)/(CLOCKS_PER_SEC );
        cout << setiosflags(ios::showpoint | ios::uppercase);
        cout << setprecision(10) << setw(20) << "Time used for the LU decomposition algorithm = " << timeused3  << endl;


    // calculating the relative errors
    for (int i = 1; i < n; i++){
        relativeError_general[i] = log10(fabs((v_general[i] - exact[i]) / exact[i]));
    }

    //relativeError_general[n+2] = 1;
    // writes values of (x, v_general, v_spes, exact, relativeError_general) to file
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << "          x,              v_generalized,  v_specialized,  exact,          relative error" << endl;
    for (int i = 0; i < n+2; i++){
        ofile << setw(20) << setprecision(8)
              << x[i] << ",     "
              << v_general[i] << ",     "
              << v_spes[i] << ",     "
              << exact[i] << ",     "
              << relativeError_general[i] << endl;
    }


    delete [] x;
    delete [] a;
    delete [] b;
    delete [] c;
    delete [] f;
    delete [] b_tilde;
    delete [] b_tilde_spes;
    delete [] v_general;
    delete [] v_spes;
    delete [] exact;
    delete [] relativeError_general;
    delete [] y;
    delete [] b_LU;
    delete [] v_LU;


    ofile.close();

    return 0;
}
