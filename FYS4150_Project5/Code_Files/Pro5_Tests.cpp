#include <iostream>
#include <armadillo>
#include <omp.h>
#include "time.h"
#include <random>
#include <chrono>

#include "Pro5_Functions.h"

using namespace std;
using namespace arma;

void test_sampling(){
    cout << "Running sampling index test" << endl;
    cout << "----------------------" << endl;
    vec M = vec(10,fill::randu);
    mat c = mat(M.n_elem,M.n_elem,fill::zeros);
    int C = 100;
    int fail = 0;
    for (int i = 0;i<100;i++){
        vector<int> index = Sampling_Rule(M,c,C,0,2);
        if (index[0] < 0 or index[0] >= M.n_elem or index[1] < 0 or index[1] >= M.n_elem or index[0] == index[1]){
            cout << "INDEX ERROR: \n ";
            cout << index[0] << "   " << index[1] << endl;
            fail = 1;
        }
    }
    if (fail == 0){
        cout << "Test passed, no index errors" << endl;
    }
    cout << "matrix for #transactions: \n ----------------------" << endl;
    cout << c << endl << accu(c)/2.0 << endl;
    cout << "----------------------" << endl;
} //end test_sampling

void test_stability(){
    cout << "Running stability test with no savings" << endl;
    cout << "----------------------" << endl;

    mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count());
    int N = 10;
    int m0 = 100;
    vec M1;
    M1 = m0*vec(N,fill::ones);
    vec M2 = M1;
    mat c = mat(M1.n_elem,M1.n_elem,fill::zeros);
    int C;
    int counter = 0;
    while (counter < 10){

    vector<int> index = Sampling_Rule(M2,c,C,0,2);
    transaction(index[0],index[1],0,M2);

    counter += 1;
    }
    cout << M1 << endl << endl << M2 << endl;

    if (sum(M1) == sum(M2)){
        cout << "Total money is conserved, and equal to " << sum(M1) << endl;
    }
    else{
        cout << "TOTAL MONEY IS NOT CONSERVED" << endl << "BEFORE: " << sum(M1) << endl << "AFTER: " << sum(M2) << endl;
    }
    cout << "----------------------" << endl;
}// end stability test

void test_stability_savings(){
    cout << "Running stability test with savings" << endl;
    cout << "----------------------" << endl;

    mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count());
    int N = 10;
    int m0 = 100;
    vec M1;
    M1 = m0*vec(N,fill::ones);
    vec M2 = M1;


    for (int i = 0; i < 10;i++){
        int i1 = static_cast<int>( (generate_canonical< double, 128 > (rng))*static_cast<double>(N/2) );
        int i2 = static_cast<int>( (generate_canonical< double, 128 > (rng))*static_cast<double>(N/2) + N/2 );
        transaction(i1,i2,0.5,M2);
    }
    cout << M1 << endl << endl << M2 << endl;

    if (sum(M1) == sum(M2)){
        cout << "Total money is conserved, and equal to " << sum(M1) << endl;
    }
    else{
        cout << "TOTAL MONEY IS NOT CONSERVED" << endl << "BEFORE: " << sum(M1) << endl << "AFTER: " << sum(M2) << endl;
    }
    cout << "----------------------" << endl;
} //end stability_savings test

void test_stability_taxes(){
    cout << "Running stability test with taxes" << endl;
    cout << "----------------------" << endl;

    mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count());
    int N = 10;
    int m0 = 100;
    vec M1;
    M1 = m0*vec(N,fill::ones);
    vec M2 = M1;


    for (int i = 0; i < 10;i++){
        int i1 = static_cast<int>( (generate_canonical< double, 128 > (rng))*static_cast<double>(N/2) );
        int i2 = static_cast<int>( (generate_canonical< double, 128 > (rng))*static_cast<double>(N/2) + N/2 );
        transaction(i1,i2,0,M2);
        wealth_tax(M2,0.25);
    }
    cout << M1 << endl << endl << M2 << endl;

    if (sum(M1) == sum(M2)){
        cout << "Total money is conserved, and equal to " << sum(M1) << endl;
    }
    else{
        cout << "TOTAL MONEY IS NOT CONSERVED" << endl << "BEFORE: " << sum(M1) << endl << "AFTER: " << sum(M2) << endl;
    }
    cout << "----------------------" << endl;
} //end stability_taxes test

/* This test only runs when Sampling_Rule is configured to accept double& P as argument
void test_P(){
    int N = 1000;
    vec M = vec(N,fill::ones);

    mat c = mat(M.n_elem,M.n_elem,fill::zeros);
    double P = 0;

    // Start MC loop
    for (int j = 0; j<100000;j++){


        vector<int> index = Sampling_Rule(M,c,P,2,2);

        transaction(index[0],index[1],0,M);

        if (j % 100 == 0){
            cout << P << endl;
        }
    }
} //end test P
*/

void test_gamma(){
    cout << "Running gamma test" << endl;
    cout << "----------------------" << endl;

    mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count());
    int N = 1000;
    int m0 = 1;
    vec A = vec(N,fill::zeros);
    mat B;
    vec M;
    double Sum = 0;

    // Start Gamma loop
    vec G = vec("0.0 1.0 2.0 3.0 4.0");

    for (int K = 0;K<5;K++){
        double g = G(K);

        for (int aa = 0;aa<N;aa++){
            A(aa) = N - 1;
        }

    M = m0*vec(N,fill::ones);
    mat c = mat(M.n_elem,M.n_elem,fill::zeros);
    mat c_ = c;
    // Start MC loop
    for (int j = 0; j<1000000;j++){

     // if (j%100 ==0){cout << j << endl; cout << c << endl << A << endl << "---------" << endl;}

    vector<int> index = Sampling_Rule_scaled(M,c,A,B,0.0,g);


    transaction(index[0],index[1],0.0,M);

    }


    vec c0 = c.col(0);
    c0 = sort(c0);
    int integer = 0;
    int I = 0;
    vector<int> counts = vector<int>{0};
    int i = 0;
    while(i < c0.n_elem){
        if ( c0(i) == integer){
            counts[I] += 1;
            i += 1;
        }
        else{
            cout << integer << " - " << counts[I] << endl;
            counts.push_back(0);
            I += 1;
            integer += 1;
        }
    }
    cout << "-------------" << endl;
    //cout << c0 << endl;
        for(int j = 0; j<M.n_elem;j++){
            if(0 != j){
            int A = 2;//cout << M(0) << "   " << c(0,j) << "   " << M(j) << "   Diff    " << abs(M(0)-M(j)) << endl;
            }
        }

    } //end gamma loop
    cout << "----------------" << "\n \n";
} //end test_gamma


void test_alpha(){
    cout << "Running alpha test" << endl;
    cout << "----------------------" << endl;

    mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count());
    int N = 10;
    int m0 = 1;
    vec A = vec(N,fill::zeros);
    vec B = vec(N,fill::zeros);
    vec M;
    double Sum = 0;

    // Start alpha loop
    vec alphas = vec("2.0");// 1.0 1.5 2.0");

    for (int K = 0;K<1;K++){
        double a = alphas(K);


        M = m0*vec(N,fill::ones);
        mat c = mat(M.n_elem,M.n_elem,fill::zeros);
        mat c_ = c;
        mat D = mat(N,N,fill::zeros);
        // Start MC loop
        for (int j = 0; j<10;j++){
            cout << endl << "-------------------" << endl;
            //if (j%100 ==0){cout << j << endl;} //cout << c << endl << A << endl << "---------" << endl;}
            for (int b=0;b<N;b++) {cout << M(b) << "  ";} cout << endl;
            vector<int> index = Sampling_Rule_scaled(M,c,A,D,a);

            transaction(index[0],index[1],0.0,M);
            for (int k = 0;k<N;k++){
                D(index[0],k) = pow(abs(M(index[0])-M(k)),-a); D(k,index[0]) = D(index[0],k);
                D(index[1],k) = pow(abs(M(index[1])-M(k)),-a); D(k,index[1]) = D(index[1],k);
            }
        D(index[0],index[0]) = 0; D(index[1],index[1]) = 0;

        for (int b=0;b<N;b++) {cout << M(b) << "  ";} cout << endl;
        cout << "\n \n";
        } // end MC loop
        cout << D << endl;
        for (int b=0;b<N;b++) {cout << M(b) << "  ";} cout << endl;
    } //end alpha loop




}// end test_alpha
