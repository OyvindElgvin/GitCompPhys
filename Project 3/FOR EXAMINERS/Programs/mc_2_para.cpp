
#include <iostream>
#include <armadillo>
#include "lib.h"
#include <random>
#include <chrono>
#include <omp.h>

using namespace std;
using namespace arma;
using namespace chrono;


// delcaration of functions
double func_importance_samp(double r1, double r2, double t1, double t2, double p1, double p2);
vector<int> readvalues(string file);

//terminal compiler
//g++-9 -o exe -std=c++11 MC_Functions.cpp mc_2_para.cpp ../Project_3/p3_functions.cpp -L/usr/local/Cellar/armadillo/9.600.6/lib/ -I/usr/local/Cellar/armadillo/9.600.6/include/ -larmadillo -fopenmp



// Monte Carlo with exponential distribution important sampling with parallelizing


int main()
{

    string save_runtimes;
    string save_results;
    cout << "do you want to save runtimes? y or n" << endl;
    cin >> save_runtimes;
    cout << "do you want to save results? y or n" << endl;
    cin >> save_results;

    vector<int> N_values = readvalues("Pro3_Nvalues.txt");
    vec runtimes(N_values.size());

    const double pi =3.141592653589793238463;
    double exact_solution = 5*pi*pi / (16*16);

for (int p = 0; p<N_values.size();p++){

    int n = N_values[p];

    double last_I = 0;
    double average_V = 0;
    double average_error = 0;
    double average_runtime = 0;

    // setting up the jaocbi
    double jacobi = 4* pow(pi,4)/16.0;

    for (int counter = 0; counter < 10; counter ++){
    double MCintIS = 0;
    //double sum_sigmaIS = 0;
    double fy = 0;
    double* f = new double[n];

    // starting clock for time keeping
    high_resolution_clock::time_point time1 = high_resolution_clock::now();


    // parallelizing 4 threads and making a seperate seed for each thread
    #pragma omp parallel num_threads(4)
    {

    // printing thread number for each thread
    //cout << "threads:" << omp_get_thread_num() << endl;

    // random number generator
    unsigned seed = system_clock::now().time_since_epoch().count();
    mt19937_64 generator (seed);

    // parallelizing the sums
    #pragma omp for reduction (+:MCintIS)
    for (int i = 1; i <= n; i++){
        //cout << "threads:" << omp_get_thread_num() << endl;
        // generating r1 and r2 with the exponential distribution
        double g;
        g = generate_canonical< double, 128 > (generator);  // random number [0,1]
        double r1 = -0.25*log(1.- g);                       // random number with exp dist
        g = generate_canonical< double, 128 > (generator);
        double r2 = -0.25*log(1.-g);

        // mapping random numbers [0,1] to theta [0,pi] and phi [0,2*pi]
        g = generate_canonical< double, 128 > (generator);
        double t1 = g*pi;
        g = generate_canonical< double, 128 > (generator);
        double t2 = g*pi;
        g = generate_canonical< double, 128 > (generator);
        double p1 = g*2*pi;
        g = generate_canonical< double, 128 > (generator);
        double p2 = g*2*pi;


        // MC integrating
        fy = func_importance_samp(r1, r2, t1, t2, p1, p2);
        f[i] = fy;
        MCintIS += fy;
        //sum_sigmaIS += fy*fy;
    }
    }

    // stops the clock
    high_resolution_clock::time_point time2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double> >(time2-time1);
    double runtime = time_span.count();


    // calculating the mean integration results and the variance
    MCintIS = MCintIS / (double (n));
    double var = 0;
    for (int i = 1; i <= n; i++){
        var +=  (f[i] - MCintIS) * (f[i] - MCintIS);
    }
    var = var*jacobi / (double (n));

    // Standard deviation
    double sigma = 0;
    sigma = sqrt(var) / sqrt(n);


    //sum_sigmaIS = sum_sigmaIS / (double (n));
    //double variance = sum_sigmaIS - MCintIS * MCintIS;


    //average_V += jacobi*variance;
    last_I = jacobi*MCintIS;
    average_V += var;
    average_error += abs(last_I-exact_solution);
    average_runtime += runtime;

    } //end of average loop

    //average_I = average_I/10.0;
    average_V = average_V/10.0;
    average_error = average_error/10.0;
    average_runtime = average_runtime/10.0;
    runtimes(p) = average_runtime;
    // printing results for Important Samplng

    cout << endl << "Parallelized Monte Carlo importance sampling used " << average_runtime
         << " seconds" << endl;

    cout << endl << "Results for MC important sampling, N = " << n << endl
         << "Variance = " << average_V << endl
         << "Integral = " << last_I << endl
         << "Error = " << average_error << endl
         << "Exact = " << exact_solution << endl;

    if (save_results == "y"){
        if(p == 0){
        ofstream output;
        output.open("Results_ISMC_para.txt",ios::out);
        output << "N = " << n << "   " << "I = " << last_I << "   " << "V = " << average_V << "   " << "Error = " << average_error << endl;
        output.close();
    }
        else{
        ofstream output;
        output.open("Results_ISMC_para.txt",ios::app);
        output << "N = " << n << "   " << "I = " << last_I << "   " << "V = " << average_V << "   " << "Error = " << average_error << endl;
        output.close();

        }
  }



} //end of N loop

if (save_runtimes == "y"){
string filenameruntimes = "ISMC_Runtimes_para.txt";
ofstream output;
output.open(filenameruntimes,ios::out);
for (int i = 0;i<N_values.size();i++){
    output << N_values[i] << endl;
}
output << endl;
output << runtimes << endl;
output.close();
}
else{}


} //end of main
