
#include <iostream>
#include <armadillo>
#include "lib.h"
#include <random>
#include <chrono>


using namespace std;
using namespace arma;
using namespace chrono;


// delcaration of function
double func_importance_samp(double r1, double r2, double t1, double t2, double p1, double p2);
vector<int> readvalues(string file);


// Monte Carlo with exponential distribution important sampling


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

    // random number generator
    unsigned seed = system_clock::now().time_since_epoch().count();
    mt19937_64 generator (seed);



    for (int i = 1; i <= n; i++){
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

    // stops the clock
    high_resolution_clock::time_point time2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double> >(time2-time1);
    double runtime = time_span.count();


    // calculating the mean integration results and the variance
    MCintIS = MCintIS / (double (n));
    double var = 0;
    for (int i = 1; i < n; i++){
        var += 1/double (n) * (f[i] - MCintIS) * (f[i] - MCintIS);
    }
    var = var*jacobi;
    double sigma = 0;
    sigma = sqrt(var) / sqrt(n);


    //sum_sigmaIS = sum_sigmaIS / (double (n));
    //double variance = sum_sigmaIS - MCintIS * MCintIS;


    //average_V += jacobi*variance;
    last_I = jacobi*MCintIS;
    average_V += var;
    average_error += abs(last_I - exact_solution);
    average_runtime += runtime;

    } //end of average loop

    //average_I = average_I/10.0;
    average_V = average_V/10.0;
    average_runtime = average_runtime/10.0;
    average_error = average_error/10.0;
    runtimes(p) = average_runtime;
    // printing results for Important Samplng

    cout << endl << "Monte Carlo importance sampling used " << average_runtime
         << " seconds" << endl;


    cout << endl << "Standard deviation = " << sqrt(average_V/(double(n))) << endl
         << "Variance = " << average_V << endl
         << "Integral = " << last_I << endl
         << "Error = " << average_error << endl
         << "Exact = " << exact_solution << endl
         <<  "N = " << n << endl;

    if (save_results == "y"){
        if(p == 0){
        ofstream output;
        output.open("Results_ISMC.txt",ios::out);
        output << "N = " << n << "   " << "I = " << last_I << "   " << "V = " << average_V << "   " << "Error = " << average_error << endl;
        output.close();
    }
        else{
        ofstream output;
        output.open("Results_ISMC.txt",ios::app);
        output << "N = " << n << "   " << "I = " << last_I << "   " << "V = " << average_V << "   " << "Error = " << average_error <<endl;
        output.close();

        }
  }



} //end of N loop

if (save_runtimes == "y"){
string filenameruntimes = "ISMC_Runtimes_non_para.txt";
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
