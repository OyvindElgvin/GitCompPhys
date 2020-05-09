#include <iostream>
#include <fstream>
#include <armadillo>
#include <omp.h>
#include <stdio.h>
#include "time.h"
#include <random>
#include <chrono>

#include "Pro5_Functions.h"

using namespace std;
using namespace arma;

int main(){

    string Tests;
    cout << "Do you want to run tests? y/n \n";
    cin >> Tests;
    if (Tests == "y"){
        //test_alpha();
        //test_gamma();
        test_sampling();
        //test_stability();
        //test_stability_savings();
        test_stability_taxes();
    }


    int Ex = 0;
    int Cycles = 0;
    cout << "Number of experiment runs? \n";
    cin >> Ex;
    cout << "Number of MC Cycles? \n";
    cin >> Cycles;


    taxes(Ex,Cycles);

    //Financial_analysis(Ex,Cycles,1000,"Money_distributions_D_7_N_1000_L_0.000000_a_2.000000_g_2.000000","Median_D_7_N_1000_L_0.000000_a_2.000000_g_2.000000",0,2.0,2.0);
    //task_d(Ex,Cycles);

    /*Median_D_7_N_1000_L_0.000000_a_2.000000_g_1.000000
     * Median_D_7_N_1000_L_0.500000_a_2.000000_g_1.000000
     */


    return 0;
} // end main
