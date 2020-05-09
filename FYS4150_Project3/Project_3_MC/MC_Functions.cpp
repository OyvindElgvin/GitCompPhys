

#include <iostream>
#include <armadillo>
#include "lib.h"
#include <random>
#include <chrono>

using namespace std;
using namespace arma;




// function for brute forcing Monte Carlo
double func_brute(double x1, double y1,double z1, double x2, double y2, double z2){
    double alpha = 2.0;
    double r1 = sqrt(x1*x1 + y1*y1 + z1*z1);
    double r2 = sqrt(x2*x2 + y2*y2 + z2*z2);
    double r1_r2 = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
    if (r1_r2 < pow(10,-10)){return 0;}
    double I = exp(-2*alpha*(r1+r2)) / r1_r2;
    return I;
}




// function for the exponentional distribution importance sampling
double func_importance_samp(double r1, double r2, double t1, double t2, double p1, double p2){
    double dr1dr2 = r1*r1 * r2*r2 *sin(t1)*sin(t2);
    double beta = cos(t1)*cos(t2) + sin(t1)*sin(t2)*cos( p1 - p2 );
    double r12 = r1*r1 + r2*r2 - 2*r1*r2 * beta;
    if (r12 < pow(10,-10)){return 0;}
    else{
    r12 = sqrt(r12);
    double I = dr1dr2 / r12;
    return I;}
}
