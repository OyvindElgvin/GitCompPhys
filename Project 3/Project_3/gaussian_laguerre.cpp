#include <iostream>
#include "lib.h"
#include <armadillo>
#include "time.h"

using namespace std;
using namespace arma;

vector<int> readvalues(string file);
double int_function_spherical(double r1,double r2,double t1,double t2,double p1,double p2);

void gauss_laguerre(double *x, double *w, int n, double alf);

int main(){


  //Decide if you want to save information or not (for easy quick-runs, or preserving old results)
  string save_runtimes;
  string save_results;
  cout << "do you want to save runtimes? y or n" << endl;
  cin >> save_runtimes;
  cout << "do you want to save results? y or n" << endl;
  cin >> save_results;

  vector<int> N_values = readvalues("Pro3_Nvalues.txt");
  vec runtimes(N_values.size());

  // Looping over all N-values in list
  for (int h = 0;h < N_values.size(); h++){

  // Setting up N and initalize arrays
  int N = N_values[h];

  double *xr = new double[N+1];
  double *wr = new double[N+1];

  double *xt = new double[N];
  double *wt = new double[N];

  double *xp = new double[N];
  double *wp = new double[N];

  // Start the clock
  time_t start, end;
  start = clock();

  // Weights and int points for Theta
  gauleg(0,M_PI,xt,wt,N);
  // Weights and int points for Phi
  gauleg(0,2*M_PI,xp,wp,N);
  // Weights and int points for radial part
  gauss_laguerre(xr,wr,N,0);

  // Perform sum over weights and function values
  double I = 0;
  for (int i=1; i<N+1;i++){
      for (int j=1; j<N+1;j++){
          for (int k=0; k<N;k++){
              for (int l=0; l<N;l++){
                  for (int m=0; m<N;m++){
                      for (int n=0; n<N;n++){
                          I += wr[i]*wr[j]*wt[k]*wt[l]*wp[m]*wp[n]*int_function_spherical(xr[i],xr[j],xt[k],xt[l],xp[m],xp[n]);
  }}}}}}

  // Stop the clock
  end = clock();

  runtimes(h) = (double)(end-start)/CLOCKS_PER_SEC;

  // Printing Results to terminal
  cout << "-----------------------------" << endl;
  cout << "N = " << N << endl;
  cout << "I = " << I << endl;
  cout << "Runtime = " << (double)(end-start)/CLOCKS_PER_SEC << endl;
  //cout << 5*M_PI*M_PI/(256) << endl;
  cout << "-----------------------------" << endl;


  // Saving results to file
  if (save_results == "y"){
      if(h == 0){
      //string filenameresults = "Results_Laguerre.txt";
      ofstream output;
      output.open("Results_Laguerre_test.txt",ios::out);
      output << "N = " << N << "   " << "I = " << I << endl;
      output.close();
  }
      else{
      //string filenameresults = "Results_Laguerre.txt";
      ofstream output;
      output.open("Results_Laguerre_test.txt",ios::app);
      output << "N = " << N << "   " << "I = " << I << endl;
      output.close();

      }
}



  } //End N loop

  // Save runtimes to file
  if (save_runtimes == "y"){
  string filenameruntimes = "Gauss_Laguerre_Runtimes.txt";
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

  return 0;} // End of main


