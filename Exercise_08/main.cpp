/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <vector>
#include <numeric>
#include "../Random/random.h"
#include "../Functions/functions.h"
#include "../Functions/analysis.h"
#include "../Metropolis/metropolis.h"

using namespace std;

double psi(double x, double mu, double sigma){
   double out;
   out = exp(-(x-mu)*(x-mu)/(2*sigma*sigma))+exp(-(x+mu)*(x+mu)/(2*sigma*sigma));
   return out;
}

double potential(double x){
   return pow(x,4)-2.5*pow(x,2);
}

double energy(double x, double mu, double sigma){

   double s2 = pow(sigma,2);
   double kinetic;
   kinetic = 1./(2*s2)*(-exp(-(pow((x-mu),2))/(2*s2))/s2*pow((x-mu),2)-exp(-(pow((x+mu),2))/(2*s2))/s2*pow((x+mu),2)+psi(x,mu,sigma));
   return  potential(x) + kinetic/psi(x,mu,sigma);
}

vector <double> unif_sampling_1D(double x, int N, double L, double mu, double sigma){


//This function is used to perform a Metropolis uniform sampling of a 1D probability distribution; x is the starting point, N the number of steps and L the length of the step; 
//func is the analytic form of the distribution
   double count = 0.;
   Random rnd;
   Analysis sys;
   Functions fun;

   int seed[4];
   int p1, p2;
   ifstream Primes("../Random/Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("../Random/seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;
   vector <double> sampling;

   sampling.push_back(x);

   for (int i = 0; i<N; i++) {

      double trial = (rnd.Rannyu()-0.5)*L+sampling[i];
      
      if(rnd.Rannyu()<pow(psi(trial, mu, sigma),2)/pow(psi(sampling[i], mu, sigma),2)){
         sampling.push_back(trial);
         count += 1;
      }
      else {
         sampling.push_back(sampling[i]);
      }
   }


   rnd.SaveSeed();

   double acceptance = count/N;

   //cout << "Acceptance with uniform sampling is: " << acceptance << endl;

   return sampling;
}







int main (int argc, char *argv[]){

   Metropolis met;
   Functions fun;
   Analysis sys;
   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("../Random/Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("../Random/seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;


//VARIATIONAL MONTE CARLO; MU AND SIGMA FIXED
   ofstream myfile;
   ofstream myfile2;

   myfile.open("sampling.txt");
   myfile2.open("integral.txt");

   double mu = 1.;
   double sigma = 0.5;
   double x_0 = 0.; //starting position
   int N = 1000000; //total number of sampling steps
   int M = 100; //number of bins
   int K = int(N/M); //number of throws in each bin
   double L = 3.0;
   vector<double> sampling = unif_sampling_1D(x_0, N, L, mu, sigma);
   fun.print_on_file(sampling,myfile);
   vector<double> integral;
   
   for(int i = 0; i<sampling.size(); i++){
      integral.push_back(energy(sampling[i], mu, sigma));
   }


   sys = fun.graph(fun.av(integral,M), fun.av2(integral,M));
   fun.print_on_file_errbar(sys.sumprog, sys.errprog, myfile2);

   myfile.close();
   myfile2.close();

//END OF VARIATIONAL MONTE CARLO

//SIMULATED ANNEALING

   double beta = 1./2;
   vector<vector<double> > mu_s_E;
   vector<double> first;
   first.push_back(mu);
   first.push_back(sigma);
   first.push_back(sys.sumprog.back());
   first.push_back(sys.errprog.back());
   mu_s_E.push_back(first);
   ofstream myfile3;
   myfile3.open("Mu_sigma.txt");
   double acc, att;
   double del_mu, del_sigma;
   for(int i = 0; i < 1000; i++){
      acc = 0.;
      att = 0.;
      if(mu_s_E.back()[2]>=0.4){
         del_mu = 0.1;
         del_sigma = 0.1;
      }
      if(mu_s_E.back()[2]>(-0.4) and mu_s_E.back()[2]<0.4){
         del_mu = 0.01;
         del_sigma = 0.005;
      }
      if(mu_s_E.back()[2]<=(-0.4)){
         del_mu = 0.005;
         del_sigma = 0.001;
         //cout << "ENERGY BELOW -0.4" << endl << endl;
      }

      for(int j = 0; j < 10 + (2*int(i/100)); j++){
         int N_1 = 100000;
         double mu_try = mu + (rnd.Rannyu()-0.5)*2*del_mu;
         double sigma_try = sigma + (rnd.Rannyu()-0.5)*2*del_sigma;
         vector<double> sampling1 = unif_sampling_1D(x_0, N_1, L, mu_try, sigma_try);
         vector<double> integral1;
         for(int i = 0; i<sampling1.size(); i++){
            integral1.push_back(energy(sampling1[i], mu_try, sigma_try));
         }
         sys = fun.graph(fun.av(integral1,M), fun.av2(integral1,M));
         double E_new = sys.sumprog.back();
         double E_old = mu_s_E.back()[2];
         double p = exp(-beta*(E_new-E_old));
         //cout << "PROBABILITY IS: " << p << endl;
         if(rnd.Rannyu() < p){
            vector<double> temp;
            temp.push_back(mu_try);
            temp.push_back(sigma_try);
            temp.push_back(E_new);
            temp.push_back(sys.errprog.back());
            mu_s_E.push_back(temp);
            mu = mu_try;
            sigma = sigma_try;
            acc ++;
            att ++;
         }
         else{
            mu_s_E.push_back(mu_s_E.back());
            att ++;
         }
      }
      beta = beta + 1;
      //std::cout << "BETA IS NOW: " << beta << endl;
      myfile3 << mu_s_E.back()[0] << "\t" << mu_s_E.back()[1] << "\t" << mu_s_E.back()[2] << "\t" << mu_s_E.back()[3]
      << "\t" <<  "ACCEPTANCE: " << double(acc/att) << endl;
   }
   myfile3.close();
   ofstream myfile4;
   ofstream myfile5;
   myfile4.open("Best_sampling.txt");
   myfile5.open("Best_integral.txt");

   sampling = unif_sampling_1D(x_0, N, L, mu, sigma);
   fun.print_on_file(sampling,myfile4);
   integral.clear();
   
   for(int i = 0; i<sampling.size(); i++){
      integral.push_back(energy(sampling[i], mu, sigma));
   }


   sys = fun.graph(fun.av(integral,M), fun.av2(integral,M));
   fun.print_on_file_errbar(sys.sumprog, sys.errprog, myfile5);
   myfile4.close();
   myfile5.close();
   
   return 0;
}


/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
