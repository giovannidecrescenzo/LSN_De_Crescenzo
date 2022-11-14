#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <vector>
#include "../Random/random.h"
#include "../Functions/functions.h"
#include "../Functions/analysis.h"

using namespace std;
 
int main (int argc, char *argv[]){

   Random rnd;
   Functions fun;

   ofstream myfile1;
   ofstream myfile2;
   ofstream myfile3;
   ofstream myfile4;

   myfile1.open("direct_sample_call.txt");
   myfile2.open("direct_sample_put.txt");
   myfile3.open("indirect_sample_call.txt");
   myfile4.open("indirect_sample_put.txt");


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



   int M = 1000000; //total trials
   int N = 100; //bins
   int K = int(M/N); //trials per bin
   double r = 0.1; //risk-free interest rate or drift of GBM
   double sigma = 0.25; //volatility
   double T = 1; //expiry date
   double s_0 = 100; //initial price
   double k = 100; //strike price
 
   vector<double> direct_call, direct_put;
   vector<double> indirect_call, indirect_put;
   Analysis sys;

//SAMPLING OF OPTION COST VIA DIRECT SAMPLE OF FINAL ASSET PRICE

   for (int i = 0; i < M; i++){
      double z = rnd.Gauss(0,1.);
      double s = s_0*exp((r-sigma*sigma/2)*T+sigma*z*sqrt(T));
      double c = exp(-r*T)*max(0.,s-k);
      direct_call.push_back(c);
   }
   sys = fun.graph(fun.av(direct_call,N), fun.av2(direct_call,N));
   fun.print_on_file_errbar(sys.sumprog, sys.errprog, myfile1);

   myfile1.close();

   for (int i = 0; i < M; i++){
      double z = rnd.Gauss(0,1.);
      double s = s_0*exp((r-sigma*sigma/2)*T+sigma*z*sqrt(T));
      double c = exp(-r*T)*max(0.,k-s);
      direct_put.push_back(c);
   }
   sys = fun.graph(fun.av(direct_put,N), fun.av2(direct_put,N));
   fun.print_on_file_errbar(sys.sumprog, sys.errprog, myfile2);

   myfile2.close();

//SAMPLING OF OPTION COST VIA INDIRECT SAMPLE OF FINAL ASSET PRICE

   double del_t =  1./100;

   for (int i = 0; i < M; i++){
      double s = s_0;
      double st = 0.;
      for(int j = 0; j < int(T/del_t); j++){
         double z = rnd.Gauss(0.,1.);
         st = s*exp((r-sigma*sigma/2)*del_t+sigma*z*sqrt(del_t));
         s = st;
      }
      double c = exp(-r*T)*max(0.,s-k);
      indirect_call.push_back(c);
   }
   sys = fun.graph(fun.av(indirect_call,N), fun.av2(indirect_call,N));
   fun.print_on_file_errbar(sys.sumprog, sys.errprog, myfile3);
   myfile3.close();

   for (int i = 0; i < M; i++){
      double s = s_0;
      double st = 0.;
      for(int j = 0; j < int(T/del_t); j++){
         double z = rnd.Gauss(0.,1.);
         st = s*exp((r-sigma*sigma/2)*del_t+sigma*z*sqrt(del_t));
         s = st;
      }
      double c = exp(-r*T)*max(0.,k-s);
      indirect_put.push_back(c);
   }
   sys = fun.graph(fun.av(indirect_put,N), fun.av2(indirect_put,N));
   fun.print_on_file_errbar(sys.sumprog, sys.errprog, myfile4);
   myfile4.close();



   rnd.SaveSeed();
   return 0;
}
