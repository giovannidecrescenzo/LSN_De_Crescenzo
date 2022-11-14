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

   ofstream myfile1;
   ofstream myfile2;
   ofstream myfile3;
   ofstream myfile4;
   ofstream myfile5;

   myfile1.open("uniform_integral.txt");
   myfile2.open("linear_sampling_integral.txt");
   myfile3.open("RW_cubic_lattice.txt");
   myfile4.open("RW_free_space.txt");
   myfile5.open("Test.txt");

   Random rnd;
   Functions fun;
   Analysis sys;


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

   //Here are defined the values for the exercise
   int M = 1000000; //total number of throws
   int N = 100; //number of blocks
   int K = int(M/N); //number of throws in each block

   vector<double> uniform;
   vector<double> importance;

   for (int i = 0; i<M; i++){
      uniform.push_back(rnd.Rannyu());
      importance.push_back(rnd.Linear(1.));
   }

   for (int i = 0; i<M; i++){
      uniform[i] = M_PI/2 * cos(M_PI/2*uniform[i]);
      importance[i] = M_PI/2 * cos(M_PI/2*importance[i]) / (-2*importance[i]+2);
   }

   sys = fun.graph(fun.av(uniform, N), fun.av2(uniform, N));
   fun.print_on_file_errbar_expected(sys.sumprog, sys.errprog, myfile1, 1.0);

   sys = fun.graph(fun.av(importance, N), fun.av2(importance, N));
   fun.print_on_file_errbar_expected(sys.sumprog, sys.errprog, myfile2, 1.0);

   myfile1.close();
   myfile2.close();

   //Now part about RW

   vector<vector<double> > pos;
   vector<double> v1;

   for (int i = 0; i<3; i++){
      v1.push_back(0.);
   }

   pos.push_back(v1);

   for (int i = 1; i<M; i++){
      pos.push_back(pos[i-1]);
      double x;
      if(rnd.Rannyu()<0.5){
         x = -1;
      }
      else{
         x = 1;
      }
      pos[i][int(rnd.Rannyu()*3)] += x;
      if(i%100 == 0){
         pos[i] = v1;
      }
   
   }
   for (int i = 0; i<M; i++){
      myfile5 << pos[i][0] << "\t" << pos[i][1] << "\t" <<pos[i][2] << endl;
   }
   //Here I have collected the 10^4 random walks


   for(int j = 0; j < 100; j++){
      vector <double> temp;
      for(int i = 0; i < M/100; i += 100){
         temp.push_back(fun.norm(pos[i+j]));
      }
      sys = fun.graph(fun.av(temp, 100), fun.av2(temp, 100));
      myfile3 << j << "\t" << sys.sumprog.back() << "\t" << sys.errprog.back() << endl;
   }



   myfile3.close();



   vector<vector<double> > pos1;
   vector<double> v2;

   for (int i = 0; i<3; i++){
      v2.push_back(0.);
   }

   pos1.push_back(v2);

   for (int i = 1; i<M; i++){
      pos1.push_back(pos1[i-1]);
      double theta = rnd.Rannyu() * M_PI;
      double phi = rnd.Rannyu() * 2 * M_PI;
      v2[0] = cos(phi) * sin(theta);
      v2[1] = sin(phi) * sin(theta);
      v2[2] = cos(theta);
      for (int j = 0; j<3; j++){
         pos1[i][j] += v2[j];
      }
      if(i%100 == 0){
         pos1[i] = v2;
      }
   }

   for(int j = 0; j < 100; j++){
      vector <double> temp;
      for(int i = 0; i < M/100; i += 100){
         temp.push_back(fun.norm(pos[i+j]));
      }
      sys = fun.graph(fun.av(temp, 100), fun.av2(temp, 100));
      myfile4 << j << "\t" << sys.sumprog.back() << "\t" << sys.errprog.back() << endl;
   }

   //sys = fun.graph(fun.av(distance1, N), fun.av2(distance1, N));


   myfile4.close();

   myfile5.close();


   rnd.SaveSeed();
   return 0;
}
