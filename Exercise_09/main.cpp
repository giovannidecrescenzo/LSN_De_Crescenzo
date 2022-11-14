#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <vector>
#include <numeric>
#include "../Random/random.h"
#include "../Functions/functions.h"
#include "population.h"

using namespace std;

int main (int argc, char *argv[]){
//FIRST PART: TSP ON A CIRCUMFERENCE
   Functions fun;
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


   int N = 34;
   //number of cities
   int M = 400;
   //number of chromosomes

   Population pop(N,M,rnd,0);

   for(int j = 0; j < pop.n_chromo; j++){
      for(int i = 0; i < pop.n_chromo; i++){
         if(pop.population[i]==pop.population[j] and i!=j){
            cout << "SAME CHROMOSOMES" << endl;
            cout << "cycling here" << endl;
         }
      }
   }
   ofstream cities;
   cities.open("circle_cities.txt");
   for(int i = 0; i < N; i ++){
      cities << pop.cities[i] << endl;
   }
   //Here i print the first try for circle
   ofstream myfile1;
   myfile1.open("first_try_circle.txt");
   for(int i = 0; i < N; i++){
      myfile1 << pop.population[0][i] << endl;
   }
   myfile1.close();
   ofstream myfile2;
   myfile2.open("Mean_best_half_circle.txt");
   for(int i = 0; i < 40000; i++){
      pop.beta = i + 1;
      pop.Mutation();
      pop.mixing();
      if ((i+1)%400 == 0) {
         double mean = 0;
         for(int j = 0; j < int(pop.fitness.size()/2); j++){
            mean += pop.fitness[j];
         }
         mean = mean/int(pop.fitness.size()/2);
         myfile2 << mean << endl;
      }
   }
   myfile2.close();
   cout << "NOW I CHECK SOME POP "  << endl;
   for(int i = 0; i < pop.n_chromo; i++){
      pop.check(i);
   }
   cout << "CHECK DONE"  << endl;
   ofstream myfile3;
   myfile3.open("best_circle.txt");
   for(int i = 0; i < N; i ++){
      myfile3 << pop.population[0][i] << endl;
   }
   myfile3.close();

//NOW SQUARE CITIES

   pop.type = 1;
   pop.Init_cities();
   pop.Init_pop();

   for(int j = 0; j < pop.n_chromo; j++){
      for(int i = 0; i < pop.n_chromo; i++){
         if(pop.population[i]==pop.population[j] and i!=j){
            cout << "SAME CHROMOSOMES" << endl;
            cout << "cycling here" << endl;
         }
      }
   }
   pop.Fill_fit();
   ofstream cities_square;
   cities_square.open("square_cities.txt");
   for(int i = 0; i < N; i ++){
      cities_square << pop.sq_cities[i][0] << "\t" << pop.sq_cities[i][1] << endl;
   }
   //Here i print the first try for circle
   ofstream myfile4;
   myfile4.open("first_try_square.txt");
   for(int i = 0; i < N; i++){
      myfile4 << pop.population[0][i] << endl;
   }
   myfile4.close();
   ofstream myfile5;
   myfile5.open("Mean_best_half_square.txt");
   for(int i = 0; i < 200000; i++){
      pop.beta = i + 1;
      pop.Mutation();
      if(rnd.Rannyu() < 0.1) pop.mixing();
      if (i%400 == 0) {
         double mean = 0;
         for(int j = 0; j < int(pop.fitness.size()/2); j++){
            mean += pop.fitness[j];
         }
         mean = mean/int(pop.fitness.size()/2);
         myfile5 << mean << endl;
      }
   }
   myfile5.close();
   cout << "NOW I CHECK SOME POP "  << endl;
   for(int i = 0; i < pop.n_chromo; i++){
      pop.check(i);
   }
   cout << "CHECK DONE"  << endl;
   ofstream myfile6;
   myfile6.open("best_square.txt");
   for(int i = 0; i < N; i ++){
      myfile6 << pop.population[0][i] << endl;
   }
   myfile6.close();



   rnd.SaveSeed();

   return 0;
}