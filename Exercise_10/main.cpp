#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <vector>
#include <numeric>
#include "../Random/random.h"
#include "../Functions/functions.h"
#include "../Exercise_09/population.h"
#include "mpi.h"

using namespace std;

int main (int argc, char *argv[]){
//MPI
   int size, rank;
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Status stat0, stat1;

   Functions fun;
   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      for (int i=0; i<rank+1; i++) {
         Primes >> p1 >> p2 ;
      }
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
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
   string input_cities = "Square/square_cities.dat";
   Population pop(N, M, rnd, 1, input_cities);


//NOW SQUARE CITIES
   ofstream cities_square;
   cities_square.open("Square/cities_" + to_string(rank) + ".txt");
   for(int i = 0; i < N; i ++){
      cities_square << pop.sq_cities[i][0] << "\t" << pop.sq_cities[i][1] << endl;
   }
   //Here i print the first try for circle
   ofstream myfile1;
   myfile1.open("Square/First_try_" + to_string(rank) + ".txt");
   for(int i = 0; i < N; i++){
      myfile1 << pop.population[0][i] << endl;
   }
   myfile1.close();

   ofstream myfile2;
   int rank_0, rank_1;
   vector<int> chrom_0(N), chrom_1(N);
   myfile2.open("Square/Mean_best_half_" + to_string(rank) + ".txt");
   for(int i = 0; i < 40000; i++){
      pop.Mutation();
      pop.mixing();
      //choosing two random cores every 100 iterations
      if((i+1) % 100 == 0){
         int itag0 = 1;
         if(rank==0){
            rank_0 = int(rnd.Rannyu(0,size));
            rank_1 = int(rnd.Rannyu(0,size));
            while (rank_0==rank_1) {
               rank_1 = int(rnd.Rannyu(0,size));
            }
         }
         MPI_Bcast(&rank_0, 4, MPI_INTEGER, 0, MPI_COMM_WORLD);
         MPI_Bcast(&rank_1, 4, MPI_INTEGER, 0, MPI_COMM_WORLD);
         //choosing best from each chosen core 
               
         if (rank==rank_0) {
            vector<int> best_0;
            best_0 = pop.population[0];
            for (int j=0; j<N; j++) chrom_0[j] = best_0[j];
         }
         if (rank==rank_1) {
            vector<int> best_1;
            best_1 = pop.population[0];
            for (int j=0; j < N; j++) chrom_1[j] = best_1[j];
         }
               
         //Switching the 2 best chromosomes; I need to send the address of the vector
         if (rank==rank_1) {
            MPI_Send(&chrom_1[0], N, MPI_INTEGER, rank_0, itag0, MPI_COMM_WORLD);
            MPI_Recv(&chrom_0[0], N, MPI_INTEGER, rank_0, itag0, MPI_COMM_WORLD, &stat0);
                  
            vector<int> best;
            for (int j=0; j < N; j++) best.push_back(chrom_0[j]);
               pop.population[0] = best;
            }
         if(rank==rank_0){
            MPI_Send(&chrom_0[0], N, MPI_INTEGER, rank_1, itag0, MPI_COMM_WORLD);
            MPI_Recv(&chrom_1[0], N, MPI_INTEGER, rank_1, itag0, MPI_COMM_WORLD, &stat0);
            vector<int> best;
            for (int j=0; j<N; j++) best.push_back(chrom_1[j]);
            pop.population[0] = best;
         }
      }
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

   ofstream myfile3;
   myfile3.open("Square/Best_individual_" + to_string(rank) + ".txt");
   for(int j = 0; j < N; j++){
      myfile3 << pop.population[0][j] << endl;
   }
   myfile3.close();
   if(rank == 0){
      cout << "=============SQUARE CITIES DONE=============" << endl;
   }
   //NOW AMERICAN CITIES
   N = 50;
   //number of cities
   M = 1000;
   //number of chromosomes
   input_cities = "America/American_coord.dat";

   Population pop_1(N, M, rnd, 1, input_cities);

   ofstream cities_america;
   cities_america.open("America/cities_" + to_string(rank) + ".txt");
   for(int i = 0; i < N; i ++){
      cities_america << pop_1.sq_cities[i][0] << "\t" << pop_1.sq_cities[i][1] << endl;
   }
   //Here i print the first try for circle
   ofstream myfile4;
   myfile4.open("America/first_try_" + to_string(rank) + ".txt");
   for(int i = 0; i < N; i++){
      myfile4 << pop_1.population[0][i] << endl;
   }
   myfile4.close();

   ofstream myfile5;
   chrom_0.clear();
   chrom_1.clear();
   for(int j = 0; j < N; j++) chrom_0.push_back(0);
   for(int j = 0; j < N; j++) chrom_1.push_back(0);
   myfile5.open("America/Mean_best_half_" + to_string(rank) + ".txt");
   for(int i = 0; i < 100000; i++){
      pop_1.Mutation();
      pop_1.mixing();
      //if ((i+1) % 100 == 0) pop_1.mixing();
      if((i+1) % 100 == 0){
         int itag0 = 1;
         if(rank==0){
            rank_0 = int(rnd.Rannyu(0,size));
            rank_1 = int(rnd.Rannyu(0,size));
            while (rank_0==rank_1) {
               rank_1 = int(rnd.Rannyu(0,size));
            }
         }
         MPI_Bcast(&rank_0, 4, MPI_INTEGER, 0, MPI_COMM_WORLD);
         MPI_Bcast(&rank_1, 4, MPI_INTEGER, 0, MPI_COMM_WORLD);
               
         if (rank==rank_0) {
            vector<int> best_0;
            best_0 = pop_1.population[0];
            for (int j=0; j<pop_1.n_cities; j++) chrom_0[j] = best_0[j];
         }
         if (rank==rank_1) {
            vector<int> best_1;
            best_1 = pop_1.population[0];
            for (int j=0; j<pop_1.n_cities; j++) chrom_1[j] = best_1[j];
         }
            
         if (rank == rank_1) {
            MPI_Send(&chrom_1[0], N, MPI_INTEGER, rank_0, itag0, MPI_COMM_WORLD);
            MPI_Recv(&chrom_0[0], N, MPI_INTEGER, rank_0, itag0, MPI_COMM_WORLD, &stat0);
                  
            vector<int> best;
            for (int j=0; j < N; j++) best.push_back(chrom_0[j]);
               pop_1.population[0] = best;
            }
         if(rank==rank_0){
            MPI_Send(&chrom_0[0], N, MPI_INTEGER, rank_1, itag0, MPI_COMM_WORLD);
            MPI_Recv(&chrom_1[0], N, MPI_INTEGER, rank_1, itag0, MPI_COMM_WORLD, &stat0);
                  
            vector<int> best;
            for (int j = 0; j < N; j++) best.push_back(chrom_1[j]);
            pop_1.population[0] = best;
         }
      }
      if ((i+1)%1000 == 0) {
         double mean = 0;
         for(int j = 0; j < int(pop_1.fitness.size()/2); j++){
            mean += pop_1.fitness[j];
         }
         mean = mean/int(pop_1.fitness.size()/2);
         myfile5 << mean << endl;
      }

   }
   myfile5.close();

   ofstream myfile6;
   myfile6.open("America/Best_individual_" + to_string(rank) + ".txt");
   for(int j = 0; j < N; j++){
      myfile6 << pop_1.population[0][j] << endl;
   }
   myfile6.close();

   pop_1.check(0);

   if(rank == 0){
      cout << "=============AMERICAN CAPITALS DONE=============" << endl;
   }

   MPI_Finalize();


   rnd.SaveSeed();

   return 0;
}