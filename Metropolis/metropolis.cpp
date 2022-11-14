#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include "../Functions/functions.h"
#include "../Random/random.h"
#include "metropolis.h"

using namespace std;



Metropolis :: Metropolis(){}

Metropolis :: ~Metropolis(){}

vector<vector<double> > Metropolis :: unif_sampling_3D(vector<double> x, int N, double L, function<double(vector<double>)> func){

//This function is used to perform a Metropolis uniform sampling of a 3D probability distribution; x is the starting point, N the number of steps and L the length of the step; 
//func is the analytic form of the distribution
    double count = 0.;
    Random rnd;

    vector<vector<double> > uniform_sampling;

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

    uniform_sampling.push_back(x);

    for (int i = 0; i<N; i++) {
        vector<double> trial1;
        for (int j = 0; j<3; j++){
        trial1.push_back((rnd.Rannyu()-0.5)*L+uniform_sampling[i][j]);
    }

    if(rnd.Rannyu()<func(trial1)/func(uniform_sampling[i])){
        uniform_sampling.push_back(trial1);
        count += 1;
    }
    else {
        uniform_sampling.push_back(uniform_sampling[i]);
    }
    }

   rnd.SaveSeed();

   acceptance = count/N;

   cout << "Acceptance with uniform sampling is: " << acceptance << endl;

   return uniform_sampling;


}



vector<vector<double> > Metropolis :: gauss_sampling_3D(vector<double> x, int N, double sigma,function<double(vector<double>)> func){
//This function is used to perform a Metropolis gaussian sampling of a 1D probability distribution; x is the starting point, N the number of steps and sigma the width of the gaussian "step"; 
//func is the analytic form of the distribution
    double count = 0.;
    Random rnd;

    vector<vector<double> > gauss_sampling;

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

    gauss_sampling.push_back(x);

    for (int i = 0; i<N; i++) {
        vector<double> trial1;
        for (int j = 0; j<3; j++){
        trial1.push_back(rnd.Gauss(0, sigma)+gauss_sampling[i][j]);
        }

        if(rnd.Rannyu()<func(trial1)/func(gauss_sampling[i])){
            gauss_sampling.push_back(trial1);
            count += 1;
        }
        else {
            gauss_sampling.push_back(gauss_sampling[i]);
        }
    }

   rnd.SaveSeed();

   acceptance = count/N;

   cout << "Acceptance with Gauss sampling is: " << acceptance << endl;

   return gauss_sampling;


}
