#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include "functions.h"
#include "analysis.h"

using namespace std;


Functions :: Functions(){}

Functions :: ~Functions(){}

double Functions :: error(vector<double> AV, vector<double> AV2, int n){
    if(n == 0) {
        return 0;
    }
    else {
        return sqrt((AV2[n] - AV[n]*AV[n])/n);
    }

}

double Functions :: GetMean (double *A, int size){
    double somma = 0;
    for(int i = 0; i<size; i++){
        somma += A[i];
    }
    return somma/size;
}

double Functions :: GetN (double *A, int n){
    return A[n];
}

double Functions :: norm (vector<double> x) {
    double sum = 0.;
    for (int i = 0; i<x.size(); i++){
        sum += x[i]*x[i];
    }

    sum = sqrt(sum);

    return sum;
}

double Functions :: average(vector<double> x){
    double sum = 0; 
    for(int i = 0; i < x.size(); i++){
        sum += x[i];
    }
    sum /= x.size();
    return sum;

}

vector<double> Functions :: av(vector<double> x, int N){

    vector<double> av;

    int K = int(x.size()/N);
    for(int i=0; i<N; i++){
      double sum = 0.;
      for (int j=0; j<K; j++){
         sum = sum + x[i*K+j];
      }
      sum = sum/K; //I obtain the mean after K throws which means the mean of the i-th block
      av.push_back(sum);
   }

   return av;
}

vector<double> Functions :: av2(vector<double> x, int N){

    vector<double> av2;

    int K = int(x.size()/N);
    for(int i=0; i<N; i++){
      double sum = 0.;
      for (int j=0; j<K; j++){
         sum = sum + x[i*K+j];
      }
      sum = sum/K; //I obtain the mean after K throws which means the mean of the i-th block
      av2.push_back(sum*sum);

   }

   return av2;
}

Analysis Functions :: graph(vector<double> x, vector<double> y){
    Functions fun;
    Analysis sum_err;
    vector<double> su2_prog;

    for (int i = 0; i<x.size(); i++){
        sum_err.errprog.push_back(0.);
        sum_err.sumprog.push_back(0.);
        su2_prog.push_back(0.);

   }
    for(int i=0; i<x.size(); i++){
      for (int j=0; j<i+1; j++){
         sum_err.sumprog[i] += x[j];
         su2_prog[i] += y[j];
         
      }
      sum_err.sumprog[i] = sum_err.sumprog[i]/(i+1);
      su2_prog[i] = (su2_prog[i]/(i+1));
      sum_err.errprog[i] = fun.error(sum_err.sumprog,su2_prog,i); //sum_prog is the progression of the mean so it should go towards the true mean
   }

   return sum_err;

}



void Functions :: print_on_file_errbar_expected(vector<double> x, vector<double> y, ostream& fout, double exp){
    for (int i = 0; i<x.size(); i++){
        fout << i << "\t" << x[i] << "\t" << y[i] << "\t" << exp << endl;
    }
}

void Functions :: print_on_file_errbar(vector<double> x, vector<double> y, ostream& fout){
    for (int i = 0; i<x.size(); i++){
        fout << i << "\t" << x[i] << "\t" << y[i] << "\t" << endl;
    }
}

void Functions :: print_on_file(vector<double> x, ostream& fout){
    for (int i = 0; i<x.size(); i++){
        fout << i << "\t" << x[i] << endl;
    }
}

vector<double> Functions :: cut_end(vector<double> x, int N){

    vector<double> v1;

    for (int i = 0; i<N ; i++){
        v1.push_back(x[i]);
    }

    return v1;
}

vector<double> Functions :: cut_begin(vector<double> x, int N){

    vector<double> v1;

    for (int i = N; i<x.size() ; i++){
        v1.push_back(x[i]);
    }

    return v1;
}