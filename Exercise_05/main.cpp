#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include "../Functions/functions.h"
#include "../Functions/analysis.h"
#include "../Random/random.h"
#include "../Metropolis/metropolis.h"
#include <vector>


using namespace std;

double probability1 (vector<double> x){

    Functions fn;


    return exp(-2*fn.norm(x));
}

double probability2 (vector<double> x){

    Functions fn;

    return x[2]*x[2]*exp(-fn.norm(x));
}



int main (int argc, char *argv[]){

   
   Functions fun;
   Analysis sys;
   Metropolis met;

   ofstream myfile1;
   ofstream myfile2;
   ofstream myfile3;
   ofstream myfile4;
   ofstream myfile5;
   ofstream myfile6;
   ofstream myfile7;
   ofstream myfile8;
   ofstream myfile9;
   ofstream myfile10;
   ofstream myfile11;
   ofstream myfile12;

   myfile1.open("Unif_sampling_1.txt");
   myfile2.open("Unif_sampling_2.txt");
   myfile3.open("Gauss_sampling_1.txt");
   myfile4.open("Gauss_sampling_2.txt");
   myfile5.open("Unif_int_eval_1.txt");
   myfile6.open("Unif_int_eval_2.txt");
   myfile7.open("Gauss_int_eval_1.txt");
   myfile8.open("Gauss_int_eval_2.txt");
   myfile9.open("Equilibration_unif_1.txt");
   myfile10.open("Equilibration_unif_2.txt");
   myfile11.open("Equilibration_gauss_1.txt");
   myfile12.open("Equilibration_gauss_2.txt");

   double L_1 = 2.5;
   double L_2 = 4.5;
   double sigma_1 = 0.7;
   double sigma_2 = 2.0;

   vector<vector<double> > unif_sampling_1;
   vector<vector<double> > unif_sampling_2;
   vector<vector<double> > gauss_sampling_1;
   vector<vector<double> > gauss_sampling_2;
   vector<double> v1;
   vector<double> v2;

   for( int i = 0; i<3; i++){
      v1.push_back(0.);
      v2.push_back(0.);
   }
   v2[2] = 1.;

   int N = 1000000;

   unif_sampling_1 = met.unif_sampling_3D(v1, N, L_1, &probability1);
   unif_sampling_2 = met.unif_sampling_3D(v2, N, L_2, &probability2);
   gauss_sampling_1 = met.gauss_sampling_3D(v1, N, sigma_1, &probability1);
   gauss_sampling_2 = met.gauss_sampling_3D(v2, N, sigma_2, &probability2);
   



   vector<double> unif_integral_1;
   vector<double> unif_integral_2;
   vector<double> gauss_integral_1;
   vector<double> gauss_integral_2;

   int time = 500;

   for (int i = 0; i<N; i++){
      unif_integral_1.push_back(fun.norm(unif_sampling_1[i]));
      unif_integral_2.push_back(fun.norm(unif_sampling_2[i]));
      gauss_integral_1.push_back(fun.norm(gauss_sampling_1[i]));
      gauss_integral_2.push_back(fun.norm(gauss_sampling_2[i]));
   }

   
   fun.print_on_file(fun.cut_end(unif_integral_1, time), myfile9);
   fun.print_on_file(fun.cut_end(unif_integral_2, time), myfile10);
   fun.print_on_file(fun.cut_end(gauss_integral_1, time), myfile11);
   fun.print_on_file(fun.cut_end(gauss_integral_2, time), myfile12);



   sys = fun.graph(fun.av(fun.cut_begin(unif_integral_1,time),100), fun.av2(fun.cut_begin(unif_integral_1, time),100));
   fun.print_on_file_errbar_expected(sys.sumprog, sys.errprog, myfile5, 1.5);

   sys = fun.graph(fun.av(fun.cut_begin(unif_integral_2,time),100), fun.av2(fun.cut_begin(unif_integral_2, time),100));
   fun.print_on_file_errbar_expected(sys.sumprog, sys.errprog, myfile6, 5.0);

   sys = fun.graph(fun.av(fun.cut_begin(gauss_integral_1,time),100), fun.av2(fun.cut_begin(gauss_integral_1, time),100));
   fun.print_on_file_errbar_expected(sys.sumprog, sys.errprog, myfile7, 1.5);

   sys = fun.graph(fun.av(fun.cut_begin(gauss_integral_2,time),100), fun.av2(fun.cut_begin(gauss_integral_2, time),100));
   fun.print_on_file_errbar_expected(sys.sumprog, sys.errprog, myfile8, 5.0);




   myfile1.close();
   myfile2.close();
   myfile3.close();
   myfile4.close();
   myfile5.close();
   myfile6.close();
   myfile7.close();
   myfile8.close();
   myfile9.close();
   myfile10.close();
   myfile11.close();
   myfile12.close();


   
   return 0;
}
