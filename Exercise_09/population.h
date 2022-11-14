#ifndef __Population__
#define __Population__
#include <vector>
#include "../Random/random.h"
using namespace std;

class Population {

private:

protected:

public:
  int n_cities, n_chromo;
  double beta;
  int type; //if type == 0 -> circle; else square
  Random class_rnd;
  vector<double> cities;
  vector<vector<double> > sq_cities;
  vector<vector<int> > population;
  vector<double> fitness;
  // constructors
  Population(int, int, Random, int);
  Population(int, int, Random, int, string);
  // destructor
  ~Population();
  // methods
  double Fit(vector<int>);
  double sq_Fit(vector<int>);
  void rand_swap(vector<int> &);
  void multi_swap(vector<int> &);
  void cont_Swap(vector<int> &);
  void inversion_swap(vector<int> &);
  void translation (vector<int> &);
  void anti_translation (vector<int> &);
  void cont_swap(vector<int> &);
  void mixing(void);
  void Init_cities(void);
  void Init_cities_input(string);
  void Init_pop(void);
  void Fill_fit(void);
  vector<size_t> sort_indices(vector<double>);
  void sort_population(void);
  void Mutation(void);
  vector<vector<int> > crossover(vector<int>, vector<int>,int,int);
  int getIndex(vector<int> v, int K);
  void check(int);
  


};

#endif // __Population__