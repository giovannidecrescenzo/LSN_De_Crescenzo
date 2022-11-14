#ifndef __Metropolis__
#define __Metropolis__
#include<vector>
#include"../Random/random.h"
#include"../Functions/functions.h"

using namespace std;

class Metropolis {

private:

protected:

public:
  // constructors
  Metropolis();
  // destructor
  ~Metropolis();
  // methods
  double acceptance;

  vector<vector<double> > unif_sampling_3D(vector<double> x, int N, double L, function<double(vector<double>)> func);
  vector<vector<double> > gauss_sampling_3D(vector<double> x, int N, double sigma, function<double(vector<double>)> func);

};

#endif // __Random__