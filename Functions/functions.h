#ifndef __Functions__
#define __Functions__
#include<vector>
#include"analysis.h"
using namespace std;
class Functions {

private:

protected:

public:
  // constructors
  Functions();
  // destructor
  ~Functions();
  // methods

  double error(vector<double> AV, vector<double> AV2, int n);
  double GetMean(double *A, int size);
  double GetN(double *A, int n);
  double norm(vector<double> x);
  double average(vector<double>);
  vector<double> av(vector<double> x, int N);
  vector<double> av2(vector<double> x, int N);
  Analysis graph(vector<double> x, vector<double> y);
  void print_on_file_errbar_expected(vector<double> x, vector<double> y, ostream& fout, double exp);
  void print_on_file_errbar(vector<double> x, vector<double> y, ostream& fout);
  void print_on_file(vector<double> x, ostream& fout);
  vector<double> cut_end(vector<double> x, int N);
  vector<double> cut_begin(vector<double> x, int N);
};

#endif // __Random__