#include <Rcpp.h>
#include <math.h>
#include "include/dist.h"

using namespace Rcpp;
double pair_dist(const Row x1, const Row x2)
{
  double dist = 0.0, tmp;
  for (int i = 0, length = x1.size(); i < length; ++i)
  {
    tmp = x1[i] - x2[i];
    dist += tmp * tmp;
  }
  return std::sqrt(dist);
}


double min_dist(const int &x, const IntegerVector y, const int num_seeds, const NumericMatrix M)
{
  double res = R_PosInf;
  double tmp_dist;
  for (int i = 0; i < num_seeds; ++i) {
    tmp_dist = pair_dist(M.row(x), M.row(y[i]));
    if (res > tmp_dist) {
      res = tmp_dist;
    }
  }
  return res;
}
