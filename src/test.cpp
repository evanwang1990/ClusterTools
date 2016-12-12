#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;
using namespace std;

void double_(double& x)
{
  x *= 2.0;
}

//[[Rcpp::export]]
NumericVector test1(NumericVector x)
{
  for (NumericVector::iterator p = x.begin(); p < x.end(); ++p)
    *p *= 2.0;
  return x;
}

//[[Rcpp::export]]
NumericVector test2(NumericVector x)
{
  for_each(x.begin(), x.end(), double_);
  return x;
}
