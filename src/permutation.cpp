#include <Rcpp.h>
#include <algorithm>

using namespace Rcpp;
using namespace std;

class Rand
{
public:
  int MaxNum;

  Rand(int max_num) {
    MaxNum = max_num;
  }

  ~Rand() {}

  int operator()(int index) {
    return rand() % MaxNum;
  }
};

//[[Rcpp::export]]
NumericMatrix permutate(NumericMatrix x) {
  NumericMatrix x_ = clone(x);
  Rand rand_(int(x_.nrow()));
  for(int i = 0; i < x_.ncol(); ++i) {
    random_shuffle(x_.column(i).begin(), x_.column(i).end(), rand_);
  }
  return x_;
}
