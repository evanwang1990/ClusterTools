#include "include/clustering.h"
#include <utility>
using namespace std;
using namespace Rcpp;

template <class ForwardIterator, class OutIterator>
OutIterator unique(ForwardIterator first, ForwardIterator last, OutIterator result_first, OutIterator result_last)
{
  OutIterator p = result_first;
  *result_first = *first;
  p++;
  while (++first != last) {
    if (p==result_last) {
      break; // early stop
    }
    bool in_ = findIt(result_first, p, *first);
    if (!in_) {
      *p = *first;
      p++;
    }
  }
  return p;
}


template <typename T, class ForwardIterator>
bool findIt(ForwardIterator first, ForwardIterator last, T x)
{
  while (*first != x && first < last) {
    first ++;
  }
  return *first == x;
}


int AsClustering(IntegerVector clustering)
{
  pair<intit, intit> min_max = minmax(clustering.begin(), clustering.end());
  if (*min_max.first != *min_max.second) {
    vector<int> unique_values(*min_max.second - *min_max.first + 1, 0);
    vector<int>::iterator it = unique(clustering.begin(), clustering.end(), unique_values.begin(), unique_values.end());
    if (it == unique_values.end()) { // if clustering is continuous, make sure clustering begins with 0
      if (*min_max.first != 0) {
        double min_ = *min_max.first;
        for (intit p = clustering.begin(); p < clustering.end(); ++p) {
          *p -= min_;
        }
      }
      return (int)unique_values.size();
    } else { // if clustering is not continuous
      map<int, int> idx;
      for (int i = 0; i < (int)(it - unique_values.begin()); ++i) {
        idx[unique_values[i]] = i;
      }

      for (intit p = clustering.begin(); p < clustering.end(); ++p) {
        *p = idx[*p];
      }
      return (int)(it - unique_values.begin());
    }
  }
}
