#include <Rcpp.h>
// #include "include/parallel_sort.hpp"
using namespace Rcpp;
using namespace std;

typedef unsigned long long ull;


NumericVector G2G3(NumericVector dist, const IntegerVector cluster, int num_dist_within = 0, bool g2 = true, bool g3 = false)
{
  NumericVector g2g3(2);
  g2g3.fill(NA_REAL);

  if (!g2 && !g3) {
    return g2g3;
  }

  const int num_cases = cluster.size();
  const int num_dist = dist.size();
  if (num_dist_within == 0) {
    for (int i = 0; i < num_cases; ++i) {
      for (int j = i+1; j < num_cases; ++j) {
        if (cluster[i] == cluster[j]) {
          num_dist_within ++;
        }
      }
    }
  }

  vector<double> dist_within(num_dist_within);
  vector<double> dist_between(num_dist - num_dist_within);
  vector<double>::iterator p1 = dist_within.begin();
  vector<double>::iterator p2 = dist_between.begin();
  int l = 0;
  for (int i = 0; i < num_cases; ++i) {
    for (int j = i+1; j < num_cases; ++j) {
      if (cluster[i] == cluster[j]) {
        *p1 = dist[l];
        p1++;
      } else {
        *p2 = dist[l];
        p2++;
      }
      l++;
    }
  }

  double s1 = 0, s2 = 0, s3 = 0;
  for (int i = 0; i < num_dist_within; i++) {
    s1 += dist_within[i];
  }
  for (int i = 0; i < num_dist; i++) {
    s2 += dist[i];
    s3 += dist[i] * dist[i];
  }

  double pearson = ((double)(num_dist) * s1 - (double)(num_dist_within) * s2) / (sqrt((double)(num_dist) * s3 - s2 * s2) * sqrt((double)(num_dist) * (double)(num_dist_within) - (double)(num_dist_within) * (double)(num_dist_within)));

  if (g2) {
    sort(dist_within.begin(), dist_within.end());
    sort(dist_between.begin(), dist_between.end());

    ull splus = 0, sminus = 0;
    int tmp_sequal = 0, tmp_sminus = 0;
    p1 = dist_within.begin();
    p2 = dist_between.begin();
    for (int i = 0; i < num_dist_within; ++i) {
      while(*p1 >= *p2) {
        if (*p1 > *p2) {
          tmp_sminus ++;
          tmp_sequal = 0;
        } else {
          tmp_sequal ++;
        }
        p2 ++;
      }
      sminus += tmp_sminus;
      splus += (num_dist - num_dist_within) - (tmp_sminus + tmp_sequal);
      p1 ++;
      if (*p1 > *p2) {
        tmp_sminus += tmp_sequal;
      }
    }
    g2g3[0] = double(splus - sminus) / double(splus + sminus);
    if (!g3) {
      return g2g3;
    }
  }

  if (g3) {
    double dmin = 0, dmax = 0;
    nth_element(dist.begin(), dist.begin()+num_dist_within, dist.end());
    for (int i = 0; i < num_dist_within; ++i) {
      dmin += dist[i];
    }
    nth_element(dist.begin(), dist.end()-num_dist_within, dist.end());
    for (int i = num_dist - num_dist_within; i < num_dist; ++i) {
      dmax += dist[i];
    }
    g2g3[1] = (s1 - dmin) / (dmax - dmin);
    if (!g2) {
      return g2g3;
    }
  }

  return g2g3;
}
