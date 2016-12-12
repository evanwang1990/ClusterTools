#include <Rcpp.h>
#include <cmath>
#include "include/setting.h"
#include "include/clustering.h"

using namespace std;
using namespace Rcpp;

ulonglong choose2(const ulonglong n);

//[[Rcpp::export]]
NumericVector ClusterCompare(const IntegerVector cluster_1, const IntegerVector cluster_2)
{
  IntegerVector cluster1 = clone(cluster_1);
  IntegerVector cluster2 = clone(cluster_2);
  const int num_cluster1 = AsClustering(cluster1);
  const int num_cluster2 = AsClustering(cluster2);

  int num_cases = cluster1.size();
  IntegerMatrix cluster_table(num_cluster1, num_cluster2);
  NumericMatrix perct_table(num_cluster1, num_cluster2);
  NumericVector num_pairs_c1(num_cluster1);
  NumericVector num_pairs_c2(num_cluster2);

  for (int i = 0; i < num_cases; ++i) {
    cluster_table(cluster1[i], cluster2[i]) ++;
  }

  int num_cell;
  ulonglong dsum = 0;
  for (int i = 0; i < num_cluster1; ++i) {
    for (int j = 0; j < num_cluster2; ++j) {
      num_cell = cluster_table(i, j);
      num_pairs_c1[i] += num_cell;
      num_pairs_c2[j] += num_cell;
      perct_table(i, j) = (double)num_cell / (double)num_cases;
      dsum += choose2((ulonglong)num_cell);
    }
  }

  ulonglong sum1 = 0;
  double pk1, h1 = 0.0;
  for (int i = 0; i < num_cluster1; ++i) {
    sum1 += choose2((ulonglong)num_pairs_c1[i]);
    pk1 = (double)num_pairs_c1[i]/(double)num_cases;
    num_pairs_c1[i] = pk1;
    h1 -= pk1 > 0 ? pk1 * std::log(pk1) : 0.0;
  }

  ulonglong sum2 = 0;
  double pk2, h2 = 0.0;
  for (int i = 0; i < num_cluster2; ++i) {
    sum2 += choose2((ulonglong)num_pairs_c2[i]);
    pk2 = (double)num_pairs_c2[i]/(double)num_cases;
    num_pairs_c2[i] = pk2;
    h2 -= pk2 > 0 ? pk2 * std::log(pk2) : 0.0;
  }

  ulonglong nsum = choose2((ulonglong)num_cases);
  // CR(corrected rand) can take values in [-1,1], where the value 1 indicates perfect agreement between the partitions,
  // whereas values near 0 (or negatives) correspond to cluster agreement found by chance. In fact, an
  // analysis by Milligan and Cooper (1986) confirmed that CR scores near 0 when presented to clusters
  // generated from random data, and showed that values lower than 0.05 indicate clusters achieved by chance.
  double adjusted_rand_idx = (double)(((ldouble)dsum - (ldouble)(sum1 * sum2) / (ldouble)nsum) / ((ldouble)(sum1 + sum2)/2 - (ldouble)(sum1*sum2)/(ldouble)nsum));
  double rand_idx = (double)((ldouble)(2 * dsum + nsum - sum1 -sum2) / (ldouble)nsum);
  double jaccard_idx = (double)((ldouble)dsum / (ldouble)(sum1 + sum2 - dsum));



  double mut_inf = 0.0;
  for (int i = 0; i < num_cluster1; ++i) {
    for (int j = 0; j < num_cluster2; ++j) {
      double& perct_cell = perct_table(i, j);
      if (perct_cell > 0) {
        mut_inf += perct_cell * std::log(perct_cell / (num_pairs_c1[i] * num_pairs_c2[j]));
      }
    }
  }

  double vi = h1 + h2 - 2 * mut_inf;
  double nmi = (h1==0 && h2==0) ? 1.0:2*mut_inf/(h1+h2);

  return NumericVector::create(jaccard_idx, rand_idx, adjusted_rand_idx, vi, nmi);
}


ulonglong choose2(const ulonglong n)
{
  return n*(n-1)/2;
}
