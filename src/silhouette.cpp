#include <Rcpp.h>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;
using namespace std;

NumericVector Silhouette (const NumericVector dist, const IntegerVector cluster, const int num_clusters)
{
  const int num_cases = cluster.size();
  NumericMatrix dists(num_cases, num_clusters);
  IntegerVector counts(num_clusters);
  int l = 0;
  int cluster_id1, cluster_id2;
  for (int i = 0; i < num_cases; ++i) {
    cluster_id1 = cluster[i];
    counts[cluster_id1] ++;
    for (int j = i+1; j < num_cases; ++j) {
      cluster_id2 = cluster[j];
      dists(i, cluster_id2) += dist[l];
      dists(j, cluster_id1) += dist[l];
      l ++;
    }
  }

  double ai = 0.0, bi, min_bi;
  IntegerVector neighbor(num_cases);
  neighbor.fill(-1);
  NumericVector si(num_cases);
  for (int i = 0; i < num_cases; ++i) {
    cluster_id1 = cluster[i];
    min_bi = R_PosInf;
    for (int j = 0; j < num_clusters; ++j) {
      if (cluster_id1 == j && counts[j] > 1) {
        ai = dists(i, j) / (counts[j] - 1);
      } else {
        bi = dists(i, j) / counts[j];
        if (bi < min_bi) {
          min_bi = bi;
          neighbor[i] = j;
        }
      }
    }
    si[i] = (ai == 0.0) ? 0.0:((min_bi - ai) / max(ai, min_bi));
  }
  return si;
}
