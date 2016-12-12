#include <Rcpp.h>
#include "include/clustering.h"


using namespace Rcpp;
using namespace std;

extern NumericVector Silhouette(const NumericVector dist, const IntegerVector cluster, const int num_clusters);

//[[Rcpp::export]]
List BasicStat(const NumericVector dist, const IntegerVector cluster_)
{
  IntegerVector cluster = clone(cluster_);
  const int num_clusters = AsClustering(cluster);
  const int num_cases = cluster.size();
  IntegerVector cluster_size(num_clusters);
  NumericMatrix dist_sep(num_clusters, num_clusters); dist_sep.fill(R_PosInf); dist_sep.fill_diag(R_NegInf); // diag:diameter,other:separation
  NumericMatrix dist_avg(num_clusters, num_clusters); // diag:average distance within cluster,other:average distance between clusters
  NumericVector dist_avg_to_other(num_clusters); // average distance between one cluster and others
  double avg_dist_within_global = 0.0, avg_dist_between_global = 0.0;
  NumericVector avg_silwidths(num_clusters);
  double avg_silwidth_global = 0.0;
  double min_separation = R_PosInf, max_diameter = R_NegInf, dunn_idx; // dunn index
  double min_dist_avg_between = R_PosInf, max_dist_avg_within = R_NegInf, dunn2_idx; // dunn2 index
  double overall_ss = 0.0, within_ss = 0.0, ch_idx; NumericVector within_ss_(num_clusters);

  int l = 0, cluster_id1, cluster_id2;
  double dist_, ss_;
  NumericVector sil = Silhouette(dist, cluster, num_clusters);
  for (int i = 0; i < num_cases; i++) {
    cluster_id1 = cluster[i];
    cluster_size[cluster_id1] ++;
    avg_silwidths[cluster_id1] += sil[i];
    for (int j = i+1; j < num_cases; j++) {
      cluster_id2 = cluster[j];
      dist_ = dist[l];
      overall_ss += dist_ * dist_;
      dist_avg(cluster_id1, cluster_id2) += dist_;
      if (cluster_id1 == cluster_id2) {
        within_ss_[cluster_id1] += dist_ * dist_;
        if (dist_sep(cluster_id1, cluster_id1) < dist_) {
          dist_sep(cluster_id1, cluster_id1) = dist_;
        }
      } else {
        if (dist_sep(cluster_id1, cluster_id2) > dist_) {
          dist_sep(cluster_id1, cluster_id2) = dist_;
        }
      }
      l++;
    }
  }
  overall_ss /= num_cases;

  double dist_avg_;
  int tot_within_pairs = 0;
  for (int i = 0; i < num_clusters; i++) {
    avg_silwidth_global += avg_silwidths[i];
    avg_silwidths[i] /= cluster_size[i];
    within_ss += within_ss_[i]/cluster_size[i];
    for (int j = i; j < num_clusters; j++) {
      if (i == j) {
        avg_dist_within_global += dist_avg(i, i);
        tot_within_pairs += cluster_size(i)*(cluster_size(i)-1) / 2;
        dist_avg(i, i) = (cluster_size[i] == 1) ? 0.0:dist_avg(i, i)/(cluster_size(i)*(cluster_size(i)-1) / 2);
        if (max_diameter < dist_sep(i, i)) { max_diameter = dist_sep(i, i); }
        if (max_dist_avg_within < dist_avg(i, i)) { max_dist_avg_within = dist_avg(i, i); }
      } else {
        dist_sep(i, j) = dist_sep(j, i) = min(dist_sep(i, j), dist_sep(j, i));
        dist_avg_ = dist_avg(i, j) + dist_avg(j, i);
        dist_avg(i, j) = dist_avg(j, i) = dist_avg_/(cluster_size[i] * cluster_size[j]);
        dist_avg_to_other[i] += dist_avg_;
        dist_avg_to_other[j] += dist_avg_;
        avg_dist_between_global += dist_avg_;
        if (min_separation > dist_sep(i, j)) { min_separation = dist_sep(i, j); }
        if (min_dist_avg_between > dist_avg(i, j)) { min_dist_avg_between = dist_avg(i, j); }
      }
    }
    dist_avg_to_other[i] /= (cluster_size[i] * (num_cases - cluster_size[i]));
  }

  avg_dist_within_global /= tot_within_pairs;
  avg_dist_between_global /= ((num_cases-1)*num_cases/2 - tot_within_pairs);
  avg_silwidth_global /= num_cases;
  dunn_idx = max_diameter > 0 ? min_separation / max_diameter : NA_REAL;
  dunn2_idx = max_dist_avg_within > 0 ? min_dist_avg_between / max_dist_avg_within : NA_REAL;
  double wb_ratio = avg_dist_within_global / avg_dist_between_global;
  ch_idx = ((overall_ss - within_ss) / (double)(num_clusters - 1)) / (within_ss / (double)(num_cases - num_clusters));
  cout<<overall_ss - within_ss<<" "<<within_ss<<endl;
  return List::create(
    _["dist_avg"] = dist_avg,
    _["dist_sep"] = dist_sep,
    _["counts"] = cluster_size,
    _["dist_avg_to_other"] = dist_avg_to_other,
    _["avg_dist_within_global"] = avg_dist_within_global,
    _["avg_dist_between_global"] = avg_dist_between_global,
    _["dunn_idx"] = dunn_idx,
    _["dunn2_idx"] = dunn2_idx,
    _["wb_ratio"] = wb_ratio,
    _["ch_idx"] = ch_idx
  );
}
