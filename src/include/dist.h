#ifndef CLUSTER_TOOLS_DIST_H_
#define CLUSTER_TOOLS_DIST_H_

#include <Rcpp.h>

typedef Rcpp::NumericMatrix::ConstRow Row;

extern double pair_dist(const Row x1, const Row x2);
extern double min_dist(const int &x, const Rcpp::IntegerVector y, const int num_seeds, const Rcpp::NumericMatrix M);

#endif
