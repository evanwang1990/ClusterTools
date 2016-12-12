#ifndef CLUSTER_TOOLS_CLUSTERING_H
#define CLUSTER_TOOLS_CLUSTERING_H

#include <Rcpp.h>
#include <map>
#include "setting.h"

template <typename T, class ForwardIterator>
bool findIt(ForwardIterator first, ForwardIterator last, T x);
template <class ForwardIterator, class OutIterator>
OutIterator unique(ForwardIterator first, ForwardIterator last, OutIterator result_begin, OutIterator result_end);
extern int AsClustering(Rcpp::IntegerVector clustering);

#endif
