#ifndef INIT_SEEDS_H_
#define INIT_SEEDS_H_

#include <Rcpp.h>
#include <vector>

extern double frunif();

extern std::vector<double> q_cpd(const int s, const Rcpp::NumericMatrix M, const int threads);
extern double q(const int x, const std::vector<double> &cpd);

extern int uniform_sample(const int n);
extern int prob_sample(std::vector<double> &cpd);


#endif
