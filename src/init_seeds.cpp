#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "include/dist.h"
#include "include/init_seeds.h"
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
IntegerVector afk_mc2(const int k, const int m, const NumericMatrix M, const int threads)
{
  IntegerVector seeds(k);
  int n = M.nrow();
  seeds[0] = uniform_sample(n);
  std::vector<double> cpd = q_cpd(seeds[0], M, threads);

  int x, y;
  double d_x, d_y;
  for (int i = 1; i < k; i++) {
    x = prob_sample(cpd);
    d_x = min_dist(x, seeds, i, M);
    for (int j = 0; j < m; ++j) {
      y = prob_sample(cpd);
      d_y = min_dist(y, seeds, i, M);
      if ((d_y * q(x, cpd)) > frunif() * (d_x * q(y, cpd))) {
        x = y;
        d_x = d_y;
      }
    }
    seeds[i] = x;
  }
  return seeds;
}


double frunif()
{
  return (double)rand() / (double)RAND_MAX;
}

int uniform_sample(const int n)
{
  return (int)((double)n * frunif());
}

std::vector<double> q_cpd(const int s, const NumericMatrix M, const int threads)
{
  int nrow = M.nrow();
  Row seed = M.row(s);
  std::vector<double> res(nrow);
  double tot_dist = 0.0, tmp;
#pragma omp parallel num_threads(threads)
{
#pragma omp for reduction(+:tot_dist) private(tmp)
  for(int i = 0; i < nrow; ++i)
  {
    tmp = pair_dist(seed, M.row(i));
    res[i] = tmp;
    tot_dist += tmp;
  }

#pragma omp for
  for(int i = 0; i < nrow; ++i)
  {
    res[i] = 0.5 * res[i] / tot_dist + 0.5 / (double)nrow;
  }

#pragma omp single
{
  for (int i = 1; i < nrow; ++i)
  {
    res[i] += res[i-1];
  }
}
}
  return res;
}

double q(const int x, const std::vector<double> &cpd)
{
  if (x == 0) {
    return cpd[0];
  }
  return cpd[x] - cpd[x-1];
}


int prob_sample(std::vector<double> &cpd)
{
  double p = frunif();
  std::vector<double>::iterator loc = std::lower_bound(cpd.begin(), cpd.end(), p);
  return std::max(0, (int)(loc - cpd.begin()) - 1);
}
