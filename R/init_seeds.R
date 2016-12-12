#' @name seeding
#' @title finding initial cluster centers
#' @description
#' A fast seeding algorithm proposing to approximate k-means++
#' using Markov chain Monte Carlo (MCMC).
#' @usage
#' seeding(X, num_seeds, m, threads)
#' @param X a numeric matrix or data frame
#' @param num_seeds number of seeds
#' @param m an integer, length of Markove chains, which is used to sample centers, 20 by default
#' @param threads an integer, number of threads to speed up computing
#' @examples
#' data(iris)
#' iris <- iris[sample(1:nrow(iris), 5000, replace = T),]
#' X <- as.matrix(iris[,1:4])
#' seeds <- seeding(X, 3, 20, 2)
#' clus1 <- kmeans(X, X[seeds,]); table(clus1$cluster, iris[,5])
#' clus2 <- kmeans(X, 3); table(clus2$cluster, iris[,5])
seeding <- function(X, num_seeds, m, threads)
{
  if (is.data.frame(X)) {
    X <- as.matrix(X)
  }
  stopifnot(is.matrix(X))

  num_cases = nrow(X)
  if (num_seeds == 1) {
    seeds <- sample(1:num_cases, 1)
  } else {
    seeds <- afk_mc2(num_seeds, m, M, threads)
    seeds <- seeds + 1
  }
  seeds
}
