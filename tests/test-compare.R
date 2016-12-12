test_that("check all the metrics to compare two clusterings are correct",{
  library(igraph)
  library(clusteval)
  data(ruspini)
  pr4 <- pam(ruspini, 4)
  pr3 <- pam(ruspini, 3)
  eval1 <- ClusterCompare(cluster1 = pr4$clustering-1, 4, cluster2 = pr3$clustering-1,3)


  as.integer(factor(x)), AsClustering(x)
})
