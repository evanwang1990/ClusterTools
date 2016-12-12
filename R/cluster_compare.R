#' @name compare
#' @title compare two clusters
#' @description
#' compare two clusters
#' @usage
#' compare(clustering, alt.clustering)
#' @param clustering an integer vector of length of the number of cases, which indicates a clustering
#' @param alt.clustering an integer vector such as for clustering, indicating an alternative clustering
#' @return jaccard
#' @return rand
#' @return ajusted.rand
#' @return vi
#' @return nmi
#' @examples
#' data(iris)
#' X <- scale(iris[,1:4])
#' c1 <- kmeans(X, 3)
#' c2 <- kmeans(X, 4)
#' compare(c1$cluster, c2$cluster)
compare <- function(clustering, alt.clustering)
{
  stopifnot(!missing(alt.clustering))
  stopifnot(length(clustering) == length(alt.clustering))
  stopifnot(length(clustering) > 1)
  if (is.character(clustering))
    clustering <- as.numeric(factor(clustering))
  if (is.character(alt.clustering))
    alt.clustering <- as.numeric(factor(alt.clustering))

  res <- ClusterCompare(clustering, alt.clustering)
  attr(res, "names") <- c("jaccard", "rand", "ajusted.rand", "vi", "nmi")
  res
}
