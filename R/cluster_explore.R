#' @name explore
#' @title explore data for clustering
#' @description
#' explore data to detect whether it can be splited into proper clusters and find which variables can be used to cluster
#' @usage
#' explore(formula, data, classifier = = c("randomForest", "ranger"), control = NULL)
#' @examples
#' data(iris)
#' iris1 <- iris[iris$Species == "versicolor",1:4]
#' e1 <- explore(~., iris1, "randomForest")
#' summary(e1) # there are no clusters in the data
#' plot(e1)
#'
#' e2 <- explore(~., iris[,1:4], "ranger")
#' summary(e2)
#' plot(e2)
explore <- function(formula, data, classifier = c("randomForest", "ranger"), control = NULL)
{
  if (missing(classifier) || is.null(classifier)) {
    stop("classifier is missing")
  }
  if (missing(classifier)) {
    classifier <- "randomForest"
  } else {
    classifier <- match.arg(classifier)
  }
  if (!exists(classifier)) {
    loadLibrary <- tryCatch({do.call(library, list(package = classifier)); 100}, error = function(e) 0)
    if (loadLibrary == 0)
      stop(classifier, " is not supported")
  }

  if (is.character(formula)) {
    formula <- as.formula(formula)
  }

  mf <- model.frame(formula, data)
  mt <- terms(mf)
  attr(mt, "intercept") <- 0
  if (!is.null(model.response(mf))) {
    stop("there should have no response variable(s)")
  }

  matrix <- model.matrix(mt, data)
  pmatrix <- permutate(matrix)
  X <- rbind(matrix, pmatrix)
  Y <- c(rep(1, nrow(matrix)), rep(0, nrow(pmatrix)))
  if (classifier == "ranger") {
    X <- cbind(Y, X)
  }

  part.formula <- switch(classifier,
              "randomForest" = "(x = X, y = factor(Y))",
              "ranger"  = "(data = X, dependent.variable.name = 'Y', classification = TRUE, importance = 'impurity')"
              )

  m <- as.list(parse(text = paste(classifier, part.formula))[[1]])

  if (!(missing(control) || is.null(control))) {
    m <- c(m, control)
  }

  fit <- eval(as.call(m))

  res <- list()
  res$classifier <- classifier
  if (classifier == "randomForest") {
    importance <- fit$importance
    res$varImp <-  data.frame(variables = rownames(importance), MeanDecreaseGini = as.numeric(importance))
    res$confusion <- fit$confusion
  } else {
    importance <- fit$variable.importance
    res$varImp <- data.frame(variables = names(importance), MeanDecreaseGini = as.numeric(importance))
    cf <- table(Y, fit$predictions)
    attr(cf, "class") <- NULL
    rowsum <- apply(cf, 1, sum)
    cf <- cbind(cf, matrix(c(cf[1,2]/rowsum[1], cf[2,1]/rowsum[2]), ncol = 1))
    colnames(cf) <- c("0", "1", "class.error")
    rownames(cf) <- c("0", "1")
    res$confusion <- cf
  }

  res$chisq.test <- chisq.test(res$confusion[,1:2])
  res$chisq.test$data.name <- "confusion matrix"
  res$chisq.test$method <- "Pearson's Chi-squared test for Confusion Matrix"
  class(res) <- "cluster.explore"
  res
}

summary.cluster.explore <- function(obj)
{
  cat("\tClusters Existence Test")
  cat("\n\n\tConfusion matrix")
  cat("\n")
  print(obj$confusion)
  print(obj$chisq.test)
}

plot.cluster.explore <- function(obj)
{
  varImp <- obj$varImp
  varImp <- varImp[order(varImp[,2], decreasing = TRUE), ]
  gplt <- ggplot(data = varImp, aes(x = MeanDecreaseGini, y = reorder(variables, MeanDecreaseGini))) +
    geom_point() +
    ggtitle("Importance of Variables") +
    labs(y = "Variables") +
    xlim(c(0, max(varImp$MeanDecreaseGini)))
  gplt
}

select_var <- function(obj, threshold = 0.1)
{
  varImp <- obj$varImp
  maxGiniDecrease <- max(varImp$MeanDecreaseGini)
  as.character(varImp$variables[varImp$MeanDecreaseGini >= maxGiniDecrease * threshold])
}
