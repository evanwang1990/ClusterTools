cluster.test <- function(formula, data, classifier = c("randomForest", "ranger"), control = NULL)
{
  if (missing(classifier) || is.null(classifier)) {
    stop("classifier is missing")
  }
  classifier <- match.arg(classifier)
  if (!exists(classifier)) {
    loadLibrary <- tryCatch(do.call(library, list(pacakge = classifier)), error = function(e) 0)
    if (loadLibrary == 0)
      stop(classifier, "is not supported")
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
              "ranger"  = "(data = X, dependent.variable.name = 1, classification = True, importance = 'impurity')"
              )

  m <- as.list(parse(text = paste(classifier, part.formula))[[1]])

  if (!(missing(control) || is.null(control))) {
    m <- c(m, control)
  }

  fit <- eval(as.call(m))

  res <- list()
  res$classifier <- classifier
  if (classifier == "randomForest") {
    res$varImp <- fit$importance
    res$confusion <- fit$confusion
  } else {
    res$varImp <- fit$variable.importance
    cf <- fit$confusion.matrix
    rowsum <- apply(cf, 1, sum)
    cf <- cbind(cf, matrix(c(cf[1,1]/rowsum[1], cf[2,2]/rowsum[2]), ncol = 1))
    colnames(cf) <- c("0", "1", "class.error")
    rownames(cf) <- c("0", "1")
  }

  res$chisq.test <- chisq.test(fit$confusion[,1:2])
  class(res) <- "cluster.test"
  res
}

summary.cluster.test <- function(obj)
{
  cat("Clusters Existence Test")
  cat("\n\nConfusion matrix\n")
  print(obj$confusion)
  cat("\nChi-squared test of confusion matrix")
  chi <- obj$chisq.test
  cat("\nX-squared = ", chi$statistic, ", df = ", chi$parameter, ", p-value = ", chi$p.value, sep = "")
}

