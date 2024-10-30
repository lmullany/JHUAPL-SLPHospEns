# This is the weighted quantile generic function from cNORM
# citation: Lenhard A, Lenhard W, Gary S (2018). “Continuous Norming (cNORM).” The Comprehensive R Archive Network. https://CRAN.R-project.org/package=cNORM.

wquantile.generic <- function (x, probs, cdf.gen, weights = NULL)
{
  if (is.null(weights)) {
    return(quantile(x, probs = probs))
  }
  n <- length(x)
  if (any(is.na(weights)))
    weights <- rep(1/n, n)
  nw <- sum(weights)/max(weights)
  indexes <- order(x)
  x <- x[indexes]
  weights <- weights[indexes]
  weights <- weights/sum(weights)
  cdf.probs <- cumsum(c(0, weights))
  sapply(probs, function(p) {
    cdf <- cdf.gen(nw, p)
    q <- cdf(cdf.probs)
    w <- tail(q, -1) - head(q, -1)
    sum(w * x)
  })
}

# this is the simple harrell davis
# again, implementation is from cNORM package

#' Generate weighted quantiles
#'
#' This function is the Harrell Davis weighted quantile estimator, and is implementing the
#' same as from cNORM package
#'
#' @param x vector of values
#' @param probs vector of probabilities
#' @param weights vector of weights (must be >0)
#' @references Lenhard A, Lenhard W, Gary S (2018). “Continuous Norming (cNORM).” The Comprehensive R Archive Network. https://CRAN.R-project.org/package=cNORM.
#' @export

weighted.quantile <- function (x, probs, weights = NULL)
{
  if (is.null(weights)) {
    return(quantile(x, probs = probs))
  }
  cdf.gen <- function(n, p) return(function(cdf.probs) {
    pbeta(cdf.probs, (n + 1) * p, (n + 1) * (1 - p))
  })
  wquantile.generic(x, probs, cdf.gen, weights)
}
