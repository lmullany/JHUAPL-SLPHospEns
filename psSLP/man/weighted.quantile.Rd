% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quantile_functions.R
\name{weighted.quantile}
\alias{weighted.quantile}
\title{Generate weighted quantiles}
\usage{
weighted.quantile(x, probs, weights = NULL)
}
\arguments{
\item{x}{vector of values}

\item{probs}{vector of probabilities}

\item{weights}{vector of weights (must be >0)}
}
\description{
This function is the Harrell Davis weighted quantile estimator, and is implementing the
same as from cNORM package
}
\references{
Lenhard A, Lenhard W, Gary S (2018). “Continuous Norming (cNORM).” The Comprehensive R Archive Network. https://CRAN.R-project.org/package=cNORM.
}
