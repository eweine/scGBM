% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/uncertainty.R
\name{get.se}
\alias{get.se}
\title{Uncertainty quantification}
\usage{
get.se(out)
}
\arguments{
\item{out}{a list that is the return value of \code{gbm.sc}}
}
\value{
Appends a two matrices \code{se_U} and \code{se_V} to
the list out containing the standard errors.
}
\description{
Perform uncertainty quantification in the low-dimensional representation
estimated by scGBM.
}
\author{
Phillip B. Nicol <philnicol740@gmail.com>
}
