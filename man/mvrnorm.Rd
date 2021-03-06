% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{mvrnorm}
\alias{mvrnorm}
\title{multivariate normal random number generator}
\usage{
mvrnorm(n, mu, sigma)
}
\arguments{
\item{n}{number of random vectors to generate}

\item{mu}{vector of means}

\item{sigma}{variance-covariance matrix}
}
\value{
a matrix with each row as a realization of multivariate normal random vector
}
\description{
multivariate normal random number generator
}
