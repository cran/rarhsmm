#' \tabular{ll}{
#' Package: \tab rarhsmm \cr
#' Type: \tab Package\cr
#' Version: \tab 1.0.7\cr
#' Date: \tab 2018-03-19\cr
#' License: \tab GPL \cr
#' LazyLoad: \tab yes\cr
#' LazyData: \tab yes\cr
#' }
#'
#' @author Zekun Xu \email{zekunxu@gmail.com}
#' @author Ye Liu \email{yliu87@ncsu.edu}
#' Maintainer: Zekun Xu \email{zekunxu@gmail.com}
#' @name package-rarhsmm
#' @aliases rarhsmm-package
#' @docType package
#' @title Regularized Autoregressive Hidden Semi Markov Models
#' @keywords hidden semi-Markov models, hidden Markov models, autoregressive, regularized
NULL

##########################################################

#' @useDynLib rarhsmm, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom graphics points
#' @importFrom stats rnorm
#' @importFrom stats lm
#' @importFrom stats coef
#' @importFrom glmnet glmnet
##coef 
##predict
#generalized logit and inverse logit function
#1 / (1+sum(exp(x[-1]))) = p1
#exp(x[k]) / (1+sum(exp(x[-1]))) = pk
frobenius <- function(A,B){
  vec <- as.vector(A-B)^2
  sqrt(sum(vec))
}

glogit <- function(p){
  k <- length(p) - 1
  if(k==0) {x <- log(p) - log(1-p)}else{
  x <- rep(NA, k)
  for(j in 1:k) x[j] <- log(p[j+1]) - log(p[1])}
  return(x)
}


ginvlogit <- function(x){
  k <- length(x) + 1
  p <- rep(NA,k)
  den <- 1+sum(exp(x))
  p[1] <- 1 / den
  for(j in 2:k) p[j] <- exp(x[j-1]) / den
  return(p)
} 






