######################################################################
#' Viterbi algorithm to decode the latent states for Gaussian hidden Markov
#' model with / without autoregressive structures
#' @param y observed series
#' @param mod list consisting the at least the following items: 
#' mod$m = scalar number of states, 
#' mod$delta = vector of initial values for prior probabilities, 
#' mod$gamma = matrix of initial values for state transition probabilies.
#' mod$mu = list of initial values for means, 
#' mod$sigma = list of initial values for covariance matrices.
#' For autoregressive hidden markov models, we also need the additional items:
#' mod$arp = scalar order of autoregressive structure
#' mod$auto = list of initial values for autoregressive coefficient matrices
#' @references Rabiner, Lawrence R. "A tutorial on hidden Markov models and 
#' selected applications in speech recognition." Proceedings of the 
#' IEEE 77.2 (1989): 257-286.
#' @return a list containing the decoded states
#' @examples 
#' set.seed(135)
#' m <- 2
#' mu <- list(c(3,4,5),c(-2,-3,-4))
#' sigma <- list(diag(1.3,3), 
#'             matrix(c(1,-0.3,0.2,-0.3,1.5,0.3,0.2,0.3,2),3,3,byrow=TRUE))
#' delta <- c(0.5,0.5)
#' gamma <- matrix(c(0.8,0.2,0.1,0.9),2,2,byrow=TRUE)
#' auto <- list(matrix(c(0.3,0.2,0.1,0.4,0.3,0.2,
#'                      -0.3,-0.2,-0.1,0.3,0.2,0.1,
#'                       0,0,0,0,0,0),3,6,byrow=TRUE),
#'             matrix(c(0.2,0,0,0.4,0,0,
#'                       0,0.2,0,0,0.4,0,
#'                      0,0,0.2,0,0,0.4),3,6,byrow=TRUE))
#' mod <- list(m=m,mu=mu,sigma=sigma,delta=delta,gamma=gamma,auto=auto,arp=2)
#' sim <- hmm.sim(2000,mod)
#' y <- sim$series
#' state <- sim$state
#' fit <- em.hmm(y=y, mod=mod, arp=2)
#' state_est <- viterbi.hmm(y=y,mod=fit)
#' sum(state_est!=state)
#' @useDynLib rarhsmm, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom graphics points
#' @importFrom stats rnorm
#' @importFrom glmnet glmnet
#' @export
viterbi.hmm <- function(y, mod){
  if(is.null(mod$auto)) result <- viterbi.mvn(y, mod)
  if(!is.null(mod$auto)) result <- viterbi.mvnarp(y, mod)
  return(result)
}

dmvn <- function(y,mu,sigma){
  p <- length(mu)
  k1 <- (2*3.1415926)^(-p/2)
  L <- t(chol(sigma))
  Lt <- t(L)
  k2 <- k1/prod(diag(L))#sqrt(det(sigma))
  diff <- y - mu
  density <- k2*exp(-0.5*t(diff)%*%
                backsolve(Lt,forwardsolve(L,diff)))
  return(density)
}

#####################################
viterbi.mvn <- function(y, mod){
  m <- mod$m
  n <- nrow(y)
  mu <- mod$mu
  sigma <- mod$sigma
  
  xi <- matrix(0, n, m)
  
  #########forward algorithm
  foo <- numeric(m)
  for(j in 1:m)
     foo[j] <- mod$delta[j] * dmvn(y[1,],mu[[j]],sigma[[j]])
  
  xi[1,] <- foo / sum(foo)
  
  for(i in 2:n){
    for(j in 1:m)
       foo[j] <- apply(xi[i-1,] * mod$gamma, 2, max)[j] * 
              dmvn(y[i,],mu[[j]],sigma[[j]])
    
    xi[i,] <- foo/sum(foo)
  }
  iv <- numeric(n)
  iv[n] <- which.max(xi[n,])
  for(i in (n-1):1)
    iv[i] <- which.max(mod$gamma[,iv[i+1]] * xi[i,])
  return(iv)
}


###################
viterbi.mvnarp <- function(y, mod){
  m <- mod$m
  n <- nrow(y)
  mu <- mod$mu
  sigma <- mod$sigma
  arp <- mod$arp
  d <- length(mu[[1]])
  
  xi <- matrix(0, n, m)
  
  #########forward algorithm
  foo <- numeric(m)
  for(j in 1:m)
    foo[j] <- mod$delta[j] * dmvn(y[1,],mu[[j]],sigma[[j]])
  
  xi[1,] <- foo / sum(foo)
  
  for(i in 2:n){
    for(j in 1:m){
      if(i<=arp)foo[j] <- apply(xi[i-1,] * mod$gamma, 2, max)[j] * 
          dmvn(y[i,],mu[[j]]+mod$auto[[j]][,(d*arp-(i-1)*d+1):(d*arp),drop=FALSE]%*%
                 as.vector(t(y[1:(i-1),])),sigma[[j]])
      
      else foo[j] <- apply(xi[i-1,] * mod$gamma, 2, max)[j] * 
          dmvn(y[i,],mu[[j]]+mod$auto[[j]]%*%
                 as.vector(t(y[(i-arp):(i-1),])),sigma[[j]])
    }
    xi[i,] <- foo/sum(foo)
  }
  iv <- numeric(n)
  iv[n] <- which.max(xi[n,])
  for(i in (n-1):1)
    iv[i] <- which.max(mod$gamma[,iv[i+1]] * xi[i,])
  return(iv)
}


