######################################################################
#' Viterbi algorithm to decode the latent states for Gaussian hidden
#' semi-Markov model with / without autoregressive structures
#' @param y observed series
#' @param mod list consisting the at least the following items: 
#' mod$m = scalar number of states, 
#' mod$delta = vector of initial values for prior probabilities, 
#' mod$gamma = matrix of initial values for state transition probabilies.
#' mod$mu = list of initial values for means, 
#' mod$sigma = list of initial values for covariance matrices.
#' mod$d = list of state duration probabilities.
#' For autoregressive hidden markov models, we also need the additional items:
#' mod$arp = scalar order of autoregressive structure
#' mod$auto = list of initial values for autoregressive coefficient matrices
#' @return a list containing the decoded states
#' @references Rabiner, Lawrence R. "A tutorial on hidden Markov models and 
#' selected applications in speech recognition." Proceedings of the 
#' IEEE 77.2 (1989): 257-286.
#' @examples 
#' set.seed(135)
#' m <- 2
#' mu <- list(c(3,4,5),c(-2,-3,-4))
#' sigma <- list(diag(1.3,3), 
#'             matrix(c(1,-0.3,0.2,-0.3,1.5,0.3,0.2,0.3,2),3,3,byrow=TRUE))
#' delta <- c(0.5,0.5)
#' gamma <- matrix(c(0,1,1,0),2,2,byrow=TRUE)
#' auto <- list(matrix(c(0.3,0.2,0.1,0.4,0.3,0.2,
#'                      -0.3,-0.2,-0.1,0.3,0.2,0.1,
#'                       0,0,0,0,0,0),3,6,byrow=TRUE),
#'             matrix(c(0.2,0,0,0.4,0,0,
#'                       0,0.2,0,0,0.4,0,
#'                      0,0,0.2,0,0,0.4),3,6,byrow=TRUE))
#' d <- list(c(0.5,0.3,0.2),c(0.6,0.4))
#' mod <- list(m=m,mu=mu,sigma=sigma,delta=delta,gamma=gamma,
#'            auto=auto,arp=2,d=d)
#' sim <- hsmm.sim(2000,mod)
#' y <- sim$series
#' state <- sim$state
#' fit <- em.semi(y=y, mod=mod, arp=2)
#' state_est <- viterbi.semi(y=y,mod=fit)
#' sum(state_est!=state)
#' @useDynLib rarhsmm, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom graphics points
#' @importFrom stats rnorm
#' @importFrom glmnet glmnet
#' @export
viterbi.semi <- function(y, mod){
  if(is.null(mod$auto)) result <- viterbi.semi.mvn(y, mod)
  if(!is.null(mod$auto)) result <- viterbi.semi.mvnarp(y, mod)
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
viterbi.semi.mvn <- function(y, mod){
  
  d <- mod$d
  ld <- sapply(d,length)
  D <- sum(ld)
  p <- length(mod$mu[[1]])
  Pi <- mod$delta
  P <- mod$gamma
  K <- mod$m
  n <- nrow(y)
  mu <- mod$mu
  sigma <- mod$sigma
  
  xi <- matrix(0, n, D)
  
  newP <- hsmm2hmm(P,d)
  newPi <- rep(NA, D)
  calc <- 0
  for(i in 1:K){
    for(j in 1:ld[i]){
      newPi[calc + j] <- Pi[i]/ld[i]
    }
    calc <- calc + ld[i]
  }
  
  B <- matrix(0,n,D)
  for(i in 1:n){
    calc <- 0
    for(j in 1:K){
      for(dur in 1:ld[j]){
        B[i,calc+dur] <- dmvn(y[i,],mu[[j]],sigma[[j]])
      }
      calc <- calc + ld[j]
    }
  }
  
  #########forward algorithm
  foo <- newPi * B[1,]
  
  xi[1,] <- foo / sum(foo)
  
  for(i in 2:n){
    foo <- apply(xi[i-1,] * newP, 2, max) * B[i,]
    xi[i,] <- foo/sum(foo)
  }
  
  iv <- numeric(n)
  iv[n] <- which.max(xi[n,])
  for(i in (n-1):1)
    iv[i] <- which.max(newP[,iv[i+1]] * xi[i,])
  
  original <- numeric(n)
  cumld <- cumsum(ld)
  for(i in 1:n){
    for(j in 1:K){
      if(iv[i]<=cumld[j]) {original[i] <- j; break}
    }
  }
  return(original)
}


###################
viterbi.semi.mvnarp <- function(y, mod){
  d <- mod$d
  ld <- sapply(d,length)
  D <- sum(ld)
  p <- length(mod$mu[[1]])
  Pi <- mod$delta
  P <- mod$gamma
  K <- mod$m
  n <- nrow(y)
  mu <- mod$mu
  auto <- mod$auto
  sigma <- mod$sigma
  arp <- mod$arp
  xi <- matrix(0, n, D)
  
  newP <- hsmm2hmm(P,d)
  newPi <- rep(NA, D)
  calc <- 0
  for(i in 1:K){
    for(j in 1:ld[i]){
      newPi[calc + j] <- Pi[i]/ld[i]
    }
    calc <- calc + ld[i]
  }
  
  B <- matrix(0,n,D)
  for(i in 1:n){
    calc <- 0
    for(j in 1:K){
      if(i==1)  diff <- mu[[j]] - y[i,]
      else if(i<=arp)diff <- mu[[j]] + auto[[j]][,(p*arp-(i-1)*p+1):(p*arp),drop=FALSE]%*%
          as.vector(t(y[(1):(i-1),])) - 
          y[i,]
      else diff <- mu[[j]] + auto[[j]]%*%
          as.vector(t(y[(i-arp):(i-1),])) - 
          y[i,]
      
      for(dur in 1:ld[j])
        B[i,calc+dur] <- dmvn(diff, rep(0,p), sigma[[j]])
          
      calc <- calc + ld[j]
    }
  }
  
  #########forward algorithm
  foo <- newPi * B[1,]
  
  xi[1,] <- foo / sum(foo)
  
  for(i in 2:n){
    foo <- apply(xi[i-1,] * newP, 2, max) * B[i,]
    xi[i,] <- foo/sum(foo)
  }
  
  iv <- numeric(n)
  iv[n] <- which.max(xi[n,])
  for(i in (n-1):1)
    iv[i] <- which.max(newP[,iv[i+1]] * xi[i,])
  
  original <- numeric(n)
  cumld <- cumsum(ld)
  for(i in 1:n){
    for(j in 1:K){
      if(iv[i]<=cumld[j]) {original[i] <- j; break}
    }
  }
  return(original)
}
