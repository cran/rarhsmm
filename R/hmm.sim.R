######################################################################
#' Simulate a Gaussian hidden Markov series with / without autoregressive
#' structures
#' @param ns length of the simulated series
#' @param mod list consisting of at least the following items: 
#' mod$m = number of states, 
#' mod$delta = vector of prior probabilities, 
#' mod$gamma = matrix of state transition probabilies.
#' mod$mu = list of means, 
#' mod$sigma = list of covariance matrices.
#' For autoregressive hidden markov models, we also need the additional items:
#' mod$auto = list of autocorrelation matrices. 
#' mod$arp = order of autoregressive.
#' @return a list containing simulated series and states
#' @references Rabiner, Lawrence R. "A tutorial on hidden Markov models and 
#' selected applications in speech recognition." Proceedings of the 
#' IEEE 77.2 (1989): 257-286.
#' @examples
#' set.seed(135)
#'#Gaussian HMM 3 hidden states (no autoregressive structure)
#' m <- 3
#' mu <- list(c(3),c(-2),c(0))
#' sigma <- list(as.matrix(1), as.matrix(0.8),as.matrix(0.3))
#' delta <- c(0.3,0.3,0.4)
#' gamma <- matrix(c(0.8,0.1,0.1,0.1,0.8,0.1,0.1,0.1,0.8),3,3,byrow=TRUE)
#' mod1 <- list(m=m,mu=mu,sigma=sigma,delta=delta,gamma=gamma)
#' sim1 <- hmm.sim(1000,mod1)
#' y1 <- sim1$series
#' fit1 <- em.hmm(y=y1, mod=mod1)
#' 
#' #AR(2) Gaussian HMM with 3 hidden states
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
#' mod2 <- list(m=m,mu=mu,sigma=sigma,delta=delta,gamma=gamma,auto=auto,arp=2)
#' sim2 <- hmm.sim(2000,mod2)
#' y2 <- sim2$series
#' fit2 <- em.hmm(y=y2, mod=mod2, arp=2)
#' @useDynLib rarhsmm, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom graphics points
#' @importFrom stats rnorm
#' @importFrom glmnet glmnet
#' @export
hmm.sim <- function(ns, mod){
  if(is.null(mod$auto)) result <- mvn.hmm.gen(ns, mod)
  if(!is.null(mod$auto)) result <- mvnarp.hmm.gen(ns, mod)
  return(list(series=result$series, state=result$state))
}

######################################
mvn.hmm.gen <- function(ns, mod){
  d <- length(mod$mu[[1]])
  mvect <- 1:mod$m
  state <- numeric(ns)
  x <- matrix(0, ns, d)
  #get the cholesky factor: sigma = L%*%L^t
  L <- lapply(1:mod$m, function(k) t(chol(mod$sigma[[k]])))
    
  state[1] <- sample(mvect, 1, prob=mod$delta)
  z <- rnorm(d)
  x[1,] <- L[[state[1]]]%*%z + mod$mu[[state[1]]] 
  
  for(i in 2:ns) {
    state[i] <- sample(mvect,1,prob=mod$gamma[state[i-1],])
    z <- rnorm(d)
    x[i,] <- L[[state[i]]]%*%z + mod$mu[[state[i]]] 
  }
  
  return(list(series=x, state=state))
}

########################################
########################################
mvnarp.hmm.gen <- function(ns, mod){
  d <- length(mod$mu[[1]])
  mvect <- 1:mod$m
  arp <- mod$arp
  
  state <- numeric(ns)
  x <- matrix(0, ns, d)
  #get the cholesky factor: sigma = L%*%L^t
  L <- lapply(1:mod$m, function(k) t(chol(mod$sigma[[k]])))
  
  state[1] <- sample(mvect, 1, prob=mod$delta)
  z <- rnorm(d)
  x[1,] <- L[[state[1]]]%*%z + mod$mu[[state[1]]] 
  
  for(i in 2:ns) {
    state[i] <- sample(mvect,1,prob=mod$gamma[state[i-1],])
    z <- rnorm(d)
    if(i > arp){
      x[i,] <- L[[state[i]]]%*%z + mod$mu[[state[i]]] + 
        mod$auto[[state[i]]]%*%as.vector(t(x[(i-arp):(i-1),]))}else{
         
          x[i,] <- L[[state[i]]]%*%z + mod$mu[[state[i]]] + 
            mod$auto[[state[i]]][,(d*arp-(i-1)*d+1):(d*arp),drop=FALSE]%*%as.vector(t(x[1:(i-1),]))
        }
  }
  
  return(list(series=x, state=state))
}



