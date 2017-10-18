######################################################################
#' Simulate a Gaussian hidden semi-Markov series with / without autoregressive
#' structures
#' @param ns length of the simulated series
#' @param mod list consisting of at least the following items: 
#' mod$m = number of states, 
#' mod$delta = vector of prior probabilities, 
#' mod$gamma = matrix of state transition probabilies.
#' mod$mu = list of means, 
#' mod$sigma = list of covariance matrices.
#' mod$d = list of state duration probabilities.
#' For autoregressive hidden markov models, we also need the additional items:
#' mod$auto = list of autocorrelation matrices. 
#' mod$arp = order of autoregressive.
#' @return a list containing simulated series and states
#' @references Rabiner, Lawrence R. "A tutorial on hidden Markov models and 
#' selected applications in speech recognition." Proceedings of the 
#' IEEE 77.2 (1989): 257-286.
#' @examples
#' set.seed(351)
#'#Gaussian HSMM 3 hidden states (no autoregressive structure)
#' m <- 3
#' mu <- list(c(3),c(-2),c(0))
#' sigma <- list(as.matrix(1), as.matrix(0.8),as.matrix(0.3))
#' delta <- c(0.3,0.3,0.4)
#' gamma <- matrix(c(0,0.5,0.5,0.5,0,0.5,0.5,0.5,0),3,3,byrow=TRUE)
#' d <- list(c(0.4,0.3,0.2,0.1), c(0.5,0.25,0.25), c(0.7,0.3))
#' mod1 <- list(m=m,mu=mu,sigma=sigma,delta=delta,gamma=gamma,d=d)
#' sim1 <- hsmm.sim(500,mod1)
#' y1 <- sim1$series
#' fit1 <- em.semi(y=y1, mod=mod1)
#' 
#' \dontrun{
#' #AR(2) Gaussian HSMM with 3 hidden states
#' m <- 2
#' mu <- list(c(3,4,5),c(-2,-3,-4))
#' sigma <- list(diag(1.3,3), 
#'             matrix(c(1,-0.3,0.2,-0.3,1.5,0.3,0.2,0.3,2),3,3,byrow=TRUE))
#' delta <- c(0.5,0.5)
#' gamma <- matrix(c(0,1,1,0),2,2,byrow=TRUE)
#' d <- list(c(0.4,0.2,0.1,0.1,0.1,0.1),c(0.5,0.3,0.2))
#' auto <- list(matrix(c(0.3,0.2,0.1,0.4,0.3,0.2,
#'                      -0.3,-0.2,-0.1,0.3,0.2,0.1,
#'                       0,0,0,0,0,0),3,6,byrow=TRUE),
#'             matrix(c(0.2,0,0,0.4,0,0,
#'                       0,0.2,0,0,0.4,0,
#'                      0,0,0.2,0,0,0.4),3,6,byrow=TRUE))
#' mod2 <- list(m=m,mu=mu,sigma=sigma,delta=delta,gamma=gamma,
#'             auto=auto,arp=2,d=d)
#' sim2 <- hsmm.sim(2000,mod2)
#' y2 <- sim2$series
#' fit2 <- em.semi(y=y2, mod=mod2, arp=2)
#' }
#' @useDynLib rarhsmm, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom graphics points
#' @importFrom stats rnorm
#' @importFrom glmnet glmnet
#' @export
hsmm.sim <- function(ns, mod){
  if(is.null(mod$auto)) result <- mvn.hsmm.gen(ns, mod)
  if(!is.null(mod$auto)) result <- mvnarp.hsmm.gen(ns, mod)
  return(list(series=result$series, state=result$state))
}


###
mvn.hsmm.gen <- function(ns, mod){
  arp <- mod$arp
  d <- mod$d
  ld <- sapply(d,length)
  p <- length(mod$mu[[1]])
  mvect <- 1:mod$m
  state <- numeric(ns)
  x <- matrix(0, ns, p)
  #get the cholesky factor: sigma = L%*%L^t
  L <- lapply(1:mod$m, function(k) t(chol(mod$sigma[[k]])))
  
  state[1] <- sample(mvect, 1, prob=mod$delta)
  
  dur <- sample(1:ld[state[1]],1,prob=d[[state[1]]])
  
  for(tt in 1:dur){
    state[tt] <- state[1]
    z <- rnorm(p)
    x[tt,] <- L[[state[tt]]]%*%z + mod$mu[[state[tt]]] 
  }
  
  total <- dur
  
  while(total<ns){
    
    state[total+1] <- sample(mvect,1,prob=mod$gamma[state[total],])
    dur <- sample(1:ld[state[total+1]],1,prob=d[[state[total+1]]])
    for(tt in 1:dur){
      if(total+tt>ns) break
      state[total+tt] <- state[total+1]
      z <- rnorm(p)
      x[total+tt,] <- L[[state[total+tt]]]%*%z + mod$mu[[state[total+tt]]] 
    }
    total <- total + dur
  }
  
  return(list(series=x, state=state))
}

##########################
mvnarp.hsmm.gen <- function(ns, mod){
  arp <- mod$arp
  d <- mod$d
  ld <- sapply(d,length)
  p <- length(mod$mu[[1]])
  mvect <- 1:mod$m
  state <- numeric(ns)
  x <- matrix(0, ns, p)
  #get the cholesky factor: sigma = L%*%L^t
  L <- lapply(1:mod$m, function(k) t(chol(mod$sigma[[k]])))
  
  state[1] <- sample(mvect, 1, prob=mod$delta)
  
  dur <- sample(1:ld[state[1]],1,prob=d[[state[1]]])
  
  for(tt in 1:dur){
    state[tt] <- state[1]
    z <- rnorm(p)
    if(tt > arp){
      x[tt,] <- L[[state[tt]]]%*%z + mod$mu[[state[tt]]] + 
        mod$auto[[state[tt]]]%*%as.vector(t(x[(tt-arp):(tt-1),]))
      }else if(tt==1){x[tt,] <- L[[state[tt]]]%*%z + mod$mu[[state[tt]]] 
      }else{
          x[tt,] <- L[[state[tt]]]%*%z + mod$mu[[state[tt]]] + 
            mod$auto[[state[tt]]][,(p*arp-(tt-1)*p+1):(p*arp),drop=FALSE]%*%
            as.vector(t(x[1:(tt-1),]))
        }
  }
  
  total <- dur
  
  while(total<ns){
    state[total+1] <- sample(mvect,1,prob=mod$gamma[state[total],])
    dur <- sample(1:ld[state[total+1]],1,prob=d[[state[total+1]]])
    for(tt in 1:dur){
      if(total+tt>ns) break
      state[total+tt] <- state[total+1]
      z <- rnorm(p)
      if(tt+total > arp){
        x[total+tt,] <- L[[state[total+tt]]]%*%z + mod$mu[[state[total+tt]]] + 
          mod$auto[[state[total+tt]]]%*%as.vector(t(x[(total+tt-arp):(total+tt-1),]))
        }else{
            x[total+tt,] <- L[[state[total+tt]]]%*%z + mod$mu[[state[total+tt]]] + 
              mod$auto[[state[total+tt]]][,(p*arp-(total+tt-1)*p+1):(p*arp),drop=FALSE]%*%
              as.vector(t(x[1:(total+tt-1),]))
          }
      
    }
    total <- total + dur
  }
   
  return(list(series=x, state=state))
}
