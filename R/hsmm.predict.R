######################################################################
#' 1-step forward prediction for (autoregressive) Gaussian hidden semi-Markov model
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
#' @return 1-step forward state probabilities and forecasts
#' @references Rabiner, Lawrence R. "A tutorial on hidden Markov models and 
#' selected applications in speech recognition." Proceedings of the 
#' IEEE 77.2 (1989): 257-286.
#' @examples
#' set.seed(15562)
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
#' mod <- list(m=m,mu=mu,sigma=sigma,delta=delta,gamma=gamma,
#'             auto=auto,arp=2,d=d)
#' sim <- hsmm.sim(2000,mod)
#' y <- sim$series
#' state <- sim$state
#' fit <- em.semi(y=y, mod=mod, arp=2,tol=1e-5)
#' forecast <- hsmm.predict(y=y,mod=fit)
#' 
#' @export
hsmm.predict <- function(y, mod){
  if(is.null(mod$auto)) result <- predict.semi.mvn(y, mod)
  if(!is.null(mod$auto)) result <- predict.semi.mvnarp(y, mod)
  return(result)
}

############################################
#em algorithm for multivariate normal
predict.semi.mvn <- function(y, mod){
  
  d <- mod$d
  ld <- sapply(d,length)
  D <- sum(ld)
  
  ns <- nrow(y)
  p <- length(mod$mu[[1]])
  Pi <- mod$delta
  P <- mod$gamma
  K <- mod$m
  mu <- mod$mu
  sigma <- mod$sigma
  
  #expanded state HMM
  newP <- hsmm2hmm(P,d)
  newPi <- rep(NA, D)
  newmu <- vector(mode="list",length=D)
  newsigma <- vector(mode="list",length=D)
  calc <- 0
  for(i in 1:K){
    for(j in 1:ld[i]){
      newPi[calc + j] <- Pi[i]/ld[i]
      newmu[[calc+j]] <- mu[[i]]
      newsigma[[calc+j]] <- sigma[[i]]
    }
    calc <- calc + ld[i]
  }
  
  #########
  alpha <- matrix(0,ns,D)
  beta <- matrix(0,ns,D)
  Scale <- rep(0, ns)
  
  #state-dependent probs
  B <- getnodeprob_nocov_mvn(y, newmu, newsigma, D, p, 0, 0)
  
  ####E-step
  #forward-backward
  
  scale <- rep(0, ns)
  
  alpha[1,] <- newPi*B[1,]
  scale[1] <- sum(alpha[1,])
  alpha[1,] <- alpha[1,]/scale[1]
  
  for(i in 2:ns){
    alpha[i,] <- (alpha[i-1,]%*%newP)*B[i,]
    scale[i] <- sum(alpha[i,])
    alpha[i,] <- alpha[i,] / scale[i]
  }
  
  foo <- alpha[ns,]%*%newP
  forwardprob <- numeric(K)
  calc <- 0
  for(i in 1:K){
    for(j in 1:ld[i]){
      forwardprob[i] <- forwardprob[i] + foo[calc+j]
    }
    calc <- calc + ld[i]
  }
  
  forecast <- lapply(1:K, function(k) mu[[k]] * forwardprob[k])
  forecast <- Reduce(`+`,forecast)
  return(list(forwardprob=forwardprob,forecast=forecast))
}

################################################
predict.semi.mvnarp <- function(y, mod){
  
  d <- mod$d
  ld <- sapply(d,length)
  D <- sum(ld)
  
  ns <- nrow(y)
  arp <- mod$arp
  p <- length(mod$mu[[1]])
  Pi <- mod$delta
  P <- mod$gamma
  K <- mod$m
  mu <- mod$mu
  sigma <- mod$sigma
  auto <- mod$auto
  
  #expanded state HMM
  #expanded state HMM
  newP <- hsmm2hmm(P,d)
  newPi <- rep(NA, D)
  newmu <- vector(mode="list",length=D)
  newsigma <- vector(mode="list",length=D)
  newauto <- vector(mode="list",length=D)
  calc <- 0
  for(i in 1:K){
    for(j in 1:ld[i]){
      newPi[calc + j] <- Pi[i]/ld[i]
      newmu[[calc+j]] <- mu[[i]]
      newsigma[[calc+j]] <- sigma[[i]]
      newauto[[calc+j]] <- auto[[i]]
    }
    calc <- calc + ld[i]
  }
  
  #initialization
  alpha <- matrix(0,ns,D)
  beta <- matrix(0,ns,D)
  Scale <- rep(0, ns)
  
  #state-dependent probs
  ycov <- rbind(rep(0,p),y[-ns,])
  autoarray <- array(as.numeric(unlist(newauto)), dim=c(p,p,D))
  muarray <- array(as.numeric(unlist(newmu)), dim=c(1,p,D))
  B <- getnodeprob_part2(y, ycov,autoarray,muarray,newsigma,D,p)
  
  ####E-step
  #forward-backward
  
  #alpha[1:arp,] <- rep(0,D) ########
  #beta[1:arp,] <- rep(0,D)
  
  scale <- rep(0, ns)
  #scale[1:arp] <- 0 ########
  
  alpha[1,] <- newPi*B[1,]
  scale[1] <- sum(alpha[1,])
  alpha[1,] <- alpha[1,]/scale[1]
  
  for(i in (2):ns){
    alpha[i,] <- (alpha[i-1,]%*%newP)*B[i,]
    scale[i] <- sum(alpha[i,])
    alpha[i,] <- alpha[i,] / scale[i]
  }
  
  foo <- alpha[ns,]%*%newP
  forwardprob <- numeric(K)
  calc <- 0
  for(i in 1:K){
    for(j in 1:ld[i]){
      forwardprob[i] <- forwardprob[i] + foo[calc+j]
    }
    calc <- calc + ld[i]
  }
  
  forecast <- lapply(1:K, function(k) (mu[[k]] + auto[[k]]%*%
                    as.vector(t(y[(ns-arp+1):(ns),]))) * forwardprob[k])
  forecast <- Reduce(`+`,forecast)
  return(list(forwardprob=forwardprob,forecast=forecast))
}
