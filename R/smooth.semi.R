######################################################################
#' Calculate the probability of being in a particular state for each observation.
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
#' @return a matrix containing the state probabilities
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
#' fit <- em.semi(y=y, mod=mod, arp=2)
#' stateprob <- smooth.semi(y=y,mod=fit)
#' head(cbind(state,stateprob),20)
#' @useDynLib rarhsmm, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom graphics points
#' @importFrom stats rnorm
#' @importFrom glmnet glmnet
#' @export
smooth.semi <- function(y, mod){
  if(is.null(mod$auto)) result <- smooth.semi.mvn(y, mod)
  if(!is.null(mod$auto)) result <- smooth.semi.mvnarp(y, mod)
  return(result)
}

############################################
#em algorithm for multivariate normal
smooth.semi.mvn <- function(y, mod){
  
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
      #state-dependent probs
  nodeprob <- getnodeprob_nocov_mvn(y, newmu, newsigma, D, p, 0, 0)
  fb <- forwardbackward(newPi,newP,nodeprob,ns,ns)
  Gamma <- fb$Gamma

    oldGamma <- t(sapply(1:ns, function(k){ 
      sapply(split(Gamma[k,],rep(1:K,ld)),sum)}))
    #oldGammasum <- colSums(oldGamma)
    
  return(oldGamma)
}

################################################
smooth.semi.mvnarp <- function(y, mod){
  
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
  
      #state-dependent probs
  ycov <- rbind(rep(0,p),y[-ns,])
  autoarray <- array(as.numeric(unlist(newauto)), dim=c(p,p,D))
  muarray <- array(as.numeric(unlist(newmu)), dim=c(1,p,D))
  nodeprob <- getnodeprob_part2(y, ycov,autoarray,muarray,newsigma,D,p)
  
  fb <- forwardbackward(newPi,newP,nodeprob,ns,ns)
  Gamma <- fb$Gamma
 
    oldGamma <- t(sapply(1:ns, function(k){ 
      sapply(split(Gamma[k,],rep(1:K,ld)),sum)}))
    #oldGammasum <- colSums(oldGamma)
    
  
  return(oldGamma)
}
