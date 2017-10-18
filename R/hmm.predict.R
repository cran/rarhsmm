######################################################################
#' 1-step forward prediction for (autoregressive) Gaussian hidden Markov model
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
#' fit <- em.hmm(y=y, mod=mod, arp=2, tol=1e-5)
#' forecast <- hmm.predict(y=y,mod=fit)
#' 
#' @export
hmm.predict <- function(y, mod){
  if(is.null(mod$auto)) result <- predict.mvn(y, mod)
  if(!is.null(mod$auto)) result <- predict.mvnarp(y, mod)
  return(result)
}

############################################
#em algorithm for multivariate normal
predict.mvn <- function(y, mod){
  
  ns <- nrow(y)
  p <- length(mod$mu[[1]])
  Pi <- mod$delta
  P <- mod$gamma
  K <- mod$m
  mu <- mod$mu
  sigma <- mod$sigma
  
  
  #initialization
  alpha <- matrix(0,ns,K)
  beta <- matrix(0,ns,K)
  Scale <- rep(0, ns)
  
  #state-dependent probs
  B <- getnodeprob_nocov_mvn(y, mu, sigma, K, p, 0, 0)
  
  ####E-step
  #forward-backward
  
  scale <- rep(0, ns)
  alpha[1,] <- Pi*B[1,]
  scale[1] <- sum(alpha[1,])
  alpha[1,] <- alpha[1,]/scale[1]
  
  for(i in 2:ns){
    alpha[i,] <- (alpha[i-1,]%*%P)*B[i,]
    scale[i] <- sum(alpha[i,])
    alpha[i,] <- alpha[i,] / scale[i]
  }
  
  forwardprob <- alpha[ns,]%*%P
  forecast <- lapply(1:K, function(k) mu[[k]] * forwardprob[k])
  forecast <- Reduce(`+`,forecast)
  return(list(forwardprob=forwardprob,forecast=forecast))
}

################################################
predict.mvnarp <- function(y, mod){
  
  ns <- nrow(y)
  p <- length(mod$mu[[1]])
  Pi <- mod$delta
  P <- mod$gamma
  K <- mod$m
  mu <- mod$mu
  sigma <- mod$sigma
  auto <- mod$auto
  arp <- mod$arp
  
  #initialization
  alpha <- matrix(0,ns,K)
  beta <- matrix(0,ns,K)
  gamma <- matrix(0,ns,K)
  B <- matrix(0,ns,K)
  
  Scale <- rep(0, ns)
  
  #state-dependent probs
  ycov <- rbind(rep(0,p),y[-ns,])
  autoarray <- array(as.numeric(unlist(auto)), dim=c(p,p,K))
  muarray <- array(as.numeric(unlist(mu)), dim=c(1,p,K))
  B <- getnodeprob_part2(y, ycov,autoarray,muarray,sigma,K,p)
  
  ####E-step
  #forward-backward
  #alpha[1:arp,] <- rep(0,K) ########
  #beta[1:arp,] <- rep(0,K)
  
  scale <- rep(0, ns)
  #scale[1:arp] <- 0 ########
  
  alpha[1,] <- Pi*B[1,]
  scale[1] <- sum(alpha[1,])
  alpha[1,] <- alpha[1,]/scale[1]
  
  for(i in 2:ns){
    alpha[i,] <- (alpha[i-1,]%*%P)*B[i,]
    scale[i] <- sum(alpha[i,])
    alpha[i,] <- alpha[i,] / scale[i]
  }
  
  forwardprob <- alpha[ns,]%*%P
  forecast <- lapply(1:K, function(k) (mu[[k]] + auto[[k]]%*%
                       as.vector(t(y[(ns-arp+1):(ns),]))) * forwardprob[k])
  forecast <- Reduce(`+`,forecast)
  return(list(forwardprob=forwardprob,forecast=forecast))
}
