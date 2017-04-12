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
#' 
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
  calc <- 0
  for(i in 1:K){
    for(j in 1:ld[i]){
      newPi[calc + j] <- Pi[i]/ld[i]
    }
    calc <- calc + ld[i]
  }
  
  #initialization
  k1 <- (2*3.1415926)^(-p/2)
  
  #recursion
  loglik <- 0

    oldlik <- loglik
    
    #get the cholesky factor: sigma = L%*%L^t
    L <- lapply(1:K, function(k) t(chol(sigma[[k]])))
    Lt <- lapply(L,t)
    k2 <- sapply(1:K,function(k) k1/prod(diag(L[[k]]))) #sqrt(det(sigma))
    
    Gamma <- NULL
    Gammasum <- matrix(0,1,D)
    oldlik <- loglik
 
    loglik <- 0 
    
      alpha <- matrix(0,ns,D)
      beta <- matrix(0,ns,D)
      gamma <- matrix(0,ns,D)
      
      B <- matrix(0,ns,D)
      
      Scale <- rep(0, ns)
      Xi <- matrix(0, ns-1, D*D)
      
      #state-dependent probs
      for(i in 1:ns){
        calc <- 0
        for(j in 1:K){
          diff <- mu[[j]] - y[i,]
          for(dur in 1:ld[j]){
            B[i,calc+dur] <- k2[[j]]*exp(-0.5*t(diff)%*%
                                           backsolve(Lt[[j]],forwardsolve(L[[j]],diff)))
          }
          calc <- calc + ld[j]
        }
      }
      
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
      
      calc <- 0
      for(i in 1:K){
        for(j in 1:ld[i]){
          beta[ns,calc+j] <- 1/K/ld[i]
        }
        calc <- calc + ld[i]
      }
      
      for(i in (ns-1):1)
        beta[i,] <- (beta[i+1,] * B[i+1,]) %*% t(newP) / scale[i]
      
      gamma <- alpha*beta
      gamma <- gamma / rowSums(gamma) 
      gammasum <- colSums(gamma) #updated
      
      
      Scale <- Scale + log(scale)
      Gamma <- rbind(Gamma, gamma)
      Gammasum <- Gammasum + gammasum

      loglik <- loglik + sum(Scale)

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
  calc <- 0
  for(i in 1:K){
    for(j in 1:ld[i]){
      newPi[calc + j] <- Pi[i]/ld[i]
    }
    calc <- calc + ld[i]
  }
  
  #initialization
  k1 <- (2*3.1415926)^(-p/2)
  
  #recursion
  loglik <- 0

    oldlik <- loglik
    
    #get the cholesky factor: sigma = L%*%L^t
    L <- lapply(1:K, function(k) t(chol(sigma[[k]])))
    Lt <- lapply(L,t)
    k2 <- sapply(1:K,function(k) k1/prod(diag(L[[k]]))) #sqrt(det(sigma))
    count <- 0
    
    Gamma <- NULL
    Gammasum <- matrix(0,1,D)
    oldlik <- loglik
 
    
    #for each subject
    loglik <- 0 
    

      alpha <- matrix(0,ns,D)
      beta <- matrix(0,ns,D)
      gamma <- matrix(0,ns,D)
      
      B <- matrix(0,ns,D)
      
      Scale <- rep(0, ns)
      
      #state-dependent probs
      for(i in 1:ns){
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
            B[i,calc+dur] <- k2[[j]]*exp(-0.5*t(diff)%*%
                                           backsolve(Lt[[j]],forwardsolve(L[[j]],diff)))
          calc <- calc + ld[j]
        }
      }
      
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
      
      calc <- 0
      for(i in 1:K){
        for(j in 1:ld[i]){
          beta[ns,calc+j] <- 1/K/ld[i]
        }
        calc <- calc + ld[i]
      }
      
      for(i in (ns-1):(1))
        beta[i,] <- (beta[i+1,] * B[i+1,]) %*% t(newP) / scale[i]
      
      gamma <- alpha*beta
      gamma <- gamma / pmax(rowSums(gamma),0.01) 
      gammasum <- colSums(gamma) #updated
      
    
      Scale <- Scale + log(scale)
      Gamma <- rbind(Gamma, gamma)
      Gammasum <- Gammasum + gammasum

  
      loglik <- loglik + sum(Scale)
 
    oldGamma <- t(sapply(1:ns, function(k){ 
      sapply(split(Gamma[k,],rep(1:K,ld)),sum)}))
    #oldGammasum <- colSums(oldGamma)
    
  
  return(oldGamma)
}
