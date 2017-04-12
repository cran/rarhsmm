######################################################################
#' EM algorithm to compute maximum likelihood estimate of Gaussian 
#' hidden Markov models with / without autoregressive structures and 
#' with / without regularization on the covariance matrices and/or
#' autoregressive structures.
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
#' @param ntimes length of each homogeneous time series. Default to NULL,
#' which means only homogeneous time series.
#' @param tol tolerance for relative change. Default to 1e-4.
#' @param maxit maximum number of iterations. Default to 100.
#' @param arp order of autoregressive. Default to 0.
#' @param cov.shrink shrinkage on the multivariate normal covariance matrix.
#' Default to 0. See references.
#' @param auto.lambda elastic net shrinkage on the autoregressive coefficients.
#' Default to 0. See references.
#' @param auto.alpha The elasticnet mixing parameter, with 0<=alpha<=1. The penalty is defined 
#' as (1-alpha)/2||.||_2^2+alpha||.||_1. 
#' alpha=1 is the lasso penalty, and alpha=0 the ridge penalty. Default to 0.
#' Same as in the glmnet package.
#' @return a list containing the fitted parameters.
#' @references Rabiner, Lawrence R. "A tutorial on hidden Markov models and 
#' selected applications in speech recognition." Proceedings of the 
#' IEEE 77.2 (1989): 257-286.
#' @references Zou, Hui, and Trevor Hastie. "Regularization and variable 
#' selection via the elastic net." Journal of the Royal Statistical 
#' Society: Series B (Statistical Methodology) 67.2 (2005): 301-320.
#' @references Friedman, Jerome, Trevor Hastie, and Robert Tibshirani. 
#' "Sparse inverse covariance estimation with the graphical lasso." 
#' Biostatistics 9.3 (2008): 432-441.
#' @examples
#' set.seed(332213)
#' data(finance)
#' x <- data.matrix(finance)
#' #log return
#' y <- x[-1,]
#' for(i in 2:nrow(x)){
#'  y[i-1,] <- log(x[i,]) - log(x[i-1,])
#' }
#' #annualize the log return
#' y <- y * 252 
#' 
#' #first, fit a Gaussian HMM without autoregressive structure
#' m <- 2
#' #initialize the list of means
#' mu <- list(apply(y,2,mean), apply(y,2,mean))
#' #initialize the list of covariance matrices
#' sigma <- list(cov(y)*1.2,cov(y)*0.8)
#' #initialize the prior probability
#' delta <- c(0.5,0.5)
#' #initialize the transition probabilities
#' gamma <- matrix(c(0.9,0.1,0.2,0.8),2,2,byrow=TRUE)
#' mod1 <- list(m=m,mu=mu,sigma=sigma,delta=delta,gamma=gamma)
#' #will not run without a shrinkage on the covariance matrices because the 
#' #series is not long enough to reliably estimate the covariance structure
#' fit1 <- em.hmm(y=y,mod=mod1,cov.shrink=0.0001)
#' st1 <- viterbi.hmm(y=y,mod=fit1)
#' sp1 <- smooth.hmm(y=y,mod=fit1)
#' 
#' \dontrun{
#' #second, fit a Gaussian HMM with 1st order autoregressive structure
#' auto <- list(matrix(rep(0,2500),50,50,byrow=TRUE),
#'              matrix(rep(0,2500),50,50,byrow=TRUE))
#' mod2 <- list(m=m,mu=mu,sigma=sigma,delta=delta,gamma=gamma,auto=auto)
#' fit2 <- em.hmm(y=y,mod=mod2,ntimes=NULL,cov.shrink=0.0001,arp=1,
#'                auto.alpha=1,auto.lambda=0.1)
#' st2 <- viterbi.hmm(y=y,mod=fit2)
#' sp2 <- smooth.hmm(y=y,mod=fit2)
#' 
#' #third, fit a Gaussian HMM with 2nd order autoregressive structure
#' auto <- list(matrix(rep(0,5000),50,100,byrow=TRUE),
#'              matrix(rep(0,5000),50,100,byrow=TRUE))
#' mod3 <- list(m=m,mu=mu,sigma=sigma,delta=delta,gamma=gamma,auto=auto)
#' fit3 <- em.hmm(y=y,mod=mod3,ntimes=NULL,cov.shrink=0.0001,arp=2,
#'                auto.alpha=1,auto.lambda=0.1)
#' st3 <- viterbi.hmm(y=y,mod=fit3)
#' sp3 <- smooth.hmm(y=y,mod=fit3)
#' }
#' @export
em.hmm <- function(y, mod, ntimes=NULL, 
                   tol=1e-4, maxit=100, arp=0, 
                   cov.shrink=0, auto.lambda=0, auto.alpha=0){
  if(arp==0) 
          result <- em.mvn(y, mod, ntimes, tol, maxit, cov.shrink)
  if(arp>0) 
    result <- em.mvnarp(y, mod, ntimes, tol, maxit, arp, cov.shrink,
                        auto.lambda, auto.alpha)
  return(result)
}

############################################
#em algorithm for multivariate normal
em.mvn <- function(y, mod, ntimes, tol, maxit,cov.shrink){
  
  ns <- nrow(y)
  if(is.null(ntimes)) ntimes <- ns
  N <- length(ntimes)
  
  p <- length(mod$mu[[1]])
  Pi <- mod$delta
  P <- mod$gamma
  K <- mod$m
  mu <- mod$mu
  sigma <- mod$sigma
  
  #initialization
  k1 <- (2*3.1415926)^(-p/2)
  
  #recursion
  loglik <- 0
  for(iter in 1:maxit){
    
    oldlik <- loglik
    
    #get the cholesky factor: sigma = L%*%L^t
    L <- lapply(1:K, function(k) t(chol(sigma[[k]])))
    Lt <- lapply(L,t)
    k2 <- sapply(1:K,function(k) k1/prod(diag(L[[k]]))) #sqrt(det(sigma))
    count <- 0
    
    Gamma <- NULL
    Gammasum <- matrix(0,1,K)
    oldlik <- loglik
    colsumXi <- rep(0, K*K)
    
    #for each subject
    loglik <- 0 
    
    for(n in 1:N){
      
      alpha <- matrix(0,ntimes[n],K)
      beta <- matrix(0,ntimes[n],K)
      gamma <- matrix(0,ntimes[n],K)
      B <- matrix(0,ntimes[n],K)
      
      Scale <- rep(0, ntimes[n])
      Xi <- matrix(0, ntimes[n]-1, K*K)
    
    #state-dependent probs
    for(i in 1:ntimes[n]){
      for(j in 1:K){
        diff <- mu[[j]] - y[count+i,]
        B[i,j] <- k2[[j]]*exp(-0.5*t(diff)%*%
                           backsolve(Lt[[j]],forwardsolve(L[[j]],diff)))
      }
    }
    
    ####E-step
    #forward-backward
    
    scale <- rep(0, ntimes[n])
    alpha[1,] <- Pi*B[1,]
    scale[1] <- sum(alpha[1,])
    alpha[1,] <- alpha[1,]/scale[1]
    
    for(i in 2:ntimes[n]){
      alpha[i,] <- (alpha[i-1,]%*%P)*B[i,]
      scale[i] <- sum(alpha[i,])
      alpha[i,] <- alpha[i,] / scale[i]
    }
    
    beta[ntimes[n],] <- rep(1/K,K)/scale[ntimes[n]]
    for(i in (ntimes[n]-1):1)
      beta[i,] <- (beta[i+1,] * B[i+1,]) %*% t(P) / scale[i]
    
    gamma <- alpha*beta
    gamma <- gamma / rowSums(gamma) 
    gammasum <- colSums(gamma) #updated
    
      xi <- matrix(0, ntimes[n]-1, K*K)
      for(i in 1:(ntimes[n]-1)){
        t <- P * ( alpha[i,] %*% t(beta[i+1,] * B[i+1,])  )
        xi[i,] <- as.vector(t) / sum(t)
      }
      
      Scale <- Scale + log(scale)
      Gamma <- rbind(Gamma, gamma)
      Gammasum <- Gammasum + gammasum
      Xi <- Xi + xi
      
      colsumXi <- colsumXi + colSums(Xi)
      loglik <- loglik + sum(Scale)
      count <- count + ntimes[n]
    }
    
    ##############
    #M-step
    Mu <- matrix(0,K,p)
    Mu <- t(Gamma)%*%y
    Mu <- Mu / as.vector(Gammasum)
    mu <- vector(mode="list", length=K)
    for(i in 1:K) mu[[i]] <- Mu[i,]
      
    sxi <- matrix(colsumXi, K,K)
    P <- sxi / rowSums(sxi)
    Pi <- rep(0, K)
    picount <- 0
    for(n in 1:N){
      Pi <- Pi + Gamma[picount+1,]
      picount <- picount + ntimes[n]
    }
    Pi <- Pi / N
    
    sigma <- vector(mode="list", length=K)
    
    for(l in 1:K){
      dd <- y - rep(1,ns) %*% t(Mu[l,])
      
      expandrow <- Gamma[,l] %*% t(rep(1,p))
      rprod <- dd * expandrow
      sigma[[l]] <- t(rprod) %*% dd 
      
      sigma[[l]] <- sigma[[l]] / Gammasum[[l]] #updated
      
      #shrinkage as weighted average
      #sigma = (1-W)sigma_ml + W*alpha*I
      alpha <- sum(diag(sigma[[l]]))/p #ensure tr(sigma_ml) = tr(I)
      W <- cov.shrink / (1+cov.shrink)
      sigma[[l]] <- sigma[[l]]*(1-W) + diag(alpha,p)*W
    }
    
    
    
    if(iter <= 1) likbase <- loglik
    
    cat("iteration ",iter,"; loglik = ", loglik, "\n")
    
    if(iter > maxit | 
       (iter>1 & (loglik - likbase)/(oldlik - likbase) < 1 + tol) ) break;
  }
  
  return(list(m=K, mu=mu, sigma=sigma, delta=Pi, gamma=P,
              loglik=loglik, iteration=iter))
}


#####################################################
#####################################################
em.mvnarp <- function(y, mod, ntimes, tol, maxit, arp,
                      cov.shrink, auto.lambda, auto.alpha){
  
  ns <- nrow(y)
  if(is.null(ntimes)) ntimes <- ns
  N <- length(ntimes)
  
  p <- length(mod$mu[[1]])
  Pi <- mod$delta
  P <- mod$gamma
  K <- mod$m
  mu <- mod$mu
  sigma <- mod$sigma
  auto <- mod$auto
  
  #initialization
  k1 <- (2*3.1415926)^(-p/2)
  
  #recursion
  loglik <- 0
  for(iter in 1:maxit){
    
    oldlik <- loglik
    
    #get the cholesky factor: sigma = L%*%L^t
    L <- lapply(1:K, function(k) t(chol(sigma[[k]])))
    Lt <- lapply(L,t)
    k2 <- sapply(1:K,function(k) k1/prod(diag(L[[k]]))) #sqrt(det(sigma))
    count <- 0
    
    Gamma <- NULL
    Gammasum <- matrix(0,1,K)
    oldlik <- loglik
    colsumXi <- rep(0, K*K)
    
    #for each subject
    loglik <- 0 
    
    for(n in 1:N){
      
      alpha <- matrix(0,ntimes[n],K)
      beta <- matrix(0,ntimes[n],K)
      gamma <- matrix(0,ntimes[n],K)
      B <- matrix(0,ntimes[n],K)
      
      Scale <- rep(0, ntimes[n])
      Xi <- matrix(0, ntimes[n]-1, K*K)
      
      #state-dependent probs
      for(i in 1:ntimes[n]){
        for(j in 1:K){
          if(i==1)  diff <- mu[[j]] - y[count+i,]
          else if(i<=arp)diff <- mu[[j]] + auto[[j]][,(p*arp-(i-1)*p+1):(p*arp),drop=FALSE]%*%
              as.vector(t(y[(count+1):(count+i-1),])) - 
               y[count+i,]
          else diff <- mu[[j]] + auto[[j]]%*%
              as.vector(t(y[(count+i-arp):(count+i-1),])) - 
              y[count+i,]
          B[i,j] <- k2[[j]]*exp(-0.5*t(diff)%*%
                                  backsolve(Lt[[j]],forwardsolve(L[[j]],diff)))
        }
      }
      
      ####E-step
      #forward-backward
      alpha[1:arp,] <- rep(0,K) ########
      beta[1:arp,] <- rep(0,K)
      
      scale <- rep(0, ntimes[n])
      scale[1:arp] <- 0 ########
      
      alpha[arp+1,] <- Pi*B[arp+1,]
      scale[arp+1] <- sum(alpha[arp+1,])
      alpha[arp+1,] <- alpha[arp+1,]/scale[arp+1]
      
      for(i in (arp+2):ntimes[n]){
        alpha[i,] <- (alpha[i-1,]%*%P)*B[i,]
        scale[i] <- sum(alpha[i,])
        alpha[i,] <- alpha[i,] / scale[i]
      }
      
      beta[ntimes[n],] <- rep(1/K,K)/scale[ntimes[n]]
      for(i in (ntimes[n]-1):(arp+1)) ##
        beta[i,] <- (beta[i+1,] * B[i+1,]) %*% t(P) / scale[i]
      
      gamma <- alpha*beta
      gamma <- gamma / pmax(rowSums(gamma),0.01) 
      gammasum <- colSums(gamma) #updated
      
      xi <- matrix(0, ntimes[n]-1, K*K)
      for(i in (1+arp):(ntimes[n]-1)){ ###############################
        t <- P * ( alpha[i,] %*% t(beta[i+1,] * B[i+1,])  )
        xi[i,] <- as.vector(t) / sum(t)
      }
      
      Scale <- Scale + log(scale)
      Gamma <- rbind(Gamma, gamma)
      Gammasum <- Gammasum + gammasum
      Xi <- Xi + xi
      
      colsumXi <- colsumXi + colSums(Xi)
      loglik <- loglik + sum(Scale[-(1:arp)])#######
      count <- count + ntimes[n]
      
    }
    
    ##############
    #M-step
    
    sxi <- matrix(colsumXi, K,K)
    P <- sxi / rowSums(sxi)
    Pi <- rep(0, K)
    picount <- 0
    for(n in 1:N){
      Pi <- Pi + Gamma[picount+arp+1,] ######
      picount <- picount + ntimes[n]
    }
    Pi <- Pi / N
    
    sigma <- vector(mode="list", length=K)
    auto <- vector(mode="list", length=K)
    
    for(l in 1:K){
      ddfor <- y[(1+arp):ns,]
      
      for(lag in 1:arp){
         if(lag==1) ddlag <- y[1:(ns-arp),]
         else ddlag <- cbind(ddlag, y[lag:(ns-arp-1+lag),])   
      }
      #model <- lm(ddfor~ddlag,weights=Gamma[-1,l])
      model <- glmnet(ddlag,ddfor,family="mgaussian",
                      lambda=auto.lambda,alpha=auto.alpha,
                      intercept = T,standardize = F,weights=Gamma[(1+arp):ns,l])
      
      tempcoef <- glmnet::coef.glmnet(model)
      coefmat <- tempcoef[[1]]
      for(i in 2:p) coefmat <- cbind(coefmat, tempcoef[[i]])
      coefmat <- as.matrix(coefmat)
      
      mu[[l]] <- coefmat[1,]
      auto[[l]] <- t(coefmat[-1,])
      colnames(auto[[l]])=rep("",arp*p)
      rownames(auto[[l]])=rep("",p)
      
      expandrow <- Gamma[(1+arp):ns,l] %*% t(rep(1,p))
      yhat <- glmnet::predict.mrelnet(model,ddlag,type="response")[,,1]
      resid <- ddfor - yhat
      rprodvar <- resid * expandrow
      sigma[[l]] <- t(rprodvar) %*% resid
      sigma[[l]] <- sigma[[l]] / Gammasum[[l]] #updated
      
      #shrinkage as weighted average
      #sigma = (1-W)sigma_ml + W*alpha*I
      alpha <- sum(diag(sigma[[l]]))/p #ensure tr(sigma_ml) = tr(I)
      W <- cov.shrink / (1+cov.shrink)
      sigma[[l]] <- sigma[[l]]*(1-W) + diag(alpha,p)*W
    }
    
    if(iter <= 1) likbase <- loglik
    
    cat("iteration ",iter,"; loglik = ", loglik, "\n")
    
    if(iter > maxit | 
       (iter>1 & (loglik - likbase)/(oldlik - likbase) < 1 + tol) ) break;
  }
  
  return(list(m=K, mu=mu, sigma=sigma, delta=Pi, gamma=P, auto=auto, arp=arp,
              loglik=loglik, iteration=iter))
}


