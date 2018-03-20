######################################################################
#' EM algorithm to compute maximum likelihood estimate of Gaussian 
#' hidden semi-Markov models with / without autoregressive structures and 
#' with / without regularization on the covariance matrices and/or
#' autoregressive structures.
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
#' @param print Default to TRUE.
#' @return a list containing the fitted parameters.
#' @references Rabiner, Lawrence R. "A tutorial on hidden Markov models and 
#' selected applications in speech recognition." Proceedings of the 
#' IEEE 77.2 (1989): 257-286.
#' @references Zou, Hui, and Trevor Hastie. "Regularization and variable 
#' selection via the elastic net." Journal of the Royal Statistical 
#' Society: Series B (Statistical Methodology) 67.2 (2005): 301-320.
#' @references Ledoit, Olivier, and Michael Wolf. "A well-conditioned estimator for
#' large-dimensional covariance matrices." Journal of multivariate analysis 88.2
#' (2004): 365-411.
#' @examples 
#' \dontrun{
#' set.seed(332213)
#' data(finance)
#' x <- data.matrix(finance)
#' #log return
#' y <- x[-1,-51]
#' for(i in 2:nrow(x)){
#'  y[i-1,] <- log(x[i,-51]) - log(x[i-1,-51])
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
#' gamma <- matrix(c(0,1,1,0),2,2,byrow=TRUE)
#' #initialize the state duration probabilities
#' d <- list(rep(0.1,10),rep(0.1,10))
#' mod1 <- list(m=m,mu=mu,sigma=sigma,delta=delta,gamma=gamma,d=d)
#' #will not run without a shrinkage on the covariance matrices because the 
#' #series is not long enough to reliably estimate the covariance structure
#' fit1 <- em.semi(y=y,mod=mod1,cov.shrink=0.0001)
#' st1 <- viterbi.semi(y=y,mod=fit1)
#' sp1 <- smooth.semi(y=y,mod=fit1)
#' 
#' #second, fit a Gaussian HSMM with 1st order autoregressive structure
#' auto <- list(matrix(rep(0,2500),50,50,byrow=TRUE),
#'              matrix(rep(0,2500),50,50,byrow=TRUE))
#' mod2 <- list(m=m,mu=mu,sigma=sigma,delta=delta,gamma=gamma,auto=auto,
#'              d=d,arp=1)
#' #increase auto.lambda to enforce stronger regularization for model to run
#' fit2 <- em.semi(y=y,mod=mod2,cov.shrink=0.001,arp=1,
#'                auto.alpha=0.8,auto.lambda=10)
#' sum(fit2$auto[[1]]==0)
#' sum(fit2$auto[[2]]==0)
#' st2 <- viterbi.semi(y=y,mod=fit2)
#' sp2 <- smooth.semi(y=y,mod=fit2)
#' 
#' #third, fit a Gaussian HSMM with 2nd order autoregressive structure
#' auto <- list(matrix(rep(0,5000),50,100,byrow=TRUE),
#'              matrix(rep(0,5000),50,100,byrow=TRUE))
#' mod3 <- list(m=m,mu=mu,sigma=sigma,delta=delta,gamma=gamma,auto=auto,
#'              d=d,arp=2)
#' #increase auto.lambda to enforce stronger regularization for model to run
#' fit3 <- em.semi(y=y,mod=mod3,ntimes=NULL,cov.shrink=0.0001,arp=2,
#'                auto.alpha=0.8,auto.lambda=30)
#' sum(fit3$auto[[1]]==0)
#' sum(fit3$auto[[2]]==0)
#' st3 <- viterbi.semi(y=y,mod=fit3)
#' sp3 <- smooth.semi(y=y,mod=fit3)
#' }
#' @useDynLib rarhsmm, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom graphics points
#' @importFrom stats rnorm
#' @importFrom glmnet glmnet
#' @export
em.semi <- function(y,mod, ntimes=NULL, 
                   tol=1e-4, maxit=100, arp=0, 
                   cov.shrink=0, auto.lambda=0, auto.alpha=0,print=TRUE){
  if(arp==0) 
    result <- em.semi.mvn(y, mod, ntimes, tol, maxit, cov.shrink,print)
  if(arp>0) 
    result <- em.semi.mvnarp(y, mod, ntimes, tol, maxit, arp, cov.shrink,
                        auto.lambda, auto.alpha,print)
  return(result)
}

#############
#convert hsmm to hmm
hsmm2hmm <- function(omega, dm, eps=1e-6){
  mv <- sapply(dm, length)
  m <- length(mv)
  G <- matrix(0,0,sum(mv))
  for(i in 1:m){
    mi <- mv[[i]]
    f <- cumsum(c(0,dm[[i]][-mi]))
    ci <- ifelse(abs(1-f)>eps, dm[[i]]/(1-f),1)
    cim <- ifelse(1-ci>0, 1-ci, 0)
    Gi <- matrix(0,mi,0)
    for(j in 1:m){
      if(i==j){ if(mi==1){
        Gi <- cbind(Gi,c(rep(0,mv[[j]]-1),cim))
      } else{
        Gi <- cbind(Gi,rbind(cbind(rep(0,mi-1),diag(cim[-mi],mi-1,mi-1)),
                             c(rep(0,mi-1),cim[[mi]])))
      }}else{ if(mi==1) {
        Gi <- cbind(Gi,matrix(c(omega[[i,j]]*ci,rep(0,mv[[j]]-1)),1))
      } else{
        Gi <- cbind(Gi,cbind(omega[[i,j]]*ci,matrix(0,mv[[i]],mv[[j]]-1)))
      }}
    }
    G <- rbind(G,Gi)
  }
  G
}


##################
em.semi.mvn <- function(y, mod, ntimes, tol, maxit,cov.shrink, print){
  
  d <- mod$d
  ld <- sapply(d,length)
  D <- sum(ld)
  
  ns <- nrow(y)
  if(is.null(ntimes)) ntimes <- ns
  N <- length(ntimes)
  
  p <- length(mod$mu[[1]])
  Pi <- mod$delta
  P <- mod$gamma
  K <- mod$m
  mu <- mod$mu
  sigma <- mod$sigma
    
  #recursion
  loglik <- 0
  for(iter in 1:maxit){
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
    
    #forward-backward
    oldlik <- loglik
    nodeprob <- getnodeprob_nocov_mvn(y, newmu, newsigma, D, p, 0, 0)
    fb <- forwardbackward(newPi,newP,nodeprob,ns,ntimes)
    colsumXi <- fb$colsumxi
    Gamma <- fb$Gamma
    Gammasum <- fb$colsumgamma
    loglik <- fb$loglik
    
    if(iter <= 1) likbase <- loglik
    
    if(iter > maxit | 
       (iter>1 & (loglik - likbase)/(oldlik - likbase) < 1 + tol) ) {
      loglik <- oldlik; break;}
    
    ##############
    #M-step
    oldGamma <- t(sapply(1:ns, function(k){ 
      sapply(split(Gamma[k,],rep(1:K,ld)),sum)}))
    oldGammasum <- colSums(oldGamma)
    
    Mu <- matrix(0,K,p)
    Mu <- t(oldGamma)%*%y
    Mu <- Mu / as.vector(oldGammasum)
    mu <- vector(mode="list", length=K)
    for(i in 1:K) mu[[i]] <- Mu[i,]
    
    cumld <- c(0,cumsum(ld))
    sxi <- matrix(colsumXi, D,D)
    nondiag <- vector(mode="list",length=K)
    for(i in 1:K){
      sublist <- matrix(0,ld[i],K-1)
      transto <- c(1:K)[-i]
      
      for(row in 1:ld[i]){
        for(col in 1:(K-1)){
          sublist[row,col] <- sxi[cumld[i]+row,cumld[transto[col]]+1]
          #state i to next smallest state in row duration
        }
      }
      nondiag[[i]] <- sublist
    }
    
    #P <- sxi / rowSums(sxi)
    if(K==2) {
      P <- matrix(c(0,1,1,0),2,2,byrow=T)}else{
      temp <- t(sapply(nondiag,colSums))
      for(i in 1:K) P[i,-i] <- temp[i,]
      P <- P/as.vector(rowSums(P))
    }
    
    d <- vector(mode="list",length=K)
    temp <- lapply(nondiag,rowSums)
    for(i in 1:K) d[[i]] <- temp[[i]]/sum(temp[[i]])
      
    Pi <- oldGamma[1,]
    
    sigma <- vector(mode="list", length=K)
    
    for(l in 1:K){
      dd <- y - rep(1,ns) %*% t(Mu[l,])
      
      expandrow <- oldGamma[,l] %*% t(rep(1,p))
      rprod <- dd * expandrow
      sigma[[l]] <- t(rprod) %*% dd 
      
      sigma[[l]] <- sigma[[l]] / oldGammasum[[l]] #updated
      
      #shrinkage as weighted average
      #sigma = (1-W)sigma_ml + W*alpha*I
      alpha <- sum(diag(sigma[[l]]))/p #ensure tr(sigma_ml) = tr(I)
      W <- cov.shrink / (1+cov.shrink)
      sigma[[l]] <- sigma[[l]]*(1-W) + diag(alpha,p)*W
    }
  
    if(print==TRUE) cat("iteration ",iter,"; loglik = ", loglik, "\n")

  }
  
  return(list(m=K, mu=mu, sigma=sigma, delta=Pi, gamma=P,d=d,
              loglik=loglik, iteration=iter))
}


#############
em.semi.mvnarp <- function(y, mod, ntimes, tol, maxit, arp,
                           cov.shrink, auto.lambda, auto.alpha,print){
  d <- mod$d
  ld <- sapply(d,length)
  D <- sum(ld)
  
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
  
  #recursion
  loglik <- 0
  
  for(iter in 1:maxit){
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
    
    oldlik <- loglik
    ################################
    #nodeprob
    if(p>1){
      ycov <- rbind(rep(0,p),y[-ntimes,])
    }else{
      ycov <- matrix(c(0,y[-ntimes,]),ncol=1)
    }
    autoarray <- array(as.numeric(unlist(newauto)), dim=c(p,p,D))
    muarray <- array(as.numeric(unlist(newmu)), dim=c(1,p,D))
    nodeprob <- getnodeprob_part2(y, ycov,autoarray,muarray,newsigma,D,p)
    
    fb <- forwardbackward(newPi,newP,nodeprob,ns,ntimes)
    colsumXi <- fb$colsumxi
    Gamma <- fb$Gamma
    Gammasum <- fb$colsumgamma
    loglik <- fb$loglik
    if(iter <= 1) likbase <- loglik
    
    if(iter > maxit | 
       (iter>1 & (loglik - likbase)/(oldlik - likbase) < 1 + tol) ) {
      loglik <- oldlik; break;}
 
    ##############
    #M-step
    oldGamma <- t(sapply(1:ns, function(k){ 
      sapply(split(Gamma[k,],rep(1:K,ld)),sum)}))
    oldGammasum <- colSums(oldGamma)
    
    cumld <- c(0,cumsum(ld))
    sxi <- matrix(colsumXi, D,D)
    nondiag <- vector(mode="list",length=K)
    for(i in 1:K){
      sublist <- matrix(0,ld[i],K-1)
      transto <- c(1:K)[-i]
      
      for(row in 1:ld[i]){
        for(col in 1:(K-1)){
          sublist[row,col] <- sxi[cumld[i]+row,cumld[transto[col]]+1]
          #state i to next smallest state in row duration
        }
      }
      nondiag[[i]] <- sublist
    }
    
    #P <- sxi / rowSums(sxi)
    if(K==2) {
      P <- matrix(c(0,1,1,0),2,2,byrow=T)
    }else{
      temp <- t(sapply(nondiag,colSums))
      for(i in 1:K) P[i,-i] <- temp[i,]
      P <- P/as.vector(rowSums(P))
    }
    
    d <- vector(mode="list",length=K)
    temp <- lapply(nondiag,rowSums)
    for(i in 1:K) d[[i]] <- temp[[i]]/sum(temp[[i]])
    
    Pi <- oldGamma[1,]
    
    sigma <- vector(mode="list", length=K)
    auto <- vector(mode="list", length=K)
    
    for(l in 1:K){
      ddfor <- y[(1+arp):ns,]
      
      for(lag in 1:arp){
        if(lag==1) ddlag <- y[1:(ns-arp),]
        else ddlag <- cbind(ddlag, y[lag:(ns-arp-1+lag),])   
      }
      
      if(p>1){
      #model <- lm(ddfor~ddlag,weights=Gamma[-1,l])
      model <- glmnet(ddlag,ddfor,family="mgaussian",
                      lambda=auto.lambda,alpha=auto.alpha,
                      intercept = T,standardize = F,weights=oldGamma[(1+arp):ns,l])
      
      tempcoef <- glmnet::coef.glmnet(model)
      coefmat <- tempcoef[[1]]
      for(i in 2:p) coefmat <- cbind(coefmat, tempcoef[[i]])
      coefmat <- as.matrix(coefmat)
      
      mu[[l]] <- coefmat[1,]
      auto[[l]] <- t(coefmat[-1,])
      colnames(auto[[l]])=rep("",arp*p)
      rownames(auto[[l]])=rep("",p)
      
      expandrow <- oldGamma[(1+arp):ns,l] %*% t(rep(1,p))
      yhat <- glmnet::predict.mrelnet(model,ddlag,type="response")[,,1]
      resid <- ddfor - yhat
      rprodvar <- resid * expandrow
      sigma[[l]] <- t(rprodvar) %*% resid
      sigma[[l]] <- sigma[[l]] / oldGammasum[[l]] #updated
      
      #shrinkage as weighted average
      #sigma = (1-W)sigma_ml + W*alpha*I
      alpha <- sum(diag(sigma[[l]]))/p #ensure tr(sigma_ml) = tr(I)
      W <- cov.shrink / (1+cov.shrink)
      sigma[[l]] <- sigma[[l]]*(1-W) + diag(alpha,p)*W
      }else{
        model <- lm(ddfor~ddlag,weights=Gamma[(1+arp):ns,l])
        tempcoef <- coef(model)
        mu[[l]] <- tempcoef[1]
        auto[[l]] <- matrix(tempcoef[2])
        sigma[[l]] <- matrix(summary(model)$sigma^2)
      }
    }
    
    if(print==TRUE) cat("iteration ",iter,"; loglik = ", loglik, "\n")
  
  }
  
  return(list(m=K, mu=mu, sigma=sigma, delta=Pi, gamma=P,d=d,auto=auto,
              arp=arp,loglik=loglik, iteration=iter))
}
