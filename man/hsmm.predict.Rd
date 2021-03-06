% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hsmm.predict.R
\name{hsmm.predict}
\alias{hsmm.predict}
\title{1-step forward prediction for (autoregressive) Gaussian hidden semi-Markov model}
\usage{
hsmm.predict(y, mod)
}
\arguments{
\item{y}{observed series}

\item{mod}{list consisting the at least the following items: 
mod$m = scalar number of states, 
mod$delta = vector of initial values for prior probabilities, 
mod$gamma = matrix of initial values for state transition probabilies.
mod$mu = list of initial values for means, 
mod$sigma = list of initial values for covariance matrices.
mod$d = list of state duration probabilities.
For autoregressive hidden markov models, we also need the additional items:
mod$arp = scalar order of autoregressive structure
mod$auto = list of initial values for autoregressive coefficient matrices}
}
\value{
1-step forward state probabilities and forecasts
}
\description{
1-step forward prediction for (autoregressive) Gaussian hidden semi-Markov model
}
\examples{
set.seed(15562)
m <- 2
mu <- list(c(3,4,5),c(-2,-3,-4))
sigma <- list(diag(1.3,3), 
            matrix(c(1,-0.3,0.2,-0.3,1.5,0.3,0.2,0.3,2),3,3,byrow=TRUE))
delta <- c(0.5,0.5)
gamma <- matrix(c(0,1,1,0),2,2,byrow=TRUE)
d <- list(c(0.4,0.2,0.1,0.1,0.1,0.1),c(0.5,0.3,0.2))
auto <- list(matrix(c(0.3,0.2,0.1,0.4,0.3,0.2,
                     -0.3,-0.2,-0.1,0.3,0.2,0.1,
                      0,0,0,0,0,0),3,6,byrow=TRUE),
            matrix(c(0.2,0,0,0.4,0,0,
                      0,0.2,0,0,0.4,0,
                     0,0,0.2,0,0,0.4),3,6,byrow=TRUE))
mod <- list(m=m,mu=mu,sigma=sigma,delta=delta,gamma=gamma,
            auto=auto,arp=2,d=d)
sim <- hsmm.sim(2000,mod)
y <- sim$series
state <- sim$state
fit <- em.semi(y=y, mod=mod, arp=2,tol=1e-5)
forecast <- hsmm.predict(y=y,mod=fit)

}
\references{
Rabiner, Lawrence R. "A tutorial on hidden Markov models and 
selected applications in speech recognition." Proceedings of the 
IEEE 77.2 (1989): 257-286.
}
