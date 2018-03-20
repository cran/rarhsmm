// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//
////////////////////

/////////////////////////
// [[Rcpp::export]]
arma::mat vectorize(arma::mat mat, int axis){
    int row = mat.n_rows;
    int col = mat.n_cols;
    int i,j;
    arma::mat result(row*col,1);
    for(i=0;i<row;i++){
        for(j=0;j<col;j++){
            if(axis==0)  result(i+j*row,0) = mat(i,j);  // by column
            else result(j+i*col,0) = mat(i,j); //by row
        }
    }
    return result;
}
///////////////////////////////////////
// [[Rcpp::export]]
arma::vec colsum(arma::mat matrix){
    long row = matrix.n_rows;
    int col = matrix.n_cols;
    long i;
    int j;
    arma::vec columnsum(col);
    
    for(j=0; j<col; j++) columnsum(j) = 0;
    for(i=0; i<row; i++){
        for(j=0; j<col; j++){
            columnsum(j) += matrix(i,j);
        }
    }
    return columnsum;
}

///////////////////////////////////////
// [[Rcpp::export]]
arma::vec rowsum(arma::mat matrix){
    long row = matrix.n_rows;
    int col = matrix.n_cols;
    long i;
    int j;
    arma::vec rowsum(row);
    
    for(i=0; i<row; i++) rowsum(i) = 0;
    for(j=0; j<col; j++){
        for(i=0; i<row; i++){
            rowsum(i) += matrix(i,j);
        }
    }
    return rowsum;
}

//////////////////////////////////////////////////////////////////////////////////
//subset a big matrix by the beginning and ending of row/col index
// [[Rcpp::export]]
arma::mat subsetmatrix(arma::mat rawmat, arma::vec rowindex, arma::vec colindex){
    int row = rowindex(1) - rowindex(0) + 1;
    int col = colindex(1) - colindex(0) + 1;
    arma::mat result(row,col);
    int i,j;
    for(i=0;i<row;i++){
        for(j=0;j<col;j++){
            result(i,j) = rawmat(rowindex(0)+i,colindex(0)+j);
        }
    }
    return result;
}

/////////////////////////////
//compute matrix power
// [[Rcpp::export]]
arma::mat matrixpower(arma::mat oldmat, int power){
    arma::mat newmat = oldmat;
    int i;
    
    for(i=1;i<power;i++) newmat *= oldmat;
    
    return newmat;
}

///////////////////////////
//compute matrix exponential
// [[Rcpp::export]]
arma::mat matrixexp(arma::mat oldmat, double t){
    arma::mat temp = t * oldmat;
    arma::mat newmat = arma::expmat(temp);
    return newmat;
}

/////////////
// [[Rcpp::export]]
double matrixsum(arma::mat mat1, arma::mat mat2){
    double tempsum=0;
    arma::mat temp = mat1 % mat2;
    int row = mat1.n_rows;
    int col = mat1.n_cols;
    int i,j;
    for(i=0;i<row;i++){
        for(j=0;j<col;j++){
            tempsum += temp(i,j);
        }
    }
    return(tempsum);
}

////////////
//integral of matrix exponential
// x and y starts from 0
// [[Rcpp::export]]
arma::mat matrixintegral(arma::mat Q, double interval, int x, int y){
    int M = Q.n_cols;
    arma::mat temp(2*M, 2*M);
    temp.zeros();
    int i,j;
    
    for(i=0;i<M;i++)
    for(j=0;j<M;j++)
    temp(i,j) = Q(i,j);
    
    for(i=M;i<2*M;i++)
    for(j=M;j<2*M;j++)
    temp(i,j) = Q(i-M,j-M);
    
    temp(x, M+y) = 1;
    
    
    arma::mat temp1 = matrixexp(temp, interval);
    //Rcpp::Rcout<<temp1<<std::endl;
    arma::vec rowindex(2);
    arma::vec colindex(2);
    rowindex(0)=0; rowindex(1) = M-1; colindex(0) = M; colindex(1) = 2*M-1;
    arma::mat result = subsetmatrix(temp1,rowindex,colindex);
    return result;
}



/////////////
//convert discrete time HSMM tpm to discrete time HMM tpm with expanded states
//[[Rcpp::export]]
arma::mat hsmm_hmm (arma::mat omega, arma::mat dm, arma::vec mv){
    //each row in dm is a dwell time pmf
    //mv is vector of the length until the first zero in each dwell time distribution
    int m = omega.n_rows;
    int dmrow = dm.n_rows;
    int dmcol = dm.n_cols;
    int dim = arma::sum(mv); // dimension of the final result
    
    int i,j,p,q,mi,rowsum,colsum;
    //double tempsum;
    arma::mat temp(dmrow,dmcol);
    arma::mat ci(dmrow,dmcol);
    arma::mat cim(dmrow,dmcol);
    arma::mat gamma(dim,dim);
    
    //Rcpp::Rcout << dim << std::endl;
    
    for(i=0;i<m;i++){
        mi = mv[i];
        
        for(j=0;j<mi;j++){
            if(j==0) temp(i,j) = 0;
            else temp(i,j) = temp(i,j-1) + dm(i,j-1);
        }
        
        for(j=0;j<mi;j++){
            if(std::abs(1-temp(i,j))>0.000000001) ci(i,j) = dm(i,j)/(1-temp(i,j));
            else ci(i,j) = 1;
            if(1-ci(i,j)>0) cim(i,j)=1-ci(i,j);
            else cim(i,j) = 0;
        }
    }
    
    rowsum = 0;
    
    
    for(i=0; i<m; i++){
        colsum = 0;
        for(j=0; j<m; j++){
            if(i==j){
                if(mv[i]==1) gamma(rowsum,colsum) = cim(i,0);
                else{
                    for(p=0; p<mv[i]; p++){
                        for(q=0; q<mv[j]; q++){
                            if((q-p)==1) gamma(rowsum+p,colsum+q)=cim(i,p);
                            else if((p==mv[i]-1) & (q==mv[j]-1)) gamma(rowsum+p,colsum+q)=cim(i,p);
                            else gamma(rowsum+p,colsum+q)=0;
                        }
                    }
                }
            }
            else{
                for(p=0; p<mv[i]; p++){
                    for(q=0; q<mv[j]; q++){
                        if(q==0) gamma(rowsum+p, colsum+q)=omega(i,j)*ci(i,p);
                        else gamma(rowsum+p, colsum+q)=0;
                    }
                }
                
            }
            colsum += mv[j];
        }
        rowsum += mv[i];
    }
    
    
    return(gamma);
    
}



/////////////////////////////////
//multinomial RNG
//' multinomial random variable generator
//' @param n number of random variables to generate
//' @param k number of categories
//' @param prob vector of probabilities summing up to 1
//' @param label vector of labels for each category
//' @return multinomial random variables
//' @export
//[[Rcpp::export]]
arma::vec rmultinomial(int n, int k, arma::vec prob, arma::vec label){
    arma::vec result(n);
    arma::vec cumprob(k);
    int i,j;
    double u;
    
    cumprob(0) = prob(0);
    for(i=1; i<k;i++){
        cumprob(i) = cumprob(i-1) + prob(i);
    }
    
    
    for(j=0;j<n;j++){
        u = Rcpp::runif(1,0,1)(0);   //to make type match
        for(i=0; i<k;i++){
            if(u<cumprob(i)) {
                result(j) = label(i);
                break;
            }
        }
    }
    return(result);
}


///////////////////////////////////////////////////////
//generate multivariate normal random variables
//' multivariate normal random number generator
//' @param n number of random vectors to generate
//' @param mu vector of means
//' @param sigma variance-covariance matrix
//' @return a matrix with each row as a realization of multivariate normal random vector
//' @export
// [[Rcpp::export]]
arma::mat mvrnorm(int n, arma::vec mu, arma::mat sigma) {
    int ncols = sigma.n_cols;
    arma::mat Y = arma::randn(n, ncols); //row vector of standard normal RVs
    return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma); // chol() returns upper by default
}

////////////////////////////////
//multivariate normal density
//' multivariate normal density
//' @param x matrix with each row as an observed vector of multivariate normal RVs
//' @param mean vector of means
//' @param sigma variance-covariance matrix
//' @param logd whether log transformation should be used
//' @return a vector of density values for each observed vector
//' @export
// [[Rcpp::export]]
arma::vec mvdnorm(arma::mat x,
                  arma::rowvec mean,
                  arma::mat sigma,
                  bool logd) {
    const double log2pi = std::log(2.0 * M_PI);
    int n = x.n_rows;
    int xdim = x.n_cols;
    arma::vec out(n);
    arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma)))); //generate upper triangular
    double rootisum = arma::sum(log(rooti.diag()));
    double constants = -(static_cast<double>(xdim)/2.0) * log2pi;
    
    for (int i=0; i < n; i++) {
        arma::vec z = rooti * arma::trans( x.row(i) - mean) ;
        out(i)      = constants - 0.5 * arma::sum(z%z) + rootisum;
    }
    
    if (logd == FALSE) {
        out = exp(out);
    }
    return(out);
}

//number of trials before failure (p=prob of failure)
//[[Rcpp::export]]
arma::vec rgeometric(int n, double pp){
    arma::vec result;
    result = Rcpp::rgeom(n,pp)+1;
    return result;
}

//[[Rcpp::export]]
double dgeometric(int x, double pp, bool loga){
    double result;
    result = R::dgeom(x-1, pp, loga);
    return result;
}




///////////////////////////////////////////////////////////////////////////
/*part 2: simulate a hidden Markov series*/


///////////////////////////
//generate multivariate normal hidden markov time series
//[[Rcpp::export]]
Rcpp::List mvnhmm_gen (long n, int M, int d, arma::vec prior, arma::mat tpm,
                       Rcpp::List meanlist, Rcpp::List sigmalist){
    
    int prev, curr, i, j;
    arma::vec label(M);
    
    arma::mat series(n,d);
    arma::vec state(n);
    
    for(i=0; i<M; i++) label(i) = i + 1;
    
    //initial state
    state(0) = rmultinomial(1, M, prior, label)(0);
    curr = state(0) - 1;
    series.row(0) = mvrnorm(1, meanlist(curr), sigmalist(curr));
    
    //recursion
    for(j=1; j<n; j++){
        prev = state(j-1) - 1;
        state(j) = rmultinomial(1, M, tpm.row(prev).t(), label)(0);
        
        curr = state(j) - 1;
        series.row(j) = mvrnorm(1, meanlist(curr), sigmalist(curr));
    }
    
    return Rcpp::List::create(Rcpp::Named("series")=series,
                              Rcpp::Named("state")=state);
    
}



/////////////////////////////////////
/*part 3: log likelihoods*/
///////////////////////////////////////////
////////////////////////////////
// [[Rcpp::export]]
arma::mat getnodeprob_nocov_mvn(arma::mat ystar, Rcpp::List mean, Rcpp::List sigma,
                                int m, int p, int arp, Rcpp::List automat){
    long ns = ystar.n_rows; //original numeric series starting from 1
    
    arma::mat nodeprob(ns, m);
    arma::mat tempmat(1,p);
    long i;
    int j;
    
    if(arp==0){
        for(i=0;i<ns;i++){
            for(j=0;j<m;j++){
                tempmat = ystar.row(i);
                nodeprob(i,j) = mvdnorm(tempmat,mean(j),sigma(j),FALSE)(0);
            }
        }
    }
    
    return nodeprob;
}

////////////////////////////////
// [[Rcpp::export]]
arma::mat getnodeprob_part2(arma::mat y, arma::mat x, arma::cube beta, arma::cube mu, Rcpp::List sigma,
                                int m, int p){
    long ns = y.n_rows; //original numeric series starting from 1
    
    arma::mat nodeprob(ns, m);
    arma::mat tempmat(1,p);
    arma::mat mean(1,p);
    long i;
    int j;
    
  
        for(i=0;i<ns;i++){
            for(j=0;j<m;j++){
                tempmat = y.row(i);
                mean = mu.slice(j) + x.row(i) * beta.slice(j).t();
                //Rcpp::Rcout<<mean<<std::endl;
                nodeprob(i,j) = mvdnorm(tempmat,mean,sigma(j),FALSE)(0);
            }
        }
    
    
    return nodeprob;
}




////////////////////////////////////////////////////
// this forward-backward can be used to build the smoothing function for prediction and interpolation
// [[Rcpp::export]]
Rcpp::List forwardbackward(arma::vec Pi, arma::mat P, arma::mat nodeprob,
                           long dim, arma::vec ntimes){
    
    
    int M = nodeprob.n_cols;
    int n = ntimes.n_rows;
    int j,m,t,i,k;
    double tempsum;
    arma::vec tempval;
    arma::mat tempmat(1,M);
    arma::mat tempmat2(M,M);
    
    arma::mat alpha(dim, M);
    arma::vec scale(dim);
    arma::mat beta(dim, M);
    arma::mat Gamma(dim, M);
    arma::mat xi(dim-1, M*M);
    
    arma::vec colsumgamma(M);
    arma::vec tempsumxi(M*M);
    arma::mat colsumxi(M,M);
    double loglik = 0;
    
    int count = 0;
    for(j=0; j<n; j++){
        
        //forward probs
        alpha.row(count) = Pi.t() % nodeprob.row(count);
        tempval = nodeprob.row(count) * Pi;
        scale(count) = tempval(0);//double
        for(m=0; m<M; m++) alpha(count, m) = alpha(count,m) / scale(count);
        
        for(t=1; t<ntimes(j); t++){
            tempmat = alpha.row(count+t-1) * P; //row matrix
            alpha.row(count+t) = tempmat % nodeprob.row(count+t);
            tempval = tempmat * nodeprob.row(count+t).t();
            scale(count+t) = tempval(0);//double
            for(m=0; m<M; m++) alpha(count+t, m) = alpha(count+t,m) / scale(count+t);
        }
        
        //backward probs and state conditional probs
        for(m=0; m<M; m++) beta(count + ntimes(j) - 1,m) = 1 / (M * scale(count + ntimes(j) - 1));
        Gamma.row(count+ntimes(j)-1) = alpha.row(count+ntimes(j)-1) % beta.row(count+ntimes(j)-1);
        tempval = alpha.row(count+ntimes(j)-1) * beta.row(count+ntimes(j)-1).t();
        for(m=0; m<M; m++) Gamma(count+ntimes(j)-1,m)=Gamma(count+ntimes(j)-1,m)/tempval(0); //double
        
        for(t=ntimes(j)-2; t>=0; t--){
            tempmat = beta.row(count+t+1) % nodeprob.row(count+t+1);
            beta.row(count+t) = tempmat * P.t();
            for(m=0; m<M; m++) beta(count+t,m) = beta(count+t,m) / scale(count + t);
            
            Gamma.row(count+t) = alpha.row(count+t) % beta.row(count+t);
            tempval = alpha.row(count+t) * beta.row(count+t).t();
            for(m=0; m<M; m++) Gamma(count+t,m)=Gamma(count+t,m)/tempval(0); //double
        }
        
        //transition probs
        for(t=0; t<ntimes(j)-1; t++){
            tempmat = beta.row(count+t+1) % nodeprob.row(count+t+1);
            tempmat2 = P % (alpha.row(count+t).t() * tempmat);
            
            tempsum = 0;
            for(i=0; i<M; i++){
                for(k=0; k<M; k++){
                    xi(count+t, i + k*M) = tempmat2(i,k);
                    tempsum += xi(count+t, i + k*M);
                }
            }
            
            for(m=0; m<M*M; m++) xi(count+t, m) = xi(count+t, m) / tempsum;
        }
        
        count += ntimes(j);
    }
    
    
    //get the column sums
    colsumgamma = colsum(Gamma);
    tempsumxi = colsum(xi);
    
    for(i=0; i<M; i++)
    for(k=0; k<M; k++)
    colsumxi(i,k) = tempsumxi(i+k*M);
    
    loglik = arma::sum(log(scale));
    
    return Rcpp::List::create(Rcpp::Named("colsumgamma")=colsumgamma,
                              Rcpp::Named("colsumxi")=colsumxi,
                              Rcpp::Named("Gamma")=Gamma,
                              Rcpp::Named("xi")=xi,
                              Rcpp::Named("loglik")=loglik);
}


//////////////////////////
// generic function for loglikelihood for HMM
// [[Rcpp::export]]
double hmmllk(arma::vec delta, arma::mat gamma, arma::mat nodeprob,
              arma::vec y, arma::vec ntimes, arma::vec timeindex,int missing){
    double loglik = 0;
    //long dim = y.n_rows;
    int M = nodeprob.n_cols;
    int n = ntimes.n_rows;
    int j,m,t;
    double tempsum;
    
    arma::vec forward(M);
    arma::vec forwardsum;  //then forwardsum(0) is double
    arma::mat tempmat(1,M);
    arma::mat transit(M,M);
    int interval=1;
    
    int count = 0;
    //Rcpp::Rcout << arma::max(timeindex) << std::endl;
    //Rcpp::Rcout << y.n_rows << std::endl;
    
    if(missing==0){ //no missing
        
        for(j=0; j<n; j++){
            
            tempsum = 0;
            //initialize the forward variable
            forward = delta % nodeprob.row(count).t(); //element wise multiplication
            forwardsum = delta.t() * nodeprob.row(count).t();
            tempsum = log(forwardsum(0));
            //Rcpp::Rcout << "loglik=" << loglik << std::endl;
            
            for(m=0; m<M; m++) forward(m) = forward(m) / forwardsum(0);
            
            //recursion
            for(t=1; t<ntimes(j); t++){
                
                tempmat = forward.t() * gamma; //row matrix
                forward = tempmat.t() % nodeprob.row(count+t).t();
                forwardsum = tempmat * nodeprob.row(count+t).t();
                tempsum += log(forwardsum(0));
                for(m=0; m<M; m++) forward(m) = forward(m) / forwardsum(0);
            }
            
            loglik += tempsum;
            count += ntimes(j);
        }
    }
    else{ //missing
        for(j=0; j<n; j++){
            
            tempsum = 0;
            //initialize the forward variable
            forward = delta % nodeprob.row(count).t(); //element wise multiplication
            forwardsum = delta.t() * nodeprob.row(count).t();
            tempsum = log(forwardsum(0));
            //Rcpp::Rcout << "loglik=" << loglik << std::endl;
            
            for(m=0; m<M; m++) forward(m) = forward(m) / forwardsum(0);
            
            //recursion
            for(t=1; t<ntimes(j); t++){
                //calculate interval from last time point
                interval = timeindex(count+t) - timeindex(count+t-1);
                transit = matrixpower(gamma, interval);
                
                tempmat = forward.t() * transit; //row matrix
                forward = tempmat.t() % nodeprob.row(count+t).t();
                forwardsum = tempmat * nodeprob.row(count+t).t();
                tempsum += log(forwardsum(0));
                for(m=0; m<M; m++) forward(m) = forward(m) / forwardsum(0);
            }
            
            loglik += tempsum;
            count += ntimes(j);
        }
        
    }
    
    return loglik;
    
}



//////////////////////////////////////////////////////////////
/*part 3: continuous time hidden markov model*/
/*missing value interpolation*/
/*smoothing, soft-cluster*/

//////////////////////////
// generic function for loglikelihood for HMM
// [[Rcpp::export]]
double hmmllk_cont(arma::vec delta, arma::mat gamma, arma::mat nodeprob,
                   arma::vec y, arma::vec ntimes, arma::vec timeindex){
    double loglik = 0;
    //long dim = y.n_rows;
    int M = nodeprob.n_cols;
    int n = ntimes.n_rows;
    int j,m,t;
    double tempsum;
    
    arma::vec forward(M);
    arma::vec forwardsum;  //then forwardsum(0) is double
    arma::mat tempmat(1,M);
    arma::mat transit(M,M);
    int interval=1;
    
    int count = 0;
    //Rcpp::Rcout << arma::max(timeindex) << std::endl;
    //Rcpp::Rcout << y.n_rows << std::endl;
    
    
    for(j=0; j<n; j++){
        
        tempsum = 0;
        //initialize the forward variable
        forward = delta % nodeprob.row(count).t(); //element wise multiplication
        forwardsum = delta.t() * nodeprob.row(count).t();
        tempsum = log(forwardsum(0));
        //Rcpp::Rcout << 0 << "," << tempsum << std::endl;
        
        for(m=0; m<M; m++) forward(m) = forward(m) / forwardsum(0);
        
        //recursion
        for(t=1; t<ntimes(j); t++){
            //calculate interval from last time point
            interval = timeindex(count+t) - timeindex(count+t-1);
            //interval = MIN(interval, 50);
            transit = matrixexp(gamma, interval);
            
            tempmat = forward.t() * transit; //row matrix
            forward = tempmat.t() % nodeprob.row(count+t).t();
            forwardsum = tempmat * nodeprob.row(count+t).t();
            tempsum += log(forwardsum(0));
            for(m=0; m<M; m++) forward(m) = forward(m) / forwardsum(0);
            //Rcpp::Rcout << t << "," << tempsum << std::endl;
        }
        
        loglik += tempsum;
        count += ntimes(j);
    }
    
    return loglik;
    
}


//////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::List fb_cont(arma::vec Pi, arma::mat P, arma::mat nodeprob,
                   long dim, arma::vec ntimes, arma::vec timeindex){
    
    
    int M = nodeprob.n_cols;
    int n = ntimes.n_rows;
    int j,m,t,i,k;
    int interval;
    arma::mat transit(M,M);
    double tempsum;
    arma::vec tempval;
    arma::mat tempmat(1,M);
    arma::mat tempmat2(M,M);
    
    arma::mat alpha(dim, M);
    arma::vec scale(dim);
    arma::mat beta(dim, M);
    arma::mat Gamma(dim, M);
    arma::mat xi(dim-1, M*M);
    
    arma::vec colsumgamma(M);
    arma::vec tempsumxi(M*M);
    arma::mat colsumxi(M,M);
    double loglik = 0;
    
    int count = 0;
    for(j=0; j<n; j++){
        
        //forward probs
        alpha.row(count) = Pi.t() % nodeprob.row(count);
        tempval = nodeprob.row(count) * Pi;
        scale(count) = tempval(0);//double
        for(m=0; m<M; m++) alpha(count, m) = alpha(count,m) / scale(count);
        
        for(t=1; t<ntimes(j); t++){
            
            interval = timeindex(count+t) - timeindex(count+t-1);
            transit = matrixexp(P,interval);
            
            tempmat = alpha.row(count+t-1) * transit; //row matrix
            alpha.row(count+t) = tempmat % nodeprob.row(count+t);
            tempval = tempmat * nodeprob.row(count+t).t();
            scale(count+t) = tempval(0);//double
            for(m=0; m<M; m++) alpha(count+t, m) = alpha(count+t,m) / scale(count+t);
        }
        
        //backward probs and state conditional probs
        for(m=0; m<M; m++) beta(count + ntimes(j) - 1,m) = 1 / (M * scale(count + ntimes(j) - 1));
        Gamma.row(count+ntimes(j)-1) = alpha.row(count+ntimes(j)-1) % beta.row(count+ntimes(j)-1);
        tempval = alpha.row(count+ntimes(j)-1) * beta.row(count+ntimes(j)-1).t();
        for(m=0; m<M; m++) Gamma(count+ntimes(j)-1,m)=Gamma(count+ntimes(j)-1,m)/tempval(0); //double
        
        for(t=ntimes(j)-2; t>=0; t--){
            interval = timeindex(count+t+1) - timeindex(count+t);
            transit = matrixexp(P, interval);
            
            tempmat = beta.row(count+t+1) % nodeprob.row(count+t+1);
            beta.row(count+t) = tempmat * transit.t();
            for(m=0; m<M; m++) beta(count+t,m) = beta(count+t,m) / scale(count + t);
            
            Gamma.row(count+t) = alpha.row(count+t) % beta.row(count+t);
            tempval = alpha.row(count+t) * beta.row(count+t).t();
            for(m=0; m<M; m++) Gamma(count+t,m)=Gamma(count+t,m)/tempval(0); //double
        }
        
        
        
        //transition probs
        for(t=0; t<ntimes(j)-1; t++){
            interval = timeindex(count+t+1) - timeindex(count+t);
            transit = matrixexp(P, interval);
            
            tempmat = beta.row(count+t+1) % nodeprob.row(count+t+1);
            tempmat2 = transit % (alpha.row(count+t).t() * tempmat);
            
            tempsum = 0;
            for(i=0; i<M; i++){
                for(k=0; k<M; k++){
                    xi(count+t, i + k*M) = tempmat2(i,k);
                    tempsum += xi(count+t, i + k*M);
                }
            }
            //Rcpp::Rcout<<count+t<<std::endl;
            for(m=0; m<M*M; m++) xi(count+t, m) = xi(count+t, m) / tempsum;
        }
        
        count += ntimes(j);
    }
    
    
    //get the column sums
    colsumgamma = colsum(Gamma);
    tempsumxi = colsum(xi);
    
    for(i=0; i<M; i++)
    for(k=0; k<M; k++)
    colsumxi(i,k) = tempsumxi(i+k*M);
    
    loglik = arma::sum(log(scale));
    
    return Rcpp::List::create(Rcpp::Named("colsumgamma")=colsumgamma,
                              Rcpp::Named("colsumxi")=colsumxi,
                              Rcpp::Named("Gamma")=Gamma,
                              Rcpp::Named("xi")=xi,
                              Rcpp::Named("loglik")=loglik);
}


