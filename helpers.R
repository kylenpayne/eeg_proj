###  ---------
## name: helpers.R
## Author: Kyle N. Payne
## Purpose: helper functions
###  ---------

cov_exp <- function(nugget, psill, decay, locations, data){
  ## let the locations be a 1 x 2 matrix of the x and y location 
  ## projection names. For the purposes of this project, this 
  ## will use the eegcoord data
  
  d <- dist(data[, locations])
  d <- as.matrix(d)
  covar <- matrix(NA, ncol=ncol(d), nrow=nrow(d))
  covar[which(d > 0), arr.ind=T] <- psill*exp(-decay*d^2)
  covar[which(d==0)] <- nugget + psill
  
  
  if(any(eigen(covar)$values < 0)){
    stop("The covariance matrix is non-positive definite")
  } else{
    ### ---- calculate the cholesky decomposition
    return(covar)
  }
}

## ---- spherical covariance

cov_sphere <- function(nugget, psill, decay, locations, data){
  d <- dist(data[, locations])
  d <- as.matrix(d)
  covar <- matrix(NA, ncol=ncol(d), nrow=nrow(d))
  covar[which(d > 1/decay), arr.ind=T] <- 0
  covar[which(d > 0 & d < 1/decay, arr.ind = T)] <-
    psill*(1- 3/2*decay*d + 1/2*(decay*d)^3)
  covar[which(d==0)] <- nugget + psill
  
  if(any(eigen(covar)$values < 0)){
    stop("The covariance matrix is non-positive definite")
  } else{
    ### ---- calculate the cholesky decomposition
    return(covar)
  }
}

cov_gauss <- function(nugget, psill, decay, locations, data){
  d <- dist(data[, locations])
  d <- as.matrix(d)
  covar <- matrix(NA, ncol=ncol(d), nrow=nrow(d))
  covar[which(d > 0), arr.ind=T] <- psill*exp(-decay^2*d^2)
  covar[which(d==0)] <- nugget + psill
  
  if(any(eigen(covar)$values < 0)){
    stop("The covariance matrix is non-positive definite")
  } else{
    ### ---- calculate the cholesky decomposition
    return(covar)
  }
}


## ---- repeats a row n many times
rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}
