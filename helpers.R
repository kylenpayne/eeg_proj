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
  t <- which(d > 0, arr.ind=T)
  covar[t] <- psill*exp(-decay*d[t]^2)
  covar[which(d==0)] <- nugget + psill
  
  
  if(any(eigen(covar)$values < 0)){
    print("The covariance matrix is non-positive definite")
    return(NA)
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
  covar[which(d > 1/decay, arr.ind=T)] <- 0
  t<-which(d > 0 & d < 1/decay, arr.ind = T)
  covar[t] <-
    psill*(1- 3/2*decay*d[t] + 1/2*(decay*d[t])^3)
  covar[which(d==0)] <- nugget + psill
  
  if(any(eigen(covar)$values < 0)){
    print("The covariance matrix is non-positive definite")
    return(NA)
  } else{
    ### ---- calculate the cholesky decomposition
    return(covar)
  }
}

cov_gauss <- function(nugget, psill, decay, locations, data){
  d <- dist(data[, locations])
  d <- as.matrix(d)
  covar <- matrix(NA, ncol=ncol(d), nrow=nrow(d))
  t <- which(d > 0,  arr.ind=T)
  covar[t] <- psill*exp(-decay^2*d[t]^2)
  covar[which(d==0, arr.ind=T)] <- nugget + psill
  
  if(any(eigen(covar)$values < 0)){
    print("The covariance matrix is non-positive definite")
    return(NA)
  } else{
    ### ---- return the covariance matrix
    return(covar)
  }
}


## ---- repeats a row n many times
rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}

## ---- returns a ggplot of the mean t value for all s,d values
means_plot <- function(t_list, num_steps, type){
  t_means_mixed <- matrix(NA, ncol=num_steps, nrow=num_steps)
  for(s in 1:num_steps){
    for(d in 1:num_steps){
      list_of_ts <-t_list[[s]][[d]][[type]]
      ts<-numeric(num_iter)
      for(i in 1:num_iter){
        print(paste(i, s, d, sep=":"))
        tmp <- list_of_ts[[i]]
        if(is.null(tmp[2])){
          ts[i] <- NA
        }else{
          ts[i] <- tmp[2]
        }
      }
      t_means[s,d] <- mean(ts)
    }  
  }
  
  row.names(t_means) <- silly
  colnames(t_means) <- decay
  t_means_melt <- melt(t_means)
  colnames(t_means_melt) <- c("psill", "decay", "t_mean")
  t_means_melt <- t_means_melt[-c(1:20),]
  library(ggplot2)
  
  p1 <- ggplot(aes(x=psill, y=decay), data=t_means_melt)
  p1 <- p1 + geom_tile(aes(fill=t_mean)) + scale_fill_gradient(low="orange", high="blue")
  return(p1)
}


## ---- returns a ggplot of the proportion t value for all s,d values

props_plot <- function(t_list, num_steps, crit = 2,  type){

  t_props_mixed <- matrix(NA, ncol=num_steps, nrow=num_steps)
  for(s in 1:num_steps){
    for(d in 1:num_steps){
      list_of_ts <-t_list[[s]][[d]][[type]]
      ts<-numeric(num_iter)
      for(i in 1:num_iter){
        tmp <- list_of_ts[[i]]
        if(is.null(tmp[2])){
          ts[i] <- NA
        }else{
          ts[i] <- tmp[2]
        }
      }
      t_props[s,d] <- mean(abs(ts) >= crit)
    }  
  }
  
  row.names(t_props) <- silly
  colnames(t_props) <- decay
  t_props_melt <- melt(t_props)
  colnames(t_props_melt) <- c("psill", "decay", "t_props")
  t_props_melt <- t_props_melt[-c(1:20),]
  p2 <- ggplot(aes(x=psill, y=decay), data=t_props_melt)
  p2 <- p2 + geom_tile(aes(fill=t_props)) + scale_fill_gradient(low="orange", high="blue")
  return(p2)
}
  

