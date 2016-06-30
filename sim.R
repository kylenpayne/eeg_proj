    ###  ---------
    ## name: mvn_sim.R
    ## Author: Kyle N. Payne
    ## Purpose: simulate sets of multivariate normal random variables in
    ## order to test the effects of autocorrelated-ness in the 
    ## random effects structure data.
    ###  ---------
    

source("helpers.R")
require(foreach)
require(lme4)
require(plyr)
require(shiny)
require(fdrtool)
require(mvtnorm)
require(reshape2)
require(eegkit)
data(eegcoord)

    
### ---- list of inputs and outputs for the program
### ---- inputs:
### ---- min_sill
### ---- max_sill
### ---- min_decay
### ---- max_decay
### ---- num_steps (decay and sill evaluation points)
### ---- region (region of the brain)
### ---- effect (effect size)
### ---- nugget
### ---- alpha 
### ---- num_subs
### ---- num_iter
### ---- samp_size
### ---- outputs: 
### ---- dot_map

### ---- server part of the program.
min_sill <- 0;
max_sill <- 10;
min_decay <- 0;
max_decay <- 10;
num_steps <- 20;
num_iter <- 50;
region <- "F"
# effect <- 0.
min_effect <- 0;
max_effect <- 10;
sd_effect <- 0.05;
alpha <- 0.05
covar_type <- "expo" # input$covar_type
### ---- number of observations
D <- 50
### ---- number of subjects
N <- 10
### ---- number of electrodes
q <- nrow(eegcoord)
nug <- 0
silly <- seq(from=min_sill, to=max_sill, length.out=num_steps)
decay <- seq(from=min_decay, to=max_decay, length.out=num_steps)
locs <- c("xproj","yproj")
    
## for each level of the effect size, run a simulation
effect <- seq(from=min_effect, to=max_effect, length.out=num_steps)

plot_list <- vector("list", num_steps)
    
i <- 0

for(eff in effect){
  i <- i + 1;
  ## ---- create a matrix of lists that contains
  ## ---- all of the different values of the psill and
  ## ---- of the decay
  
  params <- matrix(rep(list(), 
                       length(silly)*length(decay)), 
                   nrow=length(silly), ncol=length(decay))
  p_values <- matrix(rep(list(),
                         length(silly)*length(decay)),
                     nrow=length(silly), ncol=length(decay))
  
      
  
  
  ## fill the paramter matrix of lists
  for(s in 1:length(silly)){
    for(d in 1:length(decay)){
      params[s,d][[1]] <- c(silly[s], decay[d])
    }
  }  
  
  library(doParallel)
  cl <- makeCluster(3, outfile="")
  registerDoParallel(cl)
  
  opts <- list(chunkSize = 2)
  covar_list <- foreach(s=silly, .options.nws=opts) %:%
    foreach(d=decay) %dopar% {
      if(covar_type=="expo"){
        ### ---- returns the exponential covariance
        covar <- cov_exp(nugget = nug, 
                         psill = s, 
                         decay= d,
                         locations = locs,
                         data = eegcoord)
      }
      if(covar_type=="sphere"){
        ### ---- returns the spherical covariance
        covar <- cov_sphere(nugget = nug,
                            psill = s,
                            decay = d,
                            locations = locs,
                            data=eegcoord)
      }
      if(covar_type=="gauss"){
        ### ---- returns the gaussian covariance
        covar <- cov_gauss(nugget = nug,
                           psill=s,
                           decay=d,
                           locations=locs,
                           data=eegcoord)
      }
      
      ## ----- check to make sure that the 
      ## ---- covariance matrix is positive definite
      ## ---- if not, set the p-value list to 
      ## ---- NA and move onto the next iteration
      if(any(is.na(covar))){
        covar <- matrix(NA, ncol=nrow(eegcoord), nrow=nrow(eegcoord))
      }
      covar
    } ## --- end of the dopar loop
  
  ## ----- empty data frame to hold the simulations for
  ## ----- all numbers of participants
  # on second thought, this is probably not neccesary. For 
  # each iteration, there is a data_sim_final, that 
  # certainly does not need to be grown. 
  # This is why the results get so weird, 
  
  "
  data_sim_final <- data.frame(chan=character(),
                               c2=numeric(),
                               mean=numeric(),
                               xproj=numeric(),
                               yproj=numeric(),
                               subject=numeric())
  "
  ## 
  ## ----- we use this dataset for coordinates of eegcap locations
  eegcoord_edit <- cbind(eegcoord, row.names(eegcoord))
  eegcoord_edit <- eegcoord_edit[,-c(1,2,3)]
  colnames(eegcoord_edit)[ncol(eegcoord_edit)] <- "chan"
  ## ---- initialize a NULL vector that we will fill 
  ## ---- with the proportion of times that the
  ## ---- corresponding pval from the t-test
  ## ---- of the c2 effect is less than 0.05
  ### ---- create an electrode effect
  electrodes <- row.names(eegcoord)
  elect_labels <- which(substr(electrodes, 1,1) == region)
  
  
  stopCluster(cl)
  cl <- makeCluster(3, outfile="")
  registerDoParallel(cl)
  
  # N is the number of subjects
  t_list <- 
    foreach(s=1:num_steps,.verbose=TRUE,
            .packages=c('mvtnorm','reshape2',
                        'lme4','eegkit', 'foreach'
                        ,'plyr')) %:%
      foreach(d=1:num_steps) %dopar% {
      
        ## ---- create means vectors 
        mu_1 <- rmvnorm(num_iter, mean=rep(0, q), diag(q))
        mu_2 <- rmvnorm(num_iter, mean=rep(0, q), diag(q))
        
  
        
        ## ----- for all of the mu_2 rows, 
        ## ----- add in the electrode effect
        for(row in 1:num_iter){
          mu_2[row,elect_labels] <- rnorm(length(elect_labels), -1*eff, sd_effect)
        }
        
      stat_list <- vector(mode = "list", length = num_iter)
      ## ---- create a subject ID vector
      subj_id <- rep(rep(1:N, each=D),2)
      subject <- factor(subj_id) 
      
      
      
      for(i in 1:num_iter){
        ### ---- induce a pocket of electrical activity
        ### ---- that is clustered about the selected electrode sites
        ## ---- create fake EEG data with different mean shifts
        ## ---- but with identical covariances
        sig <- as.matrix(covar_list[[s]][[d]])
        if(any(is.na(sig))){
          t_stats <- NA
          next
        }
        y_1 <- rmvnorm(N*D, mean = mu_1[i,], sigma=sig) 
        y_2 <- rmvnorm(N*D, mean = mu_2[i,], sigma=sig) 
  
        ## ---- forming the response
        data_sim <- rbind(y_1, y_2)
        data_sim <- cbind(data_sim, c(rep(0, N*D), rep(1, N*D)),subject)
        data_sim_df <- data.frame(data_sim)
        colnames(data_sim_df) <- c(row.names(eegcoord), "c2", "subject")
        data_sim_df_m <- melt(data_sim_df, id=.(c2, subject))
        colnames(data_sim_df_m) <- c("c2","subject", "chan", "meanAmp")
        data_sim_df_m_merg <- merge(data_sim_df_m, eegcoord_edit, by="chan")
        data_sim_final <- data_sim_df_m_merg
        ### ---- a little bit of garbage collection
        ### ---- fit the mixed effects of model
        model_mixed <- lmer(meanAmp ~ c2 + (1|chan) + (1|subject), 
                            data=data_sim_final)
        t_stats <- summary(model_mixed)
        t_stats <- coefficients(t_stats)[,3]
        
        ## --- get all of the t-statistics for each model term
        stat_list[[i]] <- t_stats
        
        ## t_vector[i] <- t_stats[2]
          ### ----- CHANGE in strategy, as there is more considerable 
          ### ----- worry in reporting p-values, let's just return the
          ### ---- the corresponding t statistics and calculate
          ### ---- the fdr using fdr tools, that way we circumvent 
          ### ---- the issue with p-value approximations
          ## ---- normal approximation appears to be slightly anti-conservative.
          ##  p_stats <- 1-pnorm(abs(t_stats))
          ## ---- just to debug
          print(t_stats)
          ## ---- return the p-values for the "c2" (mean effect)
          ### p_stats[2] 
         
         ## ---- remove the following to try to free up some memory
         rm(y_1, y_2, 
            data_sim, 
            data_sim_df, 
            data_sim_df_m, 
            data_sim_df_m_merg, 
            data_sim_final,
            sig,
            t_stats,
            model_mixed)
            
      }    
      return(stat_list)
      gc(verbose=TRUE)
      rm(mu_1, mu_2, stat_list)
    } ## ---- sill and decay dopar loops
  
  stopCluster(cl)
  
  
  
  ## ---- 
  stopCluster(cl)
  cl <- makeCluster(8, outfile="")
  registerDoParallel(cl)
      
  t_means <- matrix(NA, ncol=num_steps, nrow=num_steps)
    for(s in 1:num_steps){
      for(d in 1:num_steps){
        list_of_ts <-t_list[[s]][[d]]
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
  p1 + geom_tile(aes(fill=t_mean)) + scale_fill_gradient(low="orange", high="blue")
  
  png(filename=paste("means_plot", eff, ".png"))
  plot(p1)
  dev.off()
  
  
  t_props <- matrix(NA, ncol=num_steps, nrow=num_steps)
      for(s in 1:num_steps){
        for(d in 1:num_steps){
          list_of_ts <-t_list[[s]][[d]]
          ts<-numeric(num_iter)
          for(i in 1:num_iter){
            tmp <- list_of_ts[[i]]
            if(is.null(tmp[2])){
              ts[i] <- NA
            }else{
              ts[i] <- tmp[2]
            }
          }
          t_props[s,d] <- mean(abs(ts) >= 2)
        }  
      }
  
  row.names(t_props) <- silly
  colnames(t_props) <- decay
  t_props_melt <- melt(t_props)
  colnames(t_props_melt) <- c("psill", "decay", "t_props")
  t_props_melt <- t_props_melt[-c(1:20),]
  p2 <- ggplot(aes(x=psill, y=decay), data=t_props_melt)
  p2 <- p2 + geom_tile(aes(fill=t_props)) + scale_fill_gradient(low="orange", high="blue")
      
  png(filename=paste("prop_plot", eff, ".png"))
  plot(p2)
  dev.off()
  
  #plot_list[[i]] <- list(p1, p2)
  
}
  