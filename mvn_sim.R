###  ---------
## name: mvn_sim.R
## Author: Kyle N. Payne
## Purpose: simulate sets of multivariate normal random variables in
## order to test the effects of autocorrelated-ness in the 
## random effects structure data.
###  ---------

source("helpers.R")

library(mvtnorm)
library(reshape2)
library(eegkit)
data(eegcoord)

### --------- Generating spatially correlated random fields
N <- 1000
q <- nrow(eegcoord)
nug <- 0.05
silly <- 0.01
phi <- 0.01
locs <- c("xproj","yproj")

### ---- returns the cholesky decomposition of the 
exp_covar <- cov_exp(nugget = nug, 
                           sill = silly, 
                           phi= phi,
                           locations = locs,
                           data = eegcoord
                           )

## create means vectors 
mu_1 <- rnorm(q, -1 ,1) 
mu_1 <- rep.row(mu_1, N)

mu_2 <- rnorm(q, -1, 1)
mu_2 <- rep.row(mu_2, N)

### ---- induce a pocket of electrical activity
### ---- that is clustered about the centroParietal sites

cp_labels <- which(substr(row.names(eegcoord),
                            start=1,
                            stop=1) == "C")

mu_2[,cp_labels] <- rnorm(length(cp_labels), -4, 2)

## ---- create fake EEG data with different mean shifts
## ---- but with identical covariances
y_1 <- rmvnorm(N, mean = rep(0, q), sigma=exp_covar) + mu_1
y_2 <- rmvnorm(N, mean = rep(0, q), sigma=exp_covar) + mu_2

## ---- fake meanAmpEEG electrode data
data_sim <- rbind(y_1, y_2)
data_sim <- cbind(data_sim, c(rep(0, N), rep(1, N)))
data_sim_df <- data.frame(data_sim)
colnames(data_sim_df) <- c(row.names(eegcoord), "c2")
data_sim_df_m <- melt(data_sim_df, .(c2))
colnames(data_sim_df_m) <- c("c2","chan", "meanAmp")

## ---- let's look at a scalp-map lattice of the meanAmpEEG
## ---- averages for each channel and across the two
## ---- conditions
library(plyr)

data_sim_p <- ddply(data_sim_df_m,
                    .(c2, chan), summarize,
                  mean=mean(meanAmp))

eegcoord_edit <- cbind(eegcoord, row.names(eegcoord))
colnames(eegcoord_edit)[ncol(eegcoord_edit)] <- "chan"

data_sim_p_mrg <- merge(data_sim_p, eegcoord_edit, by="chan")

for(i in 1:2){
  eegspace(data_sim_p_mrg[which(data_sim_p_mrg$c2 == i-1),7:8],
           voltage = data_sim_p_mrg[which(data_sim_p_mrg$c2 == i-1),"mean"], 
           vlim=c(min(data_sim_p_mrg$mean),max(data_sim_p_mrg$mean)),
           main=paste("Condition", i),cex.main=2)
}
