###  ---------
## name: mvn_data_sim.R
## Author: Kyle N. Payne
## Purpose: simulate sets of multivariate normal random variables in
## order to test the effects of autocorrelated-ness in the 
## random effects structure data.
###  ---------

### let's try eegkit first

require("eegkit")

### get the unique channel site labels from
### centroParietal_OpenClass, then cross-reference
### with that in N400Topo

chan_num <- unique(centroParietal_OpenClass$chan)
chan_labels <- N400Topo[which(N400Topo$chan %in% chan_num), "label"]


chnames <- rownames(eegcoord)
quartz(width=15,height=3)

## -----
## ----- simulate erp EEG data
## ----- 
tseq <- seq(0, 500, by=10)/1000
cps <- row.names(eegcoord)
sim_erp_data <- matrix(NA, ncol=length(cps))
colnames(sim_erp_data) <- cps

for(j in 1:length(tseq)){
 tmp <- eegsim(channel = cps, rep(tseq[j], length(cps)))
 tmp <- t(tmp)
 sim_erp_data <- rbind(sim_erp_data, tmp)
}

#### ----- highly spatially correlated
sim_erp_data <- sim_erp_data[-1,]

cond_2 <- matrix(NA, ncol=length(cps))

for(j in 1:length(tseq)){
  tmp <- eegsim(channel = cps, rep(tseq[j], length(cps)))
}


#### ----- eeg_coords

data(eegcoord)

### ----- using the gstat package
library(gstat)
library(reshape2)
sim_erp_data <- melt(sim_erp_data)

colnames(sim_erp_data) <- c("obs", "chan", "amp")
### ---- get a list of all of the channel names
chan_cp <- unique(sim_erp_data$chan)
xy <- eegcoord[which(row.names(eegcoord) %in% chan_cp), c("xproj", "yproj")]
xy$chan <- row.names(xy)

### ----- merge the location data and the erp data
sim_erp_data <- merge(sim_erp_data, xy, by="chan")

### ---- create a gstat object for the amplitude
### ---- krigged on location and time
g.sim <- gstat(formula=amp ~ 1, 
               beta=1,
               dummy=T,
               locations = ~ xproj + yproj, 
               nmax=4,
               model=vgm(psill=0.5, model="Sph", range=4))

num_sim <- 5
sim_data <- predict(g.sim, newdata=xy, nsim=num_sim)
gridded(sim_data) <- ~ xproj + yproj
spplot(sim_data[1])

for(i in 1:num_sim){
  eegspace(sim_data[,1:2], voltage = sim_data[,2+i], 
          vlim=c(-6.8,5.5),
          main=paste("sim", i),cex.main=2)
}

### ---- let's now try to re-create the distributions 
### ---- of the Payne, 2015 meanAmp data but induce 
### ---- a spatial correlational structure via the 
### ---- gstat function using centro_merged



  


