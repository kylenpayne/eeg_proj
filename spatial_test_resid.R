### --------
## name: spatial_test_resid.R
## author: Kyle N. Payne
## purpose: estimate and test spatial correlation in the residuals
## of the mixed effects models estimated for Payne, 2015. 
## These results will help frame the need for spatial correlation
## structure in the estimation of $$\Sigma_{R}$$
### --------

# install.packages("geoR") run if needed

### --------
## download and import the data

library("devtools")
# install_github("leeper/rio")
library("rio")
centroParietal_OpenClass<-import("https://raw.githubusercontent.com/payne12/singleWordEEG/master/centroParietal_OpenClass.csv")
summary(centroParietal_OpenClass) #summarize data set
### --------


### --------
## convert the polar coordinates to cartesian coordiantes
N400Topo <- read.csv("/Volumes/HAL/kyle/eeg_proj/eeg_proj/N400Topo.txt", sep="")
colnames(N400Topo)[1] <- "chan" # renaming 'cause I'm stupid.
N400Topo$x_co <- N400Topo$radius*cos(N400Topo$deg)
N400Topo$y_co <- N400Topo$radius*sin(N400Topo$deg)
### --------

### -------
## merge the cartesian coordinates
centro_merged <- merge(centroParietal_OpenClass, N400Topo[, c(1, 4,5,6)], by="chan")
### -------


### -------
library(ggplot2)
p<- ggplot(aes(x=meanAmpEEG), data=centroParietal_OpenClass)
p + geom_histogram() + facet_wrap(~SubID)
## establishes normality assumption as the meanAmpEEG
## more or less just a validation of the central limit theorem.

q <- ggplot(aes(x=meanAmpEEG), data=centroParietal_OpenClass)
q + geom_histogram(aes(fill=factor(chan)))

q2 <- q + geom_histogram(aes(color=factor(SubID)))

## looking at normality of meanAmpEEG for subject ID and channel
library(ggplot2)
p2 <- ggplot(aes(x=meanAmpEEG), data=centroParietal_OpenClass)
p2 + geom_histogram() + facet_wrap(~SubID + chan)



### -------
## fit the model used in Payne, 2015 and examine the residuals for 
## spatial correlation.
require(lme4) 
lmeModel1 <-lmer(meanAmpEEG ~ zWordOrder + C1 + C2 
                 + zWordOrder:C1 + zWordOrder:C2  
                 + (1 | SubID)
                 + (1 | chan)
                 + (1 | words)
                 + (0 + zWordOrder:C1 | SubID)
                 + (0 + zWordOrder:C2 | SubID),  
                 data = centroParietal_OpenClass)


resid_lmeModel1 <- residuals(lmeModel1)

### ------
## grab the rows that correspond to the rows that
## do not have a meanAmpEEG value of NA
centro_merged_no_NA <-
  centro_merged[which(!is.na(centro_merged$meanAmpEEG)),]
resid_lmeModel1 <- cbind(resid_lmeModel1, centro_merged_no_NA[,
                      c("label","x_co", "y_co")])


### -------
## calculate the euclidean distance between each of 

dists <- dist(N400Topo[,5:6])
summary(dists)

### ------
## based on the summary of the distances
## between electrodes, I chose
## 10 lag intervals of size .09

breaks <- seq(0, .907, 10)

### ------
## now estimate the variogram 

library(sp)
library(gstat)
coordinates(resid_lmeModel1) <- ~ x_co + y_co

vario <- variogram(resid_lmeModel1 ~ 1, data=resid_lmeModel1, cloud=TRUE)
vario.fit = fit.variogram(vario, model=vgm(10, "Sph", 300))