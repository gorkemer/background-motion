# gorkemer
# 5 August 
# analyze on the cleaned data (output of the prepareData_01.R)
# extras: #colnames(group_mean) <- c("Mean Response AR", "CuedAR", "id")
### Load libraries ###

library(reshape2); library(rlang); library(data.table); library(plyr) # data formatting
library(psych); library(pastecs) # statistics
library(lme4); library(car); library(ggpubr); library(corrplot); library(magrittr) # regression
library(ggplot2); library(Hmisc) # graphing
library(dplyr) # data wrangling
library("reshape2");
library("scatterplot3d")
library(fitdistrplus); library(MASS); library(survival); library(npsurv) ; library(lsei)
library("plot3D")
library(extrafont)
font_import("Trebuchet MS")
library(actuar)
library(statmod)
library("car")

#remove scientific notation in the entire R session
options(scipen = 100)

# set wd if needed setwd('~/Desktop/background-motion/Data Analysis'), check with getwd()
setwd('~/Desktop/background-motion/Data Analysis')
#data = read.csv('bgm_cleanData_cor.csv', header=TRUE)
data = read.csv('untouchedData.csv', header = TRUE)
#data = read.csv('bgm_cleanData_lm.csv', header = TRUE)

# define the data frame to use
data.Same <- subset(data, data$pairAR_sameCondition==1)
data.Diff <- subset(data, !(data$pairAR_sameCondition==1))
#test the local effects

#first plot the relationship between cuedAR and responseAR
testPlots <- ggplot(data.Diff, aes(x = cuedAR, y = responseAR)) + 
  geom_point(shape=19, size=1, alpha=0.1, show.legend = FALSE) +
  geom_smooth(data = data.Diff, method=lm, color = "red") +
  ylim(-.7,.7) +
  labs(x="Cued aspect ratio", y = "Reported aspect ratio") +
  ggtitle("Cleaned data") #+ facet_wrap(~id,nrow = 5)

testPlots + geom_density_2d(color="blue", alpha=0.2) + xlab((testPlots$mapping$x)) 

#get the aggregated response to cuedAR *remove id for running on the whole set
data_agg_same <- aggregate(responseAR ~ cuedAR, data = data.Same, mean)
data_agg_diff <- aggregate(responseAR ~ cuedAR, data = data.Diff, mean)
difference <- data_agg_diff - data_agg_same

t.test(difference, mu = 0) # the difference is not significantly different than 0
t.test(data_agg_diff, data_agg_same) # same and different is not significantly different

#### same-all & diff-all ####
#substracting aggregated (all participants) cuedAR & uncuedAR from identical 
#get the aggregated response to cuedAR & uncuedAR combo for all *remove id for running on the whole set
data_agg_same <- aggregate(responseAR ~ cuedAR + uncuedAR, data = data.Same, mean)
data_agg_diff <- aggregate(responseAR ~ cuedAR + uncuedAR, data = data.Diff, mean)

data.Difference <- data.frame(id=1:nrow(data_agg_diff), uncuedAR=NA, cuedAR =NA, responseAR = NA,normed_responseAR = NA, responseToIdentical = NA, selectedCuedARofIdentical =NA)
#,

for (i in 1:nrow(data_agg_same)) { # get the response to the same
  # get the id of a given trial
  #selectedID <- data_agg_same$id[i]
  selectedCuedAR <- data_agg_same$cuedAR[i]
  selectedResponseARtoSame <- data_agg_same$responseAR[i] #this is the response to the identical pair
  #look through the data and match the cuedAR with the data
  for (a in 1:nrow(data_agg_diff)){
    if (selectedCuedAR == data_agg_diff$cuedAR[a]){ #&& selectedID == data.differents$id[a]){ #find the cuedAR and id couplings and substract response to identical from that response on that trial
      #copy the results from data_agg_diff
      data.Difference$uncuedAR[a]<- data_agg_diff$uncuedAR[a] #write the uncuedAR
      data.Difference$cuedAR[a]<- data_agg_diff$cuedAR[a] #write the cuedAR
      data.Difference$responseAR[a]<-data_agg_diff$responseAR[a] 
      #write down the norm responses and cross checks
      data.Difference$normed_responseAR[a] <- data_agg_diff$responseAR[a] - selectedResponseARtoSame #positive values indicate deviation from the cuedAR towards "taller" responses
      data.Difference$responseToIdentical[a] <- selectedResponseARtoSame
      data.Difference$selectedCuedARofIdentical[a]<- selectedCuedAR #this should match with the trial's cuedAR

    }
  }
}

### 21 August after we talked with Tim about mixed approach, I came back here
# to do the population-response-version of the normalization analysis
# I do it only for coherent
data.Difference.Coherent <- subset(data.Difference, data.Difference$)
lm_motion <- lm(normed_respnseAR ~ uncuedAR * pairMot)


# run regression to see uncuedAR predicts norm_responseAR for aggregated responses to uncuedAR (all participants)
lm_main <- lm(normed_responseAR ~  cuedAR * uncuedAR , data.Difference)
summary(lm_main) # no neither of them predicts the variance in the normed_responseAR
# trying out Max's suggestion: 
lm_main <- lm(responseAR ~  cuedAR * uncuedAR + responseToIdentical , data.Difference)
summary(lm_main) # this way cuedAR still predicts the response, even after controlling for the responses to the same
#look at the 3-way interaction at a plot
# x, y and z coordinates
x <- data.Difference$uncuedAR
y <- data.Difference$normed_responseAR
z <- data.Difference$cuedAR

# plot at 3d graph
scatter3D(x, y, z, phi = 0, bty ="g")
scatterplot3d(x = data.Difference$uncuedAR, y = data.Difference$normed_responseAR, z = data.Difference$cuedAR,
              main="3D Scatter Plot",
              xlab = "uncuedAR",
              ylab = "responseAR",
              zlab = "cuedAR")

scatter3d(x = data.Difference$uncuedAR, y = data.Difference$normed_responseAR, z = data.Difference$cuedAR,
          point.col = "blue", surface=FALSE, 
          xlab = "uncuedAR", 
          ylab = "responseAR",
          zlab = "cuedAR") #groups = iris$Species # add concentration ellipses: ellipsoid = TRUE

# no apparent relationship

#### same-all (no-id) & different-id ####
#substracting invididual responses (accounting uncuedAR) by responses to identical obtained by averaging all. 
#get the aggregated response to cuedAR & uncuedAR & id
data_agg_same <- aggregate(responseAR ~ cuedAR + uncuedAR, data = data.Same, mean)
data_agg_diff <- aggregate(responseAR ~ cuedAR + uncuedAR + id, data = data.Diff, mean)

data.Difference <- data.frame(id=1:nrow(data_agg_diff), partId = NA, uncuedAR=NA, cuedAR =NA, responseAR = NA,normed_responseAR = NA, responseToIdentical = NA, selectedCuedARofIdentical =NA)
#,

for (i in 1:nrow(data_agg_same)) { # get the response to the same
  # get the id of a given trial
  selectedID <- data_agg_same$id[i]
  selectedCuedAR <- data_agg_same$cuedAR[i]
  selectedResponseARtoSame <- data_agg_same$responseAR[i] #this is the response to the identical pair
  #look through the data and match the cuedAR with the data
  for (a in 1:nrow(data_agg_diff)){
    if (selectedCuedAR == data_agg_diff$cuedAR[a]){
      #copy the results from data_agg_diff
      data.Difference$partId[a] <- data_agg_diff$id[a] #overwriting on the id
      data.Difference$uncuedAR[a]<- data_agg_diff$uncuedAR[a] #write the uncuedAR
      data.Difference$cuedAR[a]<- data_agg_diff$cuedAR[a] #write the cuedAR
      data.Difference$responseAR[a]<-data_agg_diff$responseAR[a] 
      #write down the norm responses and cross checks
      data.Difference$normed_responseAR[a] <- data_agg_diff$responseAR[a] - selectedResponseARtoSame #positive values indicate deviation from the cuedAR towards "taller" responses
      data.Difference$responseToIdentical[a] <- selectedResponseARtoSame
      data.Difference$selectedCuedARofIdentical[a]<- selectedCuedAR #this should match with the trial's cuedAR
    }
  }
}

# run regression to see uncuedAR predicts norm_responseAR for aggregated responses to uncuedAR (all participants)
lm_main <- lm(normed_responseAR ~  cuedAR * uncuedAR , subset(data.Difference,data.Difference$partId=="o1"))
summary(lm_main) # no neither of them predicts the variance in the normed_responseAR

#total id list
totalIdList <- length(unique(data.Difference$partId))

#run individual regression for every participant
data.indvRegScores <- data.frame(id=1:totalIdList, partId = NA, betaCued = NA, betaUncued = NA, betaInteraction = NA, PvalueCued = NA, PvalueUncued = NA, PvalueInteraction = NA)

for (i in 1:totalIdList){
  selectedID <- unique(data.Difference$partId)[i]
  print(selectedID)
  #run regression analysis
  lm <- lm(normed_responseAR ~  cuedAR * uncuedAR , subset(data.Difference,data.Difference$partId==selectedID))
  
  data.indvRegScores$partId[i] <- selectedID
  data.indvRegScores$betaCued[i] <- lm$coefficients[2] #cuedAR
  data.indvRegScores$betaUncued[i] <- lm$coefficients[3] #cuedAR
  data.indvRegScores$betaInteraction[i] <- lm$coefficients[3] #cuedAR
  
  data.indvRegScores$PvalueCued[i]   <- summary(lm)$coefficients[,4][2]
  data.indvRegScores$PvalueUncued[i]   <- summary(lm)$coefficients[,4][3]
  data.indvRegScores$PvalueInteraction[i]   <- summary(lm)$coefficients[,4][4]
}
##looking at the individual regression slopes and pvalues, I can see that some showed
##repulsion from uncuedAR (negative beta). However, ignoring the motion direction's role
## is not wise because on average the influence of uncuedAR might be nulled by the motion
## direction's role. 
# find number of participant with uncued beta is significant than 0.05
data.indvRegScores_subset <- subset(data.indvRegScores, data.indvRegScores$PvalueUncued<0.05)
# only 13 participants out of 77 showed signficant beta, and 8/13 showed repulsion

#### same-all & data (same&diff) -id ####
#get the aggregated response to cuedAR & uncuedAR & id
data_agg_same <- aggregate(responseAR ~ cuedAR + uncuedAR, data = data.Same, mean)
data_agg_data <- aggregate(responseAR ~ cuedAR + uncuedAR + id, data = data, mean)

data.Difference <- data.frame(id=1:nrow(data_agg_data), partId = NA, uncuedAR=NA, cuedAR =NA, responseAR = NA,normed_responseAR = NA, responseToIdentical = NA, selectedCuedARofIdentical =NA)
#,

for (i in 1:nrow(data_agg_same)) { # get the response to the same
  # get the id of a given trial
  selectedID <- data_agg_same$id[i]
  selectedCuedAR <- data_agg_same$cuedAR[i]
  selectedResponseARtoSame <- data_agg_same$responseAR[i] #this is the response to the identical pair
  #look through the data and match the cuedAR with the data
  for (a in 1:nrow(data_agg_data)){
    if (selectedCuedAR == data_agg_data$cuedAR[a]){
      #copy the results from data_agg_diff
      data.Difference$partId[a] <- data_agg_data$id[a] #overwriting on the id
      data.Difference$uncuedAR[a]<- data_agg_data$uncuedAR[a] #write the uncuedAR
      data.Difference$cuedAR[a]<- data_agg_data$cuedAR[a] #write the cuedAR
      data.Difference$responseAR[a]<-data_agg_data$responseAR[a] 
      #write down the norm responses and cross checks
      data.Difference$normed_responseAR[a] <- data_agg_data$responseAR[a] - selectedResponseARtoSame #positive values indicate deviation from the cuedAR towards "taller" responses
      data.Difference$responseToIdentical[a] <- selectedResponseARtoSame
      data.Difference$selectedCuedARofIdentical[a]<- selectedCuedAR #this should match with the trial's cuedAR
    }
  }
}

# run regression to see uncuedAR predicts norm_responseAR for aggregated responses to uncuedAR (all participants)
lm_main <- lm(normed_responseAR ~  cuedAR * uncuedAR , data.Difference)
summary(lm_main) # no neither of them predicts the variance in the normed_responseAR

#total id list
totalIdList <- length(unique(data.Difference$partId))

#run individual regression for every participant
data.indvRegScores <- data.frame(id=1:totalIdList, partId = NA, betaCued = NA, betaUncued = NA, betaInteraction = NA, PvalueCued = NA, PvalueUncued = NA, PvalueInteraction = NA)

for (i in 1:totalIdList){
  selectedID <- unique(data.Difference$partId)[i]
  print(selectedID)
  #run regression analysis
  lm <- lm(normed_responseAR ~  cuedAR * uncuedAR , subset(data.Difference,data.Difference$partId==selectedID))
  
  data.indvRegScores$partId[i] <- selectedID
  data.indvRegScores$betaCued[i] <- lm$coefficients[2] #cuedAR
  data.indvRegScores$betaUncued[i] <- lm$coefficients[3] #cuedAR
  data.indvRegScores$betaInteraction[i] <- lm$coefficients[3] #cuedAR
  
  data.indvRegScores$PvalueCued[i]   <- summary(lm)$coefficients[,4][2]
  data.indvRegScores$PvalueUncued[i]   <- summary(lm)$coefficients[,4][3]
  data.indvRegScores$PvalueInteraction[i]   <- summary(lm)$coefficients[,4][4]
}
##looking at the individual regression slopes and pvalues, I can see that some showed
##repulsion from uncuedAR (negative beta). However, ignoring the motion direction's role
## is not wise because on average the influence of uncuedAR might be nulled by the motion
## direction's role. 
# find number of participant with uncued beta is significant than 0.05
data.indvRegScores_subset <- subset(data.indvRegScores, data.indvRegScores$PvalueUncued<0.05)
# only 10 participants out of 77 showed signficant beta, and 5/13 showed repulsion


#### same-all & different-id-motionSame ####

#get the aggregated response to cuedAR & uncuedAR & id & motionSame
data_agg_same <- aggregate(responseAR ~ cuedAR + uncuedAR, data = data.Same, mean)
data_agg_diff <- aggregate(responseAR ~ cuedAR + uncuedAR + id + pairMotSame , data = data.Diff, mean)

data.Difference <- data.frame(id=1:nrow(data_agg_diff), partId = NA, uncuedAR=NA, cuedAR =NA, responseAR = NA,normed_responseAR = NA, responseToIdentical = NA, selectedCuedARofIdentical =NA, pairMotSame = NA)

for (i in 1:nrow(data_agg_same)) { # get the response to the same
  # get the id of a given trial
  #selectedID <- data_agg_same$id[i]
  selectedCuedAR <- data_agg_same$cuedAR[i]
  selectedResponseARtoSame <- data_agg_same$responseAR[i] #this is the response to the identical pair
  #look through the data and match the cuedAR with the data
  for (a in 1:nrow(data_agg_diff)){
    if (selectedCuedAR == data_agg_diff$cuedAR[a]){
      #copy the results from data_agg_diff
      data.Difference$partId[a] <- data_agg_diff$id[a] #overwriting on the id
      data.Difference$uncuedAR[a]<- data_agg_diff$uncuedAR[a] #write the uncuedAR
      data.Difference$cuedAR[a]<- data_agg_diff$cuedAR[a] #write the cuedAR
      data.Difference$responseAR[a]<-data_agg_diff$responseAR[a] 
      #write down the norm responses and cross checks
      data.Difference$normed_responseAR[a] <- data_agg_diff$responseAR[a] - selectedResponseARtoSame #positive values indicate deviation from the cuedAR towards "taller" responses
      data.Difference$responseToIdentical[a] <- selectedResponseARtoSame
      data.Difference$selectedCuedARofIdentical[a]<- selectedCuedAR #this should match with the trial's cuedAR
      data.Difference$pairMotSame[a] <- data_agg_diff$pairMotSame[a]
    }
  }
}


lm_main <- lm(normed_responseAR ~  cuedAR * uncuedAR * pairMotSame, data.Difference)
summary(lm_main) # pairMotSame interacts with the uncued's influence on the responseAR


scatter3d(x = data.Difference$uncuedAR, y = data.Difference$normed_responseAR, z = data.Difference$cuedAR,
          point.col = "blue", surface=FALSE, groups = as.factor(data.Difference$pairMotSame), grid = TRUE, fit = "linear", ellipsoid = TRUE,
          xlab = "uncuedAR", 
          ylab = "responseAR",
          zlab = "cuedAR") #groups = iris$Species # add concentration ellipses: ellipsoid = TRUE

