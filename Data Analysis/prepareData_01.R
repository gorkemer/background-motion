# Updated Analysis script for Shape Integration Experiment: Background-motion 
# 5 August 2021
# included new names following camel style in code
# ###seperate analysis for random motion and motionGrouping 

### Load libraries ###
# data formatting
library(reshape2); library(rlang); library(data.table); library(plyr)
# statistics
library(psych); library(pastecs)
# regression
library(lme4); library(car); library(ggpubr)
# graphing
library(ggplot2); library(Hmisc)
# data wrangling
library(dplyr)
library("reshape2")

# set wd if needed setwd('~/Desktop/background-motion/data'), check with getwd()
setwd('~/Desktop/background-motion/data')
d = read.csv('background-motion22June.csv', header=TRUE)
#17952, 19101, 19198, 19210 are empty participants. The reall participant number is:
print("participant number is:")
length(unique(d$participant_ID))-4
#### Organize dataframe ####
# get only the meaningful columns

data <- d[c("rt","shapeOrganizationNumber", "cueType","ellipse1_move_direction",
               "ellipse2_move_direction", "selected_ellipse_logAR", "cued_ellipse_logAR" ,"uncued_ellipse_logAR",
               "coherence_level", "cuedMotionDirection","unCuedMotionDirection", "differenceBetweenCuedAndReported",
               "trial_num","round_num","participant_ID")]

# change column names
nms <- c("rt", "shape_org", "cueType", "e1_motion_dir", "e2_motion_dir", "responseAR","cuedAR", "uncuedAR", "coherence", "cued_motion_dir", "uncued_motion_dir", "response_error", "trial_number", "round_number", "subNum")
setnames(data, nms)

# reorder columns for ease of eyeballing the data 
data <- data[,c("subNum","trial_number", "cuedAR", "uncuedAR", "responseAR","response_error","rt","shape_org", "cueType", "e1_motion_dir", "e2_motion_dir", "round_number", "coherence", "cued_motion_dir", "uncued_motion_dir")]

# remove na in r - remove rows - na.omit function / option
data <- na.omit(data)

#making characters numeric
data <- mutate_all(data, function(x) as.numeric(as.character(x)))

# what is the difference between uncued and cued?
data$arDiff <- NA
data$arDiff <- data$uncuedAR - data$cuedAR

# what is the pair average? Pair Average can be used as an index for perceptual averaging (DV)
#because it can tell how much the response were closer to the pair average
data$pairAvg <- NA
data$pairAvg <- (data$uncuedAR + data$cuedAR)/2

# what is the motion direction of the cued shape? Vertical:1 vs Horizontal:0 vs Random
data$cuedMotDir <- ifelse(data$cued_motion_dir == 90 | data$cued_motion_dir == 270, "1", ifelse(data$cued_motion_dir == 0 | data$cued_motion_dir == 180, "0", "2"))

# is the pair organized vertically (within) or horizontally (across)? [shape_org]: 1 & 2 (horiz) or 3 & 4 (verti) 
# 13.August.2021: double checked: 1 horizontal-top, 2: horizontal-bottom; 3: vertical-left, 4: vertical-right
data$pairGlobalOrg <- ifelse(data$shape_org == 3 | data$shape_org == 4, 1, 0)

# what is the shared-motion-ness? CAUTION: This takes the random motion as Same. You need to seperate out the random motion. 
data$pairMotDirSameness <- ifelse(data$cued_motion_dir == data$uncued_motion_dir, 1, 0)

# take into account the random direction (2), which is different than the same motion
data$pairMotDirSameness_wRandom <- ifelse(((as.numeric(data$cued_motion_dir) ==9999)),2,  # "RANDOM DIRECTION"
                                  ifelse(((as.numeric(data$cued_motion_dir) == as.numeric(data$uncued_motion_dir))),1, # "SAME"
                                  0) # "DIFFERENT")
)

# take into account the random direction (2), which is different than the same motion
data$pairMotSame <- ifelse(((as.numeric(data$cued_motion_dir) ==9999)),0,  # "RANDOM DIRECTION"
                                          ifelse(((as.numeric(data$cued_motion_dir) == as.numeric(data$uncued_motion_dir))),1, # "SAME"
                                                 0) # "DIFFERENT")
)

data$pairMotDiff <- ifelse(((as.numeric(data$cued_motion_dir) ==9999)),0,  # "RANDOM DIRECTION"
                           ifelse(((as.numeric(data$cued_motion_dir) == as.numeric(data$uncued_motion_dir))),0, # "SAME"
                                  1) # "DIFFERENT")
)

data$pairMotRandom <- ifelse(((as.numeric(data$cued_motion_dir) ==9999)),1,  # "RANDOM DIRECTION"
                           ifelse(((as.numeric(data$cued_motion_dir) == as.numeric(data$uncued_motion_dir))),0, # "SAME"
                                  0) # "DIFFERENT")
)

# Are the two shapes in the same category of aspect ratio? i.e., both tall?
data$pairSameARCat <- ifelse(data$cuedAR < -0.00000000  & data$uncuedAR < -0.00000000 | data$cuedAR > -0.00000000 & data$uncuedAR > -0.00000000 | data$cuedAR == data$uncuedAR, "1", "0")

# finding the trials where pairs have the same AR
data$pairAR_sameCondition <- ifelse(data$cuedAR == data$uncuedAR, 1, 0) 

#create pair motion direction composition
data$motDirComp <- 
  ifelse(data$cued_motion_dir == 90 & data$uncued_motion_dir == 90, "1", 
         ifelse(data$cued_motion_dir == 270 & data$uncued_motion_dir == 270, "1", 
                ifelse(data$cued_motion_dir == 180 & data$uncued_motion_dir == 180, "0", 
                       ifelse(data$cued_motion_dir == 0 & data$uncued_motion_dir == 0, "0", 
                              ifelse(data$cued_motion_dir == 9999 & data$uncued_motion_dir == 9999, "2", "3")))))

#define trials where motion was coherent (1) or random (0)
data$pairMotCoherent <- ifelse(!(data$cuedMotDir==2), 1, 0) 

#define trials where the motion is diagonal (1) or not (transverse; 0)
data$pairMotDiago <- ifelse(!(data$motDirComp==3), 1, 0) 

# sort the data by observer. Necessary because the next step assigns a new id when subNum changes.
data <- data[order(data$subNum),]


### Add subject IDs starting at "1" ###
# e.g., subIDs will be o1, o2, etc.

# remove na in r - remove rows - na.omit function / option
data <- na.omit(data)

sub = data$subNum[1]
number = 1

for (i in 1:nrow(data)) {
  if (data$subNum[i] != sub) {
    number = number+1
    sub = data$subNum[i]
  }
  data$id[i] = paste("o",number, sep = "")
}

singCorrs= array(data = NA,dim=c(length(unique(data$id)),3)) #empty 74x3 array
nS = length(unique(data$id)) #n iterations for the loop below. Should be 74 unique ids

for (s in 1:nS) {
  id = s
  name = paste("o",id,sep = "") 
  sing.Data <- subset(data, data$id==name) #pull up data from one observer. 
  res <- cor.test(sing.Data$cuedAR, sing.Data$responseAR, method = c("pearson")) #run correl on that person's data
  singCorrs[s,1] = name #add correlation results to array
  singCorrs[s,2] = res$estimate
  singCorrs[s,3] = res$p.value
}

# plot correlation coefficients and p-values
plot(singCorrs[,2], singCorrs[,3], main="Correlation values: Individuals",
     xlab="Corr Coeff", ylab="P-value", pch=19)

### Add individual correlation results to main data ###
for (i in 1:nrow(data)) {
  #find the row in correlation array with current number is subject
  id = data$id[i]
  r = which(singCorrs == id)
  data$corrCoeff[i] <- singCorrs[r,2]
  data$corrPval[i] <- singCorrs[r,3]
}
subset(data, data$id=="o63") #check one person's data against correlation array

#### Clean Data Creation ####
# New data set containing only Os with significant correlations
data.Clean <- subset(data, data$corrCoeff>0.2)

data.Clean_p05 <- subset(data, data$corrPval<0.05)

cleanScatter <- ggplot(data.Clean, aes(x = cuedAR, y = responseAR)) + 
  geom_point(shape=19, size=1, alpha=0.1, show.legend = FALSE, aes(color= pairAR_sameCondition) ) +
  geom_smooth(data = data.Clean, method=lm, color = "red") +
  ylim(-.7,.7) +
  labs(x="Cued aspect ratio", y = "Reported aspect ratio") +
  ggtitle("Cleaned data")

cleanScatter + geom_density_2d(color="blue", alpha=0.2)
#add x and y titles dynamically
cleanScatter + xlab((cleanScatter$mapping$x)) + ylab(cleanScatter$mapping$y)

#### Data Clean with LM ####

# find the conditions where pairs had the identical AR 
singSlopes= array(data = NA,dim=c(length(unique(data$id)),3)) #empty 74x3 array
nS = length(unique(data$id)) #n iterations for the loop below. Should be 74 unique ids
data_onlyEqualPairs <- subset(data, data$pairAR_sameCondition==1)

#then look for regression btw response AR (DV) and cued AR (IV) for each id
for (s in 1:nS) {
  id = s
  name = paste("o",id,sep = "") 
  sing.Data_slopes <- subset(data_onlyEqualPairs, data_onlyEqualPairs$id==name) #pull up data from one observer. 
  res <- lm(sing.Data_slopes$responseAR ~ sing.Data_slopes$cuedAR, data = sing.Data_slopes) #run correl on that person's data
  singSlopes[s,1] = name #add correlation results to array
  singSlopes[s,2] = res$coefficients[2]              #adding the beta coefficient (the slope) of the regression
  singSlopes[s,3] = summary(res)$coefficients[,4][2] #adding the p-values of the regression coefficients
}

# plot slope coefficients and p-values
plot(singSlopes[,2], singSlopes[,3], main="Slope Coefficient values: Individuals",
     xlab="Slope Coeff", ylab="P-value", pch=19)

# Add individual slopes results to main data
for (i in 1:nrow(data)) {
  #find the row in correlation array with current number is subject
  id = data$id[i]
  r = which(singSlopes == id)
  data$slopeCoeff[i] <- singSlopes[r,2]
  data$slopePval[i] <- singSlopes[r,3]
}
subset(data, data$id=="o63") #check one person's data against correlation array

# New data set containing only Participants with significant reg slop btw response & cuedAR
data.Clean_lm <- subset(data, data$slopeCoeff>0.2)

cleanScatter <- ggplot(data.Clean_lm, aes(x = cuedAR, y = responseAR)) + 
  geom_point(shape=19, size=1, alpha=0.1, show.legend = FALSE, aes(color= pairMotDirSameness) ) +
  geom_smooth(data = data.Clean_lm, method=lm, color = "red") +
  ylim(-.7,.7) +
  labs(x="Cued aspect ratio", y = "Reported aspect ratio") +
  ggtitle("Cleaned data for regression slope method")
cleanScatter + geom_density_2d(color="blue", alpha=0.2) + xlab((cleanScatter$mapping$x)) + ylab(cleanScatter$mapping$y)

### end of the organization of the data script ### 

# output the csv externally 

write.csv(data.Clean_p05, "data.Clean_p05.csv", row.names = FALSE)
write.csv(data, "untouchedData.csv", row.names = FALSE)
write.csv(data.Clean,"bgm_cleanData_cor.csv", row.names = FALSE)
write.csv(data.Clean_lm,"bgm_cleanData_lm.csv", row.names = FALSE)
