# gorkemer
# 5 August 
# analyze on the cleaned data (output of the prepareData_01.R)
# extras: #colnames(group_mean) <- c("Mean Response AR", "CuedAR", "id")

## figure 2B olusturabilirsin, get the pairs when pairAR_sameCondition == 1 and
## plot the cued vs responseAR and facet_wrap in global_org

### Load libraries ###

library(reshape2); library(rlang); library(data.table); library(plyr) # data formatting
library(psych); library(pastecs) # statistics
library(lme4); library(car); library(ggpubr); library(corrplot); library(magrittr) # regression
library(ggplot2); library(Hmisc) # graphing
library(dplyr) # data wrangling
library("reshape2")

#remove scientific notation in the entire R session
options(scipen = 100)

# set wd if needed setwd('~/Desktop/background-motion/Data Analysis'), check with getwd()
setwd('~/Desktop/background-motion/Data Analysis')
data = read.csv('bgm_cleanData_cor.csv', header=TRUE)

#### test before sweeny 2011 study ####

data.identicalPairs <- subset(data, data$pairAR_sameCondition == 1)

#lm1 <- lm(responseAR~ cuedAR *  pairMotCoherent * cuedMotDir, data.identicalPairs)

data.coherent <- subset(data, data$pairMotCoherent==1)
data.inCoherent <- subset(data, !(data$pairMotCoherent==1))

# overall Data
lm1 <- lm(responseAR ~ cuedAR * uncuedAR * pairAR_sameCondition * pairMotCoherent * motDirComp, data)
summary(lm1)

# only coherent
lm1 <- lm(responseAR ~ cuedAR * uncuedAR * pairAR_sameCondition * pairMotDirSameness, data.coherent)
summary(lm1)

summary(lm1)
lm1.CI = confint(lm1,  level=0.95)
plot_model(lm1, type = "slope")


# add absolute response error (without direction)
data$abs_response_error<- abs(data$response_error)

# add absolute distance of the response to the pairAvg. 
# smaller numbers mean more attraction to the average
data$abs_distTo_pairAvg <- abs(data$responseAR - data$pairAvg)

# seperate data into two defined by the motion coherence
data.random <- subset(data, data$pairMotDirSameness_wRandom==2)
data.coherent <- subset(data, !(data$pairMotDirSameness_wRandom==2))

# test plots 
testPlots <- ggplot(data.coherent, aes(x = cuedAR, y =responseAR )) + 
  geom_point(shape=19, size=1, alpha=0.1, show.legend = FALSE, aes(color= pairGlobalOrg) ) +
  geom_smooth(data = data.coherent, method=lm, color = "red") +
  ylim(-.7,.7) +
  labs(x="Cued aspect ratio", y = "Reported aspect ratio") +
  ggtitle("Cleaned data") + facet_wrap(~pairGlobalOrg)

testPlots + geom_density_2d(color="blue", alpha=0.2) + xlab((testPlots$mapping$x))+ ylab(testPlots$mapping$y) 

data_hierchical <- aggregate(responseAR  ~ cuedAR + pairGlobalOrg +pairMotDirSameness  , data = data.coherent, mean)

# grouped boxplot
abc <- ggplot(data_hierchical, aes(x=cuedAR, y=responseAR, color)) + 
  geom_boxplot() + facet_grid(~pairMotDirSameness)

abc + geom_line(aes(linetype = as.factor(pairGlobalOrg))) + 
  theme_classic(base_family = "Arial") #+ geom_jitter()


#### ANOVA ####

# Compute the analysis of variance
res.aov <- aov(responseAR ~ pairGlobalOrg + pairMotDirSameness, data = data_hierchical)
# Summary of the analysis
summary(res.aov)

#### Visually testing the assumptions of mul lin reg ####
lm1 <- lm(response_error ~arDiff + pairGlobalOrg + coherence , data)

lm2 <- lm(abs_distTo_pairAvg ~pairGlobalOrg + pairSameARCat + coherence, data)
# Assumption 1 - Linearity of Relationship
plot(lm1,1) # this should form a horizontal line

# Assumption 2 - Independence of Variables
# testing collinearity between explanatory variables
corrplot(cor(data.coherent[,c("cuedAR","uncuedAR","responseAR")]),method='circle') # smaller circles means less collinearity
# for pair
corrplot(cor(data.coherent[,c("arDiff","pairMotDirSameness","pairSameARCat", "response_error")]),method='circle') # smaller circles means less collinearity


# Assumption 3 - Normal Distribution of Residuals
#This test is to ensure that there are no other significant relationships 
#that could be explaining the variance that have not been taken into account in the linear regression
plot(lm1,2) # straight line signal that residual terms are normally distributed

# Assumption 4 - Homoscedasticity or Equal Variance of Variables
#to determine if the error terms is the same across all values of the independent variable. 
#Here we are looking for a constant spread of the residuals
plot(lm1,3) # relatively equal variance of residuals across the variable range is good.
# Residuals should evenly distributed across the range and do not appear spread or narrow at any point

#### checking and plotting Regression Model Assumptions ####
# 1- Linearity of the data
plot(fit, 1) 
# ideally, the residual plot will show no fitted pattern. 
#That is, the red line should be approximately horizontal at zero. 
#The presence of a pattern may indicate a problem with some aspect of the linear model.

# 2- Normality of residuals
plot(fit, 2)
#The normal probability plot of residuals should approximately follow a straight line.
#In our case, all the points fall approximately along this reference line, so we can assume normality.

# 3- Homogeneity of variance
plot(fit, 3)
#This plot shows if residuals are spread equally along the ranges of predictors. 
#It’s good if you see a horizontal line with equally spread points.


# Cook's distance
plot(lm1, 4)

# Residuals vs Leverage
plot(lm1, 5)

summary(lm1)


#### Subtracting responses to its identical shape from the response ####
data_check <- data


data_check$AR_repulsionIndex <- NA
data_check_random$AR_repulsionIndex <- NA
data_check_coherent$AR_repulsionIndex <- NA

# response to cuedAR is substracted by the response to that cuedAR when it was paired with its identical AR
# i.o.w when pairAR_sameCondition is = 1
# first, aggregate all by id
# second, aggregate all by cuedAR
# third, pairAR_sameCondition
# I can do this for random and coherent motion data samples seperately
# for coherent, I can add (as a fourth value) pairMotDirSameness


group_mean <- aggregate(responseAR  ~ cuedAR + subNum + pairAR_sameCondition, data = data_check_random, mean)
group_mean.sorted <- arrange(group_mean, subNum, cuedAR, desc(pairAR_sameCondition))
# sort by mpg and cyl
group_mean.sorted$AR_repulsionIndex <- NA
group_mean.sorted$cuedARCat <- NA
group_mean.sorted$assumptionTest_1_cuedNow <- NA
group_mean.sorted$assumptionTest_2_cuedPrev <- NA


# response normalization for single subject
for (i in 1:nrow(group_mean.sorted)) { # for each trial in data,

  # get the cuedAR 
  cuedAR <- group_mean.sorted$cuedAR[i];
  
  #substract the two but conditioned on the sameness
  if (group_mean.sorted$pairAR_sameCondition[i]==1){
    responseAR_same <- group_mean.sorted$responseAR[i]
    group_mean.sorted$AR_repulsionIndex[i] <- cuedAR - responseAR_same
  }
  if (group_mean.sorted$pairAR_sameCondition[i]==0){
    #substract response to each ellipse from the response when it was paired with its identical
    group_mean.sorted$AR_repulsionIndex[i] <- group_mean.sorted$responseAR[i] - responseAR_same 
    if (i != 1){
    group_mean.sorted$assumptionTest_1_cuedNow[i]  <- group_mean.sorted$cuedAR[i];
    group_mean.sorted$assumptionTest_2_cuedPrev[i] <- group_mean.sorted$cuedAR[i-1];
    }
  }
  
  # find the AR category of the cued shape
  if (group_mean.sorted$cuedAR[i] < 0) {
    group_mean.sorted$cuedARCat[i] <- "flat_ellipse"
  }
  if (group_mean.sorted$cuedAR[i] > 0) {
    group_mean.sorted$cuedARCat[i] <- "tall_ellipse"
  }
  if (group_mean.sorted$cuedAR[i] == 0) {
    group_mean.sorted$cuedARCat[i] <- "circle"
  }
}

group_mean.sorted$assumtionTest <- ifelse(group_mean.sorted$assumptionTest_1_cuedNow == group_mean.sorted$assumptionTest_2_cuedPrev, "identical", "not-identical")


# subset by same vs diff AR pairs
group_mean.sorted_onlySame <- subset(group_mean.sorted, group_mean.sorted$pairAR_sameCondition == 1)
group_mean.sorted_diff <- subset(group_mean.sorted, group_mean.sorted$pairAR_sameCondition == 0)

unique(group_mean.sorted_diff$assumtionTest)

#overall mean for the different trials: positive values indicate repulsion
m_diff <- mean(as.numeric(group_mean.sorted_diff$AR_repulsionIndex))
#overall mean for the same trials
mean(as.numeric(group_mean.sorted_onlySame$AR_repulsionIndex))

t.test(group_mean.sorted_diff$AR_repulsionIndex, mu = 0)
### overall repulsion effect is not significantly different than 0 

#get seperate mean for different ellipse shapes
# flat ellipse
group_mean.sorted_diff_flat <- subset(group_mean.sorted_diff, group_mean.sorted_diff$cuedARCat=="flat_ellipse")
mean(as.numeric(group_mean.sorted_diff_flat$AR_repulsionIndex))
# circle
group_mean.sorted_diff_circle <- subset(group_mean.sorted_diff, group_mean.sorted_diff$cuedARCat=="circle")
mean(as.numeric(group_mean.sorted_diff_circle$AR_repulsionIndex))
# tall ellipse
group_mean.sorted_diff_tall <- subset(group_mean.sorted_diff, group_mean.sorted_diff$cuedARCat=="tall_ellipse")
mean(as.numeric(group_mean.sorted_diff_tall$AR_repulsionIndex))


#### on overall data with collapsing across every trial ####
#data_collapsed <- data
#data_collapsed <- data.random
data_collapsed <- data.coherent

data_collapsed$AR_repulsionIndex <- NA

group_mean <- aggregate(responseAR  ~ cuedAR + pairAR_sameCondition + pairMotDirSameness, data = data_collapsed, mean)
group_mean.sorted <- arrange(group_mean, cuedAR, desc(pairMotDirSameness), desc(pairAR_sameCondition))

# sort by mpg and cyl
group_mean.sorted$AR_repulsionIndex <- NA

# response normalization for single subject
for (i in 1:nrow(group_mean.sorted)) { # for each trial in data,
  
  # get the cuedAR 
  cuedAR <- group_mean.sorted$cuedAR[i];
  
  #substract the two but conditioned on the sameness
  if (group_mean.sorted$pairAR_sameCondition[i]==1){
    responseAR_same <- group_mean.sorted$responseAR[i]
    group_mean.sorted$AR_repulsionIndex[i] <- cuedAR - responseAR_same
  }
  if (group_mean.sorted$pairAR_sameCondition[i]==0){
    #substract response to each ellipse from the response when it was paired with its identical
    group_mean.sorted$AR_repulsionIndex[i] <- group_mean.sorted$responseAR[i] - responseAR_same 
    if (i != 1){
      group_mean.sorted$assumptionTest_1_cuedNow[i]  <- group_mean.sorted$cuedAR[i];
      group_mean.sorted$assumptionTest_2_cuedPrev[i] <- group_mean.sorted$cuedAR[i-1];
    }
  }
  
  # find the AR category of the cued shape
  if (group_mean.sorted$cuedAR[i] < 0) {
    group_mean.sorted$cuedARCat[i] <- "flat_ellipse"
  }
  if (group_mean.sorted$cuedAR[i] > 0) {
    group_mean.sorted$cuedARCat[i] <- "tall_ellipse"
  }
  if (group_mean.sorted$cuedAR[i] == 0) {
    group_mean.sorted$cuedARCat[i] <- "circle"
  }
}

# subset by same vs diff AR pairs
group_mean.sorted_onlySame <- subset(group_mean.sorted, group_mean.sorted$pairAR_sameCondition == 1)
group_mean.sorted_diff <- subset(group_mean.sorted, group_mean.sorted$pairAR_sameCondition == 0)

lm <- lm(AR_repulsionIndex ~ cuedAR + pairMotDirSameness + pairAR_sameCondition, group_mean.sorted)
summary(lm)


#subset only trials when motion direction is same or diff 

pairMoveDiffDir <- subset(group_mean.sorted, group_mean.sorted$pairMotDirSameness == 0)
pairMoveSameDir <- subset(group_mean.sorted, group_mean.sorted$pairMotDirSameness == 1)

pairMoveDiffDir_pairsDiff <- subset(pairMoveDiffDir, pairMoveDiffDir$pairAR_sameCondition==0)
pairMoveSameDir_pairsDiff <- subset(pairMoveSameDir, pairMoveSameDir$pairAR_sameCondition==0)

mean(as.numeric(pairMoveDiffDir$AR_repulsionIndex))
mean(as.numeric(pairMoveSameDir$AR_repulsionIndex))

# grouped boxplot
abc <- ggplot(group_mean.sorted_diff, aes(x=cuedAR, y=AR_repulsionIndex)) + 
  geom_boxplot() + facet_grid(~pairMotDirSameness)

abc + geom_jitter()

# Calculate proportion of each level
proportion <- table(group_mean.sorted_diff$cuedAR)/nrow(group_mean.sorted_diff)
boxplot(group_mean.sorted_diff$AR_repulsionIndex ~ group_mean.sorted_diff$cuedAR , width=proportion , col=c("orange" , "seagreen"))

t.test(pairMoveDiffDir_pairsDiff$AR_repulsionIndex,pairMoveSameDir_pairsDiff$AR_repulsionIndex)

#overall mean for the different trials: positive values indicate repulsion
m_diff <- mean(as.numeric(group_mean.sorted_diff$AR_repulsionIndex))
m_diff

#overall mean for the same trials
mean(as.numeric(group_mean.sorted_onlySame$AR_repulsionIndex))

t.test(group_mean.sorted_diff$AR_repulsionIndex, mu = 0)

# test plots 
testPlots <- ggplot(group_mean.sorted_diff, aes(x = cuedAR, y =AR_repulsionIndex )) + 
  geom_point(shape=19, size=1, alpha=0.1, show.legend = FALSE, aes(color= pairMotDirSameness) ) + 
  geom_smooth(data = group_mean.sorted_diff, method=lm, color = "red") +
  ylim(-.7,.7) +
  labs(x="Cued aspect ratio", y = "Reported aspect ratio") +
  ggtitle("Cleaned data") + facet_wrap(as.factor(group_mean.sorted_diff$cuedARCat))

testPlots + geom_density_2d(color="blue", alpha=0.2) + xlab((testPlots$mapping$x))+ ylab(testPlots$mapping$y)





#### with global organization added ####


# adding global organization 
group_mean_wGloOrg <- aggregate(responseAR  ~ cuedAR + subNum + pairAR_sameCondition + pairGlobalOrg, data = data_check, mean)
group_mean_wGloOrg.sorted <- arrange(group_mean_wGloOrg, subNum, cuedAR, desc(pairGlobalOrg), desc(pairAR_sameCondition))


# response normalization for single subject
for (i in 1:nrow(group_mean_wGloOrg.sorted)) { # for each trial in data,
  
  # get the cuedAR 
  cuedAR <- group_mean_wGloOrg.sorted$cuedAR[i];
  
  #substract the two but conditioned on the sameness
  if (group_mean_wGloOrg.sorted$pairAR_sameCondition[i]==1){
    responseAR_same <- group_mean_wGloOrg.sorted$responseAR[i]
    group_mean_wGloOrg.sorted$AR_repulsionIndex[i] <- cuedAR - responseAR_same
  }
  if (group_mean_wGloOrg.sorted$pairAR_sameCondition[i]==0){
    responseAR_same <- group_mean_wGloOrg.sorted$responseAR[i]
    #substract response to each ellipse from the response when it was paired with its identical
    group_mean_wGloOrg.sorted$AR_repulsionIndex[i] <- responseAR_same - group_mean_wGloOrg.sorted$responseAR
  }
  
  # find the AR category of the cued shape
  if (group_mean_wGloOrg.sorted$cuedAR[i] < 0) {
    group_mean_wGloOrg.sorted$cuedARCat[i] <- "flat_ellipse"
  }
  if (group_mean_wGloOrg.sorted$cuedAR[i] > 0) {
    group_mean_wGloOrg.sorted$cuedARCat[i] <- "tall_ellipse"
  }
  if (group_mean_wGloOrg.sorted$cuedAR[i] == 0) {
    group_mean_wGloOrg.sorted$cuedARCat[i] <- "circle"
  }
}








agg <- aggregate(list(responseAR_agg=data_check$responseAR), by = list(cued_AR_agg=data_check$cuedAR, global_org_agg = data_check$global_org, motDir_Same = data_check$motDir_Same, id=data_check$id), mean)
agg_sorted <- arrange(agg, id, desc(cued_AR_agg), desc(motDir_Same))

library(dplyr)
data_check %>%
  group_by(cuedAR, id) %>% 
  summarise_each(funs(mean))

data_check <- arrange(data_check, id, desc(cuedAR), motDir_Same, desc(sameCondition))
# new dataframe with selected columns from raw dataset
data_check <- d[ c(1, 8, 17, 18, 19, 20, 21, 22, 23, 24, 25, 31, 32, 38)]
# change column names
nms <- c("rt", "shape_org", "e1_motion_dir", "e2_motion_dir", "responseAR","cuedAR", "uncuedAR", "coherence", "cued_motion_dir", "uncued_motion_dir", "response_error", "trial_number", "round_number", "subNum")
setnames(data_check, nms)
# subset random and the other motion condition into two seperate matrices
# calculate response error relative to the S.S condition
data_check <- na.omit(data_check)

data_check <- data_check %>%
  group_by(cuedAR, subNum) %>% 
  summarise_each(funs(mean))

# subset the random motion and motion cues #
data.randomOnly <- subset(data.Clean, data.Clean$motDir_Same_randomIncluded==2)
data.motionDirection <- subset(data.Clean, !(data.Clean$motDir_Same_randomIncluded==2))

#### Analysis 1: Coherent Motion Conditon ####
data_aggregated_coherent <- aggregate(data.motionDirection$responseAR, by=list(data.motionDirection$cuedAR,data.motionDirection$sameCondition, data.motionDirection$id, data.motionDirection$arDiff, data.motionDirection$motDir_Same_randomIncluded),mean)
nms <- c("cuedAR", "sameCondition", "id","arDiff","mot_dir","responseAR")
setnames(data_aggregated_coherent, nms)

data_sorted <- arrange(data_aggregated_coherent, id, desc(cuedAR), mot_dir, desc(sameCondition))

data_sorted$responseAR_rSS <- NA

# response normalization for single subject
for (i in 1:nrow(data_sorted)) { # for each trial in data,
  #print(index.Single[i,])
  
  if (data_sorted$sameCondition[i]==1){
    responseToSS = data_sorted$responseAR[i]
    data_sorted$responseAR_rSS[i] <- data_sorted$responseAR[i]
    #flag = 1
  }
  else {
    data_sorted$responseAR_rSS[i] <- data_sorted$responseAR[i] - responseToSS
    #flag == 0
  }
}

cleanScatter <- ggplot(data_sorted, aes(x = cuedAR, y = responseAR_rSS)) + 
  geom_point(shape=19, size=3, alpha=0.1, show.legend = FALSE) +
  geom_smooth(data = data_sorted, method=lm, color = "red") +
  ylim(-.7,.7) +
  labs(x="Cued aspect ratio", y = "Reported aspect ratio") +
  ggtitle("Coherent Motion Condition: Shared-motion rendered the line a little steeper didn't work")
cleanScatter + geom_density_2d(color="blue", alpha=0.2) + facet_wrap(~mot_dir)

#### Analysis 2: Random Motion Direction #### 

data_aggregated_random <- aggregate(data.randomOnly$responseAR, by=list(data.randomOnly$cuedAR,data.randomOnly$sameCondition, data.randomOnly$id, data.randomOnly$arDiff, data.randomOnly$global_org),mean)
nms <- c("cuedAR", "sameCondition", "id","arDiff","global_org","responseAR")
setnames(data_aggregated_random, nms)

data_sorted_random <- arrange(data_aggregated_random, id, desc(cuedAR), global_org,desc(sameCondition))
data_sorted_random$responseAR_rSS <- NA

# response normalization for single subject
for (i in 1:nrow(data_sorted_random)) { # for each trial in data,
  #print(index.Single[i,])
  
  if (data_sorted_random$sameCondition[i]==1){
    responseToSS = data_sorted_random$responseAR[i]
    data_sorted_random$responseAR_rSS[i] <- data_sorted_random$responseAR[i]
    #flag = 1
  }
  else {
    data_sorted_random$responseAR_rSS[i] <- data_sorted_random$responseAR[i] - responseToSS
    #flag == 0
  }
}

cleanScatter <- ggplot(data_sorted_random, aes(x = cuedAR, y = responseAR_rSS)) + 
  geom_point(shape=19, size=3, alpha=0.1, show.legend = FALSE) +
  geom_smooth(data = data_sorted_random, method=lm, color = "red") +
  ylim(-.7,.7) +
  labs(x="Cued aspect ratio", y = "Reported aspect ratio") +
  ggtitle("Cleaned data")
cleanScatter + geom_density_2d(color="blue", alpha=0.2)

#### Secondary Analysis: Response Error way ####

data_sorted$response_error <- data_sorted$responseAR - data_sorted$cuedAR
data_sorted$response_error_rSS <- data_sorted$responseAR_rSS - data_sorted$cuedAR
data_sorted_notEqual <- subset(data_sorted, data_sorted$sameCondition==0)

cleanScatter <- ggplot(data_sorted, aes(x = arDiff, y = response_error_rSS)) + 
  geom_point(shape=19, size=3, alpha=0.1, show.legend = FALSE) +
  geom_smooth(data = data_sorted_notEqual, method=lm, color = "red") +
  ylim(-.7,.7) +
  labs(x="arDiff", y = "response_error_rSS") +
  ggtitle("Motion Direction did not seem to influence the overall perceptual averaging phonemena. ")
cleanScatter + geom_density_2d(color="blue", alpha=0.2) + facet_wrap(~mot_dir)

# for random motion

data_sorted_random$response_error <- data_sorted_random$responseAR - data_sorted_random$cuedAR
data_sorted_random$response_error_rSS <- data_sorted_random$responseAR_rSS - data_sorted_random$cuedAR
data_sorted_random_notEqual <- subset(data_sorted_random, data_sorted_random$sameCondition==0)

cleanScatter <- ggplot(data_sorted_random, aes(x = arDiff, y = responseAR_rSS)) + 
  geom_point(shape=19, size=3, alpha=0.1, show.legend = FALSE) +
  geom_smooth(data = data_sorted_random_notEqual, method=lm, color = "red") +
  ylim(-.7,.7) +
  labs(x="arDiff", y = "response_error_rSS") +
  ggtitle("Motion Direction did not seem to influence the overall perceptual averaging phonemena. ")
cleanScatter + geom_density_2d(color="blue", alpha=0.2) #+ facet_wrap(~mot_dir)


#### basket oncesi, testi yapip yapmadiklari cok kesin degil cok iyi bir iliski gozukmuyor
### uncued diff iliskisine bakilabilir regression'da oda response error hesaplayip

# prepare for two different set Random & Coherent
data_random <- subset(data.Clean_lm, data.Clean_lm$motDir_Same_randomIncluded==2)
data_coherent <- subset(data.Clean_lm, (!(data.Clean_lm$motDir_Same_randomIncluded==2)))

# Run the basic model
fit = lm(responseAR ~ 1 + cuedAR*uncuedAR*global_org, data = data_coherent)
fit.stats = summary(fit)
fit.CI = confint(fit,  level=0.95)
fit.stats
#fit.stats

# Run the alternative model:  arDiff*global_org*motDir_Same
fit = lm(response_error ~ 1 + arDiff*global_org*motDir_Same, data = data_coherent)
fit.stats = summary(fit)
fit.CI = confint(fit,  level=0.95)
fit.stats
#fit.stats

#### checking and plotting Regression Model Assumptions ####
# 1- Linearity of the data
plot(fit, 1) 
# ideally, the residual plot will show no fitted pattern. 
#That is, the red line should be approximately horizontal at zero. 
#The presence of a pattern may indicate a problem with some aspect of the linear model.

# 2- Normality of residuals
plot(fit, 2)
#The normal probability plot of residuals should approximately follow a straight line.
#In our case, all the points fall approximately along this reference line, so we can assume normality.

# 3- Homogeneity of variance
plot(fit, 3)
#This plot shows if residuals are spread equally along the ranges of predictors. 
#It’s good if you see a horizontal line with equally spread points.

# plotting lm results
cleanScatter <- ggplot(data_coherent, aes(x = arDiff, y = response_error)) + 
  geom_point(shape=19, size=3, alpha=0.1, show.legend = FALSE) +
  geom_smooth(data = data_coherent, method=lm, color = "red") +
  ylim(-.7,.7) +
  labs(x="arDiff", y = "response_error") +
  ggtitle(" ")
cleanScatter + geom_density_2d(color="blue", alpha=0.2) + facet_wrap(~motDir_Same)

