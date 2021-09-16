# gorkemer
# 15 August 
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
data = read.csv('bgm_cleanData_cor.csv', header=TRUE)
#data = read.csv('untouchedData.csv', header = TRUE)
#data = read.csv('bgm_cleanData_lm.csv', header = TRUE)

# define the data frame to use
data.Same <- subset(data, data$pairAR_sameCondition==1)
data.Diff <- subset(data, !(data$pairAR_sameCondition==1))

## check for normal distribution
fit <- fitdistr(data.Same$responseAR, densfun="normal")  # we assume my_data ~ Normal(?,?)
fit
hist(data.Same$responseAR, pch=20, breaks=25, prob=TRUE, main="")
curve(dnorm(x, fit$estimate[1], fit$estimate[2]), col="red", lwd=2, add=T)

# it is not normally distributed actually. 
qqnorm(data.Same$responseAR,main="QQ plot of normal data",pch=19) #it looks okay? 
data.sample <- data.Same[sample(nrow(data.Same), 5000),]
shapiro.test(data.sample$deviationFromTheCuedAR) ## shapiro confirms that it is not normally distributed

wilcox.test(responseAR ~ pairMotSame , data = data.Same)
wilcox.test(data.Same$responseAR,data.Same$cuedAR,paired=T)

#### normalizing responseAR #### substracting from the cuedAR
data.Same$deviationFromTheCuedAR <- data.Same$responseAR - data.Same$cuedAR

fit <- fitdistr(data.Same$deviationFromTheCuedAR, densfun="normal")  # we assume my_data ~ Normal(?,?)
fit
hist(data.Same$deviationFromTheCuedAR, pch=20, breaks=25, prob=TRUE, main="")
curve(dnorm(x, fit$estimate[1], fit$estimate[2]), col="red", lwd=2, add=T)

#### test the local effects ####

# effect of pairMotSame on same responses
lm <- lm(responseAR ~ cuedAR * pairMotSame , data = data.Same)
summary(lm) # there is an overall effect of pair motion direction sameness on the same trials

lm <- lm(deviationFromTheCuedAR ~ cuedAR * pairMotSame , data = data.Same)
summary(lm) # there is an overall effect of pair motion direction sameness on the same trials

lm <- lm(deviationFromTheCuedAR ~ cuedAR * pairGlobalOrg * pairMotSame , data = data.Same)
summary(lm) # there is an overall effect of pair motion direction sameness on the same trials


