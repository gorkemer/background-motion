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
library(fitdistrplus); library(MASS); library(survival); library(npsurv) ; library(lsei)

library(extrafont)
font_import("Trebuchet MS")
library(actuar)
library(statmod)

#remove scientific notation in the entire R session
#options(scipen = 100)

# set wd if needed setwd('~/Desktop/background-motion/Data Analysis'), check with getwd()
setwd('~/Desktop/background-motion/Data Analysis')
data = read.csv('bgm_cleanData_cor.csv', header=TRUE)
#data = read.csv('untouchedData.csv', header = TRUE)
#data = read.csv('bgm_cleanData_lm.csv', header = TRUE)

# define the data frame to use
data_o1 <- data
#data_o1 <- subset(data, data$id == "o102") #| data$id == "o2" | data$id == "o4")
#data_o1 <- subset(data_o1, select = c("id","cuedAR", "pairAR_sameCondition", "responseAR", "uncuedAR", "pairMotCoherent", "pairGlobalOrg", "pairMotDirSameness") )
#data_o1 <- arrange(data_o1, id, cuedAR, desc(pairAR_sameCondition), desc(responseAR))
#data_o1 <- subset(data, data$pairMotCoherent==1) #coherent
#get the trial number of each participant data
for (i in 1:nrow(data_agg_same)) { 

}

#find the identical pairs and different pairs
data.identicals <- subset(data_o1, data_o1$pairAR_sameCondition==1)
data.identicals <- arrange(data.identicals, id, cuedAR, desc(pairAR_sameCondition), desc(responseAR))
data.differents <- subset(data_o1, data_o1$pairAR_sameCondition==0) 

#get the aggregated response to cuedAR *add for running for the whole set
data_agg_same <- aggregate(responseAR ~ cuedAR + id, data = data.identicals, mean)
data_agg_diff <- aggregate(responseAR ~ cuedAR + id, data = data.differents, mean)
data_agg_same_all <- aggregate(responseAR~cuedAR, data_agg_same,mean)

write.csv(data.identicals,"data_identicals.csv", row.names = FALSE)



#create normBySame_responseAR
data.differents$normed_responseAR <- NA
data.differents$responseToIdentical <- NA
data.differents$selectedCuedARofIdentical <- NA

# for (i in 1:nrow(data_agg_same_all)) { # get the response to the same
#   
#   # get the id of a given trial
#   #selectedID <- data_agg_same$id[i]
#   selectedCuedAR <- data_agg_same_all$cuedAR[i]
#   selectedResponseARtoSame <- data_agg_same_all$responseAR[i] #this is the response to the identical pair
# 
#   #look through the data and match the cuedAR with the data
#   for (a in 1:nrow(data.differents)){
#     if (selectedCuedAR == data.differents$cuedAR[a]){ #&& selectedID == data.differents$id[a]){ #find the cuedAR and id couplings and substract response to identical from that response on that trial
#       data.differents$normed_responseAR[a] <- data.differents$responseAR[a] - selectedResponseARtoSame #positive values indicate deviation from the cuedAR towards "taller" responses
#       data.differents$responseToIdentical[a] <- selectedResponseARtoSame
#       data.differents$selectedCuedARofIdentical[a]<- selectedCuedAR #this should match with the trial's cuedAR
#     }
#   }
# }
#### run this for individual-response-normalization approach's results
for (i in 1:nrow(data_agg_same)) { # get the response to the same
  
  # get the id of a given trial
  selectedID <- data_agg_same$id[i]
  selectedCuedAR <- data_agg_same$cuedAR[i]
  selectedResponseARtoSame <- data_agg_same$responseAR[i] #this is the response to the identical pair
  
  #look through the data and match the cuedAR with the data
  for (a in 1:nrow(data.differents)){
    if (selectedCuedAR == data.differents$cuedAR[a] && selectedID == data.differents$id[a]){ #find the cuedAR and id couplings and substract response to identical from that response on that trial
      data.differents$normed_responseAR[a] <- data.differents$responseAR[a] - selectedResponseARtoSame #positive values indicate deviation from the cuedAR towards "taller" responses
      data.differents$responseToIdentical[a] <- selectedResponseARtoSame
      data.differents$selectedCuedARofIdentical[a]<- selectedCuedAR #this should match with the trial's cuedAR
    }
  }
}

#subset to run cross checks
data.differents_Xcheck <- subset(data.differents,select = c("id","cuedAR", "pairAR_sameCondition", "uncuedAR", "responseAR","normed_responseAR","responseToIdentical", "selectedCuedARofIdentical") )
data.differents_Xcheck <- arrange(data.differents_Xcheck, id, cuedAR, desc(pairAR_sameCondition), desc(responseAR))
data.differents_Xcheck$checkOfCuedARs <- ifelse(data.differents_Xcheck$cuedAR==data.differents_Xcheck$selectedCuedARofIdentical,1,0)
unique(data.differents_Xcheck$checkOfCuedARs)
which(is.na(data.differents_Xcheck$checkOfCuedARs) == TRUE)

data.differents_randomOnly <- subset(data.differents, data.differents$pairMotCoherent==0)
data.differents_coherentOnly <- subset(data.differents, data.differents$pairMotCoherent==1)

#get only the coherent and Same
data.differents_coherentOnly_same <- subset(data.differents_coherentOnly, data.differents_coherentOnly$pairMotDirSameness==1)
data.differents_coherentOnly_diff <- subset(data.differents_coherentOnly, data.differents_coherentOnly$pairMotDirSameness==0)

#see the data to test
#data.differents_o2 <- subset(data.differents, data.differents$id=="o2")
testPlots <- ggplot(data.differents_coherentOnly, aes(x = uncuedAR, y = normed_responseAR)) + 
  geom_point(shape=19, size=1, alpha=0.1, show.legend = FALSE) +
  geom_smooth(data = data.differents_coherentOnly, method=lm, color = "red") +
  ylim(-.7,.7) +
  labs(x="Cued aspect ratio", y = "Reported aspect ratio") +
  ggtitle("Cleaned data") #+ facet_wrap(~id,nrow = 5)

testPlots + geom_density_2d(color="blue", alpha=0.2) + xlab((testPlots$mapping$x)) + theme(text=element_text(family="Trebuchet MS"))

Palette1 <- c('red','green','blue','violet','black')
#    
scatterPlot <- ggplot(data.differents, aes(x = uncuedAR, y = normed_responseAR)) + 
  geom_point(shape=19, size=3, alpha=0.03, show.legend = FALSE) +
  geom_point(aes(shape = factor(pairMotSame))) +
  geom_point(aes(color = factor(pairMotSame))) +
  #geom_smooth(data = data, method=lm, color = "red") +
  geom_smooth(method = "lm", 
              aes(color = factor(pairMotSame))) +
  labs(x="(<-flatter)   Uncued - Cued   (taller->)", y = "(<-flatter)   Response Error   (taller->)") +
  labs(title="The diff between uncued and cued AR (and motion-direction-sharedness) matters", subtitle="but is this response bias, underestimation, or a real effect?")
scatterPlot + geom_density_2d(color="blue", alpha=0.3)# + facet_wrap(~id,nrow = 4) #+ scale_colour_manual(values = Palette1)
# As the uncued shape gets taller, response errors get "taller." But is is from response bias or a real effect? 


# OVERALL ANALYSIS, pairMotDirSameness
lm_main <- lm(normed_responseAR ~  cuedAR * uncuedAR * pairGlobalOrg , data.differents)
summary(lm_main)
lm_main.stats = summary(lm_main)
lm_main.CI = confint(fit,  level=0.95)
lm_main.stats
Anova(lm_main)
plot(lm_main, las = 1)
summ(lm_main)
#### DUNNO ####
lm_main <- lm(normed_responseAR ~  cuedAR * uncuedAR * pairGlobalOrg *pairMotSame, data.differents)
summary(lm_main) #Anova(lm_main)

# SECONDARY ANALYSIS: check the result -> when pair motion moves at the same direction, increases in uncuedAR increases
# makes the responseAR increase, meaning results attract towards the uncuedAR
lm_second <- lm(normed_responseAR ~ uncuedAR * pairMotDirSameness * pairGlobalOrg, data.differents_coherentOnly)
summary(lm_second)

# third analysis checking the main analysis on the coherent data only. 
lm_third <- lm(normed_responseAR ~ uncuedAR *pairGlobalOrg, data.differents_coherentOnly)
summary(lm_third)

my_model <- lmer(normed_responseAR ~ cuedAR + uncuedAR + pairGlobalOrg : uncuedAR + (1|id), data = data.differents)
Anova(my_model) 

# forth analysis: random only
lm_forth <- lm(normed_responseAR ~ cuedAR * uncuedAR * pairGlobalOrg, data.differents_randomOnly)
summary(lm_forth)


#### checking and plotting Regression Model Assumptions ####
# 1- Linearity of the data
plot(lm_main, 1) 
# ideally, the residual plot will show no fitted pattern. 
#That is, the red line should be approximately horizontal at zero. 
#The presence of a pattern may indicate a problem with some aspect of the linear model.

# 2- Normality of residuals
plot(lm_main, 2)
#The normal probability plot of residuals should approximately follow a straight line.
#In our case, all the points fall approximately along this reference line, so we can assume normality.

# 3- Homogeneity of variance
plot(lm_main, 3)
#This plot shows if residuals are spread equally along the ranges of predictors. 
#It’s good if you see a horizontal line with equally spread points.


#' The Shapiro-Wilk test directly evaluates normality. 
#' It's in the output of stat.desc (normtest.W; test stat (W), p-val)
#' but you can also just run the test on any variable.
#' If it's significant, your data are not normally distributed.
shapiro.test(data.sample$responseAR) # clearly, our data are not normally distributed


# response error analysis
data.differents$responseError <- data.differents$responseAR - data.differents$cuedAR
lm_RE <- lm(responseError ~ uncuedAR * pairMotDirSameness * pairGlobalOrg, data.differents_coherentOnly)
summary(lm_RE)

# test plots 
testPlots <- ggplot(data.differents, aes(x = cuedAR, y =responseAR )) + 
  geom_point(shape=19, size=1, alpha=0.1, show.legend = FALSE, aes(color= id) ) +
  geom_smooth(data = data.differents, method=lm, color = "red") +
  ylim(-.7,.7) +
  labs(x="Cued aspect ratio", y = "Reported aspect ratio") +
  ggtitle("Cleaned data") + facet_wrap(~pairMotCoherent)

testPlots + geom_density_2d(color="blue", alpha=0.2) + xlab((testPlots$mapping$x))

## get the aggregated values by participants 
test_aggByPart <- aggregate(new_responseAR ~ id, data.differents, mean)
View(test_aggByPart)
test_aggByPart <- aggregate(new_responseAR ~ id + cuedAR, data.differents, mean)
test_aggByPart <- aggregate(new_responseAR ~ id + cuedAR + uncuedAR, data.differents, mean)
test_aggByPart <- aggregate(new_responseAR ~ id + cuedAR + uncuedAR + pairMotCoherent + pairGlobalOrg, data.differents, mean)

data_agg_diff$cuedARCat <- NA
### get the same value and substract it from the responses to the same 
for (i in 1:nrow(data_agg)) { # for each trial in data,
  
  # get the id of a given trial
  selectedID <- data_agg$id[i]
  selectedCuedAR <- data_agg$cuedAR[i]
  selectedResponseAR <- data_agg$responseAR[i] #this is the response to the identical pair
  
  #look through the data and match the cuedAR with the data
  for (a in 1:nrow(data_agg_diff)){
    if (data_agg_diff$cuedAR[a] == selectedCuedAR && data_agg_diff$id[a] == selectedID){
      data_agg_diff$new_responseAR[a] <- data_agg_diff$responseAR[a] - selectedResponseAR
      data_agg_diff$deductAmount[a] <- selectedResponseAR
    }
    
    # find the AR category of the cued shape
    if (data_agg_diff$cuedAR[i] < 0) {
      data_agg_diff$cuedARCat[i] <- "flat_ellipse"
    }
    if (data_agg_diff$cuedAR[i] > 0) {
      data_agg_diff$cuedARCat[i] <- "tall_ellipse"
    }
    if (data_agg_diff$cuedAR[i] == 0) {
      data_agg_diff$cuedARCat[i] <- "circle"
    }
  }
}

data_agg_diff_tall <- subset(data_agg_diff, data_agg_diff$cuedARCat=="tall_ellipse")
data_agg_diff_circle <- subset(data_agg_diff, data_agg_diff$cuedARCat=="circle")
data_agg_diff_flat_ellipse <- subset(data_agg_diff, data_agg_diff$cuedARCat=="flat_ellipse")

#aggregateById <- aggregate(new_responseAR~id,data_agg_diff, sum)
#aggregateById <- aggregate(abs(new_responseAR)~id,data_agg_diff, mean)
aggregateById_tall <- aggregate(new_responseAR~id,data_agg_diff_tall, mean)
aggregateById_flat <- aggregate(new_responseAR~id,data_agg_diff_flat_ellipse, mean)

t.test(aggregateById$new_responseAR, mu = 0)
t.test(aggregateById_tall$new_responseAR, aggregateById_flat$new_responseAR)
### overall repulsion effect is not significantly different than 0 
aggregateById$idid <- 1:nrow(aggregateById)

ggplot(aggregateById, aes(x=id, y=new_responseAR, group=idid)) +
  geom_line() +
  geom_point()

# test plots 
testPlots <- ggplot(data.differents, aes(x = cuedAR, y =new_responseAR )) + 
  geom_point(shape=19, size=1, alpha=0.1, show.legend = FALSE) +
  geom_smooth(data = data_agg_diff, method=lm, color = "red") +
  ylim(-.7,.7) +
  labs(x="Cued aspect ratio", y = "Reported aspect ratio") +
  ggtitle("Cleaned data") + facet_wrap(~id)

testPlots + geom_density_2d(color="blue", alpha=0.2) + xlab((testPlots$mapping$x))


# OVERALL ANALYSIS
lm_main <- lm(new_responseAR ~ cuedAR, data_agg_diff)
summary(lm_main)

# First make sure that none of the IVs are too correlated (r2<0.6). Manually check each pair...all are fine, except cuedAR and arDiff
res <- cor.test(data$cuedAR, data$arDiff, method = c("pearson"))
r2 <- res$estimate^2 #corr coeff ^2

#### MR: Homogeneity of variance test ####

# Is variance of one factor the same across levels of another factor?
# In our case, no.

#' Use a Levene test, which is basically a one-way anova on SDs across levels of one factor
#' Here the DV is RT and the facvtor is intensity
#' if p val is sig, variances are not homogeneous
#' Report it as a F-stat, with 2dfs, and p-val.
leveneTest(data.sample$responseAR, data.sample$cuedAR, center = mean) 


