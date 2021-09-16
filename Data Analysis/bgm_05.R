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
#### 0- Select Data ####
#data = read.csv('bgm_cleanData_cor.csv', header=TRUE)
#data$slopeCoeff <- NA # !! add slopeCoeff & slopePval columns to this data frame because the code below expects to see that.
#data$slopePval <- NA # !! Dont worry, it does not change anything. Just a cosmetic issue. 
#data = read.csv('untouchedData.csv', header = TRUE)
#data = read.csv('bgm_cleanData_lm.csv', header = TRUE)
data = read.csv('data.Clean_p05.csv', header = TRUE)
data$slopeCoeff <- NA # !! add slopeCoeff & slopePval columns to this data frame because the code below expects to see that.
data$slopePval <- NA # !! Dont worry, it does not change anything. Just a cosmetic issue. 
#doing the pooled normalization on COHERENT data only. Removing the random data from the set.
#data = subset(data, data$pairMotCoherent==1) # 1) coherent
#doing the pooled normalization on RANDOM data only. Removing the coherent data from the set.
#data = subset(data, data$pairMotCoherent==0) # 1) random
# define the data frame to use
data.Same <- subset(data, data$pairAR_sameCondition==1)
data.Diff <- subset(data, !(data$pairAR_sameCondition==1))

# look at aggregated responses 
#data_agg_same <- aggregate(responseAR ~ cuedAR, data = data.Same, mean)
#data_agg_diff <- aggregate(responseAR ~ cuedAR, data = data.Diff, mean)

# with id
data_agg_same <- aggregate(responseAR ~ cuedAR + id, data = data.Same, mean)
data_agg_diff <- aggregate(responseAR ~ cuedAR + id, data = data.Diff, mean)


fullARList <- c(-0.46332, -0.41699,-0.37066, -0.32432, -0.27799, -0.23166, -0.18533, -0.13900, -0.09266, -0.04633,  0.00000,  0.04633,  0.09266,  0.13899,  0.18533,  0.23166,  0.27799, 0.32432,  0.37065,  0.41699,  0.46332)                                      

####  1- find the number of people who have "same" responses for the full range of cuedAR ####
#pre
data_agg_same_sum <- aggregate(cuedAR ~ id, data = data_agg_same,sum)
nrow(data_agg_same_sum) #76 (104 for untouched data)
#post
data_agg_same_sum <- subset(data_agg_same_sum, data_agg_same_sum$cuedAR==data_agg_same_sum$cuedAR[1])
nrow(data_agg_same_sum) #20/76  (26/104 in untouched) of the people dont have the full range same-responses

#### 2- deviation from the population response of individual responses to identical pairs ####
popResp_same <- aggregate(responseAR ~ cuedAR, data = data.Same, mean) #population response
indResp_same<- aggregate(responseAR ~ cuedAR + id, data = data.Same, mean)#individual response 

indResp_same$popInvDiff <- NA
                         
for (i in 1:nrow(popResp_same)){
  selectedCuedAR <- popResp_same$cuedAR[i]
  popResp <- popResp_same$responseAR[i]
  for (a in 1:nrow(indResp_same)){
    if (selectedCuedAR == indResp_same$cuedAR[a]){
      indResp_same$popInvDiff[a] <-  indResp_same$responseAR[a] - popResp
      indResp_same$popResp[a] <-  popResp
    }
  }
}

# finding 2SD +- from the mean of the deviation #
mean_popResp <- mean(indResp_same$popResp)
sd_popResp <- sd(indResp_same$popResp)
indResp_same_outliers <- subset(indResp_same, indResp_same$popInvDiff <= mean_popResp - (2*sd_popResp) | indResp_same$popInvDiff >= mean_popResp + (2*sd_popResp))
nrow(indResp_same) #2119 row in total
nrow(indResp_same_outliers) # and 273 of the trial's responses are very different from the population response
unique(indResp_same_outliers$id) #almost everyone have some responseAR that deviates from the population response by a significant margin. 
# it might not be wise to go with the same-all approach. 

#### 3- for each participant, for each cuedAR, how many same trials are there? ####

data.Same_arranged <- arrange(data.Same, id, cuedAR)
#table( data.Same_arranged[ , c("Color" , "Type") ] )

selectData <- data.Same_arranged[c("cuedAR", "uncuedAR", "pairAR_sameCondition", "id" )]


only1 <- subset(selectData, selectData$id == "o1")
only1.freq <- as.data.frame(table(only1))
only1.freq$sameCheck <- ifelse(only1.freq$cuedAR == only1.freq$uncuedAR, "same", "different")
only1.freq <- subset(only1.freq, only1.freq$sameCheck=="same")

# find how many times they were faced with each cuedAR-uncuedAR identical pairs for each occurence
selectData.freq <- as.data.frame(table(selectData))
selectData.freq$sameCheck <- ifelse(selectData.freq$cuedAR == selectData.freq$uncuedAR, "same", "different")
selectData.freq <- subset(selectData.freq, selectData.freq$sameCheck=="same")
# I am thinking about sharing this table in the supplementary section

test <- aggregate(Freq~id, data = selectData.freq, sum ) #this shows the same-trial number per participant. It varies, naturally.
# we can even use this number per id, and cut at some point, and use population response for all of the responses of that particular participant. 
test2<- aggregate(pairAR_sameCondition~id, data = data.Same, sum) #these tables should be identical, and it is identical. Above method
#works. Now I will see if some AR are represented less than others. 
test3 <- aggregate(Freq~id + cuedAR, data = selectData.freq, sum ) #this shows the same-trial number per participant for each cuedAR occurences. It varies, naturally.
## test3 shows that the occurences, and it is not fairly distributed.

## more importantly, I will identify Cued-UncuedAR pairs where participants lack an appropriate number of identical-pairs to be used as a matching value 
## (used for conducting the normalization process). For those trials, population response will be used for normalization operation. 

selectData.freq$popNeededHere <- ifelse(selectData.freq$Freq<6, "here", "not-needed") # threshold "6" is arbitrarily defined 
nrow(subset(selectData.freq, selectData.freq$popNeededHere == "not-needed")) ##individual response'un yeterli oldugu trial sayisi
nrow(subset(selectData.freq, selectData.freq$popNeededHere == "here")) ## need population response for this number of trials.
# with 6 threshold. 991/1193, almost half of the data has to be removed. This shows the necessity of conducting a mixed approach. 
# we shouldn't do directly population response normalization either, see above results (" #almost everyone deviates from the population response by a significant margin."). 

# I can now add the population response for that specific TRIAL (this is important: I am not employing a "bulk" population response for a given participant)
# for that individual. This way I will conserve their individual bias as much as possible.
# technically, their biases to other trials are kept, in this way. 

#bring back the population response for each cuedAR again. 
data_agg_same <- aggregate(responseAR ~ cuedAR, data = data.Same, mean)

for (i in 1:nrow(selectData.freq)) { # get the response to the same
  # get the id of a given trial
  if (selectData.freq$popNeededHere[i]== "here"){
    cuedAR_brokenOne <- selectData.freq$cuedAR[i] #get the cuedAR to use as a index for the following line
    selectData.freq$popResponse[i] <- data_agg_same$responseAR[which(data_agg_same$cuedAR == cuedAR_brokenOne)] #find the row number and get the pop response of that value
  }
  else{
    selectData.freq$popResponse[i] <- 0.000
  }
}

# now I add individual responses to each CuedAR at the "not-needed" place
#indResp_same <- aggregate(responseAR ~ cuedAR + id, data = data.Same, mean)#individual response 

#create the attribute
selectData.freq$indvResponseValue <- NA
#go through the selectData.freq list, find the not-needed row
for (i in 1:nrow(selectData.freq)){
  if (selectData.freq$popNeededHere[i]=="not-needed"){
    idOfthatPerson <- selectData.freq$id[i] #get the id 
    corespondCuedAR <- selectData.freq$cuedAR[i] #get the cuedAR
    #now use those variables to find the corresponding individual responseAR to that cuedAR for that ID
    #selectData.freq$rowNumber[i] <- which(indResp_same$id == idOfthatPerson & indResp_same$cuedAR == corespondCuedAR)
    selectData.freq$indvResponseValue[i] <- indResp_same$responseAR[which(indResp_same$id == idOfthatPerson & indResp_same$cuedAR == corespondCuedAR)]
  }
  else if (selectData.freq$popNeededHere[i]=="here"){
    selectData.freq$indvResponseValue[i] <- 0.00
  }
}

# now I need to merge the two columns to finalize the normalization. I can change "not-needed" to 0 and then sum the two columns
# merging the two to finalize the "mixed approach" model for the new normalization method
selectData.freq$respToMatch_mixedApp <- selectData.freq$popResponse + selectData.freq$indvResponseValue

# run regression to see which pairs interact with the normed response, they all do. Not sure how to interpret this.
lm_main <- lm(respToMatch_mixedApp ~  cuedAR * uncuedAR , selectData.freq)
summary(lm_main) 

#### 4- new analysis with the mixed approach ####
#can move this to a separate script

#now with the new matching scores, calculate the normed response 
data.Diff <- subset(data, !(data$pairAR_sameCondition==1)) #get the whole data set for different trials

# aggregate to find the averaged response to CuedAR by id 
selectData.freq_agg <- aggregate(respToMatch_mixedApp ~ cuedAR + id, selectData.freq, mean)
#data_agg_same <- aggregate(responseAR ~ cuedAR + id, data = data.Same, mean) # this is the normalization with all participants (no-mixed approach, just all)

#selectData.freq_agg.singleP <- subset(selectData.freq_agg, selectData.freq_agg$id=="o1" | selectData.freq_agg$id == "o2")
#data.Diff.singleP <- subset(data.Diff, data.Diff$id=="o1" | data.Diff$id=="o2")

#find the rows where cuedAR is -0.46332
data.Diff$matchResponseAR_mix <- NA
totalData <-matrix(nrow=10,ncol=34)
totalData <- data.frame(subNum=NA, trial_number =NA, cuedAR = NA,uncuedAR = NA, responseAR = NA,
                        response_error =NA, rt=NA, shape_org =NA, cueType = NA,e1_motion_dir = NA,
                        e2_motion_dir =NA, round_number=NA, coherence =NA, cued_motion_dir = NA,uncued_motion_dir = NA,
                        arDiff =NA, pairAvg=NA, cuedMotDir =NA, pairGlobalOrg = NA,pairMotDirSameness = NA,
                        pairMotDirSameness_wRandom =NA, pairMotSame=NA, pairMotDiff =NA, pairMotRandom = NA,pairSameARCat = NA,
                        pairAR_sameCondition =NA, motDirComp=NA, pairMotCoherent =NA, pairMotDiago = NA,id = NA,
                        corrCoeff =NA, corrPval=NA, slopeCoeff =NA, slopePval = NA, matchResponseAR_mix = NA)
          
#go through the numbers and index them to add norm_response to those rows
for (i in 1:length(unique(data.Diff$id))){ #

  baseID <- unique(data.Diff$id)[i] #get the id to be worked on
  print(baseID)
  
  mainData <- subset(data.Diff, data.Diff$id ==baseID) #subset by id, do the analysis on the id-subsetted data frame only. 
  sameResponses <- subset(selectData.freq_agg, selectData.freq_agg$id== baseID) #define the data sheet with the mixed approach same-responses
  
  for (z in 1:length(fullARList)){ #go through each cuedAR for the selected baseID
    baseCuedAR <- fullARList[z]
    print(baseCuedAR)
    rowNumbers <- which(mainData$cuedAR==baseCuedAR)    #find the row numbers with that cuedAR
    print(rowNumbers)
    for (a in 1:length(rowNumbers)){ #go through each row number in the loop
      #use row numbers to find the row number with the normResponse to be included.
      mainData$matchResponseAR_mix[rowNumbers[a]] <- sameResponses$respToMatch_mixedApp[which(sameResponses$cuedAR==baseCuedAR)]#include normResponseAR_mix to those trials
      }
  }
  #totalData <- cbind(totalData, mainData)
  totalData <- rbind(totalData, mainData) #add id-subsetted data frame that was being worked on above lines to the "totalData"
}

#remove the first row, because it was generated when totalData matrix was initially created.
totalData = totalData[-1,]

#test
data.Diff.test <- totalData[,c("id","cuedAR", "uncuedAR", "responseAR", "matchResponseAR_mix")]
#okay this method works. tested with 2 participants it also worked. Okay, test for all people: normed with mixed approach is applied to the main data successfully.

#### 5- running linear analysis ####
#first I will calculate the deviation from the response to the matched.
totalData$normedResponseUsingMixApp <- totalData$responseAR - totalData$matchResponseAR_mix

# 1- main analysis
lm_main <- lm(normedResponseUsingMixApp~cuedAR * uncuedAR *pairGlobalOrg, data = totalData)
summary(lm_main) #overall attractive influence of uncuedAR, but when shapes are localized into the same hemifield, uncuedAR played a repulsive role. 
# also, global organization influenced the responses (same as reported at Sweeny2011) 

# 2- shared motion direction grouping analysis 
lm_main <- lm(normedResponseUsingMixApp~cuedAR * uncuedAR *pairGlobalOrg *pairMotSame, data = totalData)
summary(lm_main) # grouping by motion direction facilitated the attractive influence of uncuedAR, as well as the cuedAR.
# again, when shapes are localized into the same hemifield, uncuedAR played a repulsive role. 

#grouping analysis:
lm_main <- lm(normedResponseUsingMixApp~cuedAR * uncuedAR *pairGlobalOrg *pairMotDirSameness_wRandom, data = totalData)
summary(lm_main)
lm_main <- lm(normedResponseUsingMixApp~cuedAR * uncuedAR *pairMotDirSameness_wRandom, data = totalData)
summary(lm_main)

anova(lm_main)
# same Category analysis
totalData.sameARCat <- subset(totalData, totalData$pairSameARCat == 1)
lm_main <- lm(normedResponseUsingMixApp~cuedAR * uncuedAR *pairGlobalOrg, data = totalData.sameARCat)
summary(lm_main) # if they have the same cat AND located within hemifield, uncued again plays a repulsive role. 
# this is congruent with the idea that Ps might be pulling them apart to achieve a satisfactory performance at the 
# behavioral task.

#### 6- coherent only analysis #### 
#check only for coherent only responses. 
totalData.coherentOnly <- subset(totalData, totalData$pairMotCoherent==1)
totalData.randomOnly <- subset(totalData, totalData$pairMotCoherent==0)


# lm_main <- lm(normedResponseUsingMixApp~cuedAR * uncuedAR *pairGlobalOrg, data = totalData.coherentOnly)
# summary(lm_main)


lm_motion <- lm(normedResponseUsingMixApp~ uncuedAR * pairSameARCat, data = totalData)
summary(lm_motion)

lm_motionGlobalOrg <- lm(normedResponseUsingMixApp~ uncuedAR * pairMotSame * pairGlobalOrg, data = totalData.coherentOnly)
summary(lm_motionGlobalOrg)

lm_motionSameCat <- lm(normedResponseUsingMixApp~ uncuedAR * pairMotSame * pairSameARCat, data = totalData.coherentOnly)
summary(lm_motionSameCat)

lm_forth <- lm(normedResponseUsingMixApp~ uncuedAR *pairSameARCat * pairMotSame * pairGlobalOrg, data = totalData.coherentOnly)
summary(lm_forth)

lm_third <- lm(normedResponseUsingMixApp~cuedAR * uncuedAR * pairMotSame *pairGlobalOrg, data = totalData.coherentOnly)
summary(lm_third)

lm_motion <- lm(normedResponseUsingMixApp~ arDiff * pairMotSame, data = totalData.coherentOnly)
summary(lm_motion)

##### ONE BIG MODEL ####

totalData$threeLevels <- ifelse(((as.numeric(totalData$cued_motion_dir) ==9999)),0,  # "RANDOM DIRECTION"
                             ifelse(((as.numeric(totalData$cued_motion_dir) == as.numeric(totalData$uncued_motion_dir))),1, # "SAME"
                                    -1) # "DIFFERENT")
)

# seperate variables
lm_motion <- lm(normedResponseUsingMixApp~ arDiff * pairMotDiff * pairMotSame * pairMotRandom , data = totalData)
summary(lm_motion)

#single variable, -1,0,1
lm_motion <- lm(normedResponseUsingMixApp~ arDiff *threeLevels , data = totalData)
summary(lm_motion)

#### 7- Plotting ####

# New facet label names for pairSameARCat.labels variable

totalData.coherentOnly$pairSameARCat_3Levels <- ifelse(
  totalData.coherentOnly$pairSameARCat== 0 & totalData.coherentOnly$arDiff < 0, -1,
  ifelse(
    totalData.coherentOnly$pairSameARCat== 0 & totalData.coherentOnly$arDiff > 0, 0,
    1
  )
)

pairSameARCat_names <- c(
  `1` = "Same Shape Category",
  `0` = 'Different Shape Category'
)

pairSameARCat3Levels_names <- c(
  `1` = "Same Shape Category",
  `-1` = 'Different Shape Category (UC:Flat)',
  `0` = 'Different Shape Category (UC:Tall)'
)

#first plot the relationship between cuedAR and responseAR
mainPlot <- ggplot(totalData, aes(x = uncuedAR, y = normedResponseUsingMixApp)) + 
  #geom_point(shape=19, size=1, alpha=0.2, show.legend = FALSE) +
  geom_jitter(alpha=0.2, size=0.5)+
  #geom_smooth(data = totalData.coherentOnly, method=lm, color = "red") +
  #geom_point(aes(shape = factor(pairMotSame))) +
  #geom_point(aes(color = factor(pairMotSame))) +
  #geom_smooth(data = totalData, method=lm, color = "red") +
  geom_smooth(method = "lm",
              aes(color = factor(threeLevels))) +
  #geom_smooth(method = "loess", span = 0.4,
              #aes(color = factor(pairMotSame))) +
  labs(x="uncuedAR", y = "Normalized Reported AR") + #"AR Difference Between Shapes"
  ggtitle("", subtitle ="") +
  #xlab((testPlots$mapping$x)) + 
  labs(color = "Motion Sharedness") + 
  coord_cartesian(ylim=c(-0.46, 0.46))+ #, xlim= c(-0.05,0.05)) +  #+  facet_wrap(~pairGlobalOrg)+ #+ facet_wrap(~cuedMotDir) + #+ facet_wrap(~pairSameARCat) +
  #scale_color_manual(name="Motion Direction",
                   #labels=c("Unshared","Shared")
                   #,values=c("red","black")) +
                   #,values=c("#F8766D","#00BFC4")) +
  geom_segment(aes(x=-0.46,xend=0.46,y=0,yend=0), linetype="dotted")
  #scale_x_continuous(breaks = seq(-0.46, 0.46, by = 0.46)) + 
  #scale_y_continuous(breaks = seq(-0.05, 0.05, by = 0.05))


mainPlot + theme_classic() + geom_density_2d(color="blue", alpha=0.4) #+ facet_wrap(~pairSameARCat_3Levels, labeller = labeller(pairSameARCat_3Levels = pairSameARCat3Levels_names)) 
   #+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))




#first plot the relationship between cuedAR and responseAR
mainPlot_perID <- ggplot(totalData.coherentOnly, aes(x = uncuedAR, y = normedResponseUsingMixApp)) + 
  geom_point(shape=19, size=1, alpha=0.2, show.legend = FALSE) +
  #geom_smooth(data = totalData.coherentOnly, method=lm, color = "red") +
  #geom_point(aes(shape = factor(pairMotSame))) +
  #geom_point(aes(color = factor(pairMotSame))) +
  #geom_smooth(data = data, method=lm, color = "red") +
  geom_smooth(method = "lm", 
              aes(color = factor(pairMotSame))) +
  labs(x="Uncued AR", y = "Normalized Reported AR", color = "Motion Sharedness") +
  ggtitle("", subtitle ="") +
  coord_cartesian(ylim=c(-0.30, 0.30), xlim= c(-0.46,0.46)) + 
  scale_color_manual(name="Motion Direction",
                     labels=c("Unshared","Shared"),
                     values=c("red","black")) +
                     #values=c("darkorchid2","turquoise2")) +
  geom_segment(aes(x=-0.46,xend=0.46,y=0,yend=0), linetype="dotted") + 
  scale_x_continuous(breaks = seq(-0.46, 0.46, by = 0.46)) + 
  scale_y_continuous(breaks = seq(-0.30, 0.30, by = 0.30)) +
  geom_density_2d(color="blue", alpha=0.2)


mainPlot_perID + theme_classic() + theme(panel.spacing = unit(1, "lines"), axis.text.x = element_text(size=7)) + facet_wrap(~id,nrow = 6) #+ guides(x = "none") + #x = "none"

####
testPlots <- ggplot(totalData.coherentOnly, aes(x = cuedAR, y =normedResponseUsingMixApp )) + 
  geom_point(shape=19, size=1, alpha=0.7, show.legend = FALSE) +
  geom_smooth(method = "lm", 
              aes(color = factor(pairMotSame))) +
  ylim(-.7,.7) +
  labs(x="Cued aspect ratio", y = "Reported aspect ratio") +
  ggtitle("Cleaned data") + facet_wrap(~id)

testPlots + geom_density_2d(color="blue", alpha=0.2) + xlab((testPlots$mapping$x))


#### 8.a. - Population normalization approaches ####
#now with the "population responses" scores, calculate the normed response 
fullARList <- c(-0.46332, -0.41699,-0.37066, -0.32432, -0.27799, -0.23166, -0.18533, -0.13900, -0.09266, -0.04633,  0.00000,  0.04633,  0.09266,  0.13899,  0.18533,  0.23166,  0.27799, 0.32432,  0.37065,  0.41699,  0.46332)                                      

#get the whole data set for different trials
data.Diff_popResp <- subset(data, !(data$pairAR_sameCondition==1)) 

# aggregate to find the averaged response to CuedAR by id 
data_agg_same_popResp <- aggregate(responseAR ~ cuedAR, data = data.Same, mean) # this is the normalization with all participants (no-mixed approach, just all)
data_agg_same_indvResp <- aggregate(responseAR ~ cuedAR + id, data = data.Same, mean) # this is the normalization with all participants (no-mixed approach, just all)

#find the rows where cuedAR is -0.46332
data.Diff_popResp$matchResponseAR_popResp <- NA
totalData_popResp <-matrix(nrow=10,ncol=34)
totalData_popResp <- data.frame(subNum=NA, trial_number =NA, cuedAR = NA,uncuedAR = NA, responseAR = NA,
                        response_error =NA, rt=NA, shape_org =NA, cueType = NA,e1_motion_dir = NA,
                        e2_motion_dir =NA, round_number=NA, coherence =NA, cued_motion_dir = NA,uncued_motion_dir = NA,
                        arDiff =NA, pairAvg=NA, cuedMotDir =NA, pairGlobalOrg = NA,pairMotDirSameness = NA,
                        pairMotDirSameness_wRandom =NA, pairMotSame=NA, pairMotDiff =NA, pairMotRandom = NA,pairSameARCat = NA,
                        pairAR_sameCondition =NA, motDirComp=NA, pairMotCoherent =NA, pairMotDiago = NA,id = NA,
                        corrCoeff =NA, corrPval=NA, slopeCoeff =NA, slopePval = NA, matchResponseAR_popResp = NA)

#go through the numbers and index them to add norm_response to those rows
for (i in 1:length(unique(data.Diff_popResp$id))){ #
  
  baseID <- unique(data.Diff_popResp$id)[i] #get the id to be worked on
  print(baseID)
  
  mainData_popResp <- subset(data.Diff_popResp, data.Diff_popResp$id ==baseID) #subset by id, do the analysis on the id-subsetted data frame only. 
  #sameResponses <- subset(selectData.freq_agg, selectData.freq_agg$id== baseID) 
  
  for (z in 1:length(fullARList)){ #go through each cuedAR for the selected baseID
    baseCuedAR <- fullARList[z]
    print(baseCuedAR)
    rowNumbers <- which(mainData_popResp$cuedAR==baseCuedAR)    #find the row numbers with that cuedAR
    print(rowNumbers)
    for (a in 1:length(rowNumbers)){ #go through each row number in the loop
      #use row numbers to find the row number with the normResponse to be included.
      mainData_popResp$matchResponseAR_popResp[rowNumbers[a]] <- data_agg_same_popResp$responseAR[which(data_agg_same_popResp$cuedAR==baseCuedAR)]#include normResponseAR_mix to those trials
    }
  }
  #totalData <- cbind(totalData, mainData)
  totalData_popResp <- rbind(totalData_popResp, mainData_popResp) #add id-subsetted data frame that was being worked on above lines to the "totalData"
}

#remove the first row, because it was generated when totalData matrix was initially created.
totalData_popResp = totalData_popResp[-1,]

### running linear analysis ###
#first I will calculate the deviation from the response to the matched.
totalData_popResp$normedResponseUsingPopRespApp <- totalData_popResp$responseAR - totalData_popResp$matchResponseAR_popResp

### coherent only analysis ###
#check only for coherent only responses. 
totalData.coherentOnly_popResp <- subset(totalData_popResp, totalData_popResp$pairMotCoherent==1)
totalData.randomOnly_popResp <- subset(totalData_popResp, totalData_popResp$pairMotCoherent==0)

lm_motion <- lm(normedResponseUsingPopRespApp~ uncuedAR * pairMotSame, data = totalData.coherentOnly_popResp)
summary(lm_motion)

lm_motionGlobalOrg <- lm(normedResponseUsingPopRespApp~ uncuedAR * pairMotSame * pairGlobalOrg, data = totalData.coherentOnly_popResp)
summary(lm_motionGlobalOrg)

lm_motionSameCat <- lm(normedResponseUsingPopRespApp~ uncuedAR * pairMotSame * pairSameARCat, data = totalData.coherentOnly_popResp)
summary(lm_motionSameCat)

lm_forth <- lm(normedResponseUsingPopRespApp~ uncuedAR *pairSameARCat * pairMotSame * pairGlobalOrg, data = totalData.coherentOnly_popResp)
summary(lm_forth)

lm_third <- lm(normedResponseUsingPopRespApp~cuedAR * uncuedAR * pairMotSame *pairGlobalOrg, data = totalData.coherentOnly_popResp)
summary(lm_third)


####8.b. Individual Response Normalization ####

fullARList <- c(-0.46332, -0.41699,-0.37066, -0.32432, -0.27799, -0.23166, -0.18533, -0.13900, -0.09266, -0.04633,  0.00000,  0.04633,  0.09266,  0.13899,  0.18533,  0.23166,  0.27799, 0.32432,  0.37065,  0.41699,  0.46332)                                      

#get the whole data set for different trials
data.Diff_indvResp <- subset(data, !(data$pairAR_sameCondition==1)) 

# aggregate to find the averaged response to CuedAR by id 
data_agg_same_popResp <- aggregate(responseAR ~ cuedAR, data = data.Same, mean) # this is the normalization with all participants (no-mixed approach, just all)
data_agg_same_indvResp <- aggregate(responseAR ~ cuedAR + id, data = data.Same, mean) # this is the normalization with all participants (no-mixed approach, just all)

#find the rows where cuedAR is -0.46332
data.Diff_indvResp$matchResponseAR_indvResp <- NA
totalData_indvResp <-matrix(nrow=10,ncol=34)
totalData_indvResp <- data.frame(subNum=NA, trial_number =NA, cuedAR = NA,uncuedAR = NA, responseAR = NA,
                                response_error =NA, rt=NA, shape_org =NA, cueType = NA,e1_motion_dir = NA,
                                e2_motion_dir =NA, round_number=NA, coherence =NA, cued_motion_dir = NA,uncued_motion_dir = NA,
                                arDiff =NA, pairAvg=NA, cuedMotDir =NA, pairGlobalOrg = NA,pairMotDirSameness = NA,
                                pairMotDirSameness_wRandom =NA, pairMotSame=NA, pairMotDiff =NA, pairMotRandom = NA,pairSameARCat = NA,
                                pairAR_sameCondition =NA, motDirComp=NA, pairMotCoherent =NA, pairMotDiago = NA,id = NA,
                                corrCoeff =NA, corrPval=NA, slopeCoeff =NA, slopePval = NA, matchResponseAR_indvResp = NA)

#go through the numbers and index them to add norm_response to those rows
for (i in 1:length(unique(data.Diff_indvResp$id))){ #
  
  baseID <- unique(data.Diff_indvResp$id)[i] #get the id to be worked on
  print(baseID)
  
  mainData_indvResp <- subset(data.Diff_indvResp, data.Diff_indvResp$id ==baseID) #subset by id, do the analysis on the id-subsetted data frame only. 
  sameResponses <- subset(data_agg_same_indvResp, data_agg_same_indvResp$id== baseID) 
  
  for (z in 1:length(fullARList)){ #go through each cuedAR for the selected baseID
    baseCuedAR <- fullARList[z]
    print(baseID)
    print(baseCuedAR)
    rowNumbers <- which(mainData_indvResp$cuedAR==baseCuedAR)    #find the row numbers with that cuedAR
    print(rowNumbers)
    for (a in 1:length(rowNumbers)){ #go through each row number in the loop
      #print(rowNumbers[a])
      #use row numbers to find the row number with the normResponse to be included.
      if (!(length(sameResponses$responseAR[which(sameResponses$cuedAR==baseCuedAR)]) == 0)){
        # if there is a same trial for a given cuedAR, this code will run. I am checking whether there is 
        # a responseAR for a given cuedAR in the SAME list. Some people don't have same trial for some cuedAR
        # for example, "o6" doesn't have same trials with cuedAR 0.32 
        mainData_indvResp$matchResponseAR_indvResp[rowNumbers[a]] <- sameResponses$responseAR[which(sameResponses$cuedAR==baseCuedAR)]#include normResponseAR_mix to those trials
      }
      else {
        #this means there is not identical trial for a given cuedAR for a given participant
        mainData_indvResp$matchResponseAR_indvResp[rowNumbers[a]] <- NA # I just add for them      
        }
    }
  }
  totalData_indvResp <- rbind(totalData_indvResp, mainData_indvResp) #add id-subsetted data frame that was being worked on above lines to the "totalData"
}

#remove the first row, because it was generated when totalData matrix was initially created.
totalData_indvResp = totalData_indvResp[-1,]
# below is also important, I now taking those columns back again, because I need to remove trials where people had no identical pairs
# I do that by removing NA in the list. However, since, those people have columns of NA, it deletes all the data. I don't want that. 
totalData_indvResp<-totalData_indvResp[,-grep("slopeCoeff|slopePval",colnames(totalData_indvResp))]# delete columns slopeCoeff slopePval 

#remove rows with NA in them. This way I eliminate rows where participants had no matching responseAR for normalization to occur.
totalData_indvResp <- na.omit(totalData_indvResp) #around 300-400 removed. 29725/30533

### running linear analysis ###
#first I will calculate the deviation from the response to the matched.
totalData_indvResp$normedResponseUsingIndvRestApp <- totalData_indvResp$responseAR - totalData_indvResp$matchResponseAR_indvResp

### coherent only analysis ###
#check only for coherent only responses. 
totalData.coherentOnly_indvResp <- subset(totalData_indvResp, totalData_indvResp$pairMotCoherent==1)
totalData.randomOnly_indvResp <- subset(totalData_indvResp, totalData_indvResp$pairMotCoherent==0)

lm_motion <- lm(normedResponseUsingIndvRestApp~ uncuedAR * pairMotSame, data = totalData.coherentOnly_indvResp)
summary(lm_motion)

lm_motionGlobalOrg <- lm(normedResponseUsingIndvRestApp~ uncuedAR * pairMotSame * pairGlobalOrg, data = totalData.coherentOnly_indvResp)
summary(lm_motionGlobalOrg)

lm_motionSameCat <- lm(normedResponseUsingIndvRestApp~ uncuedAR * pairMotSame * pairSameARCat, data = totalData.coherentOnly_indvResp)
summary(lm_motionSameCat)

lm_forth <- lm(normedResponseUsingIndvRestApp~ uncuedAR *pairSameARCat * pairMotSame * pairGlobalOrg, data = totalData.coherentOnly_indvResp)
summary(lm_forth)

lm_third <- lm(normedResponseUsingIndvRestApp~cuedAR * uncuedAR * pairMotSame *pairGlobalOrg, data = totalData.coherentOnly_indvResp)
summary(lm_third)


#### 9- Random Only Analysis ####
totalData.randomOnly <- subset(totalData, totalData$pairMotCoherent==0)

lm_motion <- lm(normedResponseUsingMixApp~ uncuedAR * pairGlobalOrg, data = totalData.randomOnly)
summary(lm_motion)

lm_motionGlobalOrg <- lm(normedResponseUsingMixApp~ arDiff * pairGlobalOrg, data = totalData.randomOnly)
summary(lm_motionGlobalOrg)

lm_motionSameCat <- lm(normedResponseUsingMixApp~ uncuedAR * pairMotSame * pairSameARCat, data = totalData.randomOnly)
summary(lm_motionSameCat)

lm_forth <- lm(normedResponseUsingMixApp~ uncuedAR *pairSameARCat * pairMotSame * pairGlobalOrg, data = totalData.randomOnly)
summary(lm_forth)

lm_third <- lm(normedResponseUsingMixApp~cuedAR * uncuedAR * pairMotSame *pairGlobalOrg, data = totalData.randomOnly)
summary(lm_third)



pairSameARCat_names <- c(
  `1` = "Same Shape Category",
  `0` = 'Different Shape Category'
)

pairSameARCat3Levels_names <- c(
  `1` = "Same Shape Category",
  `-1` = 'Different Shape Category (UC:Flat)',
  `0` = 'Different Shape Category (UC:Tall)'
)

#first plot the relationship between cuedAR and responseAR
mainPlot <- ggplot(totalData.randomOnly, aes(x = arDiff, y = normedResponseUsingMixApp)) + 
  #geom_point(shape=19, size=1, alpha=0.2, show.legend = FALSE) +
  geom_jitter(alpha=0.2, size=0.5)+
  #geom_smooth(data = totalData.randomOnly, method=lm, color = "red") +
  #geom_point(aes(shape = factor(pairMotSame))) +
  #geom_point(aes(color = factor(pairMotSame))) +
  #geom_smooth(data = totalData.randomOnly, method=lm, color = "red") +
  geom_smooth(method = "lm",
              aes(color = factor(pairGlobalOrg))) +
  #geom_smooth(method = "loess", span = 0.4,
  #aes(color = factor(pairMotSame))) +
  labs(x="AR Difference Between Shapes", y = "Normalized Reported AR") +
  ggtitle("Random Motion Condition", subtitle ="Responses are attracted to the uncued AR") +
  #xlab((testPlots$mapping$x)) + 
  labs(color = "Motion Sharedness") + 
  coord_cartesian(ylim=c(-0.46, 0.46))+ #, xlim= c(-0.46,0.46)) +  #+  facet_wrap(~pairGlobalOrg)+ #+ facet_wrap(~cuedMotDir) + #+ facet_wrap(~pairSameARCat) +
  scale_color_manual(name="Pair Global Organization",
                     labels=c("Flat/Between","Tall/Within")
                     ,values=c("darkgray","black")) +
  #,values=c("#F8766D","#00BFC4")) +
  geom_segment(aes(x=-0.46,xend=0.46,y=0,yend=0), linetype="dotted")
#scale_x_continuous(breaks = seq(-0.46, 0.46, by = 0.46)) + 
#scale_y_continuous(breaks = seq(-0.05, 0.05, by = 0.05))


mainPlot + theme_classic() #+ facet_wrap(~pairGlobalOrg) # labeller = labeller(pairSameARCat_3Levels = pairSameARCat3Levels_names)) 
+ geom_density_2d(color="blue", alpha=0.4) + facet_wrap(~pairSameARCat_3Levels, labeller = labeller(pairSameARCat_3Levels = pairSameARCat3Levels_names)) 
#+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))



#### 10- All 3 Motion Types together ####

# pairMotDirSameness_wRandom already has 3 levels, I will use that for facet_wrap

pairMotDirSameness_wRandom_names <- c(
  `2` = "Random Motion",
  `1` = 'Coherent Motion - Same',
  `0` = 'Coherent Motion - Different'
)


mainPlot <- ggplot(totalData, aes(x = uncuedAR, y = normedResponseUsingMixApp)) + 
  #geom_point(shape=19, size=1, alpha=0.2, show.legend = FALSE) +
  geom_jitter(alpha=0.2, size=0.5)+
  #geom_smooth(data = totalData.randomOnly, method=lm, color = "red") +
  #geom_point(aes(shape = factor(pairMotSame))) +
  #geom_point(aes(color = factor(pairMotSame))) +
  #geom_smooth(data = totalData.randomOnly, method=lm, color = "red") +
  geom_smooth(method = "lm",
              aes(color = factor(pairMotDirSameness_wRandom))) +
  #geom_smooth(method = "loess", span = 0.4,
  #aes(color = factor(pairMotSame))) +
  labs(x="AR Difference Between Shapes", y = "Normalized Reported AR") +
  ggtitle("All Motion Types", subtitle ="Relationship between uncuedAR and Normed Response") +
  #xlab((testPlots$mapping$x)) + 
  labs(color = "Motion Sharedness") + 
  coord_cartesian(ylim=c(-0.46, 0.46)) +   #, xlim= c(-0.46,0.46)) +  #+  facet_wrap(~pairGlobalOrg)+ #+ facet_wrap(~cuedMotDir) + #+ facet_wrap(~pairSameARCat) +
  scale_color_manual(name="Motion Types",
                     labels=c("Coherent Motion - Different Direction", 
                              "Coherent Motion - Same Direction", 
                              "Random Motion")
                        ,values=c("red","black","blue")) + #
  #,values=c("#F8766D","#00BFC4")) +
  geom_segment(aes(x=-0.46,xend=0.46,y=0,yend=0), linetype="dotted")
#scale_x_continuous(breaks = seq(-0.46, 0.46, by = 0.46)) + 
#scale_y_continuous(breaks = seq(-0.05, 0.05, by = 0.05))


mainPlot + theme_classic() + geom_density_2d(color="blue", alpha=0.4) + facet_wrap(~pairMotDirSameness_wRandom, labeller = labeller(pairMotDirSameness_wRandom = pairMotDirSameness_wRandom_names)) 
#+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

#### 11- All Analysis Looking at Coherence ####

lm_motion <- lm(normedResponseUsingMixApp~ arDiff * coherence * pairMotSame, data = totalData)
summary(lm_motion)

#### test #### 
#loess fitting#
secondaryPlot <- ggplot(totalData.randomOnly, aes(x = uncuedAR, y = normedResponseUsingMixApp)) + 
  geom_point(shape=19, size=1, alpha=0.5, show.legend = FALSE) +
  #geom_smooth(data = totalData.coherentOnly, method=lm, color = "red") +
  #geom_point(aes(shape = factor(pairMotSame))) +
  #geom_point(aes(color = factor(pairMotSame))) +
  geom_smooth(data = totalData.randomOnly, method=lm, color = "red") +
  labs(x="arDiff", y = "Normalized Reported AR", color = "Motion Sharedness") +
  ggtitle("", subtitle ="") +
  coord_cartesian(ylim=c(-0.46, 0.46), xlim= c(-0.46,0.46)) + 
  scale_color_manual(name="Motion Direction",
                     labels=c("Unshared","Shared"),
                     values=c("red","black")) +
  geom_segment(aes(x=-0.46,xend=0.46,y=0,yend=0), linetype="dotted") + 
  scale_x_continuous(breaks = seq(-0.46, 0.46, by = 0.46)) + 
  #geom_density_2d(color="blue", alpha=0.2) + 
  theme_classic() 

secondaryPlot + facet_wrap(~id, nrow=5) # + geom_smooth(span = 0.4, method = "loess", formula = y ~ x,  aes(color = factor(pairMotSame))) # + facet_wrap(~pairMotSame)


#different formulas
secondaryPlot + geom_smooth(span = 0.8, method = lm, formula = y ~ splines::bs(x, 3),  aes(color = factor(pairMotSame))) # + facet_wrap(~pairMotSame)
#different formulas
secondaryPlot + geom_smooth(span = 0.3, method = lm, formula = y ~ poly(x, 2),  aes(color = factor(pairMotSame))) # + facet_wrap(~pairMotSame)

## per id ## 
secondaryPlot_perID <- ggplot(totalData.coherentOnly, aes(x = uncuedAR, y = normedResponseUsingMixApp)) + 
  #geom_point(shape=19, size=1, alpha=0.5, show.legend = FALSE) +
  #geom_smooth(data = totalData.coherentOnly, method=lm, color = "red") +
  #geom_point(aes(shape = factor(pairMotSame))) +
  #geom_point(aes(color = factor(pairMotSame))) +
  #geom_smooth(data = data, method=lm, color = "red") +
  geom_smooth(span = 0.8, method = "loess", 
              aes(color = factor(pairMotSame))) +
  labs(x="Uncued AR", y = "Normalized Reported AR", color = "Motion Sharedness") +
  ggtitle("", subtitle ="") +
  coord_cartesian(ylim=c(-0.30, 0.30), xlim= c(-0.46,0.46)) + 
  scale_color_manual(name="Motion Direction",
                     labels=c("Unshared","Shared"),
                     values=c("red","black")) +
  geom_segment(aes(x=-0.46,xend=0.46,y=0,yend=0), linetype="dotted") + 
  scale_x_continuous(breaks = seq(-0.46, 0.46, by = 0.46)) + 
  scale_y_continuous(breaks = seq(-0.30, 0.30, by = 0.30))


secondaryPlot_perID + theme_classic() + theme(panel.spacing = unit(1, "lines"), axis.text.x = element_text(size=7)) + facet_wrap(~id,nrow = 6) #+ guides(x = "none") + #x = "none"




  # + geom_density_2d(color="blue", alpha=0.2) 

#### below is for additional testing, will be organized later ###

nOfSame <- nrow(subset(totalData.coherentOnly, totalData.coherentOnly$pairMotSame == 1)) # 6524 same trials
nOfDiff <- nrow(subset(totalData.coherentOnly, totalData.coherentOnly$pairMotSame == 0)) # 19742 different trials
nOfRandom <- nrow(subset(totalData.randomOnly, totalData.randomOnly$pairMotSame == 0))   # 12847 random trials

ggplot(data = totalData.coherentOnly, mapping = aes(x = pairMotSame, y = normedResponseUsingMixApp, group=uncuedAR)) +
  geom_boxplot() + geom_jitter(alpha = 0.3, color = "tomato")

ggplot(data = totalData.coherentOnly, aes(x = uncuedAR, y = normedResponseUsingMixApp, color = pairMotSame)) +
  geom_line() +
  facet_wrap(facets =  vars(id))

scatter3d(x = totalData.coherentOnly$uncuedAR, y = totalData.coherentOnly$normedResponseUsingMixApp, z = totalData.coherentOnly$cuedAR,
          point.col = "blue", surface=FALSE, 
          xlab = "uncuedAR", 
          ylab = "responseAR",
          zlab = "cuedAR", groups = as.factor(totalData.coherentOnly$pairMotSame), ellipsoid = TRUE) # add concentration ellipses: ellipsoid = TRUE

library(plotly)
# in this case we use VAR4 as continuous, you can put color = ~as.factor(VAR4) to have it as factors
plot_ly(totalData.coherentOnly, x = ~uncuedAR, y = ~normedResponseUsingMixApp, z = ~cuedAR, color = ~as.factor(pairMotSame)) %>%
  add_markers(colors = c("#00FA9A34","#B22222dd","#00BFFF66"),) %>%
  layout(scene = list(xaxis = list(title = 'uncuedAR'),
                      yaxis = list(title = 'responseAR'),
                      zaxis = list(title = 'cuedAR')))

library(GGally)
ggpairs(totalData.coherentOnly, columns=c("uncuedAR","normedResponseUsingMixApp", "cuedAR"),
        upper=list(continuous="smooth_loess"),lower=list(continuous="smooth_loess"))
# no idea how to interpret this

#parallel coordinate plot 
vars <- c("uncuedAR","normedResponseUsingMixApp","cuedAR", "pairMotSame")
pcp.vars <- select(totalData.coherentOnly,vars)
ggparcoord(data = pcp.vars)

library(ggbeeswarm)
library(ggforce)
# lined up beeswarm version
# geom_sina with geom_violin
ggplot(data = totalData.coherentOnly) +
  aes(y = normedResponseUsingMixApp, x = uncuedAR, group = as.factor(pairMotSame)) +
  geom_violin() +
  geom_sina() +
  coord_flip()

#### checking overall suppression index ####

total_data_aggByID <- aggregate(normedResponseUsingMixApp ~ id, totalData, mean)
t.test(total_data_aggByID$normedResponseUsingMixApp, mu = 0) #again, not significantly different than 0.
# however, I think this is happening because the data involves both repulsive and attractive influences, and
# when averaged across, it nulls the effect. 
