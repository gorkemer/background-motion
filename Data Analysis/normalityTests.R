
# finally, making sure regression assumptions are met

#### MR: Normality test ####

#' Are the data normally distributed? It looks normal to me. 
#' Step 1 is to take a look at a histogram and compare with a normal dist
#' with the same mean and sd as the data. 
#' Here I used aes(y = ..density..) to convert count to density, so the dist matches the data:
histogram <- ggplot(totalData, aes(normedResponseUsingMixApp))
histogram + geom_histogram(aes(y = ..density..)) + #aes(y = ..density..) converts count to density. 
  stat_function(fun = dnorm, args = list(mean = mean(totalData$normedResponseUsingMixApp), sd = sd(totalData$normedResponseUsingMixApp)))

#' Then make a quantile-quantile plot, which plots cumulative probs of data against normal dist
#' If it's normal, then it should be a match, and produce a line. 
#' If not, then it'll be compressed or expanded at some points:
qqplot.RT <- qplot(sample = totalData$normedResponseUsingMixApp, stat="qq")
qqplot.RT
#As you are plotting the graph between the theoretical quantiles and observed residuals, 
#if your linear model is good enough then the distribution of these two (theoretical quantiles 
#and observed residuals) statistic variable will be very much similar, so by using QQ plot you 
#can determine this similarity between the distribution and thus verifying your linear model.
#It's important to understand that residuals only correspond to the error terms if our model is correct. 
#And because of that, a QQ-Plot can never logically VALIDATE our model - since the error terms are not observed.
#It can, however, invalidate our model and assumptions. And that's why we need to look at it.


data.sample <- totalData[sample(nrow(totalData), 5000),]
shapiro.test(data.sample$normedResponseUsingMixApp) # clearly, our data are not normally distributed

#' Use a Levene test, which is basically a one-way anova on SDs across levels of one factor
#' Here the DV is RT and the facvtor is intensity
#' if p val is sig, variances are not homogeneous
#' Report it as a F-stat, with 2dfs, and p-val.
leveneTest(data.sample$normedResponseUsingMixApp, data.sample$cuedAR, center = mean) 
# not sure how to interpret this because I don't think responseAR variance and cuedAR variance should be similar? no? Why?


#### checking and plotting Regression Model Assumptions ####
# 1- Linearity of the data
plot(lm_main, 1) 
# ideally, the residual plot will show no fitted pattern. 
#That is, the red line should be approximately horizontal at zero. 
#The presence of a pattern may indicate a problem with some aspect of the linear model.

# 2- Normality of residuals
plot(lm_main, 2)
#The normal probability plot of residuals should approximately follow a straight line.
#In our case, all the points fall approximately along this reference line, so we can assume normality. but there
#are some deviations in the data, but it is not major I think. 

# 3- Homogeneity of variance
plot(lm_main, 3)
#This plot shows if residuals are spread equally along the ranges of predictors. 
#It’s good if you see a horizontal line with equally spread points.


#' The Shapiro-Wilk test directly evaluates normality. 
#' It's in the output of stat.desc (normtest.W; test stat (W), p-val)
#' but you can also just run the test on any variable.
#' If it's significant, your data are not normally distributed.
shapiro.test(data.sample$normedResponseUsingMixApp) 
#  Take note that if the sample size is greater than 5000, you should use test statistics 
# instead of the p-value as the indicator to decide. Our sample size is 27000. I think our data is normally distributed. 

# next? https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.boxcox.html can use box-cox transformation
# in python 
# read also this: https://stats.stackexchange.com/questions/2492/is-normality-testing-essentially-useless