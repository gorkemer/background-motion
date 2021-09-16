
# fit gaussian
x <-  rgamma(50, 10)   # simulating a random sample from the gamma distribution 
ig_fit(x)              # fitting an inverse Gaussian distribution to x 

# fit a distribution
#data("groundbeef", package = "fitdistrplus")
#my_data <- groundbeef$serving
my_data <- data.frame(data.differents$uncuedAR, data.differents$normed_responseAR)     # unkonwn distribution parameters
#my_data <- data.differents$responseAR      # unkonwn distribution parameters

# remove na in r - remove rows - na.omit function / option
#my_data <- na.omit(my_data)

plot(my_data)
ggplot(my_data, aes(x = data.differents.uncuedAR, y = data.differents.normed_responseAR)) + 
  geom_point(shape=19, size=1, alpha=0.1, show.legend = FALSE) +
  #geom_smooth(data = my_data, method=lm, color = "red")
  stat_smooth(method = "nls", formula = y ~ SSasymp(x, Asym, R0, lrc), se = FALSE)

fit <- fitdistr(my_data, densfun="normal")  # we assume my_data ~ Normal(?,?)
fit
hist(my_data, pch=20, breaks=25, prob=TRUE, main="")
curve(dnorm(x, fit$estimate[1], fit$estimate[2]), col="red", lwd=2, add=T)


#Choosing which distribution to fit

plot(rnorm(100), rnorm(100), main = "Another plot.")

plot(my_data,pch=20)
plotdist(my_data, histo = TRUE, demp = TRUE) #plot the empirical density and the histogram to gain insight of the data
descdist(my_data, discrete=FALSE, boot=500) #Another tool is to show some descriptive statistics, like the moments, to help us in making a decision:
#fitting a distrubution
fit_ln <- fitdist(my_data$abs.data.differents.normed_responseAR., "lnorm")
cdfcomp(fit_ln, xlogscale = TRUE, ylogscale = TRUE)

fit_w  <- fitdist(my_data, "weibull")
fit_g  <- fitdist(my_data$abs.data.differents.normed_responseAR., "gamma")
fit_ln <- fitdist(my_data, "lnorm")

#get the start values
library(GeneralizedHyperbolic)
nvstrt <- nigFitStart(my_data, startValues = "FN")
fitinvgauss <- fitdist(my_data, "invgauss", 
                       start = list(mean=-0.02356182, dispersion=0.27256418, shape=0.5))
summary(fitinvgauss)


#The package also provides some goodness-of-it statistics:
#gofstat(list(fit_ln, fit_ll, fit_P, fit_B), fitnames = c("lnorm", "llogis", "Pareto", "Burr"))

#parameter estimates
ests <- bootdist(fit_B, niter = 1e3)


ggplotRegression <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}

### 1

st_norm <- function(x) dnorm(x,0,1)
first_deriv <- function(x) {}
second_deriv <- function(x) {}
body(first_deriv) <- D(body(st_norm), 'x')
body(second_deriv) <- D(body(first_deriv), 'x')

curve(st_norm, -4, 4, ylim = c(-0.4, 0.4), col = 'blue')
curve(first_deriv, -4, 4, add = T, col = 'red')
curve(second_deriv, -4, 4, add = T, col = 'green')
abline(v=0, h=0)

### 2 
dnorm_deriv1 <- function(x, mean = 0, sd = 1) {
  return(-((x-mean)/(sd^2))*dnorm(x, mean, sd))
} 

dnorm_deriv2 <- function(x, mean = 0, sd = 1) {
  return((((x-mean)^2) / (sd^4))*dnorm(x, mean, sd) 
         - (1/sd^2)*dnorm(x, mean, sd))
}

curve(dnorm, -4, 4, ylim = c(-0.4, 0.4), col = 'blue')
curve(dnorm_deriv1, -4, 4, add = T, col = 'red')
curve(dnorm_deriv2, -4, 4, add = T, col = ' green')
abline(v=0, h=0)

curve(dnorm(x, 2, 2), -4, 8, ylim = c(-0.1, 0.2), col = 'blue')
curve(dnorm_deriv1(x, 2, 2), -4, 8, add = T, col = 'red')
curve(dnorm_deriv2(x, 2, 2), -4, 8, add = T, col = ' green')
abline(v=2, h=0)


##### 


# required model
library(mgcv)

# make data
n <- 200
tmp <- seq(0,20*pi,,n)
x <- tmp / (2*pi)
mon <- x%%1
err <- rnorm(n, sd=0.5)
y <- sin(tmp) + err + 1
plot(x, y, t="l")
df <- data.frame(x, y, mon)

# GAM with intercept
fit1 <- gam(y ~ s(mon, bs = "cc", k = 12), data=df)
summary(fit1)
plot(fit1)

# GAM without intercept
fit2 <- gam(y ~ s(mon, bs = "cc", k = 12) - 1, data=df) # note "-1" for no intercept
summary(fit2)
plot(fit2)


#### Trying to fit non-linear data ####



x = my_data$data.differents.uncuedAR
y = my_data$data.differents.normed_responseAR

f <- nls(y~m*x+b)
f <- nls(y~x*a*w*c*exp(-(w*x^2))) 

m <- nls(y~x*a*w*c*exp(-(w*x)^2))
m <- nls(y~x*a*w*c*exp(-(w*x)^2), start=list(a=.1,w=.1))

lines(x,predict(m),col="red",lty=2,lwd=3)

#### 

#My equation search shows a good fit to a three-parameter inverse Harris yield density
#equation, "y = x / (a + b * pow(x, c))", with parameters a = 1.4956613575678071E+01, 
#b = 7.8559465184281589E-05, and c = 2.1768293119284090E+00 giving RMSE = 0.1002 and R-squared = 0.9943


# build the model
m = loess(y ~ x)

# see model summary
summary(m)

# plot points and model fit
plot(y,x)
lines(y, m$fitted, col=2)

