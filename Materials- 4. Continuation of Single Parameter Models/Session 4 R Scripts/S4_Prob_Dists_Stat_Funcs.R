#########################################
#####   Probability Distributions   #####
#####     and Statistical Models    #####
#########################################

### Probability Distributions

# R provides a comprehensive set of statistical tables. 
# Functions are provided to evaluate the cumulative 
# distribution function P(X <= x), the probability 
# density function and the quantile function 
# (given q, the smallest x such that P(X <= x) > q), 
# and to simulate from the distribution. 

#   Distribution   R name   additional arguments
#   ------------   ------   --------------------
#   beta           beta     shape1, shape2, ncp 
#   binomial       binom    size, prob 
#   Cauchy         cauchy   location, scale 
#   chi-squared    chisq    df, ncp 
#   exponential    exp      rate 
#   F              f        df1, df2, ncp 
#   gamma          gamma    shape, scale 
#   geometric      geom     prob 
#   hypergeometric hyper    m, n, k 
#   log-normal     lnorm    meanlog, sdlog 
#   logistic       logis    location, scale 
#   negative       binomial nbinom size, prob 
#   normal         norm     mean, sd 
#   Poisson        pois     lambda 
#   signed         rank     signrank n 
#   Student's t    t        df, ncp 
#   uniform        unif     min, max 
#   Weibull        weibull  shape, scale 
#   Wilcoxon       wilcox   m, n 

# Prefix the name given here by 'd' for the 
# density, 'p' for the CDF, 'q' for the quantile
# function and 'r' for simulation (random deviates). 

# For every distribution there are four commands. 
# The commands for each distribution are prepended
# with a letter to indicate the functionality:

# "d" returns the height of the probability density function 
# "p" returns the cumulative density function 
# "q" returns the inverse cumulative density function (quantiles) 
# "r" returns randomly generated numbers 

# We look at the normal distribution
help(Normal)

# we first look at dnorm. Given a set of values 
# it returns the height of the probability 
# distribution at each point. 

# If you only give the points it assumes you
# want to use a mean of zero and standard 
# deviation of one. Also can use different
# values for the mean and standard deviation:
  
dnorm(0)
# [1] 0.3989423
dnorm(0)*sqrt(2*pi)
# [1] 1
dnorm(0,mean=4)
# [1] 0.0001338302
dnorm(0,mean=4,sd=10)
# [1] 0.03682701
v <- c(0,1,2)
dnorm(v)
# [1] 0.39894228 0.24197072 0.05399097
x1 <- seq(-20,20,by=.1);x1
x2 <- seq(-4,4,by=.02);x2
# dnorm(x1) assumes mean=0 and sd=1
dnorm(x1)
y1 <- dnorm(x1);y1
# dnorm(x2) assumes mean=0 and sd=1
y2 <- dnorm(x2);y2
plot(x1,y1)
# should "look" normal with mean=0, sd=1
plot(x2,y2)
# plot with mean=4.5 and sd=5
y3 <- dnorm(x1,mean=4.5,sd=6);y3
plot(x1,y3)

# Let's look at pnorm...Given a number 
# or a list it computes the probability
# that a normally distributed random 
# number will be less than that number. 

# This function also goes by the title
# "Cumulative Distribution Function." 

# It accepts the same options as dnorm:
pnorm(0)
# [1] 0.5
pnorm(1)
# [1] 0.8413447
pnorm(0,mean=2)
# [1] 0.02275013
pnorm(0,mean=2,sd=3)
# [1] 0.2524925
v <- c(0,1,2);v
pnorm(v)
# [1] 0.5000000 0.8413447 0.9772499
x1 <- seq(-4,4,by=.02);x1
x2 <- seq(-20,20,by=.1);x2
y1 <- pnorm(x1,sd=1);y1
y2 <- pnorm(x1,sd=1.5);y2
y3 <- pnorm(x1,sd=2);y3
# successively gets "flatter"
plot(x1,y1)
plot(x1,y2)
plot(x1,y3)
y4 <- pnorm(x2,mean=-3,sd=4);y4
y5 <- pnorm(x2,mean=0,sd=4);y5
y6 <- pnorm(x2,mean=3,sd=4);y6
# successively shifts to the right
plot(x2,y4)
plot(x2,y5)
plot(x2,y6)

# Others: Converting t-values to p-values
# 2-tailed p-value for t distribution
2*pt(-2.43, df = 13)

### Examining the distribution of a set of data
# Given a (univariate) set of data we can examine 
# its distribution in a large number of ways. 

# The simplest is to examine the numbers. 
# Two slightly different summaries are given
# by summary() and fivenum() and a display 
# of the numbers by stem ("stem and leaf" plot).  

# attach faithful data
attach(faithful)
str(faithful)
dim(faithful)
head(faithful)
summary(eruptions)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  1.600   2.163   4.000   3.488   4.454   5.100

# Tukey's five number summary (minimum, lower-hinge,
# median, upper-hinge, maximum) for the input data:
fivenum(eruptions)
# [1] 1.6000 2.1585 4.0000 4.4585 5.1000

# stem-and-leaf plot
stem(eruptions)

# A stem-and-leaf plot is like a histogram:
hist(eruptions)
## make the bins smaller, make a plot of density
hist(eruptions, seq(1.6, 5.2, 0.2), prob=TRUE)
lines(density(eruptions, bw=0.1))
rug(eruptions) # show the actual data points

# We can plot the empirical cumulative distribution
# function by using the function ecdf.   
plot(ecdf(eruptions), 
     do.points=FALSE, 
     verticals=TRUE)

# This distribution is obviously far from 
# any standard distribution. How about the 
# right-hand mode, say eruptions of longer 
# than 3 minutes? Let us fit a normal 
# distribution and overlay the fitted CDF. 
long <- eruptions[eruptions > 3]
length(long)
plot(ecdf(long), 
     do.points=FALSE, 
     verticals=TRUE)
x <- seq(3, 5.4, 0.01)
lines(x, pnorm(x, mean=mean(long), 
               sd=sqrt(var(long))), lty=3)

# Quantile-quantile (Q-Q) plots can help us 
# examine this more carefully.    

# arrange for a square figure region
par(pty="s")
qqnorm(long) 
qqline(long)

# which shows a reasonable fit but a 
# shorter right tail than one would 
# expect from a normal distribution. 

# Let us compare this with some simulated 
# data from a t distribution:
x <- rt(250, df = 5)
qqnorm(x); qqline(x)

# which will usually (with a random sample) 
# show longer tails than expected for a normal.

# We can make a Q-Q plot against the generated
# distribution by 
qqplot(qt(ppoints(250), df = 5), x, 
       xlab = "Q-Q plot for t dsn")
qqline(x)

# Finally, we might want a more formal test of 
# agreement with normality (or not). R provides 
# the Shapiro-Wilk test   
shapiro.test(long)

# Shapiro-Wilk normality test

# data:  long
# W = 0.9793, p-value = 0.01052
