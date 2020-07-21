###############################
###  Multiparameter Models  ###
###############################

# So now are looking at use of R
# to summarize Bayesian models with
# several unknown parameters.

# We look at 4 examples:
# 1) unknown mean and variance of a
#    normal population;
# 2) multinomial parameters;
# 3) parameters of logistic regression
#    model; and
# 4) comparing two proportions in a
#    2 x 2 contingency table.

##########################################
# Normal Data with Both Parameters Unknown
##########################################

 ### SEE SLIDES #2 and #3

# Interested in learning about the
# distribution of completion times for
# men ages 20-29 running NY marathon.

# load LearnBayes package
library(LearnBayes)

# Construct contour plot of joint posterior
# density
data(marathontimes)

# attach data set to get at
# variable time easily
attach(marathontimes)

# use function normchi2post() which computes
# the log of a posterior density of a mean M
# and a variance S2 when a sample is taken
# from a normal density and a standard
# noninformative prior is used.

?contour

# Also use function mycontour() that adds
# to use of R contour() function and has
# four inputs: (1) name of function that
# defines log density;
d = mycontour(normchi2post,
              # (2) vector of values of rectangle
              # where density is graphed;
              c(220, 330, 500, 9000), 
              # (3) data to be used; and
              # (4) optional contour parameters
              time, xlab="mean",
              ylab="variance")

# Contour lines are 10%, 1%, and 0.1% of
# the maximum value of posterior density
# over the grid.

# We summarize this posterior distribution 
# by simulation. We simulate (mu, sigma^2)
# from joint posterior by first simulating
# sigma^2 from chi^2 distribution and mu
# from a normal distribution. Simulated
# draws of "scale times inverse chi-square"
# distribution of variance sigma^2 by
# transforming chi^2 draws.

# S is sum of mean squared error of time
S = sum((time - mean(time))^2) 
# n is number of observations
n = length(time)
# posterior variance is simulated with 1000
# draws from "scale times inverse chi^2"
sigma2 = S/rchisq(1000, n - 1)
sigma2
# mean is simulated with 1000 draws
# from normal population
mu = rnorm(1000, mean = mean(time), 
           sd = sqrt(sigma2)/sqrt(n))
mu

# NOTE: normpostsim() implements this
# simulation algorithm.
params = normpostsim(data=time, m=1000)
# one thousand simulated means
params$mu
# one thousand simulated variances
params$sigma2

# display values for mu and sigma2
points(mu, sigma2)
# repeat with LearnBayes normpostsim()
# output for points
d = mycontour(normchi2post,
              # (2) vector of values of rectangle
              # where density is graphed;
              c(220, 330, 500, 9000), 
              # (3) data to be used; and
              # (4) optional contour parameters
              time, xlab="mean",
              ylab="variance")
points(params$mu, params$sigma2)
# End up with contour plot of joint
# posterior distribution of (mu, sigma^2)
# for a normal sampling model. Points
# are simulated random sample from
# this distribution.

# Suppose you want to conduct
# inference on the population:
# 95% CI for mean
quantile(mu, c(0.025, 0.975))
# mean of the population of marathon
# running times is somewhere between
# 255.59 and 300.46 minutes in this sample
# or use normpostsim() output (different draw)
quantile(params$mu, c(0.025, 0.975))

# 95% CI for standard deviation
quantile(sqrt(sigma2), c(0.025, 0.975))
# sd is between 37.48 and 71.61 minutes
# or use normpostsim() output (different draw)
quantile(sqrt(params$sigma2), c(0.025, 0.975))

#######################
# A Multinomial Model
#######################

### SEE WIKIPEDIA

library(LearnBayes)

### SEE SLIDES 5-6
alpha = c(728, 584, 138)
alpha
theta = rdirichlet(1000, alpha)
theta

# We are interested in comparing proportions
# for Bush and Dukakis so we focus on
# difference theta[,1]-theta[,2], create
# histogram of simulated draws of differences
hist(theta[, 1] - theta[, 2], main="")
# Mass of distribution is positive showing
# evidence that proportion of Bush voters
# exceeds proportion of Dukakis voters.

### SEE SLIDES 8-9

# Based on posterior distribution of state
# proportions, can simulate from posterior
# distribution of electoral votes for
# Obama. Dataset election.2008 contains
# the percentage of voters in each state
data(election.2008)
attach(election.2008)

# User-defined function to use simulation
# from Dirichlet distribution to compute
# posterior probability that thetaOj
# exceeds thetaMj in the jth state
prob.Obama=function(j)
{
  p=rdirichlet(5000,
               500*c(M.pct[j],
                     O.pct[j],
                     100-M.pct[j]-O.pct[j])/100+1)
  mean(p[,2]>p[,1])
}

# Compute Obama win probability
# for all states
Obama.win.probs=sapply(1:51,prob.Obama)
Obama.win.probs

# With win probabilities, can simulate
# from posterior distribution of the Obama
# electoral votes by flipping set of 51
# biased coins, where coin probabilities
# correspond to the Obama state win pro-
# babilities.

# Then we compute number of Obama electoral
# votes based on results of the coin flips.

# sim.election() creates one simulation,
# repeat 1000 times with replicate()

sim.election=function()
{
  winner=rbinom(51,1,Obama.win.probs)  
  sum(EV*winner)         
}

sim.EV=replicate(1000,sim.election())
sim.EV

# construct histogram of posterior EVo
hist(sim.EV,min(sim.EV):max(sim.EV),col="blue")
# and show actual Obama electoral vote total
abline(v=365,lwd=3)  # Obama received 365 votes
text(375,30,"Actual \n Obama \n total")

############################
# A Bioassay Experiment
############################

### SEE SLIDES 11-12

# Define covariate vector x and vectors
# of sample sizes and observed success
# counts n and y
x = c(-0.86, -0.3, -0.05, 0.73)
n = c(5, 5, 5, 5)
y = c(0, 1, 3, 5)
data = cbind(x, n, y)
data

# Fit model with ML using glm()
# and summary()
glmdata = cbind(y, n - y)
glmdata
results = glm(glmdata ~ x, family = binomial)
results
summary(results)

# User has beliefs about probability
# of death at two different dose levels,
# XsubL and XsubH

# when x = -.7, median and 90th 
# percentile of p are (.2,.5)
?beta.select

a1.b1=beta.select(list(p=.5,x=.2),
                  list(p=.9,x=.5))
a1.b1

# when x = +.6, median and 90th 
# percentile of p are (.8, .98)
a2.b2=beta.select(list(p=.5,x=.8),
                  list(p=.9,x=.98))
a2.b2

### SEE SLIDE 13

# Plot conditional means prior for bioassay.
# In each bar, point corresponds to median
# and endpoints correspond to quartiles of
# prior distribution for each beta distribution
plot(c(-1,1),c(0,1),type="n",xlab="Dose",ylab="Prob(death)")
lines(-0.7*c(1,1),qbeta(c(.25,.75),a1.b1[1],a1.b1[2]),lwd=4)
lines(0.6*c(1,1),qbeta(c(.25,.75),a2.b2[1],a2.b2[2]),lwd=4)
points(c(-0.7,0.6),qbeta(.5,c(a1.b1[1],a2.b2[1]),c(a1.b1[2],a2.b2[2])),
       pch=19,cex=2)
text(-0.3,.2,"Beta(1.12, 3.56)")
text(.2,.8,"Beta(2.10, 0.74)")
response=rbind(a1.b1,a2.b2)
x=c(-0.7,0.6)
fit = glm(response ~ x, family = binomial)
curve(exp(fit$coef[1]+fit$coef[2]*x)/
        (1+exp(fit$coef[1]+fit$coef[2]*x)),add=T)

### SEE SLIDES 14, 15 and 16

# The log posterior for (Beta0, Beta1) in
# this logistic model is calculate using
# logisticpost() with data a matrix of
# dose, # successes, sample size.

# We combine data with prior data
prior=rbind(c(-0.7, 4.68, 1.12),
            c(0.6, 2.10, 0.74))
prior; data
data.new=rbind(data, prior)
data.new

# To summarize posterior distribution, we
# find rectangle that includes all of
# posterior probability. Use ML to get
# first guess

# Contour plot of posterior distribution
# of (Beta0, Beta1) for bioassay. Contour
# lines are drawn at 10%, 1%, and .1% of
# model height
mycontour(logisticpost,c(-3,3,-1,9),data.new,
          xlab="beta0", ylab="beta1")

# Use simcontour() to simulate pairs of
# (Beta0, Beta1) from posterior density
# computed on this rectangular grid
s=simcontour(logisticpost,c(-2,3,-1,11),data.new,1000)
s
# We put points on contour grid...have 
# contour plot of posterior distribution
# of (Beta0, Beta1) with simulated
# random sample
points(s)

# Display density estimate of simulated
# values of slope parameter Beta1. All
# of mass of density of Beta 1 is positive,
# indicating evidence that increasing
# levels of the dose does increase the
# probability of death
plot(density(s$y),xlab="beta1")

# Parameter of interest is LD-50, value
# of dose X such that prob of death is
# one half. It is equal to -Beta0/Beta1
theta=-s$x/s$y
hist(theta,xlab="LD-50",breaks=20)

# We compute 95% CI
quantile(theta,c(.025,.975))

# So probability that theta is contained
# in (-0.3467, 0.4424) is 0.95

#############################
# Comparing Two Proportions
#############################

# Howard considers general problem of
# comparing proportions from two
# independent binomial distributions.

# Suppose we observe y1 distributed as
# binomial(n1,p1) and y2 distributed as
# binomial(n2,p2). We want to know
# if data favor H1: p1 > p2 or H2: p1 < p2

# We also want a measure of strength of
# the evidence in support of one hypothesis.

# Important task is construction of
# appropriate prior distribution. But
# issue of independence is questionable.

# Suppose we are give info that one
# proportion is equal to a particular
# value, p1 = .8. This can influence
# a user's prior beliefs about p2 by
# believing p2 is also close to .8.

### SEE SLIDE 21

# Suppose alpha = beta = lambda = gamma = 1,
# reflecting prior beliefs about each
# individual parameter. The logarithm of
# the dependent prior is defined by
# howardprior()

# Using mycontour() we create contour
# plots of dependent prior for values
# of association parameters 
# sigma = 2, 1, .5, and .25.
sigma=c(2,1,.5,.25)
plo=.0001;phi=.9999
par(mfrow=c(2,2))
for (i in 1:4)
  mycontour(howardprior,c(plo,phi,plo,phi),c(1,1,1,1,sigma[i]),
            main=paste("sigma=",as.character(sigma[i])),
            xlab="p1",ylab="p2")

# Note that as sigma goes to zero, prior
# is placing more of mass along line where
# two proportions are equal.

### SEE SLIDE 22

# Since posterior is same functional
# form as prior, we use howardprior()
# for posterior calculations and create
# contour plots of posterior for four
# values of association parameter sigma

sigma=c(2,1,.5,.25)
par(mfrow=c(2,2))
for (i in 1:4)
{
  mycontour(howardprior,c(plo,phi,plo,phi),
            c(1+3,1+15,1+7,1+5,sigma[i]),
            main=paste("sigma=",as.character(sigma[i])),
            xlab="p1",ylab="p2")
  lines(c(0,1),c(0,1))
}

# Now we can test hypothesis H1: p1 > p2
# by computing prob of this region of
# the parameter space. Use simcontour()
# to get simulated sample from posterior
# distribution of (p1, p2), and then find
# proportion of simulated pairs where
# p1 > p2 for sigma = 2, 1, .5 and .25

# sigma = 2
s=simcontour(howardprior,c(plo,phi,plo,phi),
             c(1+3,1+15,1+7,1+5,2),1000)
s
# P(p1 > p2)
sum(s$x>s$y)/1000

# sigma = 1
s=simcontour(howardprior,c(plo,phi,plo,phi),
             c(1+3,1+15,1+7,1+5,1),1000)
# P(p1 > p2)
sum(s$x>s$y)/1000

# sigma = .5
s=simcontour(howardprior,c(plo,phi,plo,phi),
             c(1+3,1+15,1+7,1+5,.5),1000)
# P(p1 > p2)
sum(s$x>s$y)/1000

# sigma = .25
s=simcontour(howardprior,c(plo,phi,plo,phi),
             c(1+3,1+15,1+7,1+5,.25),1000)
# P(p1 > p2)
sum(s$x>s$y)/1000
