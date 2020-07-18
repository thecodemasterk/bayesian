#################################
# Single Parameter Models
#################################

# We introduce use of R in summarizing 
# posterior distributions for several 
# single-parameter models.

# Begin by describing Bayesian inference
# for a variance for a normal popu-
# lation and inference for a Poisson
# mean when informative prior infor-
# mation is available.

# We simulate posterior distributions
# from the exponential family.

# In estimating a normal mean, we il-
# lustrate the use of two distinct
# priors in modeling beliefs and show
# that inferences may or may not be
# sensitive to the choice of a prior.

#################################
# Normal Distribution with 
# Known Mean but Unknown Variance
#################################
# Gelman et al (2003) consider a problem
# of estimating an unknown variance with
# American football scores. Looks at
# difference d between a game outcome
# (winning score minus losing score)
# and a published point spread. We
# observe d1, ..., dn, the observed
# difference between game outcomes
# and point spreads for n football
# games. 

### SEE SLIDES 2 and 3

# Load LearnBayes package
library(LearnBayes)

# get the data for 672 games
data(footballscores)

# so can manipulate variables directly
attach(footballscores)

# for each games, have actual scores
# on favorite and underdog teams and
# the spread, the published point spread:
d = favorite - underdog - spread;d

# sample size
n = length(d);n

# sum of the squares of differences:
v = sum(d^2);v

# We simulate 1000 values from posterior
# distribution of sd sigma in two steps:

# 1) Simulate values of precision parameter
# P = 1/sigma-squared from scaled chi-square
# distribution
P = rchisq(1000, n)/v;P

# 2) Perform transformation of sigma (s)
# to get values from the posterior distri-
# bution of standard deviation sigma
s = sqrt(1/P)

# Histogram of sigma:
hist(s)

# Extract 2.5%, 50%, and 97.5% percentiles
# of this simulated posterior distribution. 
# A point estimate for sigma is the posterior
# median 13.86 (this estimate will vary
# with each resample) and the extreme
# percentiles (13.18, 14.62) represent a
# 95% probability interval for sigma.
quantile(s, probs = c(0.025, 0.5, 0.975))

######################################
# Estimating a Heart Transplant 
# Mortality Rate - Poisson mean
######################################

# Consider the problem of learning about
# the rate of success of heart transplant
# surgery at a particular hospital in the
# United States. We observe the number of
# transplant surgeries n, and the number
# of deaths within 30 days of surgery y
# is recorded.

# Also, one can predict the probability of
# death for an individual patient based on
# various factors. 

# Based on these predicted probabilities,
# we can obtain an expected number of deaths
# denoted as e.

# We assume number of deaths y follows a
# Poisson distribution with mean e*lambda
# (lambda is exposure) and the objective 
# is to estimate the mortality rate per 
# unit exposure lambda.

# Standard estimate of lambda is MLE
# estimate lambda-hat = y/e but this estimate
# can be poor when number of deaths y is
# close to zero. Better to use a Bayesian
# estimate with small death counts that uses
# prior knowledge about the size of the
# mortality rate. We choose a prior distri-
# bution which is a member of gamma(a, b)
# density of the form

# SEE SLIDES 5 thru 9

# We consider inference about the heart
# transplant death rate for two hospitals,
# one that has experienced a small number
# of surgeries and a second that has had
# many surgeries.

# Hospital A had 1 death so y-obs = 1
# with an exposure of e = 66. Standard
# estimate of Hospital A's rate, 1/66
# falls into suspect category due to small
# number of deaths.

# So we perform bayesian calculations in R.
# We define gamma prior parameters alpha and
# beta and exposure ex, then find the pre-
# dictive density of the values y = 0, 1,
# ..., 10 using dpois() and dgamma().

# Values of f(y) are are computed for the
# prior mean value lambda = a/b. Note that
# most of the probability of the predictive
# density is concentrated on the two values
# y = 0 and y = 1. We have no reason to
# doubt our Bayesian model since the observed
# number of deaths (y-obs = 1) is in the
# middle of this predictive distribution.

alpha=16;beta=15174
yobs=1; ex=66
y=0:10; y
lam=alpha/beta;lam

# Prior predictive density
py=dpois(y, lam*ex)*dgamma(lam, shape = alpha,
                           rate = beta)/dgamma(lam, shape= alpha + y,
                                               rate = beta + ex)
py

# most of predictive density is between
# two values 0 and 1
cbind(y, round(py, 3))

# We summarize posterior density lambda by
# simulating 1000 values from gamma density:
lambdaA = rgamma(1000, 
                 shape = alpha + yobs, 
                 rate = beta + ex)
lambdaA

# Now we consider the estimation of Hospital B
# which had y-obs = 4 deaths with an exposure
# of e = 1767.
ex = 1767; yobs=4
y = 0:10; y

# Prior predictive density:
py = dpois(y, lam * ex) * dgamma(lam, shape = alpha, 
                                 rate = beta)/dgamma(lam, shape = alpha + y,
                                                     rate = beta + ex)
py

# We again see that observed number of deaths
# is consistent with this model as y-obs = 4
# is not in the extreme tails of this distri-
# bution.
cbind(y, round(py, 3))

# We again summarize posterior density lambda
# by simulating 1000 values from gamma density:
lambdaB = rgamma(1000, 
                 shape = alpha + yobs, 
                 rate = beta + ex)
lambdaB

# So we have two prior densities and we look
# at impact on inference. We display prior
# and posterior density distributions for 
# heart transplant death rate for the two
# hospitals on the same graph.

# puts both on same plot:
par(mfrow = c(2, 1))
# Hospital A posterior:
plot(density(lambdaA), 
     main="HOSPITAL A", 
     xlab="lambdaA", lwd=3)
# Hospital A prior:
curve(dgamma(x, shape = alpha, 
             rate = beta), add=TRUE)
legend("topright",
       legend=c("prior","posterior"),
       lwd=c(1,3))
# Hospital B posterior:
plot(density(lambdaB),
     main="HOSPITAL B", 
     xlab="lambdaB", lwd=3)
# Hospital B prior:
curve(dgamma(x, shape = alpha, 
             rate = beta), add=TRUE)
legend("topright",
       legend=c("prior","posterior"),
       lwd=c(1,3))

##########################################
# An Illustration of Bayesian Robustness
##########################################

# Is possible to have incomplete prior infor-
# mation about a parameter such that a number
# of different priors match the given info.

# Say you believe a priori that median of parameter
# theta is 30 and its 80th percentile is 50,
# there might be many priors that match. Here
# it is desirable that inferences from the pos-
# terior NOT be dependent on some exact func-
# tional form of the prior.

# A Bayesian analysis is said to be 'robust'
# to the choice of prior if the inference is
# insensitive to different priors that may
# match the user's beliefs.

# Suppose you want to estimate the true IQ theta
# for Joe. You believe Joe has average intelligence
# with a median of your prior distribution of 100.
# Also, you are 90% confident that Joe's IQ falls
# between 80 and 120.

# Can use normal.select() to find values of mean
# and standard deviation of normal density that
# match the beliefs that the median is 100 and
# the 95th percentile is 120.

library(LearnBayes)

quantile1=list(p=.5,x=100); quantile1
quantile2=list(p=.95,x=120); quantile2

# We see that normal density with mu = 100
# and r = 12.16 matches this prior information.
normal.select(quantile1, quantile2)

# Then Joe takes 4 IQ tests and his scores
# are y1, y2, y3, y4. If we assume that an
# individual score y is distributed N(theta,sigma)
# with a known sd sigma = 15, the observed
# mean score y-bar N(theta, sigma/4)

# SEE SLIDES 11 and 12

# We illustrate posterior calculations for
# three hypothetical test results for Joe
# with observed mean test scores of y-bars
# equal to 110, 125 or 140. In each case, 
# we compute the posterior mean (mu1) and 
# sd (tau1) of Joe's true IQ theta

# From before:
mu = 100
tau = 12.16
sigma = 15
n = 4
se = sigma/sqrt(4); se
# Joe's test scores:
ybar = c(110, 125, 140); ybar
# posterior sd:
tau1 = 1/sqrt(1/se^2 + 1/tau^2); tau1
# posterior mean:
mu1 = (ybar/se^2 + mu/tau^2) * tau1^2; mu1
summ1=cbind(ybar, mu1, tau1)
summ1

# SEE SLIDE 13 about t density

# So scale parameter lambda
tscale = 20/qt(0.95, 2)
tscale

# We display the normal and t priors in a
# single graph. They have the same basic shape
# but the t density has flatter tails which
# will impact the posterior density for
# "extreme" test scores

par(mfrow=c(1,1))

# Normal and t priors for representing prior
# opinion about a person's true IQ score
curve(1/tscale*dt((x-mu)/tscale,2),
      from=60, to=140, xlab="theta", 
      ylab="Prior Density")
curve(dnorm(x,mean=mu,sd=tau), 
      add=TRUE, lwd=3)
legend("topright",
       legend=c("t density","normal density"),
       lwd=c(1,3))

# SEE SLIDE 15

# We construct a grid of theta values that
# cover the posterior density, compute the
# product of the normal likelihood and the
# t prior on the grid, and convert to prob-
# abilities by dividing by the sum.

# Basically, we are approximating a continuous
# posterior density by a discrete distribution
# on this grid.

# Then we use this discrete distribution to
# compute the posterior mean and sd. The function
# norm.t.compute() implements this computation
# for a single value of y-bar. Using sapply()
# we apply it to 3 values of y-bar.

norm.t.compute=function(ybar) {
  theta = seq(60, 180, length = 500)
  like = dnorm(theta,mean=ybar,sd=sigma/sqrt(n))
  prior = dt((theta - mu)/tscale, 2)
  post = prior * like
  post = post/sum(post)
  m = sum(theta * post)
  s = sqrt(sum(theta^2 * post) - m^2)
  c(ybar, m, s) }

summ2=t(sapply(c(110, 125, 140),
               norm.t.compute))
summ2
dimnames(summ2)[[2]]=c("ybar","mu1 t","tau1 t")

# 2nd and 3rd columns represent mean and sd
# of posteriors
summ2

# We compare posteriors of theta using
# the two priors
cbind(summ1,summ2)

# When y-bar is 110, means and sd's are similar
# using normal and t prior. However, they can
# be different when the observed mean score
# is inconsistent with the prior mean.

# We graph posterior densities for the two
# priors (normal and t) in the extreme case
# where y-bar = 140.
theta=seq(60, 180, length=500)
normpost = dnorm(theta, mu1[3], tau1)
normpost = normpost/sum(normpost)
plot(theta,normpost,type="l",
     lwd=3,ylab="Posterior Density")
like = dnorm(theta,mean=140,
             sd=sigma/sqrt(n))
prior = dt((theta - mu)/tscale, 2)
tpost = prior * like / sum(prior * like)  
lines(theta,tpost)
legend("topright",
       legend=c("t prior","normal prior"),
       lwd=c(1,3))

# Inference about the mean is robust to
# choice of prior (normal or t) when observed
# mean IQ score is consistent with prior
# beliefs. But when an extreme IQ score
# is observed, the inference is not 
# robust to the choice of prior density.
