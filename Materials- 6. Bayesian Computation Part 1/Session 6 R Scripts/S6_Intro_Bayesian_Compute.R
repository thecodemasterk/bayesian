########################################
# Introduction to Bayesian Computation #
########################################

# Previously, we used two types of strategies
# for the summarization of posterior distributions:

# (1) If the sampling density has a familiar
# functional form, such as exponential family, 
# and a conjugate prior is used for the parameter, 
# then the posterior distribution often is
# expressible in terms of familiar probability
# distributions. In this case, we can
# simulate parameters directly by using the 
# R collection of random variate functions
# (such as rnorm, rbeta, and rgamma), and 
# we can summarize the posterior using computations
# on this simulated sample. 

# (2) The "brute-force" method is a second
# type of computing strategy. Where the
# posterior distribution is not a familiar 
# functional form, one can simply compute
# values of the posterior over a grid of 
# points and then approximate the continuous
# posterior by a discrete posterior concentrated
# on the values of the grid. The brute-force 
# method can be applied for one and for
# two-parameter problems.

# We now describe the Bayesian computational
# problem and introduce some more sophisticated
# computational methods. 

# One approach is based on the behavior of
# the posterior distribution about its mode. 
# This gives a multivariate normal approximation
# to the posterior that serves as a good
# first approximation in the development
# of more exact methods. 

# We then provide an introduction to the use
# of simulation in computing summaries of 
# the posterior distribution.

# When one can directly simulate samples
# from the posterior distribution, then
# the Monte Carlo algorithm gives an 
# estimate and associated standard error
# for the posterior mean for any function 
# of the parameters of interest. 

# Where the posterior distribution is not
# a standard functional form, we can use
# rejection sampling with a suitable choice
# of proposal density as an alternative
# method for producing posterior draws. 

# Importance sampling and sampling importance
# resampling (SIR) algorithms are alternative 
# general methods for computing integrals 
# and simulating from a general posterior
# distribution. The SIR algorithm is useful 
# when one wishes to investigate the sensi-
# tivity of a posterior distribution with 
# respect to changes in the prior and
# likelihood functions.

#####################
# Computing Integrals
#####################

# SEE SLIDES 2-3

##########################
# Setting up Problem in R
##########################

# Before we describe some general 
# summarization methods, we describe
# setting up a Bayesian problem in R. 

# Suppose one is able to write an explicit
# expression for the joint posterior density. 

# In writing this expression, one doesn't
# need to include normalizing constants 
# that don't involve the parameters. 

# For the algorithms described here, we might
# reparameterize all parameters so that 
# they are all real-valued. If one has a
# positive parameter such as a variance, 
# then transform using a log function. If
# one has a proportion parameter p, then 
# it can be transformed to the real line
# by the logit function logit(p) = log(p/(1 - p)).

# After the posterior density has been
# expressed in terms of transformed
# parameters, the first step in summarizing
# this density is to write an R function
# defining the logarithm of the joint
# posterior density.

# To apply the functions described here, 
# theta is assumed to be a vector with 
# k components: theta = (theta1, ..., thetak). 

# The input data is a vector of observed
# values or a list of data values and 
# other model specifications such as
# the values of prior hyperparameters. 

# The function returns a single value of
# the log posterior density.

# SEE SLIDE 4

# We write the function mylogposterior()

# theta contains joint parameters
mylogposterior=function(theta,data)
{
  # how many observations
  n=length(data)
  # extract mean and standard deviation
  mu=theta[1]; sigma=exp(theta[2])
  # create function logf to compute 
  # posterior real value log
  logf = function(y, mu, sigma)
    # get the density with probability
    # p returned as the log(p)
    dnorm(y,mean=mu,sd=sigma,log=TRUE)
  # add to it the sum of log likelihood terms
  val=dnorm(mu, mean=10, sd=20,
            log=TRUE)+sum(logf(data,mu,sigma))
  return(val)
}

##########################################
# A Beta-Binomial Model for Overdispersion
##########################################

# SEE SLIDES 5-6-7

# We write function betabinexch0() to compute
# logarithm of the posterior density. Inputs
# are theta, a vector containing the values
# of ?? and K, and data, a matrix having as
# columns the vector of counts {yj} and 
# the vector of sample sizes {nj}.

betabinexch0=function (theta, data)
{
  eta = theta[1]
  K = theta[2]
  y = data[, 1]
  n = data[, 2]
  N = length(y)
  logf = function(y, n, K, eta) lbeta(K * eta + y, K * (1 - eta) + n - y) - lbeta(K * eta, K * (1 - eta))
  val = sum(logf(y, n, K, eta))
  val = val - 2 * log(1 + K) - log(eta) - log(1 - eta)
  return(val)
}

# load LearnBayes package
library(LearnBayes)

# read in dataset
data(cancermortality)

# We use function mycontour() with
# the log density function above to 
# display a contour plot of the
# posterior density of (eta, K)
mycontour(betabinexch0,
          c(.0001,.003,1,20000),
          cancermortality,
          xlab="eta",ylab="K")

# Note the strong skewness in the
# density, especially with large values
# of the precision parameter K. This
# right-skewness is common with likeli-
# hood functions of a precision or
# variance parameter.

# SEE SLIDE 8

# We write another function betabinexch()
# which computes the logarithm of the 
# posterior for the parameters (logit
# mean and log precision) in a beta-
# binomial model

betabinexch=function (theta, data)
{
  eta = exp(theta[1])/(1 + exp(theta[1]))
  K = exp(theta[2])
  y = data[, 1]
  n = data[, 2]
  N = length(y)
  logf = function(y, n, K, eta) lbeta(K * eta + y, K * (1 - eta) + n - y) - lbeta(K * eta, K * (1 - eta))
  val = sum(logf(y, n, K, eta))
  val = val + theta[2] - 2 * log(1 + exp(theta[2]))
  return(val)
}

# Now we create a contour plot of the 
# posterior of (??1, ??2) using the
# mycontour() function.
mycontour(betabinexch,
          c(-8,-4.5,3,16.5),
          cancermortality,
          xlab="logit eta",
          ylab="log K")

# Although the density has an unusual 
# shape, the strong skewness has been
# reduced and the distribution is more 
# amenable to the computational methods
# described later in this session.

#########################################
# Approximations Based on Posterior Modes
#########################################

# SEE SLIDES 11-12-13

# After one writes an R function to 
# evaluate the log posterior density, the
# function laplace() in the LearnBayes 
# package finds the joint posterior mode 
# by using optim and the default Nelder-Mead algorithm. 

# The inputs to laplace() are the function
# defining the joint posterior, an 
# intelligent guess at the posterior
# mode, and data and parameters used in 
# the definition of the log posterior.

# The choice of "intelligent guess" can 
# be important since the algorithm may
# fail to converge with a poor choice of 
# starting value. Suppose that a suitable
# starting value is used and laplace is 
# successful in finding the posterior mode.

# The output of laplace() is a list with 
# four components. The component mode
# gives the value of the posterior mode ^??, 
# the component var is the associated
# variance-covariance matrix V , the 
# component int is the approximation to
# the logarithm of the prior predictive 
# density, and converge indicates if the
# algorithm converged.

#################
# The Example
#################

# library(LearnBayes)

# data(cancermortality)

# We illustrate the use of the function
# laplace() for beta-binomial modeling
# example. Based on our contour plot, 
# we start the Nelder-Mead method with
# the initial guess (logit(??), logK) = (???7, 6).

fit=laplace(betabinexch,c(-7,6),cancermortality)
fit

# We find the posterior mode to be (???6.82, 7.58). 
# From the output of laplace, we have the 
# approximation that (logit(??), logK) is 
# approximately bivariate normal with mean
# vector fit$mode and variance-covariance matrix fit$var.

# Using mycontour() function with the log 
# bivariate normal function lbinorm(),
# we display the contours of the approximate 
# normal density. Comparing the previous
# plot, we see significant differences 
# between the exact and approximate normal posteriors.

npar=list(m=fit$mode,v=fit$var)
mycontour(lbinorm,c(-8,-4.5,3,16.5),npar,
          xlab="logit eta", ylab="log K")

# An advantage of this algorithm is 
# that one obtains quick summaries of
# the parameters by using the multivariate
# normal approximation. By using the 
# diagonal elements of the variance-covariance
# matrix, one can construct approximate 
# probability intervals for logit(??) and logK. 

# The following code constructs 90% probability
# intervals for the parameters:
se=sqrt(diag(fit$var))
fit$mode-1.645*se
fit$mode+1.645*se

############################################
# Monte Carlo Method for Computing Integrals
############################################

# SEE SLIDE 15

# The Monte Carlo approach is an effective
# method for summarizing a posterior distribution
# when simulated samples are available from 
# the exact posterior distribution. 

# For an illustration of the Monte Carlo method, 
# return to where we were interested in the
# proportion of heavy sleepers p at a college.

# With the use of a beta prior, the posterior
# distribution for p was beta(14.26, 23.19). 

# Suppose we are interested in the 
# posterior mean of p2. (This is the predictive
# probability that two students in a future
# sample will be heavy sleepers.) 

# We simulate 1000 draws from the beta
# posterior distribution.

# If {pj} represent the simulated sample, 
# the Monte Carlo estimate at this
# posterior mean will be the mean of 
# the {(pj)2}, and the simulated standard
# error is the standard deviation of the
# {(pj)2} divided by the square root of
# the simulation sample size.

p=rbeta(1000, 14.26, 23.19)
est=mean(p^2)
se=sd(p^2)/sqrt(1000)
c(est,se)

# The Monte Carlo estimate at E(p2|data) 
# is 0.150, with an associated simulation
# standard error of 0.002.

######################
# Rejection Sampling
######################

# library(LearnBayes)

# data(cancermortality)

# Previously, we were able to produce 
# simulated samples directly from the 
# posterior distribution since the 
# distributions were familiar functional forms. 

# So we could obtain Monte Carlo estimates
# of the posterior mean for any function 
# of the parameters of interest.

# But in other situations, such as the 
# beta-binomial example today, the
# posterior does not have a familiar form 
# and so we need to use an alternative
# algorithm for producing a simulated sample.

# A general-purpose algorithm for 
# simulating random draws from a given
# probability distribution is rejection sampling. 

# Suppose we wish to produce an independent
# sample from a posterior density g(??|y) 
# where the normalizing constant may not be known. 

# The first step in rejection sampling
# is to find another probability density 
# p(??) such that:
# . It is easy to simulate draws from p.
# . The density p resembles the posterior 
#   density of interest g in terms of
#   location and spread.
# . For all ?? and a constant c, g(??|y) ??? cp(??).

# Suppose we are able to find a density p 
# with these properties. Then one obtains
# draws from g using the following 
# accept/reject algorithm:
#  1. Independently simulate ?? from p and 
#     a uniform random variable U on the
#     unit interval.
# 2.  If U ??? g(??|y)/(cp(??)), then accept ?? 
#     as a draw from the density g; 
#     otherwise reject ??.
# 3.  Continue steps 1 and 2 of the algorithm
#     until one has collected a sufficient
#     number of "accepted" ??.

# We return to the beta-binomial example. 
# We want to find a proposal density of
# a simple functional form that, when 
# multiplied by an appropriate constant,
# covers the posterior density of interest. 

# One choice for p would be a multivariate
# t density with mean and scale matrix
# chosen to match the posterior density.

# We use a multivariate t density with
# location fit$mode, scale matrix
# 2 fit$var, and 4 dof.

# We write a new function betabinT()
# with two inputs, the parameter
# theta and a list datapar with
# components data, the data matrix,
# and par, a list with the parameters
# of the t proposal density (mean,
# scale matrix, and degrees of freedom)

fit=laplace(betabinexch,c(-7,6),cancermortality)

betabinT=function(theta,datapar)
{
  data=datapar$data
  tpar=datapar$par
  d=betabinexch(theta,data)-dmt(theta,mean=c(tpar$m),S=tpar$var,df=tpar$df,log=TRUE)
  return(d)
}

# We define parameters of t proposal density
# and the list datapar:
tpar=list(m=fit$mode,var=2*fit$var,df=4)
tpar
datapar=list(data=cancermortality,par=tpar)
datapar

# We run laplace() with the above
# function and an "intelligent"
# starting value
start=c(-6.9,12.4)
fit1=laplace(betabinT,start,datapar)
fit1$mode

# We find that the maximum value d 
# occurs at the value ?? = (???6.889, 12.422).

# The value of d is found by evaluating
# the function at the modal value.
betabinT(fit1$mode,datapar)

# We use rejectsampling() using
# constant value of d and simulate
# 10,000 draws from proposal density
theta=rejectsampling(betabinexch,tpar,-569.2813,10000,cancermortality)
dim(theta)
# theta has 2406 rows so acceptance
# rate is 2406/10000 = .24

# We plot simulated draws from rejection 
# sampling on contour plot of previous
# log posterior density plot. Most 
# of the draws are within the inner
# contour of the exact density
mycontour(betabinexch,c(-8,-4.5,3,16.5),cancermortality, xlab="logit eta",ylab="log K")
points(theta[,1],theta[,2])

