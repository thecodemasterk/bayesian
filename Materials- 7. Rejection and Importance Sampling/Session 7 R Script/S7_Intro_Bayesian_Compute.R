########################################
# Introduction to Bayesian Computation #
########################################

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
# sample from a posterior density g(theta|y) 
# where the normalizing constant may not be known. 

# The first step in rejection sampling
# is to find another probability density 
# p(theta) such that:
# - It is easy to simulate draws from p.
# - The density p resembles the posterior 
#   density of interest g in terms of
#   location and spread.
# - For all theta and a constant c, 
#   g(theta|y) l.t.e.t. cp(theta).

# Suppose we are able to find a density p 
# with these properties. Then one obtains
# draws from g using the following 
# accept/reject algorithm:
#  1. Independently simulate theta from p and 
#     a uniform random variable U on the
#     unit interval.
#  2. If U l.t.e.t. g(theta|y)/(cp(theta)), then accept theta
#     as a draw from the density g; 
#     otherwise reject theta.
#  3. Continue steps 1 and 2 of the algorithm
#     until one has collected a sufficient
#     number of "accepted" theta's.

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

fit=laplace(betabinexch,
            c(-7,6),
            cancermortality)
fit

betabinT=function(theta,datapar)
{
  data=datapar$data
  tpar=datapar$par
  d=betabinexch(theta,data)-dmt(theta,
                                mean=c(tpar$m),
                                S=tpar$var,
                                df=tpar$df,
                                log=TRUE)
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
# occurs at the value theta = (-6.889, 12.422).

# The value of d is found by evaluating
# the function at the modal value.
betabinT(fit1$mode,datapar)

# We use rejectsampling() using
# constant value of d and simulate
# 10,000 draws from proposal density
theta=rejectsampling(betabinexch,tpar,
                     -569.2813,10000,
                     cancermortality)
dim(theta)
# theta has 2406 rows so acceptance
# rate is 2406/10000 = .24

# We plot simulated draws from rejection 
# sampling on contour plot of previous
# log posterior density plot. Most 
# of the draws are within the inner
# contour of the exact density
mycontour(betabinexch,c(-8,-4.5,3,16.5),
          cancermortality, xlab="logit eta",
          ylab="log K")
points(theta[,1],theta[,2])

###########################
### Importance Sampling ###
###########################

### see slides 6-9 for mathematical details ###

# As in rejection sampling, the issue in 
# designing a good importance sampling 
# estimate is finding a suitable sampling 
# density p. This density should be of a
# familiar functional form so simulated 
# draws are available. 

# The density should mimic the posterior 
# density g and have relatively flat tails 
# so that the weight function w(theta) is 
# bounded from above. One can monitor the 
# choice of p by inspecting the values of 
# the simulated weights w(thetaj). If there
# are no unusually large weights, then it 
# is likely that the weight function is bounded
# and the importance sampler is providing a 
# suitable estimate.

# library(LearnBayes)

data(cancermortality)

fit=laplace(betabinexch,c(-7,6),
            cancermortality)
fit

# We write betabinexch.cond() this posterior
# density conditional on the value of theta1
# -6.819. The function allows the input of
# the vector of values theta2 = log K. The
# function returns the value of the density
# (not log density).

betabinexch.cond=function (log.K, data)
{
  eta = exp(-6.818793)/(1 + exp(-6.818793))
  K = exp(log.K)
  y = data[, 1]; n = data[, 2]; N = length(y)
  logf=0*log.K
  for (j in 1:length(y))
    logf = logf + lbeta(K * eta + y[j], K * (1 - eta) + n[j] - y[j]) - lbeta(K * eta, K * (1 - eta))
  val = logf + log.K - 2 * log(1 + K)
  return(exp(val-max(val)))
}

# To compute the mean of logK for the cancer 
# mortality data, suppose we let the proposal
# density p be normal with mean 8 and standard
# deviation 2.

# In this R code, we use the integrate() 
# to find the normalizing constant of the 
# posterior density of logK. 

# Then, using the curve function, we display
# the conditional posterior density of logK 
# and the normal proposal density in the top 
# left graph. The top right graph displays the
# weight function, the ratio of the posterior
# density to the proposal density.

I=integrate(betabinexch.cond,2,16,cancermortality)
I
par(mfrow=c(2,2))
curve(betabinexch.cond(x,cancermortality)/I$value,
      from=3,to=16, ylab="Density", 
      xlab="log K",lwd=3, main="Densities")
curve(dnorm(x,8,2),add=TRUE)
legend("topright",legend=c("Exact","Normal"),lwd=c(3,1))
curve(betabinexch.cond(x,cancermortality)/I$value/dnorm(x,8,2),
      from=3,to=16, ylab="Weight",xlab="log K",
      main="Weight = g/p")

# Although the normal proposal density 
# resembles the posterior density with 
# respect to location and spread, the 
# posterior density has a flatter right 
# tail than the proposal and the weight 
# function is unbounded for large logK. 

# Suppose instead that we let the proposal
# density have the t functional form with 
# location 8, scale 2, and 2 degrees of 
# freedom. Using more R commands, the
# bottom graphs display the posterior and 
# proposal densities and the weight function. 

# Here the t proposal density has flatter
# tails than the posterior density and the 
# weight function is bounded. Here the t 
# functional form is a better proposal 
# for importance sampling.

curve(betabinexch.cond(x,cancermortality)/I$value,
      from=3,to=16, ylab="Density", xlab="log K",
      lwd=3, main="Densities")
curve(1/2*dt(x-8,df=2),add=TRUE)
legend("topright",legend=c("Exact","T(2)"),lwd=c(3,1))
curve(betabinexch.cond(x,cancermortality)/I$value/
        (1/2*dt(x-8,df=2)),from=3,to=16, 
      ylab="Weight",xlab="log K",
      main="Weight = g/p")

##############################
### Using a Multivariate t ###
### as a Proposal Density  ###
##############################

# For a posterior density of a vector of 
# real-valued parameters, a convenient
# choice of sampler p is a multivariate 
# t density. The R function impsampling
# implements importance sampling for an 
# arbitrary posterior density when p is 
# a t density. 

# There are five inputs to this function: 
# - logf() is the function defining the 
#   logarithm of the posterior, 
# - tpar() is a list of parameter values of
#   the t density, 
# - h() is a function defining the function 
#   h(theta) of interest, 
# - n is the size of the simulated sample, and 
# - data is the vector or list used in the 
# definition of logf(). 

# In the function impsampling(), the functions
# rmt() and dmt() from the mnormt library 
# are used to simulate and compute values 
# of the t density. 

# In this R code from impsampling(), we 
# simulate draws from the sampling density, 
# compute values of the log sampling density
# and the log posterior density at the 
# simulated draws, and compute the weights
# and importance sampler estimate.

# theta = rmt(n, mean = c(tpar$m), S = tpar$var, df = tpar$df)
# lf = matrix(0, c(dim(theta)[1], 1))
# lp = dmt(theta, mean = c(tpar$m), S = tpar$var, df = tpar$df,
#         log = TRUE)
# md = max(lf - lp)
# wt = exp(lf - lp - md)
# est = sum(wt * H)/sum(wt)

# Note that the value md is the maximum 
# value of the difference of logs of the
# posterior and proposal density - this value 
# is used in the computation of the weights
# to prevent possible overflow. 

# The output of impsampling() is a list with
# four components: 

# - est is the importance sampling estimate, 
# - se is the corresponding simulation standard
#   error, 
# - theta is a matrix of simulated draws from
#   the proposal density p, and 
# - wt is a vector of the corresponding
#   weights.

# To illustrate importance sampling, we 
# return to our beta-binomial example and
# consider the problem of estimating the 
# posterior mean of logK.

# The proposal density used in the development
# of a rejection algorithm seems to be a 
# good choice for importance sampling. 

# We choose a t density where the location
# is the posterior mode (found from laplace), 
# the scale matrix is twice the estimated 
# variance-covariance matrix, and the number
# of degrees of freedom is 4. This choice for
# p will resemble the posterior density and
# have flat tails that we hope will result 
# in bounded weights. 

# We define myfunc() to compute function h(). 
# Since we are interested in the posterior
# mean of logK, we define the function to
# be the second component of the vector theta.

tpar=list(m=fit$mode,var=2*fit$var,df=4)
tpar
myfunc=function(theta)
  return(theta[2])
s=impsampling(betabinexch,tpar,
              myfunc,10000,
              cancermortality)
s
cbind(s$est,s$se)

# We see from the output that the importance
# sampling estimate of the mean of logK 
# is 7.918 with an associated standard 
# error of 0.0184.

##############################################
# Section 5.10 Sampling Importance Resampling
##############################################

library(LearnBayes)
data(cancermortality)
fit=laplace(betabinexch,c(-7,6),cancermortality)

tpar=list(m=fit$mode,var=2*fit$var,df=4)

theta.s=sir(betabinexch,tpar,10000,cancermortality)

S=bayes.influence(theta.s,cancermortality)

plot(c(0,0,0),S$summary,type="b",
     lwd=3,xlim=c(-1,21),
     ylim=c(5,11), 
     xlab="Observation removed",
     ylab="log K")
for (i in 1:20)
  lines(c(i,i,i),S$summary.obs[i,],type="b")
