###########################################
###  Introduction to Bayesian Thinking  ###
###########################################

# Here we introduce basic bayesian inferential
# approach with problem of learning about
# a population proportion.

# Before taking data, one has beliefs about
# the value of the proportion and one
# models his or her beliefs in terms
# of a prior distribution.

# This 'prior' can take different functional
# forms.

# After data is observed, one updates ones
# beliefs about the proportion by computing
# the posterior distribution.

# One may also want to predict likely outcomes
# of a new sample taken from the population.

### LEARNING ABOUT PROPORTION
### OF HEAVY SLEEPERS

# Sleeping habits of American college students.
# Doctors recommend 8 hours. What proportion
# get at least 8 hours.

# Population consists of all American
# college students and p represents
# the proportion who sleep (on a typical
# night during week) an average of 8 hours.

# We want to know about the location of
# p which is unknown. In Bayesian terms,
# a person's beliefs about the uncertainty
# in this proportion can be represented
# by a probability distribution placed
# on this parameter. The distribution
# reflects the person's subjective prior
# opinion about plausible values of p.

# Person reads in newspaper that most
# students only get 6 hours. Second article
# says based on sample of 100, about 70%
# reported getting 5-6 hours sleep on
# weekdays, 28% reported 7-8 hours, and
# only 2% getting "healthy nine hours".

# So: p, proportion that sleep < 8 hours
# is likely smaller than 50% or 0.5

# Best guess is that value of p is 0.3

# But is plausible that p could be
# anywhere between 0 and 0.5.

# We take a sample of 27 students, 11
# say they had at least 8 hours the
# previous night. Now we want to predict
# proportion p if we draw a new sample
# of 20.

### SEE SLIDE

####################################
### Using a Discrete Prior
####################################

# One approach to get prior for p is to
# write down posible proportion values
# and to assign weights to the values.

# .05, .15, .25, .35, .45, .55, .65, .75, .85, .95
# with weights:
# 1, 5.2, 8, 7.2, 4.6, 2.1, 0.7, 0.1, 0, 0

# convert to prior probabilities by
# dividing each weight by the sum.

library(LearnBayes)

# proportion p vector:
p = seq(0.05, 0.95, by = 0.1); p

# corresponding weights
prior = c(1, 5.2, 8, 7.2, 4.6, 2.1, 0.7, 0.1, 0, 0); prior

# prior probabilities of each p proportion
prior = prior/sum(prior); prior

# here we graph this discrete prior 
# distribution for porportions p:
plot(p, prior, type = "h", ylab="Prior Probability")

# we determine likelihood function for 
# 11 of 27 students that sleep a sufficient
# number of hours (SEE NEXT SLIDE AFTER PLOT)

## PLEASE SEE WIKIPEDIA FOR BETA DISTRIBUTION

# beta density parameters:
data = c(11, 16); data

# posterior densities
post = pdisc(p, prior, data)
post

# show p, priors and posteriors
round(cbind(p, prior, post),2)

# Here we make a fancy plot
library(lattice)

# add a new column up front,
# fill with "prior"
PRIOR=data.frame("prior",p,prior);PRIOR

# same with "posterior" for previous
# posterior numbers
POST=data.frame("posterior",p,post);POST

# names for prior column is "Type"
names(PRIOR)=c("Type","P","Probability");PRIOR

# names for posterior column is "Type"
names(POST)=c("Type","P","Probability");POST

# rowbind them (stack them up)
data=rbind(PRIOR,POST);data

# Do not need windows() command
# it just opens a new window for
# the xyplot
# windows()
xyplot(Probability~P|Type,
       data=data,layout=c(1,2),
       type="h",lwd=3,col="black")

# Note that most of posterior probability
# is concentrated on values of p = .35 
# and p = .45. If we combine the probabilities
# for the three most likely values, we can
# conclude that posterior prob is that p
# falls in the set {.25, .35, .45} = 0.94

#############################
### Using a Beta Prior
#############################

# A proportion is a continuous parameter,
# so an alternative approach is to construct
# a density g(p) on the interval (0, 1) that
# represents the person's initial beliefs.

# Suppose belief is that proportion is equally
# likely to be smaller or larger than 0.3.
# Also, she is 90% confident that p < 0.5.

#### SEE NEXT SLIDE

library(LearnBayes)

quantile2=list(p=.9,x=.5)
quantile2
quantile1=list(p=.5,x=.3)
quantile1
ab=beta.select(quantile1,quantile2)
ab

### SEE NEXT SLIDE

a = ab[1];a
b = ab[2];b
s = 11
f = 16

# We display the three densities: the prior
# density (g(p)), the likelihood function
# L(p), and the posterior density (g(p|data)
# for learning about the proportion p. 

# We see that (in this case) the posterior 
# density is acompromise between the
# initial prior beliefs and the info
# provided by the data.
curve(dbeta(x,a+s,b+f), from=0, to=1, 
      xlab="p",ylab="Density",lty=1,lwd=4)
curve(dbeta(x,s+1,f+1),add=TRUE,lty=2,lwd=4)
curve(dbeta(x,a,b),add=TRUE,lty=3,lwd=4)
legend(.7,4,c("Prior","Likelihood","Posterior"),
       lty=c(3,2,1),lwd=c(3,3,3))

# So there are different ways of sum-
# marizing the beta posterior distribution
# to make inferences about the proportion
# of heavy sleepers p.

# The beta cum density function (cdf) and
# inverse cdf functions pbeta and qbeta are
# useful in computing probabilities and 
# creating interval estimates for p.

# Is it likely that the proportion of heavy
# sleepers > 0.5? We compute the posterior
# prob P(p >= .5|data)

1 - pbeta(0.5, a + s, b + f)

# It is small so it is unlikely that more 
# than half of the students are heavy sleepers.

# We can find a 90% interval estimate for 
# p by computing the 5th and 95th percen-
# tiles of the beta density:

qbeta(c(0.05, 0.95), a + s, b + f)

# So a 90% posterior credible interval 
# for the proportion p is (0.256, 0.513)

# These summaries are exact because are
# based on R functions for the beta
# posterior density.

### SIMULATION AS ALTERNATIVE
# Can also use simulation to summarize
# a posterior density....we simulate a
# large number of values from the beta
# posterior bensity and summarize the
# simulated output.

# We use pseudo-random beta generator
# rbeta() to simulate 1000 random pro-
# portion values from the beta (a+s, b+f)
# posterior with:
ps = rbeta(1000, a + s, b + f)
ps

# and then display the proportion as
# a histogram of the simulated values.

# We estimate the probability that the
# proportion is larger than 0.5 using
# the proportion of simulated values
# in the range of the simulated sample
# from the beta posterior distribution
# of p:
hist(ps,xlab="p")
sum(ps >= 0.5)/1000

# Note these summaries of the posterior
# density for p based on simulation are
# very close to the exact values based
# on the calculations from the beta
# distribution:

# Here we get 90% interval estimate from
# the 5th and 95th sample quantiles:
quantile(ps, c(0.05, 0.95))

#####################################
### Using a Histogram Prior
#####################################

library(LearnBayes)

# There are computational advantages to
# using a beta prior, can also perform
# posterior computations for any choice
# of prior.

### BRUTE FORCE TECHNIQUE: NEXT SLIDE

# We illustrate for a 'histogram' prior
# that may better reflect the person's prior
# opinion about proportion p.

# They state prior beliefs about propor-
# tion of heavy sleepers by dividing the
# range of p into ten subintervals:
# (0, .1), (.1, .2) . . . (.9, 1) and
# then assigning probabilities to these
# intervals.

# Person assigns weights: 1, 5.2, 8, 7.2,
# 4.6, 2.1, 0.7, 0.1, 0, 0

# We view this as a continuous version
# of the discrete prior used earlier

# So histogram prior has
# midpoints of intervals:
midpt = seq(0.05, 0.95, by = 0.1)
midpt

# the associate prior weights:
prior = c(1, 5.2, 8, 7.2, 4.6,
          2.1, 0.7, 0.1, 0, 0)
prior

# We convert prior weights to probs
# by dividing each weight by the sum:
prior = prior/sum(prior); prior

# We graph with curve() and with
# histprior() from LearnBayes package
# to yield a histogram prior for pro-
# portion p:
curve(histprior(x,midpt,prior), from=0, to=1,
      ylab="Prior density",ylim=c(0,.3))

# from previously:
s = 11
f = 16

# We compute posterior density by multiplying
# the histogram prior by the likelihood function
# (likelihood function for a binomial
# density is given by beta(s+1,f+1)) density
# which we can call with dbeta()

# We again graph with curve() and with
# histprior() from LearnBayes package. So
# here we have the posterior density for
# proportion p using a histogram prior:
curve(histprior(x,midpt,prior) * dbeta(x,s+1,f+1), 
      from=0, to=1, ylab="Posterior density")

# To obtain a simulated sample from posterior
# density we construct an equally spaced
# grid of values of proportion p and
# compute the product of the prior and
# the likelihood on this grid. Then we 
# convert the products on the grid to
# probabilities:
p = seq(0, 1, length=500)
p
post = histprior(p, midpt, prior)*dbeta(p, s+1, f+1)
post
post = post/sum(post)
post

# Here we take a sample with replacement
# from the grid:
?sample
ps = sample(p, replace = TRUE, prob = post)
ps

# Here we create a histogram of the simulated
# draws from the posterior distribution 
# of p with the use of a histogram prior.
hist(ps, xlab="p", main="")

########################
### Prediction
########################

# Suppose someone wants to predict the
# number of heavy sleepers y-tilde in a
# future sample of m=20 students.

### SEE NEXT TWO SLIDES (15 and 16)

# pdiscp() computes the predictive probabilities
# when p is given a discrete distribution. 
# Like before, p is a vector of proportion
# values and prior is a vector of current
# probabilities.

# Other arguments are future sample size m,
# vector ys of number of successes of interest.

# Output is a vector of the corresponding
# predictive probabilities

p=seq(0.05, 0.95, by=.1);p
prior = c(1, 5.2, 8, 7.2, 4.6, 2.1, 0.7, 0.1, 0, 0)
prior
prior=prior/sum(prior)
prior
m=20; ys=0:20; m; ys
pred=pdiscp(p, prior, m, ys); pred
# predictive probabilities: we see
# most likely number of successes in
# future sample is y-tilde = 5 and
# y-tilde = 6
round(cbind(0:20,pred),3)

### SEE NEXT SLIDE (#17)

# We illustrate this computation using
# the beta(3.26, 7.19) used before to
# reflect the person's beliefs about
# the proportion of heavy sleepers:
ab=c(3.26, 7.19);ab
m=20; ys=0:20; m; ys
# This is the predictive density:
pred=pbetap(ab, m, ys); pred

# Note one can always compute a predictive
# density for any prior with simulation.

# To get y-tilde we first simulate p-star
# from g(p) and then simulate y-tilde from
# the binomial distribution

# Here we do it for beta(3.26, 7.19) prior.
# We simulate 1000 draws from the prior
# and put the simulated values in p:
p=rbeta(1000,3.26, 7.19);p

# Now we simulate values for y-tilde for
# these random p's using the rbinom()
# function
y = rbinom(1000, 20, p);y

# We tabulate the distinct values:
table(y)

# We save these frequencies and then convert
# to probabilities by dividing each freq
# by sum and use plot() to graph predictive
# distribution
freq=table(y)
freq
ys=as.integer(names(freq))
ys
predprob=freq/sum(freq)
predprob

# Predictive probabilities of the number
# of sleepers y-tilde in a future sample
# of size 20 when proportion is assigned
# a beta (3.26, 7.19) prior:
plot(ys,predprob,type="h",xlab="y",
     ylab="Predictive Probability")

# Now we want to summarize this discrete
# predictive distribution by an interval
# that covers at least 90% of the prob
# with discint()

# Vector of ys contain values of y-tilde 
# and predprob contains associated probs
# from table output.

# Matrix dist contain this prob dist
dist=cbind(ys,predprob)
dist

# discint() has inputs matrix dist and a
# given coverage prob covprob:
covprob=.9

# Output is a list where component 'set'
# gives the credible set and 'prob' gives
# the exact coverage probability

# So prob y-tilde falls in interval
# {1,2,3,4,5,6,7,8,9,10,11} is 91.1%
discint(dist,covprob)

# OR: if y-tilde/20 is the proportion of
# heavy sleepers in the future sample, the
# probability that this sample proportion
# is in this interval is 91.1%. Note that
# this interval is much wider than 91.1%
# probability interval for the population
# proportion p.

# When predicting a future sample proportion
# there are two sources of uncertainty:

# uncertainty in the value of p

# AND

# binomial uncertainty in the value of y-tilde

# so the predictive interval is longer
# than the population interval
