###################################
###  Bayesian Model Comparison  ###
###################################

# Comparing hypotheses

# We use the notion where we judge two hypotheses
# concerning a parameter using prior and
# posterior odds ratios to derive a Bayes
# factor which is the ratio of the posterior
# odds to the prior odds 

# BF = posterior odds / prior odds

# A One-Sided Test of a Normal Mean

# Example where one is interested in determining
# his true weight from a variable bathroom scale.

# We assume measurements are normally distributed
# with mean mu and standard deviation sigma.

# The person weighs himself 10 times at 182,
# 172, 173, 176, 176, 180, 173, 174, 179, and
# 175.

# We assume accuracy of scale has sigma (standard
# deviation) of 3 pounds.

# If mu is person's true weight, let's say he
# wants to assess if his true weight is more
# than 175 pounds.

# So are testing the two hypotheses:
# (Null) Ho: mu =< 175
# (Alt.) H1: mu > 175

# Person assumes prior mu is 170 and sd is 5
# or mu is distributed as N(170, 5)

# Prior odds of null is 
# P(mu =< 175)/P(mu > 175)

# load library
library(LearnBayes)

# prior mean is 170, 
# prior variance is 5^2=25
pmean=170; pvar=25 

# We compute prior odds (prob) that weight is
# less than or equal to 175 from N(170,5) 
# density w/ pnorm(). 
?pnorm
pmean
pvar
probH=pnorm(175,pmean,sqrt(pvar))
# is 84% based on prior knowledge
probH
# prior probability of alternative 
# (that weight is greater than 175)
# is 1 - probH
probA=1-probH
# or 16%
probA
# So the prior odds ratio of
# null compared to alternative is
prior.odds=probH/probA 
prior.odds 

# So a priori, the null hypothesis is 5 times
# more likely than the alternative hypothesis

# We enter the 10 weight measurements
weights=c(182, 172, 173, 176, 
          176, 180, 173, 174, 
          179, 175)
# is a vector with 10 elements
weights
# We compute the sample mean
ybar=mean(weights)
# is 176 pounds
ybar
# and associated sample variance
sigma2=3^2/length(weights)
# is 0.9 pounds
sigma2

# Use normal density/normal prior
# updating formula and so we estimate
# posterior precision (inverse of 
# variance) of mu as sum of precisions of
# the observed data:
1/sigma2
# and the prior:
1/pvar
# posterior precision is:
post.precision=1/sigma2+1/pvar
post.precision
# posterior variance is inverse
post.var=1/post.precision
post.var

# posterior mean of mu is weighted average
# of the sample mean and the prior mean,
# where weights are proportional to the
# respective precisions
post.mean=((ybar/sigma2)+(pmean/pvar))/post.precision 
post.mean
# so posterior mean and sd is
c(post.mean,sqrt(post.var)) 

# So the posterior density of mu is
# N(175.79, 0.93)

# Using this normal posterior density,
# we calculate odds of null hypothesis
post.odds=pnorm(175,post.mean,sqrt(post.var))/(1-pnorm(175,post.mean,sqrt(post.var))) 
post.odds 

# So Bayes factor in support of null
# hypothesis (weight =< 175)  is
BF = post.odds/prior.odds 
BF 

# From prior probabilities and Bayes factor,
# we compute posterior probability of null
# hypothesis weight =< 175 as
postH=probH*BF/(probH*BF+probA) 
postH 

# So it is unlikely that his weight is
# at most 175 pounds

# If we just used the Frequentist approach
# to compute the p-value for same evidence
z=sqrt(length(weights))*(mean(weights)-175)/3
# z is standard normal variate
z
#p-value is 1 - pnorm(z) or 15%
1-pnorm(z) 

# Here we repeat this Bayesian analysis
# using a very flat prior with mean of 170
# and standard deviation of 1000

# same weights
weights=c(182, 172, 173, 176, 
          176, 180, 173, 174, 
          179, 175)
# data for input into mnormt.onesided()
# LearnBayes function
?mnormt.onesided
#data
data=c(mean(weights),length(weights),3)
data
# this is mean and sd
prior.par=c(170,1000)
prior.par
# use mnormt.onesided()
mnormt.onesided(m0=175,
                normpar=prior.par,
                data=data)

# Note that probability of null hypothesis postH
# is approximately equal to p-value. This is
# the general result: that a Bayesian
# probability of a hypothesis is equal to
# the p-value for one-sided testing problems
# when a vague (too broad) prior distribution
# is placed on the parameter.

###############################################
# A Two-Sided Test of a Normal Mean

# Suppose person from previous example knows
# that his weight last year was 170 pounds
# and he wonders if he still weighs 170 this
# year so
# Ho: true current weight mu = 170
# H1: weight is now larger or smaller than 170

# These assumptions impose a unique construction
# on the prior, the person basically assigns
# the statement that mu = 170 with a P(.5)

# What about mu if Ho is not true?

# If weight did change from last year, he may
# think that mu is closer to last year's weight
# than farther from it.

# So can model mu as N(170, unknown sd).

# But how do we arrive at the standard deviation?
# One way is to think of possible alternative
# values for the standard deviation of this
# weight change from last year is to think
# of possible alternative values for mu
# and then solve for this standard
# deviation by setting the 95% range of the
# normal distribution to this range.

# same weights (10 obs with pop sd of 3)
weights=c(182, 172, 173, 176, 
          176, 180, 173, 174, 
          179, 175)
# input for mnormt.twosided()
data=c(mean(weights),length(weights),3)
data
# t is possible values for the sd
# of the weight change
t=c(.5,1,2,4,8)
t
?mnormt.twosided
mnormt.twosided(170,.5,t,data)

# For each value of the likely standard 
# deviation t we get the Bayes factor in
# support of the hypothesis that mu takes
# on the specific value and the posterior
# probability that hypothesis H is true
# (weight is 170).

# For example, if t=2 so are using N(170,2)
# density to reflect alternative values for
# weight mu, then the Bayes factor in support
# of the hypothesis mu = 170 is equal to 3rd one
# 0.0000002, which is much smaller than 
# person's prior probability of 0.5. He
# should conclude that his current weight
# is not 170.

#############################################
# Is Baseball Hitter Really Streaky?

## Patterns of Dependence in a Sequence

# Sports fans love athletes who exhibit 
# patterns of extreme performance. 
# Records are kept of extreme winning 
# or losing patterns, such a getting a 
# base hit in baseball....
# ... a "hitting streak"

# In the 2006 baseball season, Chase Utley 
# of the Philadelphia Phillies had a hitting
# streak of 35 games, one of the best
# in baseball history.

# But how "significant" was this streak? 
# Utley was a good hitter anyway, so it 
# might be expected.

# As another example, long runs of heads
# and tails can be observed flipping a coin.

# We investigate with Monte Carlo simulation:

## Writing a Function to Compute Streaks

# We represent hitting success or failure 
# with 0's and 1's

y = c(1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1)
y

# What are the lengths of the hitting streaks?

# Add 0's to front and end to bracket it
# 'where' is vector of TRUE's (if 0) 
# and FALSE's (if 1)
where = (c(0, y, 0) == 0)
where
n = length(y)

# are 15 elements in y:
n
# loc.zeros records locations in sequence
# where 0's occur at 0  8  9  11  16
loc.zeros = (0:(n+1))[where]
loc.zeros

# So no base hit in games 0, 8, 9, 11 and 16

# We compute length of streaks:
streak.lengths = diff(loc.zeros) - 1
streak.lengths
# We take out the zeros:
streak.lengths = streak.lengths[streak.lengths > 0]
streak.lengths

# [1] 7 1 4

# Write function longest.streak by taking the max()
# out of the above vector:
longest.streak=function(y){
  where = c(0, y, 0) == 0
  n = length(y)
  loc.zeros = (0:(n+1))[where]
  streak.lengths = diff(loc.zeros) - 1
  streak.lengths = streak.lengths[streak.lengths > 0]
  return(max(streak.lengths))
}

# Make sure "utley2006.txt" is in path of default directory
# file.exists("utley2006.txt")

# read it in
dat = read.table(file.choose(), 
                 header=TRUE, sep="\t") # tab sep
head(dat)
# are counting 2's, 3's etc as 1's:
utley = as.numeric(dat$H > 0) # Make sure H is numeric
utley
# apply longest.streak() function to this vector
longest.streak(utley)

# So Utley's longest streak in 2006 was 35 games.

## Writing a Function to Simulate Hitting Data

# We apply Monte Carlo to understand "significance"
# Utley played in 160 games and hit in 116 of them

# What if 116 "hit" games were randomly distributed?
# The what would be length of longest sequence?

# Write function random.streak() with binary sequence y as
# argument
random.streak=function(y){
  # first it randomly permutes y and stores in mixed.up.y
  mixed.up.y = sample(y)
  # then find longest streak of 1's in vector mixed.up.y
  longest.streak(mixed.up.y)
}

# replicate random.streak 100000, store in L
L = replicate(100000, random.streak(utley))
L
# tabulate values in L and plot
plot(table(L))
# superimpose line for 2006 season:
abline(v=35, lwd=3)
# add text label in plot at 38, 10000
text(38, 10000, "Utley")

# We see long hitting streaks from 8 to 18 are quite common
# With a hitting streak of 11 most likely

# Utley's 35 streak is in far right tail, indicating it is 
# very rare. Can measure "how extreme" by estimating probability
# that a random sequence would be like this (35 or more)

# Probability....only 7 out of 100000 = 0.00007...VERY RARE !!

# Streakiness in Baseball for Bayesian

# Again, we look at measuring the success of
# a hitter as the batting average or the
# proportion of base hits. A hitter may have
# both hot and cold streaks during a season.

# Interesting question is: What do these
# 'streaky' data say about the ability of 
# a player to be streaky?

# Are two possible outcomes for each "at-bat":
# a 'hit' (success) or an 'out' (failure).

# Suppose we divide all 'at-bats' in a season
# into N periods. We let pi denote the prob
# that a player gets a hit in a single at-bat
# during the ith period, i = 1, . . . , N.

# If a player is consistent, or non-streaky,
# then the probability of a hit stays constant
# across all periods. This is non-streaky model
# Mo: p1 = ... = pN = p. To specify, we assign
# the common probability value p a uniform
# prior on (0, 1).

# On the other hand, if a player is truly
# 'streaky' then prob of a hit pi will change
# across the season.

# We can model this variation in the probs
# by assuming that p1, ..., pN are a random
# sample from a beta density g, where eta is
# the mean and K is a precision parameter.

# The we can represent the streaky model by
# Mk: pi, ..., pN as iid beta(K(eta), K(1-eta)).

# For non-streaky model we place a uniform
# prior on mean parameter eta, reflecting
# little knowledge about location of random
# effects distribution.

# We compare the models Mo and Mk. We need
# to compute the associated marginal densities.

# Under model Mo, numbers of hits y1, ..., yN
# are independent, where yi is binomial(ni,p)
# with assumption that p is uniform(0,1), we
# obtain a different marginal density

# The Bayes factor in support of "streaky" model
# Hk compared with "nonstreaky" model Ho is
# Bk = mk(y)/mo(y)

# We use laplace() to compute integral in the
# Bayes factor Bk.

# Complete hitting data for Jeter. He played
# in 154 games that season. We group the data
# into 5 game intervals to try to detect 
# streakiness
data(jeter2004) 
attach(jeter2004)
# gather up hits and at-bats
data=cbind(H,AB) 
data
# group data into periods of 5 games
data1=regroup(data,5)
# data1 contains grouped hitting data (yi,ni),
# i = 1, ..., 30i, where yi is number of hits
# by Jeter in ni at-bats in ith interval of
# games
data1
# We compute Bayes factor for sequence of values
# of log K using laplace() and definition of
# log integral defined in bfexch().

# Wrapper function that computes the log
# Bayes factor for a single value of log K.
log.marg=function(logK) 
  laplace(bfexch,0,list(data=data1,K=exp(logK)))$int 
# Vector log.K contains values log(K) =
# 2, 3, 4, 5, and 6.
log.K=seq(2,6) 
K=exp(log.K) 
# By using sapply() we store corresponding 
# values of the log Bayes factor log Bk in
# variable log.BF
log.BF=sapply(log.K,log.marg) 
BF=exp(log.BF)
# We display values of log K, values of K,
# values of log Bk, and values of Bayes
# factor Bk (BF column)
round(data.frame(log.K,K,log.BF,BF),2)

# We see that value log K = 4 is most supported
# by the data with a corresponding Bayes factor
# of Bk = 2.51.

# This particular streaky model is approximately
# two and one half times as likely as the
# consistent model. This indicates that Jeter
# did indeed display some true streakiness
# in his hitting behavior for this season.