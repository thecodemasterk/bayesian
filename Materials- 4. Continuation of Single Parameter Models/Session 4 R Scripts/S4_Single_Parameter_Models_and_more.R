#################################
# Single Parameter Models
#################################

####################################
# Mixtures of Conjugate Priors
####################################

library(LearnBayes)

# So far we have shown how with the binomial, 
# Poisson, and normal sampling models, the 
# use of a conjugate prior where the prior and 
# posterior distributions have same functional form.

# You can extend the family of conjugate priors 
# by using discrete mixtures. We illustrate the 
# use of  a mixture of beta densities to learn
# about the probability that a biased coin lands
# heads.

### SEE SLIDES #2 and #4

# Mixture of beta densities prior distribution that
# reflects belief that a coin is biased.
curve(.5*dbeta(x, 6, 14) + .5*dbeta(x, 14, 6), 
      from=0, to=1, xlab="P", ylab="Density")

# The R function binomial.beta.mix computes the 
# posterior distribution when the proportion p has 
# a mixture of betas prior distribution. The inputs
# to this function are probs, the vector of mixing 
# probabilities; betapar, a matrix of beta shape 
# parameters where each row corresponds to a component 
# of the prior; and data, the vector of the number of 
# successes and number of failures in the sample.

# The output of the function is a list with two components -
# probs is a vector of posterior mixing probabilities and
# betapar is a matrix containing the shape parameters of 
# the updated beta posterior densities.

probs=c(.5,.5)
beta.par1=c(6, 14)
beta.par2=c(14, 6)
betapar=rbind(beta.par1, beta.par2)
data=c(7,3)
post=binomial.beta.mix(probs,betapar,data)
post

## SEE SLIDE #5

# Can open up separate graphics window:
windows()

# Prior and posterior densities of a proportion
# for a biased coin example:

# Here is posterior curve based on observed data:
curve(post$probs[1]*dbeta(x,13,17)+post$probs[2]*dbeta(x,21,9),
      from=0, to=1, lwd=3, xlab="P", ylab="DENSITY")

# Here is the prior curve again:
curve(.5*dbeta(x,6,12)+.5*dbeta(x,12,6),0,1,add=TRUE)
# legend of what is what:
legend("topleft",legend=c("Prior","Posterior"),lwd=c(1,3))

##############################################
# A Bayesian Test of the Fairness of a Coin
##############################################

### SEE SLIDE #7

pbinom(5, 20, 0.5)

# The p-value here is 2 x .021 = .042. Since
# this value is smaller than the common significance
# level of .05, you would decide to reject
# the hypothesis H and conclude that the coin
# is not fair.

### SEE SLIDES #8, 9 and 10

n = 20
y = 5
a = 10
p = 0.5
m1 = dbinom(y, n, p) * dbeta(p, a, a)/dbeta(p, a + y, a + n - y)
lambda = dbinom(y, n, p)/(dbinom(y, n, p) + m1)
lambda

# We get the surprising result that the 
# posterior probability of the hypothesis
# of fairness H is .28, which is less evidence 
# against fairness than is implied by the p-value
# calculation above.

# The function pbetat in the LearnBayes package 
# performs a test of a binomial proportion.

# The inputs to the function are the value of 
# p to be tested, the prior probability of that 
# value, a vector of parameters of the beta prior
# when the hypothesis is not true, and a vector of 
# numbers of successes and failures. 

# In this example, the syntax would be

pbetat(p,.5,c(a,a),c(y,n-y))

# The output variable post is the posterior 
# probability that p = .5, which agrees with the
# calculation. The output variable bf is the Bayes
# factor in support of the null hypothesis, which
# we discuss later.

# Since the choice of the prior parameter a = 10
# in this analysis seems arbitrary, it is natural
# to ask about the sensitivity of this posterior
# calculation to the choice of this parameter. 

# To answer this question, we write a short
# function prob.fair() that computes the 
# probability of a fair coin as a function
# of log a.
prob.fair=function(log.a)
{
  a = exp(log.a)
  m2 = dbinom(y, n, p) * dbeta(p, a, a)/
    dbeta(p, a + y, a + n - y)
  dbinom(y, n, p)/(dbinom(y, n, p) + m2)
}

# Using curve() function we graph the posterior
# probability for a range of values of log a.
n = 20; y = 5; p = 0.5

# Posterior probability that a coin is fair graphed
# against values of the prior parameter log a:
curve(prob.fair(x), from=-4, to=5, xlab="log a", 
      ylab="Prob(coin is fair)", lwd=2)

# We see from this graph that the probability 
# of fairness appears to be greater than .2 
# for all choices of a. We must remember that 
# the p-value is not interpretable as a probability
# of fairness, although it is sometimes mistakenly
# viewed as this probability. 

# But the Bayesian posterior probability of .2 
# is larger than the p-value calculation of .042, 
# suggesting that the p-value is overstating
# the evidence against the hypothesis that the 
# coin is fair.

### SEE SLIDE #12

n=20
y=5
a=10
p=.5
m2=0
for (k in 0:y)
  m2=m2+dbinom(k,n,p)*dbeta(p,a,a)/dbeta(p,a+k,a+n-k)
lambda=pbinom(y,n,p)/(pbinom(y,n,p)+m2)
lambda

# Note that the posterior probability of fairness 
# is .218 based on the data "5 heads or fewer." 

# This posterior probability is smaller than the 
# value of .280 found earlier based on y = 5. 

# This is a reasonable result since observing
# 5 heads or fewer" is stronger evidence against
# fairness than the result "5 heads."

#######################################
# In Class Exercise
# Bayesian robustness
#######################################

# Suppose you are about to flip a coin that you
# believe is fair. If p denotes the probability
# of flipping a head, then your "best guess" 
# at p is 0.5. Moreover, you believe that it
# is highly likely that the coin is close to fair,
# which you quantify P(.44 < p < .56) = 0.9.

# Consider the following two priors for p:

# P1 is a non-mixture prior
# P1: p distributed as beta(100,100)

# P2 is a mixture prior
# P2: p distributed according to the
#     mixture prior:
# g(p) = .9fB(p;500,500)+.1fB(p;1,1)
# where fB(p;a,b) is the beta density 
# with parameters a and b.

# a) Simulate 1000 values from each prior density
#    P1 and P2. By summarizing the simulated samples,
#    show that both priors match the given prior
#    beliefs about the coin flipping probability p.

# simulate 1000 draws from the beta(100, 100) prior
p.prior1=rbeta(1000,100,100)
p.prior1

# simulate 1000 draws from the mixture prior
# .9 beta(500, 500) + .1 beta(1, 1)
u=runif(1000)
p.prior2=(u<.9)*rbeta(1000,500,500)+(u>=.9)*rbeta(1000,1,1)
p.prior2

# confirm that each prior matches belief
# P(.44 < p < .56) = .9
mean((p.prior2>.44)&(p.prior2<.56))
mean((p.prior1>.44)&(p.prior1<.56))

# b) Suppose you flip the coin 100 times and obtain
#    45 heads. Simulate 1000 values from the posteriors
#    from priors P1 and P2, and compute 90% probability
#    intervals.

# observe 45 heads and 55 tails 

# simulate 1000 draws from beta posterior (beta prior)
s=45; f=55
p.post1=rbeta(10000,100+s,100+f)
p.post1

# compute 90% interval estimate
interval.1=quantile(p.post1,c(.05,.95))
interval.1

# c) Suppose that you only observe 30 heads out
#    of 100 flips. Again simulate 1000 values
#    from the two posteriors and compute 90%
#    probability intervals.

# for the mixture prior
# set up a grid of values of p
p=seq(.001,.999,by=.001);p

# compute the posterior density of this grid
# using prior x like recipe

# Mixture prior:
prior=.9*dbeta(p,500,500)+.1*dbeta(p,1,1)
prior

# Compute post from this mixture prior
post=prior*dbeta(p,s+1,f+1)
post=post/sum(post)
post

# for the mixture prior
# set up a grid of values of p
p=seq(.001,.999,by=.001)
p

# simulates from the discrete distribution on the grid
p.post2=sample(p,size=10000,prob=post,replace=TRUE)
p.post2

# Confidence interval
interval.2=quantile(p.post2,c(.05,.95))
interval.2

# compare interval.1 and interval.2 to check 
# the robustness of the inference with respect to the prior

# repeats the computations for the 
# data 30 heads and 70 tails
s=30; f=70
p.post1=rbeta(1000,100+s,100+f)
p.post1
interval.1=quantile(p.post1,c(.05,.95))
interval.1

# Set up grid of values
p=seq(.001,.999,by=.001)
# Use mixture prior
prior=.9*dbeta(p,500,500)+.1*dbeta(p,1,1)
# Compute posterior
post=prior*dbeta(p,s+1,f+1)
post=post/sum(post)
post
p.post2=sample(p,size=1000,prob=post,replace=TRUE)
p.post2
interval.2=quantile(p.post2,c(.05,.95))
interval.2

# d) Looking at results from (b) and (c), comment
#    on the robustness of the inference with respect
#    to the choice of prior density in each case.
