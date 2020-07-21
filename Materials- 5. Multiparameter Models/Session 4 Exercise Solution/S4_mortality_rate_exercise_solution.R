##########################################
# Exercise Solution: Learning about a
# Mortality Rate using a Mixture Prior
##########################################

# load LearnBayes package
library(LearnBayes)

# look up poisson.gamma.mix() function
?poisson.gamma.mix

# a) graph the prior density
curve(.5*dgamma(x,shape=1.5,rate=1000)+.5*dgamma(x,shape=7,rate=1000),
  from=0, to=.02)

# b) compute the posterior distribution of lambda
yobs=4; e=1767
data=list(y=yobs, t=e)
probs=c(.5,.5)
gammapar=rbind(c(1.5,1000),c(7,1000))
post=poisson.gamma.mix(probs,gammapar,data)
post

# c) show prior and posterior on same graph
# this is posterior
curve(.7597*dgamma(x,shape=5.5,rate=2767)+
      .2403*dgamma(x,shape=11,rate=2767),
       from=0, to=.02, col="red")
# this is prior
curve(.5*dgamma(x,shape=1.5,rate=1000)+.5*dgamma(x,shape=7,rate=1000),
       from=0, to=.02, add=TRUE)

# d) compute probability mortality rate 
#    exceeds .005
.7597*(1-pgamma(.005,shape=5.5,rate=2767))+
      .2403*(1-pgamma(.005,shape=11,rate=2767))

# e) based on mixing probabilities, were the data more
#    consistent with beliefs of first or second expert?

# First expert...the posterior is more consistent
# with curve of gamma() function of second expert:

# graph the prior density of first expert
curve(dgamma(x,shape=1.5,rate=1000),
      from=0, to=.02, col="blue")
# this is posterior
curve(.7597*dgamma(x,shape=5.5,rate=2767)+
        .2403*dgamma(x,shape=11,rate=2767),
      from=0, to=.02, col="red", add=TRUE)
# graph the prior density of second expert
curve(dgamma(x,shape=7,rate=1000),
      from=0, to=.02, col="green", add=TRUE)
