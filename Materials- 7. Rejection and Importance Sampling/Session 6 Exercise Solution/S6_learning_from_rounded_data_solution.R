###############################
# Exercise Solution: Inference
# about a Normal Population
###############################

# load LearnBayes package
library(LearnBayes)

# (a) Assuming observations represent a
# random sample from a normal population
# with mean mu and variance sigma^2 and
# the usual noninformative prior is placed
# on (mu, sigma^2), simulate a sample of
# 1000 draws from the joint posterior
# distribution

# sleeping times as vector y
y=c(9.0, 8.5, 7.0, 8.5, 6.0, 
    12.5, 6.0, 9.0, 8.5, 7.5,
    8.0, 6.0, 9.0, 8.0, 7.0, 
    10.0, 9.0, 7.5, 5.0, 6.5)

# sum of mean squared error
S=sum((y-mean(y))^2)
S

# sample size
n=length(y)
n

# simulate variance using inverse chi
# simulation 1000 times
sigma2=S/rchisq(1000, n - 1)
sigma2

# simulate mu with rnorm() 1000 times
mu = rnorm(1000, mean = mean(y), sd = sqrt(sigma2)/sqrt(n))
mu

# (b) Use the simulated sample to
# find 90% interval estimates for
# mean mu and standard deviation sigma

# 90% interval estimate for mu
int.est.mu=quantile(mu,c(.05, .95))
int.est.mu

# 90% interval estimate for sigma
int.est.sigma=quantile(sqrt(sigma2),c(0.05, .95))
int.est.sigma

# (c) Suppose we are interested in
# estimating upper quartile p75 of
# normal population. Using p75= mu + 0.674(sigma)
# find the posterior mean and posterior
# standard deviation of p75

# compute 1000 values of p75
p.75=mu+0.674*sqrt(sigma2)
p.75

# derive the mean
post.mean=mean(p.75)
post.mean

#derive the standard deviation
post.sd=sd(p.75)
post.sd
