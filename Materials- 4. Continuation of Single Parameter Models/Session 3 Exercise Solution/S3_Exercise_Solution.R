############################################
### Exercise Session 3: How Many Taxis ? ###
############################################

# a) Use R to compute the posterior
# probabilities of N on a grid of values.

# Observed taxi numbers:
data=c(43, 24, 100, 35, 85); data

# Highest N and number of taxis observed:
yn = max(data); n = length(data)
yn; n

# Belief that there cannot be more
# than 200 taxis
B=200; B

# grid of values from highest observed
# to highest possible number 200
N=yn:B; N

# b) Compute the posterior mean and 
# posterior standard deviation of N:

# posteriors:
post=1/N^n; post
hist(post)

# posterior as probabilities
post=post/sum(post); post
hist(post)

# posterior mean
m=sum(N*post);m

# posterior standard deviation
s=sqrt(sum((N-m)^2*post));s

# Find the probability that there 
# are more than 150 taxis in the city.

# mean of N being greater than
# 150 times the posterior probabilities
post.prob=mean((N>150)*post)
post.prob
