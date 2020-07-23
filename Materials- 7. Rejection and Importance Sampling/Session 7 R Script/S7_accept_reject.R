####################################################
#####                                          #####
#####       RANDOM  VARIABLE  GENERATION       #####
#####                                          #####
####################################################

#####   ACCEPT-REJECT METHODS

# So have many distributions where inverse
# or general transforms fail. What to do?

# Use indirect method: Generate 'candidate'
# random variable and accept it subject to
# passing a test.

## Must know:

# Functional form of the density f (target
# density) up to a multiplicative constant M

# We simulate with a simpler density g
# (candidate density) to generate rv.

# Constraints on candidate density g:

# 1) f and g have compatible supports
# (e.g. g(x) > 0 when f(x) > 0)

# 2) There is a constant M with 
# f(x)/g(x) l.t.e.t. M for all x.

# If these constraints are met, can test with
# an equality (SEE SLIDES). If inequality is 
# met, accept Y; otherwise, discard both
# Y and U and start over

# R implementation is of accept/reject 
# algorithm is straightforward. If randg()
# is a function that randomly generates
# from g() (like rnorm() or rpois())

# this code produces a single generation 
# y from f:
u=runif(1)*M
y=randg(1)
while (u>f(y)/g(y)){
  u=runif(1)*M
  y=randg(1)}

# Can prove that cdf of accepted random
# variable Y is exactly the cdf of X

# Accept-Reject Example: Beta Random Variable
# SEE WIKIPEDIA BETA DISTRIBUTION

# Cannot generate Beta rvs with general transforms
# But can do it with accept-reject method
# Use as instrumental distribution U[0,1]
# where both alpha and beta > 1 (note generic
# rbeta() function in R does not impose this
# restriction).

# Upper bound M is maximum of beta density...
# can obtain with optimize() function

?optimize
?dbeta # distribution function for beta

# This gives us M 
# optimize(f=function(x){dbeta(x,2.7,6.3)},
#         + interval=c(0,1),max=T)$objective
# [1] 2.669744

# Since the candidate density g is equal to one
# the proposed value Y is accepted if
# M x U < f(Y), that is, if M x U is under
# the beta density f at that realization.

# So alternative R code for Accept-Reject 
# with alpha=2.7 and beta=6.3 is:

Nsim=2500
a=2.7;b=6.3
M=2.67
u=runif(Nsim,max=M) #uniform over (0,M)
y=runif(Nsim) #generation from g
x=y[u<dbeta(y,a,b)] #accepted subsample

### REJECTION METHOD FOR CONTINUOUS RVs
### WITH A DISCRETE INTERVAL USING UNIFORM
## Alternative to Inversion
## Only works with intervals

# SEE SLIDE #5 - continuous random variable X
# with pdf fx in interval(0,4) see slide
# We "sprinkle" points uniformly at random
# under the density function.

# Small target square under pdf has same chance
# of being hit wherever it is located. But must
# 'sprinkle' in rectangle [0,4] x [0,0.5] 
# because 'inversion' has failed. 

# Then reject those above the function line....
# that serves as your simulated sample

## EXAMPLE: Triangular Density
# Triangular pdf Fx(x) =   x   if (0 < x < 1);
# Triangular pdf Fx(x) = (2-x) if (1 <= x <= 2);
# Triangular pdf Fx(x) =   0   otherwise;

# program spuRs/resources/scripts/rejecttriangle.r
?runif
rejectionK <- function(fx, a, b, K) {
  # simulates from the pdf fx 
  # using the rejection algorithm
  # assumes fx is 0 outside [a, b] 
  # and bounded by K
  # note that we exit the infinite loop 
  # using the return statement
  while (TRUE) {
    # 'sprinkles' within [a,b]
    x <- runif(1, a, b)
    # sprinkles in interval[0,K]
    y <- runif(1, 0, K)
    # keep 'sprinkles' under fx(x)
    if (y < fx(x)) return(x) # forces us out
  }
}

# here is fx
fx<-function(x){
  # triangular density
  if ((0<x) && (x<1)) {
    return(x) # pdf in this interval is x
  } else if ((1<x) && (x<2)) {
    return(2-x) # here pdf is 2-x
  } else {
    return(0) # is beyond [0,2]
  }
}

# generate a sample
set.seed(21)
nreps <- 3000 # going for 3000
Observations <- rep(0, nreps) # initialize variable
for(i in 1:nreps)   {
  # fill observations with non-rejected candidates
  Observations[i] <- rejectionK(fx, 0, 2, 1)
}

# plot a scaled histogram of
# the sample and the density on top
hist(Observations, 
     breaks = seq(0, 2, by=0.1), 
     freq = FALSE,
     ylim=c(0, 1.05), 
     main="")
lines(c(0, 1, 2), c(0, 1, 0))