####################################
# An Example of Rejection Sampling #
####################################

# We are first going to look at a simple 
# example of rejection sampling of the random
# variable Z which has pdf f(z)=6z(1-z) on [0,1]. 

# Note that Z has a Beta(2,2) distribution. 
# Also note that f(z) has a maximum of 3/2 
# at 1/2. 

# So to simulate Z we can simulate X from 
# Unif(0,1) and Y from Unif(0,3/2) and 
# accept X if Y < f(X). 

# We can simulate a single observation
# of Z as follows. 

# runif generates a random deviate from min to max

rz1<-function () 
{
  repeat {
    x <- runif(1, 0, 1)
    y <- runif(1, 0, 3/2)
    fx <- 6 * x * (1 - x)
    if (y < fx) 
      return(x)
  }
}

# try it out:
rz1()

# what is happening in function rz1()?:
rz1.explained<-function () 
  { 
  repeat 
    { 
      # we simulate z from x with a Unif(min=0, max=1)
      x <- runif(1, min=0, max=1); print("x is: ", quote=FALSE)
      print(x)
      # we simulate y with a Unif(min=0, max=3/2)
      y <- runif(1, min=0, max=3/2); print("y is: ", quote=FALSE)
      print(y)
      # fx is 6 times x times (1-x)
      fx <- 6 * x * (1 - x); print("fx is: ", quote=FALSE)
      print(fx)
      # keeps repeating until 
      # y < fx is true
      if (y < fx) 
        # and then returns x
        return(x) 
      } 
  } 

# Try it out:
rz1.explained()

# Note that we repeat the simulation block 
# indefinitely, until the acceptance criterion 
# is satisfied, at which point we break the 
# loop by returning the accepted value. 

# We can now use this to simulate a vector 
# of n such quantities as follows. 

rz2<-function(n) 
  { 
  # creates vector to store values
  zvector <- vector("numeric", n) 
  for (i in 1:n) { 
    # fills up the vector
    zvector[i] <- rz1() 
  } 
  zvector 
  } 

# try it out
rz2(10)

# Of course, this uses explicit looping, 
# and hence will be slow in R. However, we can 
# vectorize this procedure. 

# We do this by simulating vectors of 
# X and Y and keep only those which 
# satisfy the acceptance criterion.

# Consider the following function: 
rz3<-function(n) 
  { 
  # vector of n x's
  xvec <- runif(n, 0, 1)
  # vector of n y's
  yvec <- runif(n, 0, 3/2)
  # vector of fx if criterion is satisfied
  fvec <- 6 * xvec * (1 - xvec)
  # only returns those in vector of 10
  # which satisfy the criterion
  xvec[yvec < fvec] 
  }

# Try it out:
rz3(10)

# The last line of the function body 
# selects and returns those elements 
# of the X vector satisfying the 
# acceptance criterion. 

# So the length of the returned vector
# will be considerably less than n. 

# How much less will depend on the 
# acceptance probability. 

# Note that the length of the returned 
# vector divided by n is an estimate 
# of the acceptance probability. 

# We can use this estimate in order 
# to get an idea of how many simulations 
# we should start off with in order 
# to end up with a vector of about 
# the length we want. 

# Consider the following function. 
rz4<-function (n) 
  { 
  x <- rz3(n) 
  len <- length(x)
  # calculate estimate of
  # acceptance probability
  aprob <- len/n
  # how many more obs we would like:
  shortby <- n - len 
  # how many obs we need to start with
  # to end up with shortby:
  n2 <- round(shortby/aprob) 
  # runs rz3 again with that number
  x2 <- rz3(n2)
  # and adds that output to the vector
  x3 <- c(x, x2) 
  x3 
  } 

# This function calls rz3 to generate 
# a vector of Zs. aprob is an estimate 
# of the acceptance probability. shortby
# is how many more observations we would 
# like, and n2 is an estimate of how 
# many observations we need to start with
# to end up with shortby. The two 
# resulting vectors are joined together
# and returned. 

# So, this function will return a vector
# whose length is approximately n. 

# This is good enough for most purposes, 
# and because it is vectorized, 
# is much faster than rz2. 

?hist

hist(rz4(10000),breaks=30) 
