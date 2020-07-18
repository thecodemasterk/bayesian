###########################################
###    Session 2 Exercise Solutions     ###
###########################################

# Discrete Prior
p <- seq(0, 1, by = 0.01)
p
prior <- 1 / 101 + 0 * p
prior
plot(p, prior,
     type="h",
     main="Prior Distribution")
# The prior distribution looks uniform

# Posterior
library(LearnBayes)
post <- pdisc(p, prior, c(20, 12))
post
plot(p, post,
     type="h",
     main="Posterior Distribution")
# The posterior distribution looks normal

# Highest Probabilities
?discint
discint(cbind(p, post), 0.90)
# The probability that p falls in the interval 
# {0.49, 0.75} is approximately 0.90.
discint(cbind(p, post), 0.95)
# The probability that p falls in the interval 
# {0.46, 0.77} is approximately 0.95.
discint(cbind(p, post), 0.99)
# The probability that p falls in the interval 
# {0.40, 0.81} is approximately 0.99.

