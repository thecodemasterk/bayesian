###############################
# Exercise Solution:
# Comparing Two Proportions
###############################

# load LearnBayes package
library(LearnBayes)

# yN and yS have independent beta distributions with
# yN distributed beta(1602, 162528)
# yS distributed beta(511, 412369)

# (b) Use rbeta() to simulate 1000 values
#     joint posterior distribution of pN and pS
?rbeta
pN=rbeta(1000, 1602, 162528)
pN
pS=rbeta(1000, 511, 412369)
pS

# (c) Using your sample, construct a histogram
#     of relative risk pN/pS. Find a 95%
#     interval estimate of this risk.
rel.risk=pN/pS
rel.risk
hist(rel.risk)
# Here is 95% CI
int.estimate=quantile(rel.risk,c(.025,.975))
int.estimate

# (d) Construct a histogram of the difference
#     in risks pN - pS.
diff.risks=pN-pS
diff.risks
hist(diff.risks)
int.estimate=quantile(diff.risks,c(.025,.975))
int.estimate

# (e) Compute the posterior probability that
#     the difference in risks exceeds 0.
prob=mean(diff.risks>0)
prob
