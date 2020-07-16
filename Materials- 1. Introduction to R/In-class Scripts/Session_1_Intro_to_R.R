########################################
# Chapter 1:  An Introduction to R
########################################

# Two main chapter goals:

# 1) Introduce basic commands for summarizing
#    and graphing data; and
# 2) Introduce R as an environment for running
#    Monte Carlo simulations

### EXPLORING STUDENT DATABASE
# Bowling Green State University
# Students in Introductory Statistics class

# Reading data into R

# install.packages(LearnBayes) # only do this once

# to load package into current R session (environment):
library(LearnBayes)

# make data action
data(studentdata)

# read first row or record
studentdata[1,]

# attach it (don't always need to do this)
attach(studentdata)

# R commands to summarize and graph

# look at whole dataset
str(studentdata)

# look at categorical Drink variable
# indicates student's drinking preference
Drink

# tabulate command for 'Drink' column values
table(Drink) # more than half prefer water

# draw a barlot of tabulated results
# we graph the frequencies of the levels
# we add labels to horizontal and vertical axes
barplot(table(Drink),xlab="Drink",ylab="Count")

# what is the barplot function?
# call up help
?barplot # also help(barplot)

# or we could save output of table(Drink) to t
t <- table(Drink) # same as t=table(Drink)

# and then plot t: shows barplot of
# the drinking preferences of statistics students
barplot(t, xlab="Drink",ylab="Count")

# This command is unused, is for batching this file in
# to cause execution to stop at this point
# S=readline(prompt="Type  <Return>   to continue : ")
# S

# This command is also unused, is for opening a separate
# graphics window which is not necessary for us.
# windows()

# Want to compute how many hours students slept:
# We compute the number of sleeping hours:
hours.of.sleep = WakeUp - ToSleep

# then to retrieve descriptive stats 
# on this variable:
summary(hours.of.sleep)

# So, on average students slept 7.5 hours and
# one half slept between 6.5 and 8.5 hours

# draws a histogram:
hist(hours.of.sleep, main="Histogram of Sleep Hours")

# Histogram is symmetric about 7.5 average

# This command is unused, is for batching this file in
# to cause execution to stop at this point
# S=readline(prompt="Type  <Return>   to continue : ")
# S

# This command is also unused, is for opening a separate
# graphics window which is not necessary for us.
# windows()

# Section 1.2.4 R Commands to Compare Batches

# Want to make comparisons of men vs. women

# draw parallel boxplots of hours of sleep 
# by each level (male, female) of gender:
boxplot(hours.of.sleep~Gender,
  ylab="Hours of Sleep")

# Looks like men and women are similar with
# respect to their sleeping times

# For other variables, there are substantial
# differences between the two genders

# We divide the haircut prices into
# the two groups, '==' is a logical operator

# Price of female haircuts:
female.Haircut=Haircut[Gender=="female"]
female.Haircut

# Price of male haircuts:
male.Haircut=Haircut[Gender=="male"]
male.Haircut

# Get quantitative summaries of each:
summary(female.Haircut)
summary(male.Haircut)

# There are large differences between the men
# and women. The men average about $10.54 per
# haircut and the women about $34.08

# Section 1.2.5 Commands for Studying Relationships

# Is length of sleep related to time one
# goes to bed? Use a scatterplot but it is
# difficult to see a pattern since many of
# the points are identical
plot(ToSleep,hours.of.sleep)

# Use jitter() function which adds small amount
# of 'noise' . . . moves the points a bit
plot(jitter(ToSleep),jitter(hours.of.sleep))

# We fit a linear regression line to see the
# decreasing pattern. lm() is least squares
fit=lm(hours.of.sleep~ToSleep)

# We store fitted object in object 'fit'
# which shows the intercept and slope
fit

# slope is approximately -0.5 which means that
# a student loses about a half hour of sleep

# there is much more information in the
# lm()-fitted object
str(fit)

# draw an abline of fitted object 
# on the existing plot
abline(fit)

# Exploring the Robustness of t statistic
# Section 1.3.2 Write a Function to compute t-stat

# We generate some random data, normally distributed:
# x and y have 10 obs with mean of 50 and sd of 10
x=rnorm(10,mean=50,sd=10)
mean(x)
var(x)
sd(x)

# x is vector of 10 observations (values)
x
y=rnorm(10,mean=50,sd=10)
mean(y)
var(y)
sd(y)

# y is a vector of 10 points
y
m=length(x);m

# 10 observations in x and y
n=length(y);n

# Here we compute the pooled standard deviation
sp=sqrt(((m-1)*sd(x)^2+(n-1)*sd(y)^2)/(m+n-2))
sp

# with m, n and sp defined, we can manually
# compute the t-stat
t.stat=(mean(x)-mean(y))/(sp*sqrt(1/m+1/n))
t.stat

# But this is cumbersome, let's write a function
# to do these calculations 'all at once':
tstatistic=function(x,y)
{
  m=length(x)
  n=length(y)
  sp=sqrt(((m-1)*sd(x)^2+(n-1)*sd(y)^2)/(m+n-2))
  t.stat=(mean(x)-mean(y))/(sp*sqrt(1/m+1/n))
  return(t.stat)
}

# make up some data for x and y
data.x=c(1,4,3,6,5)
data.y=c(5,4,7,6,10)

# and compute the t-stat for the difference
# in the means between x and y:
tstatistic(data.x, data.y)

# Programming a Monte Carlo Simulation
# Section 1.3.3

# We want to learn about the true significance level 
# for this t-statistic of mean differences when the
# populations do not follow the standard assumptions
# of normality and equal variances

# The true significance level depends on:
# 1) stated level of significance alpha
# 2) the shape of the populations (normal, 
#    skewed, heavy-tailed)
# 3) the spreads of the two populations as
#    measured by two standard deviations
# 4) sample sizes m and n

# That is, given a particular alpha level, shape,
# spreads, and sample sizes, we want to estimate
# the true significance level given by

# alphaT = P(|T|>= tn+m-2,alpha/2)

# Steps of the algorithm are:
# 1) get random samples from x and y pops
# 2) Compute the t-stat that the means of
#    the two samples are different
# 3) Decide if |T| exceeds the critical point
#    and reject null hypothesis Ho that the
#    means are not different

# We repeat these 3 steps N times and then we
# estimate the true significance level by:

# alphaT = number of rejections of Ho / N

# simulation algorithm for normal populations
# (means pops have mean of 0 and sd of 1):
# sets alpha, m, n
alpha=.1; m=10; n=10

# sets the number of simulations
N=10000

# counter of num. of rejections, set to 0
n.reject=0

# Here we implement the simulation algorithm:
# loop through it 10000 times:
for (i in 1:N)
{
  # simulates xs from population 1:
  x=rnorm(m,mean=0,sd=1)
  # simulates ys from population 2:
  y=rnorm(n,mean=0,sd=1)
  # computes the t statistic:
  t.stat=tstatistic(x,y)
  # reject if |t| exceeds critical pt:
  # qt(p,df) is the pth quantile of a 
  # t distribution with df degrees of freedom
  if (abs(t.stat)>qt(1-alpha/2,n+m-2))
    n.reject=n.reject+1
}
# est. is proportion of rejections
# true.sig.level is observed sig level
true.sig.level=n.reject/N
true.sig.level

# But what if normality assumptions are not met?
# What then is the behavior of the
# true significance level?

# simulation algorithm for normal and exponential populations
# storing the values of the t statistic in vector tstat

m=10; n=10 
my.tsimulation=function()
  # rnorm() is normal; rexp is exponential
  tstatistic(rnorm(m,mean=10,sd=2), rexp(n,rate=1/10))
# we run the function 10000 times
# tstat.vector stores the simulated values of the
# t-statistic
tstat.vector=replicate(10000, my.tsimulation())
# take a look:
tstat.vector

# the R command density() constructs a nonparametric
# density estimate of the exact sampling distribution
# of the t-statistic.
plot(density(tstat.vector),xlim=c(-5,8),ylim=c(0,.4),lwd=3)
# curve() plots the t density with 18 dof on top
curve(dt(x,df=18),add=TRUE)

# Note that the actual sampling distribution of the
# t statistic is right-skewed, which accounts for
# the large true significance level.

