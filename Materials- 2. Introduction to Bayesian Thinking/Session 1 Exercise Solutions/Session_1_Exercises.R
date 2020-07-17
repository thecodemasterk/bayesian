########################################
# Session 1 Introduction to R Exercises
########################################

# The variable Dvds in the student dataset 
# contains the number of movie DVDs owned 
# by students in the class.
# a) Construct a histogram of this variable 
#    using the hist() command.

# to load package into current R session (environment):
library(LearnBayes)

# make student data available
data(studentdata)

# attach dataset
attach(studentdata)

### MOVIE DVDs OWNED BY STUDENTS

# a) Construct a histogram of this
#    variable using the hist() command.
hist(Dvds)

# b) Summarize this variable using 
#    the summary() command.
summary(Dvds)

# c) Use the table() command to construct
#    a frequency table of the individual 
#    values of Dvds that were observed. 
table(Dvds)

# If one constructs a barplot of these 
# tabled values using the command 
# barplot(table(Dvds)) one will see 
# that particular response values are 
# very popular. 
barplot(table(Dvds))

# Is there any explanation for these
# popular values for the number of DVDs
# owned? It looks like multiples of 5's
# are very popular on the low end (say
# less than 50 total). Could this perhaps
# be to marketing programs that sell
# DVDs to students in groups or sets
# of 5?

### STUDENT HEIGHTS

# The variable Height contains the 
# height (in inches) of each student 
# in the class.

# a) Construct parallel boxplots of 
#    the heights using the Gender variable.
boxplot(Height~Gender, ylab="Height by Gender")

# b) If one assigns the boxplot output 
# to a variable
output=boxplot(Height~Gender)
# then output is a list that contains 
# statistics used in constructing the 
# boxplots. Print output to see the 
# statistics that are stored.
output

# c) On average, how much taller are 
#    male students than female students?

# average male height:
male.height=mean(Height[Gender=="male"], na.rm=T)
male.height

# average female height:
female.height=mean(Height[Gender=="female"], na.rm=T)
female.height

# mean difference is:
male.height-female.height
