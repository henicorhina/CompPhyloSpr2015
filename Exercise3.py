#! /usr/bin/env python
# -*- coding: utf-8 -*-


from __future__ import division
from scipy.stats import binom
import matplotlib.pyplot as plt
import scipy as sp
import scipy.stats
import random
import numpy
import math


def binomDraws(n=5, p=0.5):
    """ 
    draws from a binomial disribution using the built in scipy 
    binomial function
    """
    # binom.rvs draws a k from the binomial distribution
    data = binom.rvs(n,p)
    return data


# testing the binomDraws function 
# using the defaults of n = 5 and p = 0.5
test = binomDraws()
#print test

numTrials = input("# of trials? ") # supply the number of trials here (this is n)
data = input("# of observed successes? ")  # Supply observed number of successes here (this is k).


# Set up a list with all relevant values of p

#pList = numpy.arange(0, 1.05, 0.05) # this was returning an array...
pList = [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0]
#print pList


# Calculate the likelihood scores for these values of p, in light of the data you've collected

# empty list to store likelihood values
likeScores = []
"""
Just for reference:
n = numTrials
k = data
p = pList
"""

def binomPMF(k, n, p):
    """
    Returns the Binomial PMF given inputs of n, k, and p
    """
    pmf = binom.pmf(k, n, p)
    return pmf

# testing the function:
#x = binomPMF(4, 5, 0.5)
#print x 

# Takes each value in pList and runs it through the Binomial PMF function with the 
# user inputs of n and p and appends them to the likeScores list
for i in pList:
    x = binomPMF(data, numTrials, i)
    likeScores.append(x)
    

# Calculates the maximum likelihood value of p (at least, the max in this set)
print "\n", "the maximum likelihood value of p is: ", max(likeScores)


# What is the strength of evidence against the most extreme values of p (0 and 1)?

"""
Because their likelihood values are both 0.0
"""

# Calculates the likelihood ratios comparing each value (in the numerator) 
# to the max value (in the denominator)
ratios = []
for x in likeScores:
    y = x / max(likeScores)
    ratios.append(y)

# pulls the probability associated with the maximum ratio value
maxRatio = ratios.index(1.0)
print "\n", "the probability is: ", pList[maxRatio]    


# plot the likelihood ratios
plt.plot(ratios, 'bo', ms=8)
plt.xlabel("probability index")
plt.ylabel("likelihood ratios")
#plt.xlim(xmax=numTrials) # adjust the max leaving min unchanged
plt.show()

# When is the ratio small enough to reject some values of p?


def mean_confidence_interval(data, confidence=0.95):
    """
    Calculutes a 95% confidence interval, assuming a normal distribution
    for a given data set (data)
    Returns a sequence of the mean, mean minus 2 SE and mean plus 2 SE
    """
    b = scipy.stats.sem(data)
    n = len(data)
    h = b * sp.stats.t._ppf((1+confidence)/2., n-1)
    mean = 1.0
    return mean - h

mci = mean_confidence_interval(ratios)
#print mci

# this is the interval of probabilities
newProbList = []

# takes each value in ratios and sees if it is within 2 SE of the mean
# if so, it returns the associated pList value and appends it to a list
# if it is below 2 SE of the mean, it does nothing.
for x in ratios:
    if x > mci:
        y = ratios.index(x)
        g = pList[y]
        newProbList.append(g)
    else:
        ratios = ratios

print "\n", "the interval of probabilities within 2 SE is: ", newProbList

# **** EVERYTHING ABOVE HERE TO BE POSTED TO GITHUB BY TUESDAY, FEB. 3RD. ****
# **** CODE BELOW TO BE POSTED TO GITHUB BY THURSDAY, FEB. 5TH ****

"""
Sometimes it will not be feasible or efficient to calculate the likelihoods for every
value of a parameter in which we're interested. Also, that approach can lead to large
gaps between relevant values of the parameter. Instead, we'd like to have a 'hill
climbing' function that starts with some arbitrary value of the parameter and finds
values with progressively better likelihood scores. This is an ML optimization
function. There has been a lot of work on the best way to do this. We're going to try
a fairly simple approach that should still work pretty well, as long as our likelihood 
surface is unimodal (has just one peak). Our algorithm will be:
(1) Calculate the likelihood for our starting parameter value (we'll call this pCurr)
(2) Calculate likelihoods for the two parameter values above (pUp) and below (pDown)
our current value by some amount (diff). So, pUp=pCurr+diff and pDown=pCurr-diff. To
start, set diff=0.1, although it would be nice to allow this initial value to be set
as an argument of our optimization function.
(3) If either pUp or pDown has a better likelihood than pCurr, change pCurr to this
value. Then repeat (1)-(3) until pCurr has a higher likelihood than both pUp and
pDown.
(4) Once L(pCurr) > L(pUp) and L(pCurr) > L(pDown), reduce diff by 1/2. Then repeat
(1)-(3).
(5) Repeat (1)-(4) until diff is less than some threshold (say, 0.001).
(6) Return the final optimized parameter value.
Write a function that takes some starting p value and observed data (k,n) for a
binomial as its arguments and returns the ML value for p.
To write this function, you will probably want to use while loops. The structure of
these loops is
while (someCondition):
    code line 1 inside loop
    code line 2 inside loop
    
As long as the condition remains True, the loop will continue executing. If the
condition isn't met (someCondition=False) when the loop is first encountered, the 
code inside will never execute.
If you understand recursion, you can use it to save some lines in this code, but it's
not necessary to create a working function.
"""

# Write a function that finds the ML value of p for a binomial, given k and n.

#pCurr = 0.9
#diff = input("what is the difference between the parameter values? (is set to 0.1 initially) ")

k, n, p  = data, numTrials, random.random()


def optimize(k, n, pCurr, diff=0.1):
    pCurr = binom.pmf(k, n, pCurr)
    pUp = pCurr + diff
    pDown = pCurr - diff
    binpCurr = binom.pmf(k, n, pCurr)
    binpUp = binom.pmf(k, n, pUp)
    binpDown = binom.pmf(k, n, pDown)
    
   # def insideloop(binpCurr = binpCurr, binpUp = binpUp, binpDown = binpDown, pCurr = pCurr, pUp = pUp, pDown = pDown, diff = diff):
    while diff > 0.001:
        while diff > 0.001:  
                
            if binpCurr < binpUp:
                pCurr = pUp
                binpCurr = binom.pmf(k, n, pCurr)
                pUp = pCurr + diff
                binpUp = binom.pmf(k, n, pUp)
                                
            elif binpCurr > binpDown:
                pCurr = pDown
                binpCurr = binom.pmf(k, n, pCurr)
                pDown = pCurr - diff
                binpDown = binom.pmf(k, n, pDown)

            else:
                diff *= 0.5
                #return insideloop()
        return pCurr

    #return insideloop()

ML = optimize(k, n, p) #input("what is the difference between the parameter values? (defaults to 0.1) "))
print "the maximum likelihood value is: ", ML

"""
pScores = []
for i in numpy.arange(0,1, 0.0001):
    x = binomPMF(data, numTrials, i)
    pScores.append(x)
"""

"""
In the exercise above, you tried to find an intuitive cutoff for likelihood ratio
scores that would give you a reasonable interval in which to find the true value of 
p. Now, we will empirically determine one way to construct such an interval. To do 
so, we will ask how far away from the true value of a parameter the ML estimate 
might stray. Use this procedure: (1) start with a known value for p, (2) simulate
a bunch of datasets, (3) find ML parameter estimates for each simulation, and then 
(4) calculate the likelihood ratios comparing the true parameter values and the ML
estimates. When you do this, you will be constructing a null distribution of
likelihood ratios that might be expected if the value of p you picked in (1)
was true. Note that the ML values for these replicates are very often greater than
L(true value of P), because the ML value can only ever be >= L(true value). Once 
you have this distribution, find the likelihood ratio cutoff you need to ensure 
that the probability of seeing an LR score that big or greater is <= 5%. 
"""

# Set a starting, true value for p

trueP = 0.5
#trueP = input("enter the true value for P ")

# Simulate 1,000 datasets of 200 trials from a binomial with this p
# If you haven't already done so, you'll want to import the binom class from scipy:
# from scipy.stats import binom
# binom.rvs(n,p) will then produce a draw from the corresponding binomial.

# this is the matrix of 1000 datasets
matrix = []

# list for ML parameter estimates for each trial
mlList = []

# matrix of ML ratios
mlRatios = []

# pulls 1000 datasets and appends them to matrix
for i in range(10):
    inter = [] # holds each dataset from the for loop below before appending = k
    mlList2 = [] # list for ML ratios from within the bottom for loop
    listRatios = [] # list of ratios - l(trueP):ML
    for n in range(200):
        inter.append(binom.rvs(n, trueP)) # calls the binomial distribution for each value of n)
    for x in inter:
        mlList2.append(optimize(x, n, trueP)) # calls the ML function for each value in the inter datasets
    for val in mlList2:    
        listRatios.append(binom.pmf(x, n, trueP) / val)
    mlRatios.append(listRatios)
    mlList.append(mlList2) # appends the list of MLs to the mlList
    matrix.append(inter)     # appends the newly generated datasets to the matrix 


#for cols in zip(*matrix):
#    print cols
''' 
for i in matrix:  
    rangeList = []
    for y in range(200):
        
    mlList.append(rangeList)
'''
# Calculate likelihood ratios comparing L(trueP) in the numerator to the maximum
# likelihood (ML) in the denominator. Sort the results and find the value
# corresponding to the 95th percentile.

percentile = []

for i in mlRatios:
    perc = numpy.percentile(i, 95)
    percentile.append(perc)
    

# Now, convert the likelihood ratios (LRs) to -2ln(LRs) values.
# Find the 95th percentile of these values. Compare these values to this table:
# https://people.richland.edu/james/lecture/m170/tbl-chi.html. In particular, look
# at the 0.05 column. Do any of these values seem similar to the one you calculated?
# Any idea why that particular cell would be meaningful?

lnPerc = []

for i in percentile:
    lnPerc.append(-2 * math.log1p(i))

# Based on your results (and the values in the table), what LR statistic value 
# [-2ln(LR)] indicates that a null value of p is far enough away from the ML value
# that an LR of that size is <=5% probable if that value of p was true?

plt.plot(percentile, "b-", lnPerc, "r-")
plt.show()



# Using this cutoff, what interval might you report for the 5- and 20-trial data
# sets above?



# We've talked in previous classes about two ways to interpret probabilities. Which
# interpretation are we using here to define these intervals?

