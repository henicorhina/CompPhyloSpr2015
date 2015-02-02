#! /usr/bin/env python
# -*- coding: utf-8 -*-


from __future__ import division
from scipy.stats import binom
import matplotlib.pyplot as plt
import scipy as sp
import scipy.stats


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

data = input("# of observed successes? ")  # Supply observed number of successes here.
numTrials = input("# of trials? ") # supply the number of trials here


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


"""
Now let's try this all again, but with more data. 
This time, we'll use 20 draws from our cup of marbles.
Run the file again with k (data) = 4 and n (numTrials) = 20
"""

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
    return mean, mean - h

mci = mean_confidence_interval(ratios)
#print mci

neg2SE = mci[1]
#print neg2SE

# this is the interval of probabilities
newProbList = []

# takes each value in ratios and sees if it is within 2 SE of the mean
# if so, it returns the associated pList value and appends it to a list
# if it is below 2 SE of the mean, it does nothing.
for x in ratios:
    if x > neg2SE:    
        y = ratios.index(x)        
        g = pList[y]
        newProbList.append(g)

    else:
        ratios = ratios

print "\n", "the interval of probabilities within 2 SE is: ", newProbList
