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

k = input("# of observed successes? ")  # Supply observed number of successes here (this is k).
n = input("# of trials? ") # supply the number of trials here (this is n)


# Set up a list with all relevant values of p

#pList = numpy.arange(0, 1.05, 0.05) # this was returning an array...
pList = [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0]
#print pList


# Calculate the likelihood scores for these values of p, in light of the data you've collected

# empty list to store likelihood values
likeScores = []

def functionMult(y, x):
    """
    multiplies all consecutively decreasing numbers 
    between a maximum and a minimum supplied as arguments
    """
    # increases the maximum value by one, since a range does not include the max value
    x += 1
    
    # creates a list of numbers between the min and max values
    list = range(y,x) 
    
    # uses reduce and lambda to multiply all the values of a list together
    return(reduce(lambda x, y: x*y, list))
    
# calls the function with user input of max and min values
#print(functionMult(input("What is the minimum value? "), input("What is the maximum value? ")))


# Question 2a

# I think this is the way that you wanted us to do it, but I'm really not sure
def binomialCoefficientA(n,k):
    """
    a function that calculates the binomial coefficient
    
    uses the formula: n(n-1)...(n-k+1) / k!
    calling the factorial function. Note that you need to specify 1 as the minimum value
    the first argument is the numerator and the second two arguments is the denominator
    """    
    return (functionMult((n-k+1), n) / (functionMult(1, k)))

# calls the function with the user inputs of n and k
bino = binomialCoefficientA(n, k)
#print  "\n", "the binomial coefficient is: ", bino

# Question 4

def bernoulli(p):
    """
    function that returns the binomial(n,p) distribution of k successes 
    in n bernoulli trials, given a user-input of probability (prob)
    note that "bino" is the return from the binomial coefficient
    """
    return bino * (pow(p,k) * pow((1 - p),(n-k)))


# testing the function:
#x = binomPMF(4, 5, 0.5)
#print x 

# Takes each value in pList and runs it through the Binomial PMF function with the 
# user inputs of n and p and appends them to the likeScores list
for i in pList:
    x = bernoulli(i)
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
