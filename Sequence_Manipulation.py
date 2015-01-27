#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sun 25 Jan 2015

@author: Oscar Johnson github.com/henicorhina

Assignment 2 for Computational Phylogenetics at LSU
"""
# Practice with using Discrete Sampling

import random
from scipy.stats import rv_discrete
import matplotlib.pyplot as plt

# Question 1

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
    #print(reduce(lambda x, y: x*y, list)) 
    return(reduce(lambda x, y: x*y, list))
    
# calls the function with user input of max and min values
print(functionMult(input("What is the minimum value? "), input("What is the maximum value? ")))


# Question 2a
n = input("What is n? ")
k = input("What is k? ")

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
print  "\n", "the binomial coefficient is: ", bino

# Question 2b

def binomialCoefficientB(n,k):
    """
    a function that calculates the binomial coefficient
    
    uses the formula: n! / (n-k)! k!
    calling the factorial function. Note that you need to specify 1 as the minimum value
    the first argument is the numerator
    the second two arguments are the denominator
    """    
    return (functionMult(1, n)) / ((functionMult(1, (n-k)) * (functionMult(1, k))))
    
# calling the function for n and k
#slow = binomialCoefficientB(n, k)
#print  "\n", "the binomial coefficient is: ", slow
    

# Question 3

"""
I ran the two equations with n = 100,000 and k = 500
The first equation (binomialCoefficientA) was significantly faster than B
A took less than 1 second, B took 7 seconds
This time increased even more when n = 1,000,000 and k = 1,000
My computer almost crashed when I tried to run B, but A still ran in ~1 second
"""


# Question 4

def bernoulli(p):
    """
    function that returns the binomial(n,p) distribution of k successes 
    in n bernoulli trials, given a user-input of probability (prob)
    note that "bino" is the return from the binomial coefficient
    """
    return bino * (pow(p,k) * pow((1 - p),(n-k)))

prob = input("what is the probability of the bernoulli distribution (enter a decimal value): ")
print "\n", "the binomial(n,p) distribution is: ", bernoulli(prob)

size = input("how many trials do you want to run? ")

# Plot a Probability Mass Function (PMF) distribution	

# newList to append bernoulli probabilities
newList = []	
for x in range(size):
    # appends bernoulli probabilities based on the user input number of trials
    newVal = random.random() # creates random probabilities
    newList.append(bernoulli(newVal)) # runs through the bernoulli function

# plot to histogram using matplotlib

plt.hist(newList)
plt.xlabel("k (Success!!!)")
plt.ylabel("Probability Mass Function")
plt.show()	


# Question 5

def discreteSamp(happyCat, sadCat, s):
    """
    function to sample from a discrete distribution. Most of this was copied off
    the internet, so i'm not entirely sure how it works (ask Glaucia). 
    takes three inputs, happyCat = list of events, sadCat = list of probabilities,
    and s = number of trials
    """
    trials = rv_discrete(name='trials', values=(happyCat, sadCat))
    sample = trials.rvs(size = s)
    x = []
    x.append(sample)
    return x

# just a list of six events (die rolls) and associated probabilties to test the function
xk = [1, 2, 3, 4, 5, 6]
pk = [0.1, 0.2, 0.3, 0.1, 0.1, 0.1]

# calls the function
discrete = discreteSamp(happyCat=xk, sadCat=pk, s=size)

print "\n", " the list of results is: ", discrete.pop(0)


# Question 6









    
