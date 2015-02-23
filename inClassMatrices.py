#! /usr/bin/env python

"""
In-Class Markov Chain Exercise
2.10.15
@author: henicorhina
"""


from __future__ import division
import random
import numpy as np
import matplotlib.pyplot as plt

"""
Recall from your reading that any irreducible and aperiodic Markov chain has a 
stationary distribution. To convince ourselves that things will converge for 
such a chain with arbitrary transition probabilities, let's give it a try.
Work in pairs for this. It's more fun to be social.
"""

# Paste your Markov chain simulation function below, where the starting state
# is drawn with uniform probability from all possible states. Remember to also
# copy any import statements or other functions on which your simulator is
# dependent.

tup = ('A', 'B')

matrix = [[0.9, 0.1], [0.3, 0.7]]


def marcov(i = random.choice(tup), step = 10, matrix = [[0.9, 0.1], [0.3, 0.7]]):
    """
    Marcov chain simulator with a state space of 2
    i = the initial state, which is set to draw a random state to initiate the chain
    step = how many steps to run the simulation for    
    """
    step -= 1 # to accomodate the initial probability choice from tup of 0.5
    probs = [] # List to hold the probabilities of the Marcov chain
    sims = [] # List to hold the results of the Marcov chain
    sims.append(i) # state to start the chain
    probs.append(0.5)
    for x in range(step):        
        if sims[-1] == 'A': # looks at the last value in the chain and decides what to do
            x = np.random.random() # Random number generator
            if matrix[0][0] > x:
                sims.append('A')
                probs.append(matrix[0][0]) 
            else:
                sims.append('B') 
                probs.append(matrix[0][1])
        else:
            y = np.random.random()
            if matrix[1][0] > y:
                sims.append('A')
                probs.append(matrix[1][0])                
            else:
                sims.append('B')
                probs.append(matrix[1][1])
    #return sims # if you only want to the return the chain and not the probabilities
    return sims, probs
            
run = marcov(random.choice(tup), input("how many steps would you like to run the simulator for? "), matrix)
print "\n", run

probs = run[1] # the probabilities associated with each transition in the chain
chain = run[0] # the marcov chain

# sums all the probabilities for the transition states in the marcov chain
sum = 1 
for i in probs:
    sum *= i

print "\n", "the sum of probabilities is: ", sum

"""
# just a little recursive function for a factorial
# has no real purpose here...
# I just needed a place to save it

def factorial (num):
	if num != 1:
		return num * factorial(num - 1)
	elif num == 1:
		return 1
"""

# runs the marcov chain simulator for 100 chains each of length 100
# appends the last value of each chain to the endVal list
endVal = []
for x in range(100):
    run = marcov(random.choice(tup), 100, matrix)
    endVal.append(run[-1])

print "\n", "there are", endVal.count("A"), "A end values and", endVal.count("B"), "B end values"


# Calculate the probability of observing the state in step 3, given the initial
# state in step 1 (i.e., as if you didn't know the state in step 2).

if chain[-1] == "A":
    p = (matrix[0][0] * matrix[0][1]) + (matrix[0][0] * matrix[1][0])
    print "the probability of observing A in step 3 is: ", p
else: 
    p = (matrix[1][0] * matrix[1][1]) + (matrix[0][1] * matrix[1][1])
    print "the probability of observing B in step 3 is: ", p


# calculates the standing frequencies of a matrix using a numpy array
# to calculate the probability of observing A (or B) over n time steps
n = input("how many time steps would you like to run the matrix for? ")
array = np.matrix(matrix)
future = array ** n

#print "\n", "the standing state frequencies are: ", future

print "the probability of observing A is: ", future[0, 0]
print "the probability of observing B is: ", future[1, 0]


# Now think of the chain progressing in the opposite direction. What is the
# probability of the progression through all 3 states in this direction? How
# does this compare to the original direction?

"""
It would depend on your matrix of probabilities. If the transition from A -> B = B -> A, then 
it would be identical. Otherwise, it will be different.
"""

# reverse the chain and the array
# note that reversing the chain is not particularly necessary here
chainReversed = chain[::-1]
reversed_array = array[::-1]
reversed_future = reversed_array ** n

# print the reversed array over n time steps
print "\n", "the reversed array over n time steps: ", reversed_future

# Try the same "forward" and "reverse" calculations as above, but with this
# transition matrix:
revMat = [[0.77,0.23],
          [0.39,0.61]]
# and these starting frequencies for "a" and "b"
freqA = 0.63   
freqB = 0.37

# chooses the starting state i from the starting frequencies
choice = random.random()
if choice > freqA:
    i = "A"
else: 
    i = "B"

# runs the marcov simulator for 3 times steps and the above starting frequencies and matrix
revRun = marcov(i, 3, revMat)
#print revRun

# reverses the matrix
revRevMat = revMat[::-1]

# runs the marcov simulator for 3 time steps and using the reversed matrix
revRevRun = marcov(i, 3, revRevMat)
#print revRevRun


# What is (roughly) true about these probabilities?
"""
they appear to actually be close to the standing state frequencies, surprisingly.
forward: (['B', 'B', 'A'], [0.37, 0.61, 0.39])
reversed: (['B', 'A', 'B'], [0.37, 0.77, 0.61])
"""

# Simulate 1,000 replicates  (or 10K if your computer is fast enough) of 25 
# steps. What are the frequencies of the 2 states across replicates through time?

def mcStateFreqSum(sims,state="A"):
    """
    Pass this function a list of lists. Each individual list should be the
    states of a discrete-state Markov chain through time (and all the same 
    length). It will return a list containing the frequency of one state 
    ("A" by default) across all simulations through time.
    """
    freqs = []
    for i in range(len(sims[0])):  # Iterate across time steps
        stateCount = 0
        for j in range(len(sims)): # Iterate across simulations
            if sims[j][i] == state:
                stateCount += 1
        freqs.extend([float(stateCount)/float(len(sims))])
    return freqs

# Run replicate simulations 
yetAnotherMatrix = [[0.9, 0.1], 
                 [0.1, 0.9]]

list10000 = [] # blank list for the 10000 runs
for x in range(10000):
    run10000 = marcov(random.choice(tup), 25, yetAnotherMatrix)
    list10000.append(run10000[0])

    
# Summarize the frequency of one state through time
summaryA = mcStateFreqSum(list10000, "A")
summaryB = mcStateFreqSum(list10000, "B")

plt.plot(summaryA)
plt.plot(summaryB)
plt.show()


# What do you notice about the state frequencies through time? Try another round
# of simulations with a different transition matrix. How do the state freq.
# values change?


"""
After the first few runs, the frequncies reach the standing state frequencies.
matrix = [[0.77, 0.23],
          [0.39, 0.61]]
using 10,000 runs: 
A = 0.63012380952380964
B = 0.36987619047619047
The vector probabilities for the focal state: 
A = 0.62903226
B = 0.37096774
Pretty close, but the values seem to tend towards each other.



Using the matrix 
Matrix = [[0.9, 0.1], 
          [0.1, 0.9]]
using 10,000 runs: 
A = 0.50035714285714283
B = 0.49964285714285694
the vector probabilities are all 0.5
Really close!
"""


# Now, calculate a vector of probabilities for the focal state (e.g., 'a')
# based on the transition matrix directly (not by simulation). How do these
# values compare to the simulated frequencies?

# I did this above

