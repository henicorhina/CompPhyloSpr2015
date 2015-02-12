"""
In-Class Markov Chain Exercise
2.10.15
@author: jembrown
"""


from __future__ import division
import random
import numpy as np

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


def marcov (i = random.choice(tup), step = 10, matrix):
    """
    Marcov chain simulator with a state space of 2
    i = the initial state, which is set to draw a random state to initiate the chain
    step = how many steps to run the simulation for
    
    """
    #step += 1
    probs = []
    sims = [] # List to hold the results of the Marcov chain
    sims.append(i)
    probs.append(0.5)
    for x in range(step):        
        if sims[-1] == 'A':
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
    return sims, probs
#    return sims[1:-1]
            

# draws a random state from the tuple of state space
#state = random.choice(tup)

run = marcov(random.choice(tup), input("how many steps would you like to run the simulator for? "), matrix)

print "\n", run


probs = run[1]


sum = 1 
for i in probs:
    sum *= i

#print sum

"""

def factorial (num):
	if num != 1:
		return num * factorial(num - 1)
	elif num == 1:
		return 1
"""


endVal = []

for x in range(100):
    run = marcov(random.choice(tup), 100)
    endVal.append(run[-1])


#print "\n", "there are", endVal.count("A"), "A end values and", endVal.count("B"), "B end values"



# Calculate the probability of observing the state in step 3, given the initial
# state in step 1 (i.e., as if you didn't know the state in step 2).

n = input("how many time steps would you like to run the matrix for? ")

array = np.matrix(matrix)
future = array ** n

print future



# Now think of the chain progressing in the opposite direction. What is the
# probability of the progression through all 3 states in this direction? How
# does this compare to the original direction?

"""
It would depend on your matrix of probabilities. If the transition from A -> B = B -> A, then 
it would be identical. Otherwise, it will be different.
"""




# Try the same "forward" and "reverse" calculations as above, but with this
# transition matrix:
# revMat = [[0.77,0.23],
#           [0.39,0.61]]
# and these starting frequencies for "a" and "b"
# freq(a) = 0.63   freq(b) = 0.37




# What is (roughly) true about these probabilities?





# Simulate 1,000 replicates  (or 10K if your computer is fast enough) of 25 
# steps. What are the frequencies of the 2 states across replicates through time?

# NOTE: Here is a function that reports the frequencies of a state through time 
# for replicate simulations. You'll need to do this several times during this exercise.

def mcStateFreqSum(sims,state="a"):
    """
    Pass this function a list of lists. Each individual list should be the
    states of a discrete-state Markov chain through time (and all the same 
    length). It will return a list containing the frequency of one state 
    ("a" by default) across all simulations through time.
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



    
# Summarize the frequency of one state through time




# What do you notice about the state frequencies through time? Try another round
# of simulations with a different transition matrix. How do the state freq.
# values change?






# Now, calculate a vector of probabilities for the focal state (e.g., 'a')
# based on the transition matrix directly (not by simulation). How do these
# values compare to the simulated frequencies?









