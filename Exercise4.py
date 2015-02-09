#! /usr/bin/env python
# -*- coding: utf-8 -*-

# author: henicorhina

# Discrete-time Markov chains


from __future__ import division
import random
import numpy as np

"""
In this exercise, we will explore Markov chains that have discrete state spaces
and occur in discrete time steps. To set up a Markov chain, we first need to 
define the states that the chain can take over time, known as its state space.
To start, let's restrict ourselves to the case where our chain takes only two
states. We'll call them A and B.
"""

# Create a tuple that contains the names of the chain's states

# tuple of state space
tup = ('A', 'B');


"""
The behavior of the chain with respect to these states will be determined by 
the probabilities of taking state A or B, given that the chain is currently in 
A and B. Remember that these are called conditional probabilities (e.g., the 
probability of going to B, given that the chain is currently in state A is 
P(B|A).)
We record all of these probabilities in a transition matrix. Each row
of the matrix records the conditional probabilities of moving to the other
states, given that we're in the state associated with that row. In our example
row 1 will be A and row 2 will be B. So, row 1, column 1 is P(A|A); row 1, 
column 2 is P(B|A); row 2, column 1 is P(A|B); and row 2, column 2 is P(B|B). 
All of the probabilities in a ROW need to sum to 1 (i.e., the total probability
associated with all possibilities for the next step must sum to 1, conditional
on the chain's current state).
In Python, we often store matrices as "lists of lists". So, one list will be 
the container for the whole matrix and each element of that list will be 
another list corresponding to a row, like this: mat = [[r1c1,r1c2],[r2c1,r2c2]]. 
We can then access individual elements use two indices in a row. For instance,
mat[0][0] would return r1c1. Using just one index returns the whole row, like
this: mat[0] would return [r1c1,r1c2].
Define a transition matrix for your chain below. For now, keep the probabilties
moderate (between 0.2 and 0.8).
"""

# Define a transition probability matrix for the chain with states A and B

matrix = [[0.8, 0.2], [0.3, 0.7]]

# Try accessing a individual element or an individual row 
# Element

# Accessing the first element of the second row
print "\n", "the first element of the second row is: ", matrix[1][0]

# Row

# Accessing the first row
print "\n", "the first row of the matrix is: ", matrix[0]


"""
Now, write a function that simulates the behavior of this chain over n time
steps. To do this, you'll need to return to our earlier exercise on drawing 
values from a discrete distribution. You'll need to be able to draw a random
number between 0 and 1 (built in to scipy), then use your discrete sampling 
function to draw one of your states based on this random number.
"""

def marcov (i = random.choice(tup), step = 10):
    """
    Marcov chain simulator with a state space of 2
    i = the initial state, which is set to draw a random state to initiate the chain
    step = how many steps to run the simulation for
    
    """
    sims = [] # List to hold the results of the Marcov chain
    for x in range(step):
        
        if i == 'A':
            x = np.random.random() # Random number generator
            if matrix[0][0] > x:
                sims.append('A')
            else:
                sims.append('B')
        else:
            y = np.random.random()
            if matrix[1][0] > y:
                sims.append('A')
            else:
                sims.append('B')
    return sims
            

# draws a random state from the tuple of state space
#state = random.choice(tup)

run = marcov(random.choice(tup), input("how many steps would you like to run the simulator for? "))

print "\n", run



# Run a simulation of 10 steps and print the output.

# the output: ['B', 'A', 'B', 'A', 'B', 'B', 'A', 'B', 'A', 'B']



# ----> Try to finish the above lines before Tues, Feb. 10th <----

# Now try running 100 simulations of 100 steps each. How often does the chain
# end in each state? How does this change as you change the transition matrix?

endVal = []

for x in range(100):
    run = marcov(random.choice(tup), 100)
    endVal.append(run[-1])


print "\n", "there are", endVal.count("A"), "A end values and", endVal.count("B"), "B end values"

"""
with a transition matrix of: [[0.8, 0.2], [0.3, 0.7]]
there are 60 A end values and 40 B end values


with a transition matrix of: [[0.5, 0.5], [0.5, 0.5]]
there are 60 A end values and 40 B end values
"""



# Try defining a state space for nucleotides: A, C, G, and T. Now define a 
# transition matrix with equal probabilities of change between states.


stateSpace = ('A', 'C', 'G', 'T');

def marcovNuc (i = random.choice(stateSpace), step = 100):
    """
    Marcov chain simulator with a state space of 4 (e.g. nucleotides)
    i = the initial state, which is set to draw a random state to initiate the chain
    step = how many steps to run the simulation for
    
    """
    # matrix of transition probabilities
    # matrix = [[0.25, 0.25, 0.25, 0.25], [0.25, 0.25, 0.25, 0.25], [0.25, 0.25, 0.25, 0.25], [0.25, 0.25, 0.25, 0.25]] 
    matrix = [[0.4, 0.3, 0.2, 0.1], [0.4, 0.3, 0.2, 0.1], [0.4, 0.3, 0.2, 0.1], [0.4, 0.3, 0.2, 0.1]] 
    sims = [] # List to hold the results of the Marcov chain
    for x in range(step):
        
        if i == 'A':
            w = np.random.random() # Random number generator
            if matrix[0][0] > w:
                sims.append('A')
            elif matrix[0][1] + matrix[0][0] > w:
                sims.append('C')
            elif matrix[0][2] + matrix[0][1] + matrix[0][0] > w:
                sims.append('G')
            else:
                sims.append('T')
        elif i == 'C':
            x = np.random.random()
            if matrix[1][0] > x:
                sims.append('A')
            elif matrix[1][1] + matrix[1][0] > x:
                sims.append('C')
            elif matrix[1][2] + matrix[1][1] + matrix[1][0] > x:
                sims.append('G')
            else:
                sims.append('T')
    
        elif i == 'G':
            y = np.random.random()
            if matrix[2][0] > y:
                sims.append('A')
            elif matrix[2][1] + matrix[2][0] > y:
                sims.append('C')
            elif matrix[2][2] + matrix[2][1] + matrix[2][0] > y:
                sims.append('G')
            else:
                sims.append('T')

        else:
            z = np.random.random()
            if matrix[3][0] > z:
                sims.append('A')
            elif matrix[3][1] + matrix[3][0] > z:
                sims.append('C')
            elif matrix[3][2] + matrix[3][1] + matrix[3][0] > z:
                sims.append('G')
            else:
                sims.append('T')

    return sims
            

         
# Again, run 100 simulations of 100 steps and look at the ending states. Then
# try changing the transition matrix.

simsNuc = marcovNuc()

print simsNuc



"""
here are the ending states for a matrix of equal probabilities and 100 steps:

['C', 'G', 'T', 'C', 'T', 'A', 'G', 'C', 'T', 'C', 'C', 'G', 'T', 'C', 'A', 
'G', 'A', 'G', 'T', 'G', 'T', 'C', 'C', 'G', 'C', 'A', 'C', 'G', 'C', 'T', 
'G', 'A', 'C', 'A', 'G', 'A', 'G', 'A', 'G', 'A', 'C', 'G', 'G', 'C', 'C', 
'C', 'A', 'C', 'T', 'G', 'T', 'G', 'A', 'A', 'G', 'T', 'C', 'G', 'G', 'C', 
'T', 'G', 'A', 'T', 'C', 'C', 'G', 'G', 'C', 'C', 'G', 'G', 'G', 'C', 'G', 
'G', 'A', 'T', 'C', 'C', 'C', 'A', 'G', 'G', 'T', 'T', 'A', 'G', 'G', 'T', 
'G', 'C', 'C', 'C', 'G', 'A', 'T', 'T', 'C', 'C']

Counts:
A = 17
C = 32
G = 33
T = 18

Why is the GC content 2x???

"""

"""
here are the ending states for 100 steps and the matrix:

matrix = [[0.4, 0.3, 0.2, 0.1], [0.4, 0.3, 0.2, 0.1], [0.4, 0.3, 0.2, 0.1], [0.4, 0.3, 0.2, 0.1]] 


['A', 'A', 'C', 'G', 'C', 'C', 'A', 'T', 'C', 'G', 'G', 'G', 'G', 'G', 'A', 
'G', 'C', 'A', 'A', 'A', 'A', 'G', 'C', 'G', 'C', 'T', 'C', 'C', 'A', 'G', 
'A', 'A', 'G', 'A', 'C', 'A', 'C', 'A', 'A', 'A', 'C', 'G', 'T', 'G', 'G', 
'A', 'A', 'A', 'G', 'A', 'A', 'G', 'G', 'A', 'G', 'A', 'G', 'A', 'C', 'C', 
'C', 'G', 'C', 'C', 'A', 'A', 'C', 'A', 'A', 'A', 'C', 'A', 'C', 'G', 'G', 
'C', 'C', 'C', 'T', 'A', 'C', 'C', 'A', 'A', 'A', 'C', 'A', 'A', 'G', 'A', 
'C', 'G', 'C', 'G', 'C', 'A', 'C', 'T', 'A', 'A']

Counts:
A = 40
C = 30
G = 25
T = 5

This matches the 4:3:2:1 ratio of the matrix very well.

"""





