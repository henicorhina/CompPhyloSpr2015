#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Markov Chain Exercise using object-oriented Python
2.24.15
@author: henicorhina
"""

from __future__ import division
import random
import numpy as np
import scipy as sp
from scipy.stats import binom
from scipy.linalg import expm


class markov(object):
    """base class for markov chain simulators
    
       Attributes:
       stateSpace = character state space of the simulation
       numSims = how many simulations to run the simulator for
       chainLength = desired length of the branches for a ct chain
                   = number of time steps (length of chain) for a dt chain
       matrix = transition probability matrix
       Q = Q-matrix
    """    
    
    def __init__(self, stateSpace = ('A', 'C', 'G', 'T'), numSims = 1, chainLength = 10, 
                matrix = [[0.8, 0.2], [0.3, 0.7]],
                Q = np.array([[-1.916, 0.541, 0.787, 0.588], 
                              [0.148, -1.069, 0.415, 0.506], 
                              [0.286, 0.170, -0.591, 0.135],
                              [0.525, 0.236, 0.594, -1.355]])):                       
        self.stateSpace = stateSpace
        self.numSims = numSims
        self.chainLength = chainLength
        self.matrix = matrix 
        self.Q = Q # A Q-matrix
        self.chain = [] # List to hold the results of the Marcov chain
        self.times = [] # List to hold the probabilities of the Marcov chain
        self.margQ = () # holds the array for the marginal probabilities of a Q-matrix
        self.tProb = () # holds the array for the transitional probabilities of a Q-matrix
        
    def discSamp(self, events, probs):
        """
        This function samples from a list of discrete events provided in the events
        argument, using the event probabilities provided in the probs argument. 
        These lists must:
            - Be the same length
            - Be in corresponding orders
        Also, the probabilities in probs must sum to 1.
        """
        ranNum = sp.random.random()
        cumulProbs = []
        cumulProbs.extend([probs[0]])
        for i in range(1,len(probs)):
            cumulProbs.extend([probs[i]+cumulProbs[-1]])
        for i in range(0,len(probs)):
            if ranNum < cumulProbs[i]:
               return events[i]
        return None

    def standFreq(self):
        """calculates the standing frequencies of a matrix using a numpy array
           to calculate the probability of observing A (or B) over n time steps
        """
        n = input("how many time steps would you like to multiply the matrix for? ")
        future = np.matrix(self.matrix) ** n
        return future

class dtmarkov(markov):
    """ Subclass for running discrete-time markov chains and associated functions
       
      Contains simulators written by Oscar Johnson (only works for a character space
      of two, and rather clumsy) and Jeremy Brown (a more general simulator that
      works a lot better).
    """    
    
    def simulator_oscar(self):
        """
        Discrete-time Markov chain simulator with a state space of 2
        i = the initial state, which is set to draw a random state to initiate the chain
        step = how many steps to run the simulation for    
        """
        i = random.choice(self.stateSpace)
        self.chainLength -= 1 # to accomodate the initial probability choice from tup of 0.5
        self.chain.append(i) # state to start the chain
        self.times.append(0.5)
        for x in range(self.chainLength):        
            if self.chain[-1] == 'A': # looks at the last value in the chain and decides what to do
                x = np.random.random() # Random number generator
                if self.matrix[0][0] > x:
                    self.chain.append('A')
                    self.times.append(self.matrix[0][0]) 
                else:
                    self.chain.append('B') 
                    self.times.append(self.matrix[0][1])
            else:
                y = np.random.random()
                if self.matrix[1][0] > y:
                    self.chain.append('A')
                    self.times.append(self.matrix[1][0])                
                else:
                    self.chain.append('B')
                    self.times.append(self.matrix[1][1])
        return self.chain, self.times
    
    
    def dmcSim(self):
        """
        @author: JeremyBrown
        
        This function simulates the progression of a discrete-time, discrete-state
        Markov chain. It uses 3 arguments: (1) The number of steps (chainLength), (2) the 
        state space, and (3) the transition matrix (matrix). It returns a list containing 
        the progression of states through time. This list should have length chainLength.
        
        The chain will be initiated with a randomly drawn state.
        """
        
        # Draw a state to initiate the chain
        currState = self.discSamp(self.stateSpace,[1.0/len(self.stateSpace) for x in self.stateSpace])
        self.chain.extend(currState)
    
        # Simulate the chain over chainLength-1 steps following the initial state
        for step in range(1,self.chainLength):
            probs = self.matrix[self.stateSpace.index(currState)] # Grabbing row associated with currState
            currState = self.discSamp(self.stateSpace,probs) # Sample new state
            self.chain.extend(currState)        
        return self.chain


    def multSimulate(self):
        """ runs multiple simulations of the discrete time markov simulator
    
            Attributes:
            numSims = number of simulations to be run. defined in base class
        """
        for x in range(self.numSims): 
            run = self.dmcSim()
            self.chain.append(run)
        return self.chain


    def endVals(self):
        """ Calculates the frequency of the ending values of a chain 
            over a user-defined number of simulations
        """
        
        endVal = [] # list to hold ending values of the simulations 
        self.chainLength = input("how long do you want each chain in you simulation to be? ")
        for x in range(input("how many times would you like to run the simulator? ")): 
            run = self.dmcSim()
            endVal.append(run[-1])

        print "\n", "the frequency of A is: ", float(endVal.count("A") / len(endVal))
        print "\n", "the frequency of A is: ", float(endVal.count("B") / len(endVal))
        return endVal
        
        
    def sumProb(self):
        """ Sums all of the probabilities for the transition states in the marcov chain
        """        
        sum = 1 
        for i in self.times:
            sum *= i
        
        print "\n", "the sum of probabilities is: ", sum
        return sum


    def mcStateFreqSum(self,state="A"):
        """
        Pass this function a list of lists. Each individual list should be the
        states of a discrete-state Markov chain through time (and all the same 
        length). It will return a list (freqs) containing the frequency of one state 
        ("A" by default) across all simulations through time.
        """
        freqs = []
        for i in range(len(self.chain[0])):  # Iterate across time steps
            stateCount = 0
            for j in range(len(self.chain)): # Iterate across simulations
                if self.chain[j][i] == state:
                    stateCount += 1
            freqs.extend([float(stateCount)/float(len(self.chain))])
        return freqs


    def likelihood(self):
        """Calculates likelihood scores given discrete time marcov chains
        """
        likeScoresA = []
        likeScoresC = []
        likeScoresG = []
        likeScoresT = []
        pList = [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0]
        
        for i in pList:
            likeScoresA.append(binom.pmf(self.chain.count("A"), self.chainLength, i))
            likeScoresC.append(binom.pmf(self.chain.count("C"), self.chainLength, i))
            likeScoresG.append(binom.pmf(self.chain.count("G"), self.chainLength, i))
            likeScoresT.append(binom.pmf(self.chain.count("T"), self.chainLength, i))
            
        return likeScoresA, likeScoresC, likeScoresG, likeScoresT
        
            
class ctmarkov(markov):
    """ Continuous-time Markov chain simulator
        
        Attributes: 
        stateSpace = list of characters of the state space. Designed for 4 characters (eg nucleotides)
        chainLength = branch length (time)
        RateMatrix = Q-Matrix of transition probabilities for a continuous markov chain        
    """
    
    
    def ctmcSim(self):
        """ A continuous-time markov chain simulator. 
        
            Takes the 4-state nucleotide State Space and 
            nucleotide Rate Matrix defined above
        """
        
        # current state if sampling from the stationary distribution
        currState = self.discSamp(self.stateSpace, self.margProb(chainLength=100)[0])
        
        # current state if sampling from a uniform distribution
        #currState = self.discSamp(self.stateSpace,[1.0/len(self.stateSpace) for x in self.stateSpace])
        
        self.chain.extend(currState)
        
        # samples the waiting time from the exponential distribution using the appropriate value from the rateMatrix
        while sum(self.times) < self.chainLength: # stops the chain when the branch length exeeds the sum of the times list
            timeSample = (random.expovariate(-(self.Q[self.stateSpace.index(currState), self.stateSpace.index(currState)])))
            self.times.append(timeSample)        
            
            x = 1
            row = self.Q[self.stateSpace.index(currState)] # row from appropriate index from the stateSpace
            # converts that row to a one-dimensional array of transition probabilities
            for i in row:
                if i < 0:
                    x = i
                val = (row / x) * -1
            x = val.tolist() # converts array to list
            x.pop(self.stateSpace.index(currState)) # removes the negative value (was the diagonal)
            x.insert(self.stateSpace.index(currState), 0) # replaces the diagonal with prob of "0"
            #print x

            currState = self.discSamp(self.stateSpace, x)         
            self.chain.extend(currState)
            
        return self.chain, self.times
        
        
    def margProb(self, chainLength = 10):
        """ Calculates the marginal probabilities of a Q-matrix
            There is a built in function to do this. In Scipy.linalg.exm(Q*v)
            pass it Q, which has to be an array
            Q = scipy.array(Q)
        """
        return sp.linalg.expm(self.Q*chainLength)

        # failed attempts:
        #self.margQ = sp.linalg.expm(self.Q*chainLength)
        #return self.margQ

    def multSimsCTMC(self):
        """ runs multiple simulations of the continuous time markov simulator
    
            Attributes:
            numSims = number of simulations to be run. defined in base class
        """
        for x in range(self.numSims): 
            run = self.ctmcSim()
            self.chain.append(run[0])
        return self.chain

    def stateProb(self):
        """a work in progress
        """
        sum = 1
        for i in self.times:
            sum *= i
        for i in self.chain:
            x = self.stateSpace.index[i]
            
            if i > 0:
                sum *= i
        
            


# stuff to pass the ctmarkov class in an example
simulation1 = ctmarkov(stateSpace = ("A", "C", "G", "T"))
print simulation1

""" IT'S WORKING!!!!
    This is from a continuous time marcov chain simulator
(['C', 'G', 'A', 'G', 'A', 'T', 'A', 'C', 'T', 'G', 'T'], 
[3.6766402766957738, 1.0117011882308808, 0.28369755057830681, 0.78290088271085601, 
0.0084820317566074652, 0.28676516298665916, 0.13703870620384057, 5.4660886053329039])
"""
