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
import matplotlib.pyplot as plt
from scipy.stats import binom
from scipy.stats import norm


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
                matrix = [[0.25, 0.25, 0.25, 0.25], 
                          [0.1, 0.1, 0.1, 0.7], 
                          [0.9, 0.04, 0.02, 0.02], 
                          [0.5, 0.2, 0.1, 0.2]],
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
        self.tProblist = [] # holds the array for the transitional probabilities of a Q-matrix
        self.tProb() # creates a transitional probability matrix
        self.margProb() # creates a marginal probability
    
        
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


    def optimize(self, k = 4, n = 5, pCurr  = random.random(), diff=0.01):
        """ Optimization hill climbing function
        """
        #pCurr = binom.pmf(k, n, pCurr)
        pUp = pCurr + diff
        pDown = pCurr - diff
        binpCurr = binom.pmf(k, n, pCurr)
        binpUp = binom.pmf(k, n, pUp)
        binpDown = binom.pmf(k, n, pDown)
        
        while diff > 0.001:  
                
            if binpCurr < binpUp:
                while binpCurr < binpUp:
                    pCurr = pUp
                    binpCurr = binom.pmf(k, n, pCurr)
                    pUp = pCurr + diff
                    binpUp = binom.pmf(k, n, pUp)
                
            elif binpCurr > binpDown:
                while binpCurr > binpDown:
                    pCurr = pDown
                    binpCurr = binom.pmf(k, n, pCurr)
                    pDown = pCurr - diff
                    binpDown = binom.pmf(k, n, pDown)

            else:
                diff *= 0.5
            
            return pCurr



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
        
        if len(self.stateSpace) != 2.0:
            return "error with state space"
        
        if len(self.matrix) != 2.0:
            return "error with matrix"
            
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
        if len(self.chain) != 0:
            x = raw_input("chain is not empty. are you sure that you want to run the simulator? y/n ").lower()
            if x == "n":
                return "the simulation has been terminated"
            elif x == "y":
                x = x
                        
        endVal = [] # list to hold ending values of the simulations 
        self.chainLength = input("how long do you want each chain in your simulation to be? ")
        for x in range(input("how many times would you like to run the simulator? ")): 
            run = self.dmcSim()
            endVal.append(run[-1])

        print "\n", "the frequency of A is: ", float(endVal.count("A") / len(endVal))
        print "\n", "the frequency of C is: ", float(endVal.count("C") / len(endVal))
        print "\n", "the frequency of G is: ", float(endVal.count("G") / len(endVal))
        print "\n", "the frequency of T is: ", float(endVal.count("T") / len(endVal)), "\n"
    
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
            nucleotide Rate Matrix defined in the base class
        """        
        # current state if sampling from the stationary distribution
        currState = self.discSamp(self.stateSpace, self.margProb(chainLength=100)[0])
        # current state if sampling from a uniform distribution
        self.chain.extend(currState)        
        while sum(self.times) < self.chainLength: # stops the chain when the branch length exeeds the sum of the wait times list
            # samples the waiting time from the exponential distribution using the appropriate value from the QMatrix
            timeSample = (random.expovariate(-(self.Q[self.stateSpace.index(currState), self.stateSpace.index(currState)])))
            self.times.append(timeSample)
            row = self.Q[self.stateSpace.index(currState)] # takes the row from appropriate index from the stateSpace
            x = 1
            # converts that row variable to a one-dimensional array of transition probabilities
            for i in row:
                if i < 0:
                    x = i
                val = (row / x) * -1
            probs = val.tolist() # converts array to list
            probs.pop(self.stateSpace.index(currState)) # removes the negative value (was the diagonal)
            probs.insert(self.stateSpace.index(currState), 0) # replaces the diagonal with probability of 0
            
            if sum(probs) > 1.0: # error check
                print "error. you've made a dumb. your probabilities are not probabilities"
            #print probs

            currState = self.discSamp(self.stateSpace, probs) # sample a new current state and append
            self.chain.extend(currState)
            
            # these are all failed attemp ts at getting the ct chain to work:
            #currState = self.discSamp(self.stateSpace, self.margProb.margQ[0])                     
            #probs = self.Q[self.stateSpace.index(currState)] # Grabbing row associated with currState
            #currState = self.margProb[0][self.discSamp(self.stateSpace,probs)]
        return self.chain, self.times


    def tProb(self):
        """ Creates a transition probability matrix from the Q-Matrix
        """
        
        for val in self.Q:
            # converts that row to a one-dimensional array of transition probabilities
            for i in val:
                if i < 0.0:
                    sums = val / i # find the diagonal and divide by it
                    sums2 = sums * -1.0
                    probs = sums2.tolist() # converts array to list
                    list = val.tolist()
                    probs.pop(list.index(i)) # removes the negative value (was the diagonal)
                    probs.insert(list.index(i), 0) # replaces the diagonal with prob of "0"
                    self.tProblist.append(probs)
        
        # error check
        for i in self.tProblist:
            if sum(i) != 1.0:
                return "error. you've made a dumb. there is probably an error in your Q matrix"
        
        # convert self.tProb to a numpy array
        self.tProblist = np.array(self.tProblist)

        return self.tProblist


    def margProb(self, chainLength = 100):
        """ Calculates the marginal probabilities of a Q-matrix
            There is a built in function to do this. In Scipy.linalg.exm(Q*v)
            pass it Q, which has to be an array
            Q = scipy.array(Q)
        """
        values = sp.linalg.expm(self.Q*chainLength)
        self.margQ = values
        
        return self.margQ

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


    def nucStateProbs(self):
        """ Calculates the probabilities across all waiting 
            times and nucleotide states in the chain
            using the format: 
            P(nuc1)*P(t1)*P(nuc2|nuc1)...*P(1 - cdf(tlast))
            This is useful for stochastic character mapping
        """
        if len(self.chain) == 0:
            return "chain is empty"
        if len(self.tProblist) == 0:
            return "run the tProb function, the transition probability matrix is empty"
        
        list = [] # blank list for probabilities
        # multiply all of the waiting times together, except the last value
        # and convert waiting times to a probability using the pdf
        sum = 1
        for i in self.times[:-1]:
            x = norm.pdf(i)
            sum *= x
            list.append(x)
            
        # tally all of the chain transition probabilities except for the first value
        for a, b in zip(self.chain, self.chain[1:]):            
            curr = self.stateSpace.index(b)
            prev = self.stateSpace.index(a)
            #print curr, prev            
            val = self.tProblist[prev, curr]
            sum *= val
            list.append(val)

        # calculates the probability of the first state in the chain
        # using the stationary frequency
        row = self.margProb(chainLength=100)[0]
        firstProb = row[self.stateSpace.index(self.chain[0])]
        sum *= firstProb
        list.append(firstProb)

        # calculate the last waiting time: 1 - cdf(t last)
        y = 1.0 - (norm.cdf(self.times[-1]))
        sum *= y        
        list.append(y)

        return sum, list



    def contLike(self):
        """ Continuous-time markov chain likelihood calculator. 
            A work in progress
            2 * difference in log(likelihood) = chi-squared
            degrees of freedom
        """
        list = []
        for x in self.nucStateProbs()[1]:
            
            list.append()


    def optimize(self, diff=0.01):
        """ Function to find the maximum likelihood for a continuous time markov chain
            An optimization hill climbing function
            Attributes:
                diff = value steps by which the function goes up or down to find the maximum
        """
        if len(self.chain) == 0: # error check
            return "chain is empty, run ctmcSim"

        pCurr = self.chainLength
        values = [] # list to hold the intermediate values while optimizing
        # first and last are the stateSpace index values for the first and last states of the chain
        first = self.stateSpace.index(self.chain[0])
        last = self.stateSpace.index(self.chain[-1])

        pUp = pCurr + diff
        pDown = pCurr - diff
        # sample the appropriate transition value from the marginal probability matrix
        bpCurr = self.margProb(chainLength=pCurr)[first, last]
        bpUp = self.margProb(chainLength=pUp)[first, last]
        bpDown = self.margProb(chainLength=pDown)[first, last]
        values.append(bpCurr)

        while diff > 0.001: # runs the function until the diff drops below 0.001. Adjust if needed

            if bpCurr < bpUp: # going up!
                while bpCurr < bpUp: # runs the function as long as you're still going up.
                    pCurr = pUp
                    pUp = pCurr + diff
                    bpCurr = self.margProb(chainLength=pCurr)[first, last]
                    bpUp = self.margProb(chainLength=pUp)[first, last]
                    values.append(bpCurr)
                    if len(values) > 1000:
                        return bpCurr, values

            elif bpCurr > bpDown: # going down!
                while bpCurr > bpDown:
                    pCurr = pDown
                    bpCurr = self.margProb(chainLength=pCurr)[first, last]
                    values.append(bpCurr)
                    pDown = pCurr - diff
                    bpDown = self.margProb(chainLength=pDown)[first, last]
                    if len(values) > 1000:
                        return bpCurr, values

            else:
                diff *= 0.5

            plt.plot(values)
            plt.show

            return bpCurr, values


# stuff to pass the ctmarkov class in an example
simulation1 = ctmarkov(stateSpace = ("A", "C", "G", "T"))
#print simulation1.ctmcSim()

""" IT'S WORKING!!!!
    This is from a continuous time marcov chain simulator
(['C', 'G', 'A', 'G', 'A', 'T', 'A', 'C', 'T', 'G', 'T'], 
[3.6766402766957738, 1.0117011882308808, 0.28369755057830681, 0.78290088271085601, 
0.0084820317566074652, 0.28676516298665916, 0.13703870620384057, 5.4660886053329039])
"""
