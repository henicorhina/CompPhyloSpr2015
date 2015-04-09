"""
Exercise 6 - Creating and Using Node and Tree Classes
@author: henicorhina
Below is the beginning of a Node class definition and a simple example of how
to link nodes to form a tree. Use this as a springboard to start thinking about:
- What other attributes of a Node might we like to store?
- How do we define a Tree class? What attributes should it have?
- Can you write a function to print out a parenthetical tree string 
   (e.g., ((spA,spB),spC)) if the only argument passed to the function is a
   root node? This will require recursion.
"""
from marcovObjects import ctmarkov
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp

# ---> Defining Node and Tree classes <---

class Node:
    def __init__(self,name="",parent=None,children=None):
        self.name = name
        self.parent = None
        if children is None:
            self.children = []
        else:
            self.children = children
        
        
        
# ---> Creating and linking nodes <---
 
# Creating nodes to build this simple three-taxon tree: ((spA,spB),spC)
       
#  spA     spB  spC
#    \    /     /
#     \  /     /
#      \/     /
#       \    /
#        \  /
#         \/
#         |
"""
# Define the root node to start. It currently has no parents or children.
root = Node("root") 

# Define a node for species C. It is a direct descendant of the root.
spC = Node("SpeciesC",parent=root)
root.children.append(spC)   # Adds spC as a child of the root

# Define a node for the ancestor of species A and B, descending from the root.
ancAB = Node("ancAB",parent=root)
root.children.append(ancAB)
spA = Node("SpeciesA",parent=ancAB) # Creates spA with ancAB as its parent.
spB = Node("SpeciesB",parent=ancAB) # Creates spB with ancAB as its parent.
ancAB.children.append(spA)
ancAB.children.append(spB)


print("ancAB's children: ")
for child in ancAB.children:
    print child.name
    
print("")
print("root's children: ")
for child in root.children:
    print child.name
"""
# Play around with nodes and see if you can build more complicated trees!


# Eventually, we will want to create a Tree class, where a parenthetical tree
# string is passed as an argument to the constructor and it automatically creates
# all the nodes and links them together. Start thinking about how to do that.


# Let's go ahead and define a Tree object that houses all these nodes and 
# organizes methods associated with them.



class Tree(object):
    """
    Defines a class of phylogenetic tree, consisting of linked Node objects.
    """
    
    def __init__(self):
        """
        The constructor really needs to be more flexible, but for now we're 
        going to define the whole tree structure by hand. This just uses
        the same statements we used above. By next Thurs (3/19), see if you can
        write a constructor that takes a parenthetical tree as its argument and 
        builds the corresponding tree in memory. 
        """
        self.root = Node("root") 
        self.spC = Node("SpeciesC",parent=self.root)
        self.root.children.append(self.spC)
        self.ancAB = Node("ancAB",parent=self.root)
        self.root.children.append(self.ancAB)
        self.spA = Node("SpeciesA",parent=self.ancAB)
        self.spB = Node("SpeciesB",parent=self.ancAB)
        self.ancAB.children.append(self.spA)
        self.ancAB.children.append(self.spB)
        # Now, let's add branch lengths to our Node objects (remember, these fields
        # can be added arbitrarily in Python). In the future, we should probably include
        # branch lengths in the Node constructor.
        self.spA.brl = 0.1
        self.spB.brl = 0.1
        self.spC.brl = 0.2
        self.ancAB.brl = 0.1
        self.root.brl = 0
        # We're also going to add lists to each node that will hold simulated
        # sequences.
        self.spA.seq = []
        self.spB.seq = []
        self.spC.seq = []
        self.ancAB.seq = []
        self.root.seq = []
        self.totalTreeLength = 0 # counter for total tree length calculator
        self.setModels(self.root)

    # Write a recursive function that takes the root node as its only argument and
    # prints out all the names of the terminal nodes in the tree. Due next Tues (3/17).
 
    def printNames(self,node):
        """
        A method of a Tree object that will print out the names of its
        terminal nodes.
        """
        if len(node.children) == 0: # checks if the node has no children
            print node.name # if no children, print your name
        else:        
            for child in node.children: # otherwise, call printNames again
                self.printNames(child)
        # I successfully printed out all of the child names calling the root as node.                
                
    """ A failed attempt at getting the printNames function to work
    def printNames2(self, node):
            
            A method of a Tree object that will print out the names of its
            terminal nodes.
            
            # creates a list to hold the names of the tip nodes        
            tip_nodes = []
            if len(node.children) > 0: # checks to see if the node given is not a tip node
                for child in node.children: # for every child of that internal node
                    if len(child.children) > 0: # checks if the child is not a tip node
                        tip_nodes.append(self.printNames2(child)) # use the function in each child of the argument node
                        # appends the name of the child node if it is a tip node
                    else:
                        tip_nodes.append(child.name)
            else:
                tip_nodes.append(node.name) # appends the name of the argument node if it is a tip node
            return tip_nodes
    """ 
        
    
    # Write a recursive function to calculate the total tree length (the sum of
    # all the branch lengths). Again, the root node of a tree should be the only 
    # argument the first time this function is called.
    
    def treeLength(self,node):
        """
        A method to calculate and return total tree length.
        Pass it the root of the tree. 
        """
        if len(node.children) > 0: # checks if the node has children
            self.totalTreeLength += node.brl
            for child in node.children: 
                self.treeLength(child) # call treeLength again
        else:        
            self.totalTreeLength += node.brl
        
        return self.totalTreeLength
        

    # Write a recursive function that takes the root node as one of its arguments
    # and prints out a parenthetical (Newick) tree string. Due next Tues (3/17).
    
    def newick(self,node):
        """
        A method of a Tree object that will print out the Tree as a 
        parenthetical string (Newick format).
        There is still a random comma after tips that have no sister, which I 
        am working on getting rid of...
        """
        newickString = "(" # the start of the string
        
        if len(node.children) == 0: # checks if the node is a tip
            return node.name + ":" + str(node.brl)# returns name : branch length

            #if newickString[-1] is not "(":
            #    return "," +node.name + ":" + str(node.brl) # returns name : branch length
            #else:
            #    return node.name + ":" + str(node.brl) # returns name : branch length
                
        else:
            for child in node.children: 
                newickString += self.newick(child) # runs the function for all the children
                if node.children[-1] == child: # supposedly checks if the previous child is a child
                    newickString = newickString 
                else:
                    newickString += "," 
            if node.brl == 0: # checks if the node is the root
                newickString += ")" 
            else:
                newickString += ")"+ str(node.brl) # otherwise, adds the brl for the ancestor
            
            return newickString
        
        
    # Now, let's write a recursive function to simulate sequence evolution along a
    # tree. This amounts to simply simulating evolution along each branch 
    # from the root towards the tips. We'll need to use our ctmc class for setting the 
    # conditions of our simulation, which is why we imported it above our tree 
    # class definition. In this case, we've stored the definition of our ctmc 
    # class in a separate file (ctmc.py) to keep our tree code compact.
    # Now, let's add a ctmc object to each internal node in our tree (except the
    # root). Again, it would be best to add the ctmcs as part of the Node
    # constructor, if we know that we'll be simulating data.
        
    def setModels(self,node):
        """
        This method of a Tree object defines a ctmc object associated with all
        nodes that have a branch length (i.e., all but the root).
        """


    def simulate(self,node):
        """
        This method simulates evolution along the branches of a tree, taking
        the root node as its initial argument.
        """
             
    def printSeqs(self,node):
        """
        This method prints out the names of the tips and their associated
        sequences as an alignment (matrix).
        """

        
