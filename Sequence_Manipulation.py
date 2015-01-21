# -*- coding: utf-8 -*-
"""
Created on Tue Jan 20 11:27:04 2015

@author: Oscar Johnson github.com/henicorhina

Assignment 1 for Computational Phylogenetics at LSU
"""
#! /usr/bin/env python

# Will manipulate a sequence of DNA nucleotides. Converts to the RNA equivalent, the reverse complement of the RNA strand, and the amino acid sequence


# this string is a nucleotide sequence
# The sequence lenth was not divisible by 3, so i deleted the last two codons.
mySequence = "aaaagctatcgggcccataccccaaacatgttggttaaaccccttcctttgctaattaatccttacgctatctccatcattatctccagcttagccctgggaactattactaccctatcaagctaccattgaatgttagcctgaatcggccttgaaattaacactctagcaattattcctctaataactaaaacacctcaccctcgagcaattgaagccgcaactaaatacttcttaacacaagcagcagcatctgccttaattctatttgcaagcacaatgaatgcttgactactaggagaatgagccattaatacccacattagttatattccatctatcctcctctccatcgccctagcgataaaactgggaattgccccctttcacttctgacttcctgaagtcctacaaggattaaccttacaaaccgggttaatcttatcaacatgacaaaaaatcgccccaatagttttacttattcaactatcccaatctgtagaccttaatctaatattattcctcggcttactttctacagttattggcggatgaggaggtattaaccaaacccaaattcgtaaagtcctagcattttcatcaatcgcccacctaggc"

# this prints the length of the DNA sequence
print "the length of the sequence is", (len(mySequence)), "nucleotides" 

# replacing thymine with uracil to convert the DNA sequence to an RNA sequence
print"\n", "this is the RNA equivalent of the DNA sequence: "
print(mySequence.replace("t","u")) # the replace function

# note that to actually save the replaced codons, you
# need to redefine mySequence each time (so, mySequence = mySequence.append)
# also need to add an intermediate placeholder so that it doesn't just end up with 
# a long string of g's or c's
mySequence = mySequence.replace("c", "x") # x is the placeholder for c
mySequence = mySequence.replace("t", "y") # y is a placeholder for u
mySequence = mySequence.replace("a", "u") # converts a to u
mySequence = mySequence.replace("g", "c") # converts g to c
mySequence = mySequence.replace("x", "g") # gets rid of that placeholder and replaces it with a g
mySequence = mySequence.replace("y", "a") # gets rid of that placeholder and replaces it with an a

print"\n", "this is the reverse complement of the RNA sequence: "
print(mySequence) # prints the new version of mySequence


print"\n", "the 13th codon is: ", mySequence[36:39] # prints the 13th codon. 
# note that this is number 14 in the index.
print"the 14th codon is: ", mySequence[39:42], "\n" # prints the 14th codon.

#Convert to upper case to work with the dictionary
mySequence = mySequence.upper()
#print mySequence

# dictionary to store the Vertebrate mitochondrial genetic code
codonDict = {"UUU":"F","UCU":"S","UAU":"Y","UGU":"C",
"UUC":"F","UCC":"S","UAC":"Y","UGC":"C",
"UUA":"L","UCA":"S","UAA":"*","UGA":"W",
"UUG":"L","UCG":"S","UAG":"*","UGG":"W",
"CUU":"L","CCU":"P","CAU":"H","CGU":"R",
"CUC":"L","CCC":"P","CAC":"H","CGC":"R", 
"CUA":"L","CCA":"P","CAA":"Q","CGA":"R", 
"CUG":"L","CCG":"P","CAG":"Q","CGG":"R",
"AUU":"I","ACU":"T","AAU":"N","AGU":"S",
"AUC":"I","ACC":"T","AAC":"N","AGC":"S",
"AUA":"M","ACA":"T","AAA":"K","AGA":"*",
"AUG":"M","ACG":"T","AAG":"K","AGG":"*",
"GUU":"V","GCU":"A","GAU":"D","GGU":"G",  
"GUC":"V","GCC":"A","GAC":"D","GGC":"G", 
"GUA":"V","GCA":"A","GAA":"E","GGA":"G", 
"GUG":"V","GCG":"A","GAG":"E","GGG":"G"}
#codonDict["CAG"] # testing that the dictionary works
#print codonDict

codonSeq = [] # codon sequence goes in this list
aminoSeq = [] # amino acid sequence goes in this list

# This is supposed to be the function to convert the codons to amino acids, 
# if I can figure out what i'm doing
def seqConvert(sequence): # here I'm defining a function "seqConvert" to convert a nucleotide sequence
    y = 0    # these are the starting index values for the nucleotide sequence string
    z = 3
    codonSeq.append(sequence[y:z]) # appending the first codon to the codonSeq list
    for x in sequence: # Using a for loop to do this same thing for every set of three nucleotides
        y += 3 # increases y and z by 3 to move down the list
        z += 3
        codonSeq.append(sequence[y:z]) # appends the new codon to the codonSeq list
    #aminoSeq = [codonDict("") for item in codonSeq]
    for item in codonSeq:
  #      x  = codonDict.get(item)
        aminoSeq.append(codonDict.get(item))

     
seqConvert(mySequence) # calling the function seqConvert with mySequence

# for some reason, the seqConvert function is appending a bunch of extraneous stuff
# to the end of the lists, so this is to remove that
aminoSeq = aminoSeq[0:(len(mySequence)/3)] 

# I had this line of code because the dictionary was initially lower case, but I left it here just in case.
# codonSeq = [item.lower() for item in codonSeq] # convert list to lower case to work with dictionary

#print aminoSeq # just for testing purposes

# this converts the list aminoSeq to a string and then prints it
string = ''.join(aminoSeq) # converts the aminoSeq list to a string
print "the amino acid sequence is:", "\n", string # prints the new string of amino acids

