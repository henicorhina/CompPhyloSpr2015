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

# use this sequence for converting directly from DNA to AA
mySequence2 = "aaaagctatcgggcccataccccaaacatgttggttaaaccccttcctttgctaattaatccttacgctatctccatcattatctccagcttagccctgggaactattactaccctatcaagctaccattgaatgttagcctgaatcggccttgaaattaacactctagcaattattcctctaataactaaaacacctcaccctcgagcaattgaagccgcaactaaatacttcttaacacaagcagcagcatctgccttaattctatttgcaagcacaatgaatgcttgactactaggagaatgagccattaatacccacattagttatattccatctatcctcctctccatcgccctagcgataaaactgggaattgccccctttcacttctgacttcctgaagtcctacaaggattaaccttacaaaccgggttaatcttatcaacatgacaaaaaatcgccccaatagttttacttattcaactatcccaatctgtagaccttaatctaatattattcctcggcttactttctacagttattggcggatgaggaggtattaaccaaacccaaattcgtaaagtcctagcattttcatcaatcgcccacctaggc"

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

mySequence = mySequence[::-1] # reverses the string

print"\n", "this is the reverse complement of the RNA sequence: "
print(mySequence) # prints the new version of mySequence


print"\n", "the 13th codon is: ", mySequence[36:39] # prints the 13th codon. 
# note that this is number 14 in the index.
print"the 14th codon is: ", mySequence[39:42], "\n" # prints the 14th codon.

#Convert to upper case to work with the dictionary
mySequence = mySequence.upper()
#print mySequence

mySequence2 = mySequence2.upper()

# dictionary to store the Vertebrate mitochondrial genetic code
codonDictRNA = {"UUU":"F","UCU":"S","UAU":"Y","UGU":"C",
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
#codonDictRNA["CAG"] # testing that the dictionary works
#print codonDictRNA

codonDictDNA = {"TTT":"F","TCT":"S","TAT":"Y","TGT":"C",
"TTC":"F","TCC":"S","TAC":"Y","TGC":"C",
"TTA":"L","TCA":"S","TAA":"*","TGA":"W",
"TTG":"L","TCG":"S","TAG":"*","TGG":"W",
"CTT":"L","CCT":"P","CAT":"H","CGT":"R",
"CTC":"L","CCC":"P","CAC":"H","CGC":"R", 
"CTA":"L","CCA":"P","CAA":"Q","CGA":"R", 
"CTG":"L","CCG":"P","CAG":"Q","CGG":"R",
"ATT":"I","ACT":"T","AAT":"N","AGT":"S",
"ATC":"I","ACC":"T","AAC":"N","AGC":"S",
"ATA":"M","ACA":"T","AAA":"K","AGA":"*",
"ATG":"M","ACG":"T","AAG":"K","AGG":"*",
"GTT":"V","GCT":"A","GAT":"D","GGT":"G",  
"GTC":"V","GCC":"A","GAC":"D","GGC":"G", 
"GTA":"V","GCA":"A","GAA":"E","GGA":"G", 
"GTG":"V","GCG":"A","GAG":"E","GGG":"G"}



codonSeq = [] # codon sequence goes in this list
aminoSeq = [] # amino acid sequence goes in this list
    
# This is supposed to be the function to convert the codons to amino acids, 
# if I can figure out what i'm doing
def seqConvert(sequence):
    """
    seqConvert will take a nucleotide sequence and convert it to an amino acid sequence
    and prints it to the screen
    user inputs whether the sequence is a DNA sequence ("y") or an RNA sequence
    """    
    aminoSeq = []
    y = 0    # these are the starting index values for the nucleotide sequence string
    z = 3
    codonSeq.append(sequence[y:z]) # appending the first codon to the codonSeq list
    for x in sequence: # Using a for loop to do this same thing for every set of three nucleotides
        y += 3 # increases y and z by 3 to move down the list
        z += 3
        codonSeq.append(sequence[y:z]) # appends the new codon to the codonSeq list
    if nuc == "dna":
        for item in codonSeq:
            aminoSeq.append(codonDictDNA.get(item)) # "gets" the item from the dictionary and appends to the amino acid sequence
    else:
        for item in codonSeq:
            aminoSeq.append(codonDictRNA.get(item)) # "gets" the item from the dictionary and appends to the amino acid sequence
  
    # for some reason, the seqConvert function is appending a bunch of extraneous stuff
    # to the end of the lists, so this is to remove that
    
    aminoSeq = aminoSeq[0:(len(mySequence)/3)] 

    # removes the * from the amino acid sequence
    while "*" in aminoSeq:
        aminoSeq.remove("*")
    while "None" in aminoSeq:
        aminoSeq.remove("None")
        
    # this converts the list aminoSeq to a string and then prints it
    AAstring = ''.join(aminoSeq) # converts the aminoSeq list to a string
    
    if nuc == "dna":        
        print "the amino acid sequence from the DNA sequence is:", "\n", AAstring # prints the new string of amino acids

    else:
        print "the amino acid sequence from the reverse complement RNA sequence is:", "\n", AAstring # prints the new string of amino acids

nuc = raw_input("do you want to convert the DNA or RNA sequence? (enter DNA or RNA): ").lower()    
print "\n"
if nuc == "dna":
    seqConvert(mySequence2) # calling the function seqConvert with the DNA sequence

else:
    seqConvert(mySequence) # calling the function seqConvert with the RNA sequence

# Genbank BLAST came back as Pseudacris: http://blast.ncbi.nlm.nih.gov/Blast.cgi#alnHdr_594593098
