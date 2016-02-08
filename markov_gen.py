#!/usr/bin/python2.7
############################################################################
# simulation.py
# Code for generating a simulated chromosome based on a k-degree markov chain and a model sequences
# Written by: Jiajun Wang, John Karro
# Date: March, 2014
# Email: karroje@miamiOH.edu

import re
import sys
import argparse
import random
import unittest
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pickle

fn={
    'A':0, 'a':0,
    'C':1, 'c':1,
    'G':2, 'g':2,
    'T':3, 't':3,
    'N':-1, 'n':-1
}

bases = "ACGTacgt"
baseSet = set(bases)
base2index = {x:y for x,y in zip("acgtnACGTN", [0,1,2,3,-1,0,1,2,3,-1])}

def ProbTuple(Mark,kmo):
    total = 0
    tlist=[0]*4
    for t in range(4):
        total += Mark[(kmo<<2)|t]
    if total == 0 :
        return [0,0,0,0]     
    probsum = 0
    for t in range(4):
        tlist[t] = probsum + Mark[(kmo<<2)|t]/float(total)
        probsum = tlist[t]
    return tlist
    
def pickFromProbTuple(probTuple):
    """Pick a random index based on a probability tuple.  Assumes last element of tuple is 1."""
    r = random.random()
    i = 0
    while r > probTuple[i]:
        i += 1
    return i

#return the index of the first base without 'N' in following it for k length 
def NextIndex(i,k,seq): # Should be having inSeq as a parameter
    """Return the smallest j >= i such that inSeq[j:j+k] contains only bases"""
    j = i
    while j + k < len(seq) and seq[j] not in baseSet:
        j += 1
    if k == 0:
        return j
    while j + k < len(seq):
        p = 1
        while p < k and seq[j+p] in baseSet:
            p += 1
        if p == k:
            return j
        j += p + 1
    return len(seq)

   
def Markov(k,inSeq): # inSeq as parameter?
    """
    Create a k-th order Markov chain based on the distribution
    of inSeq
    """
    Mark= [0]*(4**(k+1))   
    mask=(1<<((k+1)*2))-1
    t=0
    #get first k+1 subsequence
    i=NextIndex(0,k,inSeq)
    for j in range(k+1):
        t=(t<<2) | base2index[inSeq[i+j]] 

    Mark[t] = 1
    #go through the whole seq
    i += k+1     
    while i<len(inSeq) : 
        if inSeq[i] not in bases:
             i=NextIndex(i,k,inSeq)
             if (i+k)>=len(inSeq):
                 break
             t=0
             for j in range(k+1):
                 t = (t<<2) | base2index[inSeq[i+j]]
             i += k+1 
        else:
             t = ((t<<2) | base2index[inSeq[i]]) & mask     
             i += 1             
        Mark[t] += 1   
        
    
    Klist=[[0,0,0,0]]*(4**k)  
    for KMOne in range(4**k):         
        Klist[KMOne] = ProbTuple(Mark,KMOne)
    return Klist       

def MarkovArray(k, inSeq): 
    """
    Create an array such that M[i] is the ith-order for inSeq
    (0 <= i <= k).
    """
    return [Markov(i,inSeq) for i in range(0,k+1)]


def generate_sequence(markov_list, seq_len):
    """Generate a random sequence of length length using the list of k-order Markov chain modeled on sequence seq"""
    k = len(markov_list)-1
    mask = (1 << (k*2)) - 1

    j = 0
    seq = ['']*seq_len
    key = 0
    for i in range(0, seq_len):
        next_tuple = markov_list[j][key]
        if next_tuple[-1] < 0.9999999:
            j = 0
            key = 0
            next_tuple = markov_list[j][key]
        base = pickFromProbTuple(next_tuple)
        key = ((key << 2) | base) & mask
        seq[i] = bases[base]
        j = min(j+1,k)

    return "".join(seq)

def pickle_markov_list(markov_list, file):
    pickle.dump(markov_list, open(file, "wb"))
    

def create_pmck(seqFile, k, output = None):
    s = "N".join([str(r.seq) for r in SeqIO.parse(seqFile, 'fasta')])
    M = MarkovArray(k,s)
    if not ouput:
        output = re.sub("\.(fa|fasta)$", ".pmc%d" % (k), seqFile)
    pickle_markov_list(M, output)

def read_pmck(file):
    return pickle.load(open(file, "rb"))

    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Generate simulated sequence using a k-order Markov chain")
    parser.add_argument('-k', type = int, help = "Order of Markov chain", default = 5)
    parser.add_argument('-p', '--pickle_file', help = "Name of pickled file (default: <seq file>.pmc<k>)")
    parser.add_argument('seq_file', help = "Model sequence file (fasta format)")
    args = parser.parse_args()

    create_pmck(args.seq_file, args.k, args.output)
                
