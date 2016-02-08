#!/usr/bin/python2.7
############################################################################
# simulation.py
# program for generating simulated sequence using a kth order markov chain
# Written by: Jiajun Wang
# Date: Aug. 21, 2013
# Email: wangj29@miamioh.edu

import argparse
import random
import unittest
import sys
from Bio import SeqIO

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
    
#return the index of the first base without 'N' in following it for k length 
def NextIndex(i,k,inSeq): # Should be having inSeq as a parameter
    while i<len(inSeq) and inSeq[i]=='N':
        i+=1
    if i+k>=len(inSeq):
        return i
    
    if 'N' in inSeq[i:i+k+1]:
        while inSeq[i]!='N':
            i+=1
        return NextIndex(i,k,inSeq)
    return i  
   
def Markov(k,inSeq): # inSeq as parameter?
    """Create Markov chain
    Input: k is the degree of the Markov chain
    """
    Mark= [0]*(4**(k+1))   
    mask=(1<<((k+1)*2))-1
    t=0
    #get first k+1 subsequence
    i=NextIndex(0,k,inSeq)
    for j in range(k+1):
        t=(t<<2) | fn[inSeq[i+j]]
    Mark[t] = 1
    #go through the whole seq
    i += k+1     
    while i<len(inSeq) : 
        if inSeq[i]=='N':
             i=NextIndex(i,k,inSeq)
             if (i+k)>=len(inSeq):
                 break
             t=0
             for j in range(k+1):
                 t=(t<<2) | fn[inSeq[i+j]]
             i += k+1 
        else:
             t = ((t<<2) | fn[inSeq[i]])& mask     
             i += 1             
        Mark[t] += 1   
        
    
    Klist=[[0,0,0,0]]*(4**k)  
    for KMOne in range(4**k):         
        Klist[KMOne] = ProbTuple(Mark,KMOne)
    return Klist       

fn={
    'A':0, 'a':0,
    'C':1, 'c':1,
    'G':2, 'g':2,
    'T':3, 't':3,
    'N':-1, 'n':-1
}
gn={
   0:'A',
   1:'C',
   2:'G',
   3:'T'
}


def create_markov_chain(chr_file, k):
    """
    Input: 
    * A .fa file containing a single sequences
    * An integer k specifying the degree of the Markov chain
    Output: A markov chain for use in the simulation
    """
   
    
    pass


def create_simulated_sequence(M, l):
    """
    Input:
    * M: A markov chain (as created by create_markov_chain)
    * l: Length of the sequence
    Output:
    * A simulated sequence of length l
    """

    
    pass

def go(file_name, name, k, length = None):
    """
    Place the main code (everything except for argparser code) in here
    """
    handle = open(file_name)
    seq_record = SeqIO.read(handle,"fasta")
  
    inSeq=str(seq_record.seq).upper()
    if length:
        n = length
    else:
        n = len(inSeq)
    #Zero Markov probabilities list
    total=0
    base_count={}
    for i in range(4):   
        base_count[i]=inSeq.count(gn[i])
        total += base_count[i]        
    plist_zero=ProbTuple(base_count,0)
    
    #generate #0 to #k-1 bases AND kth markov prob list
    sublist = plist_zero
    plist=[]
    key=0
    simSeqStr=''
    for i in range(1,k+1):          
          r,newBase = random.random(),0
          if sublist == [0,0,0,0]:
              sublist = plist_zero        
          while r > sublist[newBase]:
              newBase += 1      
          simSeqStr += gn[newBase]
          key = (key<<2) | newBase     
          plist = Markov(i, inSeq)
          sublist=plist[key]

    #generate the rest sequence using kth Markov prob         
    mask=(1<<(k*2))-1
    for j in range(k,n):
        sublist=plist[key]
        r,newBase = random.random(),0
        if sublist == [0,0,0,0]:
            sublist = plist_zero        
        while r > sublist[newBase]:
            newBase+=1         
        simSeqStr += gn[newBase]       
        key = ((key<<2) | newBase)& mask 

    return simSeqStr

    
    if args.name:
        seqName = name
    else :
        seqName = "SimSeq.fa" 

    return simSeqStr

    #fo=open(seqName,"w+");
    #fo.write(">"+seqName[:seqName.find('.')]+"\n")        
    #fo.write(simSeqStr)
    #fo.close()


class Tester(unittest.TestCase):
    def testProbTuple(self):
        
        """
        1. Create the inputs that will be passed into ProbTuple
        2. Plug inputs into ProbTuple and run ProbTuple
        3. Obtain output from ProbTuple
        4. Compare output from ProbTuple with expected output (look at self.assert)
        """    
        # mark = [0]*4
        # kmo =1 
        # result = ProbTuple(mark,kmo)
        #expectedOutput = 
        # self.assertEquals(result,expectedOutput)
        pass

    def testNextIndex(self):
       i = 1
       k = 0
       inSeq = "ACGT"
       result = NextIndex(i,k,inSeq)
       expectedOutput = len(inSeq) 
       self.assertEquals(result,expectedOutput)
       pass

    def testMarkov(self):
       # k = 5
       # result = Markov(k)
       # expectedOutput = 
       # self.assertEquals(result,expectedOutput)
        pass

if __name__=='__main__':
    parser = argparse.ArgumentParser(description="Generate Simulated Sequence")
    parser.add_argument(\
        '--seq',action='store',required=False,
        help="sequence fasta file")     
    parser.add_argument(\
        '-k',action='store',dest='k',type=int,default=5,
        help="order of Markov chain, default is 5")     
    parser.add_argument(\
        '-l','--length',action='store',dest='length',type=int,
        help="length of sequence, same as the orginal sequence by default") 
    parser.add_argument(\
        '-n','--name',action='store',dest='name',
        help="name of simulated sequence")
    parser.add_argument(\
        '--test',action='store_const',const=True,
        help="name of simulated sequence")
    
    args = parser.parse_args()

    if (args.test == True):
        del sys.argv[1:]
        unittest.main()
    #else:
    #    go(args.seq)
