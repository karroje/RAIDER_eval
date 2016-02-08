import argparse
import simulation
import Bio
import re
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import random


def main(seq,k,r,name,negative):
    # 1) Read in sequence using biopython
    directory = "../genomes"
    genomeFile = directory+"/"+seq
    record = SeqIO.read(genomeFile, "fasta")

    # 2) Generate a simulated sequence s using go()
    length = len(record)
    newSeq = simulation.go(genomeFile, name, k, length)# Doesn't seem to making a long enough sequence
    newList = list(newSeq)
    print "Kth markov sequence: " + Seq("".join(newSeq))

    
    #Read in chr22.fa.out  (using re.split("\s+", line))
    repeatFile =directory+"/test.fa.out"
    repeatRec = open(repeatFile,'r')
    
    lines = repeatRec.readlines()
    for i in range(3,len(lines)):
        line = lines[i]

        # Reading in first repeat
        line = line.rstrip()
        toks = re.split("\s+", line)

        # Checking if including repeats on negative strand

        if negative is False and toks[9] == 'C':
            System.out.println("negative not included")
        else:
            #Replacing portions of s with portions of chr22.fa as dictated by lines of chr22.fa.out
            start,end = int(toks[6]), int(toks[7])
            newList[start:end] = list(record)[start:end]

        newSeq = Seq("".join(newList))
        print "New sequence: " + newSeq


    #Write simulated sequence to new file using biopython
    output_handle = open(name, "w")
    simSeq = SeqRecord(newSeq)
    SeqIO.write(simSeq, output_handle, "fasta")

if __name__ == "__main__":
    # main()

    parser = argparse.ArgumentParser(description="Generate Simulated Sequence")
    parser.add_argument(\
        '-s', '--seq',action='store',required=False,type = str,default = "test.fa",
        help="sequence fasta file")
    parser.add_argument(\
        '-k',action='store',dest='k',type=int,default=5,
        help="order of Markov chain, default is 5")
    # Will be more applicable later, at the moment, the test repeat file only has two repeats
    parser.add_argument(\
        '-r','--repeat',action='store',dest='repeat',type=str,
        help="The repeat that you are choosing to include in the simulated sequence")
    parser.add_argument(\
        '-n','--name',action='store',dest='name',default = "result.fa",
        help="name of simulated sequence")
    parser.add_argument(\
        '-c','--negative',action='store',dest='negative',type = bool, default = False,
        help = "boolean value for if the repeats on the negative strandare included")

    args = parser.parse_args()

    main(args.seq,args.k,args.repeat,args.name,args.negative)
