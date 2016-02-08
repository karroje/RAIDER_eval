"""Given a seq.fa file and the corresponding seq.fa.out file, create a BLAST database of the repeat sequences."""
import subprocess
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import time
import os
import sys
import re

filter_words = {'Simple', 'Satellite', 'complexity', 'telo'}

def createDatabase(db_file, force = True):
    """Create BLAST db out of elements file"""
    if force or not os.path.isfile(db_file + ".nhr"):
        cmd = "makeblastdb -in %s -dbtype nucl" % (db_file)
        p = subprocess.Popen(re.split("\s+", cmd), stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        output, err = p.communicate()
        if err:
            sys.stderr("madkeblastdb error:\n" + err.decode() + "\n")
            sys.stderr('-------------------------')

def create_target(seq_file, rm_file, out_file):
    """from seq_file and rm_file create an .fa file of repeat sequences in out_file"""
    assert rm_file.endswith(".fa.out")

    D = {}
    chr_seq = SeqIO.read(seq_file, 'fasta').seq

    L = []
    fp = open(rm_file)
    fp.readline()
    fp.readline()
    for line in fp:
        if line.rstrip():
            A = re.split("\s+", line.strip())
            start = int(A[5])
            stop = int(A[6]) - 1
            location = "{chr}:{start}-{stop}".format(chr=A[4], start=start, stop=stop)
            #if any([x in A[10] for x in filter_words]):
            #    continue
            id = "{class_name}|{family}|{location}".format(class_name=A[9], family=A[10], location=location)
            D[id] = 0 if id not in D else D[id] + 1
            id += " " + str(D[id])
            description = ""
            R = SeqRecord(seq = chr_seq[start:stop], id = id, description=description)
            L.append(R)

    SeqIO.write(L, out_file, "fasta")

if __name__ == "__main__":
    seq_file = sys.argv[1]
    rm_file  = sys.argv[2]
    output_file = sys.argv[3]

    create_target(seq_file, rm_file, output_file)

    createDatabase(output_file, True)
    
