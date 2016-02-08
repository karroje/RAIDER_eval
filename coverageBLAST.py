import re
import sys
import subprocess
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SeqIO
import time
import os

def createDatabase(db_file, force = True):
    """Create BLAST db out of elements file"""
    if force or not os.path.isfile(db_file + ".nhr"):
        cmd = "makeblastdb -in %s -dbtype nucl" % (db_file)
        p = subprocess.Popen(re.split("\s+", cmd), stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        output, err = p.communicate()
        if err:
            sys.stderr("madkeblastdb error:\n" + err.decode() + "\n")
            sys.stderr('-------------------------')


def invokeBlast(query_file, db_file, out_file = None, e_val = 0.0001):
    """Blast each ancestor against the elements database; return name of output file"""
    print("HERE1")
    output = out_file if out_file else query_file + ".xml"
    blastn_cline = NcbiblastnCommandline(query=query_file, db=db_file, outfmt=5, out=out_file if out_file else query_file + ".xml", evalue = e_val)
    stdout, stderr = blastn_cline()

    #assert not stderr, "BLAST: " + stderr

    return output

def coverageBLAST(query_file, target_file, e_val = 0.0001, isRerun = False):
    """
    INPUT: Two fasta files: query and target
           isRerun: Assume the files already exist and skip rerunning blast.  
                    (Intended for debugging.)
    OUTPUT: Two dictionaries
    * First dictionary: map each sequence id in the query file to x, where
      x is the fraction of bases in the sequence that are part of at least one
      significant local alignment with some sequence in the target file.
    * Second dictionary:  map each sequence id in the target file to x, where
      x is the fraction of bases in the sequence that are part of at least one
      significant local alignment with some sequence in the query file.
    SIDE EFFECTS: Will create several new files for the BLAST database
    REQUIRES: biopython, BLAST 
    """
    query_map = {}
    target_map = {}


    
    # 1) Create blast database
    if not isRerun:
        createDatabase(target_file)

    # 2) blast query sequence against database
    output_file = invokeBlast(query_file, target_file) if not isRerun else (query_file + ".xml")

    # 3) Check covereage of each sequence
    fp = open(output_file)
    blast_records = NCBIXML.parse(fp)

    for blast_obj in blast_records:
        query_id = re.match("(\w+)", blast_obj.query).group(1)
        if not query_id in query_map:
            query_map[query_id] = set()

        for align_obj in blast_obj.alignments:
            target_id = re.match("\S+\s+(\S+)", align_obj.title).group(1)

            if not target_id in target_map:
                target_map[target_id] = set()

            for hsp_obj in align_obj.hsps:
                query_start = hsp_obj.query_start-1
                query_finish = query_start + len([1 for x in hsp_obj.query if x != '-'])
                query_map[query_id] |= set(range(query_start, query_finish))

                subject_start = hsp_obj.sbjct_start-1
                subject_finish = subject_start + len([1 for x in hsp_obj.sbjct if x != '-'])
                target_map[target_id] |= set(range(subject_start,subject_finish))

    # 4) Calculate converage of target sequences
    target_coverage = {}
    query_coverage = {}

    for seq_record in SeqIO.parse(target_file, "fasta"):
        try:
            target_coverage[seq_record.id] = len(target_map[seq_record.id]), float(len(seq_record))
        except:
            target_coverage[seq_record.id] = 0, float(len(seq_record))

    for seq_record in SeqIO.parse(query_file, "fasta"):
        try:
            query_coverage[seq_record.id] = len(query_map[seq_record.id]), float(len(seq_record))
        except:
            query_coverage[seq_record.id] = 0, float(len(seq_record))

    return query_coverage,target_coverage

