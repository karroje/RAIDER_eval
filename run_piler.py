import sys
import subprocess
import re
import os
import glob

"""Helper script for running piler"""

def file_base(file):
    return os.path.basename(file)

def file_dir(file):
    return file.rstrip(file_base(file)).rstrip("/")

if __name__ == "__main__":
    seq_file = sys.argv[1];
    piler_dir = sys.argv[2]
    output_file = sys.argv[3]

    seq_base = file_base(seq_file);
    seq_dir = file_dir(seq_file);


    gff_file = re.sub("\.fa$|\.fasta$", ".gff", seq_base)
    trs_file = re.sub("\.fa$|\.fasta$", ".trs.gff", seq_base)
    fam_dir =  re.sub("\.fa$|\.fasta$", "", seq_base).upper()

    try:
        os.mkdir(piler_dir + "/" + fam_dir)
    except FileExistsError:
        pass

    cmd1 = "./pals -self {seq_file} -out {piler_dir}/{gff_file}"
    cmd2 = "./piler2 -trs {piler_dir}/{gff_file} -out {piler_dir}/{fam_dir}/{trs_file}"
    cmd3 = "./piler2 -trs2fasta {piler_dir}/{fam_dir}/{trs_file} -seq {seq_file} -path {piler_dir}/{fam_dir}"
    cmd = "; ".join([cmd1, cmd2, cmd3]).format(seq_file=seq_file, gff_file=gff_file, trs_file = trs_file, fam_dir=fam_dir, piler_dir=piler_dir)

    print("\n".join(cmd.split("; ")))
    exit(1)
    subprocess.call(cmd, shell = True)

    cmd = ""
    for fam in glob.glob(piler_dir + "/" + fam_dir + "/*"):
        cmd += "./muscle -in {fam} -out {fam}.aligned -maxiters 1 -diags1; ".format(fam=fam)
        cmd += "./piler2 -cons {fam}.aligned -out {fam}.cons -label {fam}; ".format(fam=fam)
    cmd += "cat {piler_dir}/{fam_dir}/*.cons > {piler_dir}/{output_file}".format(fam_dir = fam_dir, piler_dir = piler_dir, output_file = output_file)
    #print("\n".join(cmd.split("; ")))
    subprocess.call(cmd, shell = True)

    

    
