import argparse
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import random
import sys
import re
import os.path

import markov_gen


def parse_params(args):
    parser = argparse.ArgumentParser(description = "Generate simulated chromosome")
    # parser.add_argument('-c', '--cutoff', type = int, help = "Limit model to first c non-N bases")
    parser.add_argument('-k', type = int, help = "Order of Markov chain", default = 5)
    parser.add_argument('-s', '--seed', '-rng_seed', dest = 'seed', type = int, help = "RNG seed", default = None)
    parser.add_argument('-n', '--negative_strand', action = "store_true", help = "Use repeats on negative string", default = False)
    parser.add_argument('--family_file', help = "List of repeat families to use", default = None)
    parser.add_argument('-m', '--mask', action = "store_true", help = "Turn masking on (all repeats printed as lower case).", default = False)
    parser.add_argument('--mc', '--mc_file', dest = 'mc_file', help = "MC File (by default -- look in local directory; generates if not found).", default = None)
    parser.add_argument('-S', '--suppress_pmck', action = "store_true", help = "Suppress the generation of a .pmc<k> file to store the markov chain for re-use")
    parser.add_argument('--mi', '--max_interval', dest = "max_interval", type = int, help = "Maximum allowed length of interval between repeats; -1 value (default) means no maximum", default = -1)
    parser.add_argument('--mi2', '--min_interval', dest = "min_interval", type = int, help = "Min allowed length of interval between repeats", default = 0)
    
    parser.add_argument('--rn', '--retain_n', dest = "retain_n", action = 'store_true', help = "If used, will use the whole chromosome.  Otherwise, cuts of Ns at either end.", default = False)
    parser.add_argument('--nr', '--num_repeats', dest = 'num_repeats', type = int, help = "Specify the number of repeats.  Simulation will terminate either 1000 bases or max interval bases past the nth instance of a repeat (excluding any other repeats in that range).", default = None)
    parser.add_argument('-l', '--max_length', dest = 'max_length', type = int, help = "Maximum allowed length of simulated sequence.", default = None)
    parser.add_argument('--lc', '--low_complexity', dest = 'low_complexity', action = 'store_true', help = "Keep low complexity and simple repeats (kept by default)", default = False)
    parser.add_argument('--rb', '--rep_base', dest = 'rep_base', help = "Replace each TE with a ful copy of its ancestral seqeunce in the specified RepBase file", default = None)
    parser.add_argument('-f', '--family_min', dest = "family_min", type = int, help = "Number of elements per family", default = 2)
    parser.add_argument('--nf', '--num_family', dest = 'num_family', type = int, help = "Number of families", default = None)

    #parser.add_argument('-o', '--output', help = "Output file (Default: replace chomosome file \".fa\" with \".sim.fa\")")
    parser.add_argument("seq_file", help = "Sequence file (must be .fa)")
    parser.add_argument("repeat_file", help = "RepeatMasker file (.fa.out)")
    parser.add_argument("output", help = "Output file")
    return parser.parse_args(args)


def nextRepeat(rpt_file, use_negative = True, S = {}, E = {}, I = {}):
    """Generator: each invokation returns the chromosome, start, finish, starand, 
    and family for the next repeat of the repeatmasker .fa.out files.  S, if not empty,
    is a filter for which repeats to use."""
    fp = open(rpt_file)
    fp.readline()
    fp.readline()
    for line in fp:
        if line.rstrip():
            A = re.split("\s+", line.strip())
            chr, start, finish, strand, family, rpt_class, rpt_id = A[4], int(A[5])-1, int(A[6]), A[8], A[9], A[10], A[14]
            if strand == '-' and not use_negative:
                continue
            if S and any([s in family for s in S]):
                continue
            if E and any([e in rpt_class for e in E]):
                continue
            if I and not family in I:
                continue
            if (strand == '+' or use_negative) and ((family in S) or not S) and not (rpt_class in E):
                yield chr, start, finish, strand, family, rpt_class, int(rpt_id)

# fa_out_header: The fixed header lines for the .fa.out file
fa_out_header = "\tSW\tperc\tperc\tperc\tquery\tposition in query\tmatching\trepeat\tposition in  repeat\n\tscore\tdiv.\tdel.\tins.\tsequence\tbegin\tend\t(left)\trepeat\tclass/family\tbegin\tend (left)\tID\n"
# fa_out_template: A template for creating lines for the .fa.out file.
fa_out_template = "\t0\t0\t0\t0\t{chr}\t{start}\t{finish}\t({left})\t{strand}\t{family}\t{rpt_class}\t0\t0\t(0)\t{rpt_id}\n"
def generate_chromosome(seq, markov_list, chr_start, chr_finish, rpt_gen, mask = False, max_interval = -1, min_interval = 0,num_repeats = None, max_length = None, limiting_chr = None, rep_base_hash = None):
    """
    Generate a syntehtic sequence with real repeats:
    * seq: A sequence (as a string).
    * markov_list: List of the k+1 i-th order markov chains (from the markov_gen module).
    * start/finish: Defined the coordinates of our actual template sequence.  (We are ignoring anything that occurs before/faster.
    *               Allows us to cut of a prefix and/or suffix.
    * rpt_gen: A generating function returning the repeat information (created by nextRepeat)
    * mask: If true, all repeats will be lower-case.  Otherwise, upper case.)
    * max_interval: Maximum inter-repeat length.
    * min_interval: Minimum allowed length of a sequence between repeats.  If two repeats are closer than this,
    *             extend the length.
    * max_interval: Minimum allowed length of a sequence between repeats.  If two repeats are closer than this,
    *             cut the length.
    """
    last_end = chr_start
    if max_interval == -1:
        max_interval = len(seq)

    sim_seq = ""          # Simulated sequence
    fa_out = []           # Hold the new .fa.out file contents (by line)

    rpt_count = 0         # Count of repeats (so we can quit when we reach num_repeats, if applicable)

    for chr, start, finish, strand, family, rpt_class, rpt_id in rpt_gen:
        if limiting_chr and chr not in limiting_chr:    # Skip if we are on the wrong chromsome
            continue

        if start >= chr_finish:     # Quit if we have gone past the allowed range (repeats are assumed to be sorted by start)
            break
        
        if start < chr_start or finish > chr_finish:   # Skip if we are outside the allowed range
            continue

        if start < last_end:      # Skip if this repeat overlapped the last one
            continue

        rpt_count += 1

        # Add the next inter-TE sequence
        inter_seq_len = max(min_interval, min(start - last_end, max_interval))
        sim_seq += markov_gen.generate_sequence(markov_list, inter_seq_len)
        
        # Add the next sequence
        if rep_base_hash:
            rpt_seq = rep_base_hash[family]
        else:
            rpt_seq = seq[start:finish]

        fa_out.append([chr, len(sim_seq)+1, len(sim_seq) + len(rpt_seq), strand, family, rpt_class, rpt_id])   # Coords adjusted for biologist notation
        sim_seq += rpt_seq.lower() if mask else rpt_seq.upper()

        
        if rpt_count == num_repeats:
            break

        last_end = max(last_end, finish)

    # Add final sequence on
    final_seq_len = max(min_interval, min(chr_finish - last_end, max_interval))
    sim_seq += markov_gen.generate_sequence(markov_list, inter_seq_len)
    
    sim_seq_len = len(sim_seq)
    fa_out_str = fa_out_header
    for chr, start, finish, strand, family, rpt_class, rpt_id in fa_out:
        fa_out_str += fa_out_template.format(chr=chr, start=start, finish=finish, left = sim_seq_len - finish, strand=strand, family=family, rpt_class=rpt_class, rpt_id=rpt_id)


    return sim_seq, fa_out_str

bases = set("ACGTacgt")
def loadSeqAndChain(seq_file, k, suppress_save = False, mc_file = None, retain_n = False):
    """Load the sequence and the Markov Chain List.
    Load the MC list from a file if it exists.  If not, create the chain
    and save it to the file for the next use (skip the save if suppressed).
    Parameters:
    * seq_file: The sequence file.
    * k: The order of the markov chain.
    * suppress_save: Boolean.  If true, don't save the generated MC file.  (Can't imagine why we would want this.)
    * mc_file: The name of the mc_file to use.  (Derive from seq_file if not provided.)
    * retrain_n: If false, we will be cutting of the largest possible N* prefix and suffix.
    Return: A tuple:
    1. The chromosome sequence.
    2. The markov chain
    3. Where we will start in the template sequence (in case a prefix has been removed).
    4. Where we will end in the templace sequence (in case a suffix has been removed).
    """

    template_seq = str(SeqIO.read(seq_file, 'fasta').seq)

    # Cut out all the maximul prefix and suffix of ambiguity codes -- which will have no effect on the Markov chain construction.
    start, finish = 0, len(template_seq)
    if not retain_n:   # Cut down the chromsome to the first real base at each end -- eliminate trailing Ns.
        while template_seq[start] not in bases: start += 1
        while template_seq[finish-1] not in bases: finish -= 1
    
    mc_file = re.sub("\.(fa|fasta)$", ".pmc%d" % (k), seq_file) if mc_file is None else mc_file
    if os.path.exists(mc_file):
        markov_list = markov_gen.read_pmck(mc_file)
    else:
        markov_list = markov_gen.MarkovArray(k, template_seq)
        if not suppress_save:
            markov_gen.pickle_markov_list(markov_list, mc_file)

    return template_seq, markov_list, start, finish
                        

def readRepBase(file):
    return {R.id:"".join([x for x in str(R.seq) if x.upper() in {'A', 'C', 'G', 'T'}]) for R in SeqIO.parse(file, 'fasta')}


low_complexity = {'Low_complexity', 'Simple', 'Satellite'}
def select_families(repeat_file, f, num_fams, use_3prime, filter_set, toss_low, rep_base_hash):
    """Used to select those families that have at least f members on the template chromosome.
    Parameters:
    * repeat_file: the .fa.out file
    * f minimum number of allowed instances in a family.
    * num_fams: Number of families to be choosen
    * use_3prime: if false, ignore instances on the 3' strand
    * filter_set: families that should be ignored
    * toss_low: if true, ignore low-complexity families
    * rep_base_hash: a hash table mapping family names to their rep_base sequences
    Returns: 
    * List of families chosen
    """
    C = {}   # Family name -> count
    for T in nextRepeat(repeat_file, use_3prime, filter_set, E = low_complexity if toss_low else {}):
        if rep_base_hash and  not T[4] in rep_base_hash:
            continue
        if T[4] in C:
            C[T[4]] += 1
        else:
            C[T[4]] = 1

    L = [k for k in C if C[k] >= f]
    if num_fams == None:
        return L

    if num_fams > len(L):
        sys.stderr.write("Not enough families for f\n")
        exit(1);
        
    return L[:num_fams]

    
def create_chromosome_file(seq_file, repeat_file, output_file, k = 5, use_3prime = True, filter_file = "rpt_list.txt", mask = False, seed = None, suppress = False, max_interval = -1, min_interval = 0, retain_n = False, num_repeats = None, max_length = None, toss_low = False, rep_base = None, f = 1, num_fams = None):
    """
    Create a simualted chrosome with real repeat sequences from a chromsoe file.
    Parameters:
    * seq_file: fasta <seq>.fa, file containing the template sequence.
      -- Assumed to exist a file <seq>.fa.out containing the repeatmasker annotations.
    * k: Use a k-order markov chain.  There must exists a markov chain file <seq>.pmc<k>.
    * output_file: Fasta file to print sequence to.
    * use_3prime: If false, only sequence on the 5' strand will be used.  Default: True
    * filter_file: A list of the repeats that should be used.  If empty: all repeats.  Default: "rpt_list.txt"
    * mask: If true: copied repeats will be in lower case.  Default: False
    * seed: RNG seed
    """
    if not output_file.endswith(".fa"):
        output_file += ".fa"
    random.seed(args.seed)

    # First: load in the template sequence, markov chain, and the start/end coords of what we are using.
    template_seq, markov_list, chr_start, chr_finish  = loadSeqAndChain(args.seq_file, args.k, suppress, args.mc_file, args.retain_n)

    # Read in the set of families to be ignored
    filter_set = {y.strip() for line in open(filter_file) for y in re.split("\s+", line.rstrip())} if filter_file else {}

    # Read in the RepBase sequence: maps name -> RepBaseSequence
    rep_base_hash = readRepBase(rep_base) if rep_base else None   # Hash of repeats ID -> sequences)

    # Pick which families we are using.
    selected = select_families(repeat_file, f, num_fams, use_3prime, filter_set, toss_low, rep_base_hash)

    # Create a sequence generator
    rpt_gen = nextRepeat(repeat_file, use_3prime, filter_set, E = low_complexity if toss_low else {}, I = selected)

    # Create the simulated sequence
    simulated_sequence, fa_out = generate_chromosome(seq = template_seq, markov_list = markov_list, chr_start = chr_start, chr_finish = chr_finish, rpt_gen = rpt_gen, mask = mask, max_interval = max_interval, min_interval = min_interval, num_repeats = num_repeats, max_length = max_length, rep_base_hash = rep_base_hash)
    
    # Write output to file        
    SeqIO.write([SeqRecord(seq = Seq(simulated_sequence), id = "seq_file", description = "Simulated sequence from %s using order %d markov chain" % (seq_file, len(markov_list)-1))], output_file, 'fasta')
    open(output_file + ".out", "w").write(fa_out)



if __name__ == "__main__":
    args = parse_params(sys.argv[1:])
    create_chromosome_file(seq_file = args.seq_file, k = args.k, output_file = args.output, 
                           repeat_file = args.repeat_file, use_3prime = args.negative_strand, 
                           filter_file = args.family_file, mask = args.mask, seed = args.seed,
                           max_interval = args.max_interval, min_interval = args.min_interval, num_repeats = args.num_repeats,
                           max_length = args.max_length, toss_low = not args.low_complexity,
                           rep_base = args.rep_base, f = args.family_min, num_fams = args.num_family)
