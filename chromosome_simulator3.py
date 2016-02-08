"""Generate a chromosome using the shuffle method instead of a Markov chain.
Written only for debugging purposes -- not well constructed."""
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
    # parser.add_argumemnt('-c', '--cutoff', type = int, help = "Limit model to first c non-N bases")
    parser.add_argument('-k', type = int, help = "Order of Markov chain", default = 5)
    parser.add_argument('-s', '--seed', '-rng_seed', dest = 'seed', type = int, help = "RNG seed", default = None)
    parser.add_argument('-n', '--negative_strand', action = "store_true", help = "Use repeats on negative string", default = False)
    parser.add_argument('-f', '--family_file', help = "List of repeat families to use", default = None)
    parser.add_argument('-m', '--mask', action = "store_true", help = "Turn masking on (all repeats printed as lower case).", default = False)
    parser.add_argument('--mc', '--mc_file', dest = 'mc_file', help = "MC File (by default -- look in local directory; generates if not found).", default = None)
    parser.add_argument('-S', '--suppress_pmck', action = "store_true", help = "Suppress the generation of a .pmc<k> file to store the markov chain for re-use")
    parser.add_argument('--mi', '--max_interval', dest = "max_interval", type = int, help = "Maximum allowed length of interval between repeats; -1 value (default) means no maximum", default = -1)
    parser.add_argument('--rn', '--retain_n', dest = "retain_n", action = 'store_true', help = "If used, will use the whole chromosome.  Otherwise, cuts of Ns at either end.", default = False)
    parser.add_argument('--nr', '--num_repeats', dest = 'num_repeats', type = int, help = "Specify the number of repeats.  Simulation will terminate either 1000 bases or max interval bases past the nth instance of a repeat (excluding any other repeats in that range).", default = None)
    parser.add_argument('-l', '--max_length', dest = 'max_length', type = int, help = "Maximum allowed length of simulated sequence.", default = None)
    parser.add_argument('--lc', '--low_complexity', dest = 'low_complexity', action = 'store_false', help = "Toss low complexity and simple repeats (tossed by default)", default = True)
    #parser.add_argument('-o', '--output', help = "Output file (Default: replace chomosome file \".fa\" with \".sim.fa\")")
    parser.add_argument("seq_file", help = "Sequence file (must be .fa)")
    parser.add_argument("repeat_file", help = "RepeatMasker file (.fa.out)")
    parser.add_argument("output", help = "Output file")
    return parser.parse_args(args)


def nextRepeat(rpt_file, use_negative = True, S = {}, E = {}):
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
            if (strand == '+' or use_negative) and (family in S or not S) and not (rpt_class in E):
                yield chr, start, finish, strand, family, rpt_class, int(rpt_id)

# fa_out_header: The fixed header lines for the .fa.out file
fa_out_header = "\tSW\tperc\tperc\tperc\tquery\tposition in query\tmatching\trepeat\tposition in  repeat\n\tscore\tdiv.\tdel.\tins.\tsequence\tbegin\tend\t(left)\trepeat\tclass/family\tbegin\tend (left)\tID\n"
# fa_out_template: A template for creating lines for the .fa.out file.
fa_out_template = "\t0\t0\t0\t0\t{chr}\t{start}\t{finish}\t({left})\t{strand}\t{family}\t{rpt_class}\t0\t0\t(0)\t{rpt_id}\n"
def generate_chromosome(seq, markov_list, coord_adjust, rpt_gen, mask = False, max_interval = None, num_repeats = None, max_length = None, limiting_chr = None):
    """
    Generate a syntehtic sequence with real repeats:
    * seq: A sequence (as a string).
    * markov_list: List of the k+1 i-th order markov chains (from the markov_gen module).
    * coord_adjust: Size of the prefix that has been cut off the template sequence (requiring that
    *               that .fa.out coordinate be adjusted).
    * rpt_gen: A generating function returning the repeat information (created by nextRepeat)
    * mask: If true, all repeats will be lower-case.  Otherwise, upper case.)
    * max_interval: Maximum inter-repeat length.
    """
    current_coord = coord_adjust
    if max_interval == -1:
        max_interval = len(seq)

    s = []                # Hold the sequence (in chunks)
    fa_out = []           # Hold the new .fa.out file contents (by line)

    rpt_count = 0
    length = min(len(seq), max_length) if max_length else len(seq)
    debug_sim_len = 0
    
    for chr, start, finish, strand, family, rpt_class, rpt_id in rpt_gen:
        if limiting_chr and chr not in limiting_chr:
            continue

        if start >= current_coord:
            
            rpt_count += 1
            inter_seq_len = min(start-current_coord, max_interval)
            #inter_seq = markov_gen.generate_sequence(markov_list, inter_seq_len)
            inter_seq = "".join([random.choice("ACGT") for i in range(inter_seq_len)])
            assert len(inter_seq) == inter_seq_len
            s.append(inter_seq)
            debug_sim_len += len(inter_seq)
            coord_adjust += max(0, start-current_coord-max_interval)

            rpt_seq = seq[start:finish]
            s.append(rpt_seq.lower() if mask else rpt_seq.upper())
            debug_sim_len += len(rpt_seq)
            
            fa_out.append([chr, start+1-coord_adjust, finish-coord_adjust, strand, family, rpt_class, rpt_id])
            #fa_out.append(fa_out_template.format(chr=chr, start=start+1-coord_adjust, finish=finish-coord_adjust, strand=strand, family=family, rpt_class=rpt_class, rpt_id=rpt_id))
            
            if num_repeats and rpt_count == num_repeats:
                break

            current_coord = finish
    
    if num_repeats:
        max_interval = min(1000, max_interval)

    tail_length = min(max_interval, length-current_coord)
    if tail_length > 0:
        s.append("".join([random.choice("ACGT") for i in range(tail_length)]))
    

    sim_seq = "".join(s)
    sim_seq_len = len(sim_seq)
    fa_out_str = fa_out_header
    for chr, start, finish, strand, family, rpt_class, rpt_id in fa_out:
        fa_out_str += fa_out_template.format(chr=chr, start=start, finish=finish, left = sim_seq_len - finish, strand=strand, family=family, rpt_class=rpt_class, rpt_id=rpt_id)


    return sim_seq, fa_out_str

bases = set("ACGTacgt")
def loadSeqAndChain(seq_file, k, suppress_save = False, mc_file = None, retain_n = False):
    """Load the sequence and the Markov Chain List.
    Load the MC list from a file if it exists.  If not, create the chain
    and save it to the file for the next use (skip the save if suppressed)."""
    template_seq = str(SeqIO.read(seq_file, 'fasta').seq)

    # Cut out all the maximul prefix and suffix of ambiguity codes -- which will have no effect on the Markov chain construction.
    if not retain_n:
        start = 0
        while template_seq[start] not in bases: start += 1
        finish = len(template_seq)
        while template_seq[finish-1] not in bases: finish -= 1
        coord_adjust = start
        template_seq = template_seq[start:finish]
    else:
        coord_adjust = 0

    
    #mc_file = re.sub("\.(fa|fasta)$", ".pmc%d" % (k), seq_file) if mc_file is None else mc_file
    #if os.path.exists(mc_file):
    #    markov_list = markov_gen.read_pmck(mc_file)
    #else:
    #    markov_list = markov_gen.MarkovArray(k, template_seq)
    #    if not suppress_save:
    #        markov_gen.pickle_markov_list(markov_list, mc_file)

    markov_list = None
    return template_seq, markov_list, coord_adjust
                        

    

def create_chromosome_file(seq_file, repeat_file, output_file, k = 5, use_3prime = True, filter_file = "rpt_list.txt", mask = False, seed = None, suppress = False, max_interval = -1, retain_n = False, num_repeats = None, max_length = None, toss_low = False):
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
    random.seed(args.seed)
    template_seq, markov_list, coord_adjust  = loadSeqAndChain(args.seq_file, args.k, suppress, args.mc_file, args.retain_n)

    filter_set = {y.strip() for line in open(filter_file) for y in re.split("\s+", line.rstrip())} if filter_file else {}
    rpt_gen = nextRepeat(repeat_file, use_3prime, filter_set, E = {'Low_complexity', 'Simple_repeat'} if toss_low else {})

    simulated_sequence, fa_out = generate_chromosome(seq = template_seq, markov_list = markov_list, coord_adjust = coord_adjust, rpt_gen = rpt_gen, mask = mask, max_interval = max_interval, num_repeats = num_repeats, max_length = max_length)
            
    SeqIO.write([SeqRecord(seq = Seq(simulated_sequence), id = "seq_file", description = "Simulated sequence from %s using uniform random base distribution" % (seq_file))], output_file, 'fasta')
    open(output_file + ".out", "w").write(fa_out)

if __name__ == "__main__":
    args = parse_params(sys.argv[1:])
    create_chromosome_file(seq_file = args.seq_file, k = args.k, output_file = args.output, 
                           repeat_file = args.repeat_file, use_3prime = args.negative_strand, 
                           filter_file = args.family_file, mask = args.mask, seed = args.seed,
                           max_interval = args.max_interval, num_repeats = args.num_repeats,
                           max_length = args.max_length, toss_low = not args.low_complexity)
