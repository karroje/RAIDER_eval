#!/software/python/3.3.3/bin/python3.3
import sys
import subprocess
import os
import os.path
import shutil
import argparse
from redhawk import *
import tempfile
import re
import perform_stats

#################################################################
# The following global variables are related to debugging issues.
show_progress = False
job_index = {}
time_limit = "4:00:00"

#################################################################

#################################################################
# These global variables have to do with executable locations.
MacLocations = {'build_lmer_table':'/usr/local/RepeatScout/build_lmer_table',
                'RptScout':'/usr/local/RepeatScout/RepeatScout',
                'filter_stage-1':'/usr/local/RepeatScout/filter-stage-1.prl',
                'filter_stage-2':'/usr/local/RepeatScout/filter-stage-2.prl',
                'raider':'./raider',
                'bigfoot':'./bigfoot',
                'python':'python3.4',
                'araider':'./araider'}
RedhawkLocations = {'build_lmer_table':'./build_lmer_table',
                    'RptScout':'./RepeatScout',
                    'filter_stage-1':'./filter-stage-1.prl',
                    'filter_stage-2':'./filter-stage-2.prl',
                    'raider':'./raider',
                    'bigfoot':'./bigfoot',
                    'python':'python3.3',
                    'araider':'./araider'}
Locations = None;    # This will be set to one of the above two, and references to find exectuable locations.

#########
# Utility functions
def sum_resources(T1, T2):
    if T1[0] == -1 or T2[0] == -1:
        return [-1]*4
    return [T1[0] + T2[0], T1[1] + T2[1], max(T1[2], T2[2]), max(T1[3], T2[3])]

def file_base(file):
    return os.path.basename(file)

def get_job_index(s):
    global job_index
    if s not in job_index:
        job_index[s] = 0
    v = job_index[s]
    job_index[s] += 1
    return v

def file_dir(file):
    return file.rstrip(file_base(file)).rstrip("/")




def parse_params(args):
    """Parse command line arguments using the argparse library"""
    parser = argparse.ArgumentParser(description = "Evaluate RAIDER against RepeatScout")
    
    # GENERAL ARGUMENTS
    #parser2 = parser.add_mutually_exclusive_group()
    #parser2.add_argument('--organize', action = "store_true", help = "Create directory for all Raider Eval output", default = False)
    #parser2.add_argument('--no', '--named_organize', dest = "named_organize", help = "Organize under a named directory", default = None)


    # TOOL SELECTION
    parser_tools = parser.add_argument_group("tool selection (all on by default)")
    parser_tools.add_argument('-R', '--raider_on', dest = 'run_raider', action = 'store_true', help = 'Turn RAINDER on', default = False)
    parser_tools.add_argument('--AR', '--araider_on', dest = 'run_araider', action = 'store_true', help = 'Turn ARAIDER on', default = False)
    parser_tools.add_argument('--RS', '--repscout_on', dest = 'run_repscout', action = 'store_true', help = 'Turn RAINDER on', default = False)
    parser_tools.add_argument('-B', '--bigfoot_on', dest = 'run_bigfoot', action = 'store_true', help = 'Turn BIGFOOT on', default = False)
    parser_tools.add_argument('-P', '--piler_on', dest = 'run_piler', action = 'store_true', help = 'Turn PILER on', default = False)
    parser_tools.add_argument('-A', '--all_tools', dest = 'all_tools', action = 'store_true', help = 'Turn all tolls on (overide all other tool arguments)', default = False)
    parser_tools.add_argument('--tl', '--time_limit', dest = 'time_limit', help = 'Redhawk time limit (max: 400:00:00 default: 4:00:00)', default = "4:00:00")
    # Will later add: RepeatModeler, RECON, PILER (other?)


    # I/O ARGUMENTs
    parser_io = parser.add_argument_group("i/o arguments")
    parser_io.add_argument('-r', '--results_dir', dest = "results_dir", help = "Directory containing all results", default = "EVAL")
    parser_io.add_argument('--nuke', dest ='nuke', action = "store_true", help = "Nuke the results directory", default = False)
    parser_io.add_argument('--rd', '--raider_dir', dest = "raider_dir", help = "Subdirectory containing raider results", default = "RAIDER")
    parser_io.add_argument('--ard', '--araider_dir', dest = "araider_dir", help = "Subdirectory containing araider results", default = "ARAIDER")
    parser_io.add_argument('--rsd', '--rptscout_dir', dest = 'rptscout_dir', help = "Subdirectory containing rpt scout results", default = "RPT_SCT")
    parser_io.add_argument('--bfd', '--bigfoot_dir', dest = 'bigfoot_dir', help = "Subdirectory containing bigfoot results", default = "BIGFOOT")
    parser_io.add_argument('--pd', '--pilder_dir', dest = 'piler_dir', help = "Subdirectory containing piler results", default = "PILER")
    parser_io.add_argument('--dd', '--data_dir', dest = 'data_dir', help = "Directory containing the resulting simulated chromosome", default = "SOURCE_DATA")

    
    # RAIDER ARGUMENTS
    raider_argument = parser.add_argument_group("RAIDER parameters")
    raider_argument.add_argument('-f', type = int, help = "E.R. occurrence threshold", default = 2)
    raider_argument.add_argument('-d', '--output_dir', help = "Raider output directory", default = None)
    raider_argument.add_argument('-e', '--output_ext', help = "Output Extension", default = None)
    raider_argument.add_argument('-C', '--cleanup_off', dest = "cleanup", action = "store_false", help = "Turn off file cleanup", default = True)
    raider_argument.add_argument('--raider_min', '--raider_min', type = int, help = "Minimum repeat length. Defaults to pattern length.", default = None)
    seed_group = raider_argument.add_mutually_exclusive_group(required = False)     
    seed_group.add_argument('-s', '--seed', dest = "seed", help = "Spaced seed string", default = "111111111111111111111111111111")    
    seed_group.add_argument('--sf', '--seed_file', dest = 'seed_file', help = 'File containing raider seeds', default = None)
    
    
    # REPSCOUT ARGUMENTS
    repscout_argument = parser.add_argument_group("REPSCOUT parameters")
    raider_argument.add_argument('--repscout_min', type = int, help = "Minimum repeat length for repscout.", default = 10)
    raider_argument.add_argument('--uff', '--use_first_filter', dest = "use_first_filter", action = "store_true", help = "Use the first RepScout filter", default = False)
    raider_argument.add_argument('--usf', '--use_second_filter', dest = "use_second_filter", action = "store_true", help = "Use the second RepScout filter", default = False)

    # BIGFOOT ARGUMENTS
    bigfoot_arguments = parser.add_argument_group("BIGFOOT parameters")
    bigfoot_arguments.add_argument('-L', '--bigfoot_L', type = int, help = "Minimum repeat length. Defaults to 20.", default = 20)
    bigfoot_arguments.add_argument('-min', '--bigfoot_min', type = int, help = "E.R. occurrence threshold", default = 2)
    bigfoot_arguments.add_argument('-I', '--bigfoot_I', type = float, help = "Minimum percent frequency of more frequent Lmer a less frequent Lmer must have to be part of the same family", default = 0.75)
    bigfoot_arguments.add_argument('-T', '--bigfoot_T', type = float, help = "Minimum percent of time a base must occur after an Lmer to be considered significant", default = 0.75)

    # REPEAT MASKER ARGUMENTS
    repeatmasker_arguments = parser.add_argument_group("RepeatMasker parameters")
    repeatmasker_arguments.add_argument('--masker_dir', help = "Repeat masker output directory", default = None)
    repeatmasker_arguments.add_argument('-p', '--pa', type = int, help = "Number of processors will be using", default = 1)

    # STATISTICS ARGUMENT
    stats_group = parser.add_argument_group(title = "Statistics argument")
    stats_group.add_argument('--stats_file', dest = 'stats_file', help = "Statistics output directory", default = "stats.txt")
    #stats_group.add_argument('--print_reps', action = "store_true", help = "Print out repeats in statistics file", default = False)

    # DEBUGGING ARGUMENTS
    debug_group = parser.add_argument_group(title = "debugging")
    debug_group.add_argument('--sp', '--show_progress', dest = 'show_progress', action = 'store_true', help = "Print reports on program progress to stderr", default = False)
    debug_group.add_argument('--so', '--simulate_only', dest = 'simulate_only', action = 'store_true', help = "Quit after creating simulated file", default = False)

    subparsers = parser.add_subparsers(dest="subparser_name")
    
    # SEQUENCE FILE OPTION ARGUMENTS
    parser_seqs = subparsers.add_parser("seq_files")
    parser_seqs.add_argument('seq_files', nargs = '+', help = "Use files directly (no simulation)", default = None)
    
    # CHROMOSOME SIMULATION OPTION ARGUMENTS
    parser_chrom = subparsers.add_parser("chrom_sim")
    parser_chrom.add_argument('-k', type = int, help = "Order of markov chain", default = 5)  # KARRO: Added this
    parser_chrom.add_argument('--rng_seed', type = int, help = "RNG seed", default = None) 
    parser_chrom.add_argument('-n', '--negative_strand', action = "store_true", help = "Use repeats on negative string", default = False)
    parser_chrom.add_argument('--family_file', help = "List of repeat families to use", default = None)
    parser_chrom.add_argument('--mc', '--mc_file', dest = 'mc_file', help = "Markov Chain file", default = False)
    parser_chrom.add_argument('--mi', '--max_interval', dest = "max_interval", type = int, 
                              help = "Maximum allowed length of interval between repeats; -1 value (default) means no maximum", default = None) 
    parser_chrom.add_argument('--rn', '--retain_n', dest = "retain_n", action = 'store_true', 
                              help = "If used, will use the whole chromosome.  Otherwise, cuts of Ns at either end.", default = False)
    parser_chrom.add_argument('--nr', '--num_repeats', dest = 'num_repeats', type = int, 
                              help = "Specify the number of repeats.  Simulation will terminate either 1000 bases or max interval bases past the nth instance of a repeat (excluding any other repeats in that range).", default = None)
    parser_chrom.add_argument('-l', '--length', type = int, help = "Simulated sequence length", default = None)
    parser_chrom.add_argument('-o', '--output', help = "Output file (Default: replace chromosome file \".fa\" with \".sim.fa\")")
    parser_chrom.add_argument('-t', '--num_sims', type = int, dest = "num_sims", help ="Number of simulations", default = 1)
    parser_chrom.add_argument('--lc', '--low_complexity', dest = 'low_complexity', action = 'store_false', help = "Toss low complexity and simple repeats (tossed by default)", default = True)

    parser_chrom.add_argument('chromosome', help = "Template chromosome file")


    arg_return =  parser.parse_args(args)

    global time_limit
    time_limit = arg_return.time_limit

    global show_progress
    show_progress = arg_return.show_progress

    if arg_return.all_tools:
        arg_return.run_raider = True
        arg_return.run_araider = True
        arg_return.run_repscout = True
        arg_return.run_bigfoot = True
        arg_return.run_piler = True

    #### The following is to set the global debugging variables 
    if arg_return.simulate_only:    # Set to supress all tools
        arg_return.run_raider = False
        arg_return.run_araider = False
        arg_return.run_repscout = False
        arg_return.run_bigfoot = False
        arg_return.run_piler = False



    return arg_return


############################################################
# Main functions 
def simulate_chromosome(chromosome, repeat, rng_seed, length, neg_strand, fam_file, data_dir, output_file, file_index, k, mc_file, mi, retain_n, num_repeats, low_complexity):
    """Given chromosome file and repeat file and rng_seed, runs chromosome 
    simulator and then passes raider params (including path to new simulated chromosome 
    file) into run_raider"""

    # Output file is either specified or replace .fa with .sim.#.fa
    length_arg = "-l %d" % (length) if length else ""
    k_arg = "-k %d" % (k)
    seed_arg = "-s %d" % (rng_seed) if rng_seed else ""
    neg_arg = "-n" if neg_strand else ""
    fam_arg = "-f %s" % (fam_file) if fam_file else ""
    mi = ("--mi %d" % (mi)) if mi else ""
    retain_n = "--rn" if retain_n else ""
    num_repeats = ("--nr %d" % (num_repeats)) if num_repeats else ""
    low_complexity = "--lc" if low_complexity else ""
    seq_arg = chromosome
    repeat_arg = repeat

    output_file = (output_file if output_file else re.sub(".fa$", ".sim.%d.fa" % (file_index), file_base(chromosome)))
    output_path = "%s/%s" % (data_dir, output_file)

    mc = "--mc %s" % mc_file if mc_file else ""
    cmd = "{python} chromosome_simulator.py {mi} {length} {mc} {k} {seed} {neg} {fam} {seq} {retain_n} {num_repeats} {lc} {repeat} {output}".format(python = Locations['python'], mi=mi, mc=mc, length=length_arg, k=k_arg, seed=seed_arg, neg=neg_arg, fam=fam_arg, seq=seq_arg, retain_n=retain_n, num_repeats=num_repeats, lc=low_complexity, repeat=repeat_arg, output=output_path)

    if show_progress:
        sys.stderr.write("Creating simulation:\n%s\n" % (cmd))
        sys.stderr.flush()
    progress_fp.write("Creating simulation:\n%s\n" % (cmd))
    progress_fp.flush()

    batch_name = data_dir + "/" + output_file + ".sim.batch"
    job_name = "simulation.%d" % (get_job_index("simulation"))
    
    #print("Sim batch: %s\n" % (batch_name))
    p = pbsJobHandler(batch_file = batch_name, executable = cmd, job_name = job_name,
                      stdout_file = output_file + ".stdout", stderr_file = output_file + ".stderr", 
                      output_location = data_dir, walltime = time_limit)
    p.submit()

    p.output_file = output_file
    p.seq_file = file_base(output_file)
    p.sim_output = output_path
    p.index = file_index
    return p

# def run_raider_chrom(p, seed, f, m, output_dir):
#    """Given the pbs object used to start simulated chromosome job as well as 
#    raider parameters, wait until the job is done and then call run_raider 
#    on the newly generated chromosome output file"""
#    p.wait()
#    return run_raider(seed, f, m, p.chrom_output, output_dir, p.curr_dir)
def run_raider(seed, seed_num, f, m, input_file, raider_dir):
    """Given raider parameters and an input file, run RAIDER and put the output into
    the directory specified in output_dir (creating a random name is none is
    specified."""

    input_base = file_base(input_file).rstrip(".fa")
    output_dir = raider_dir + "/" + input_base.upper() + ".s" + str(seed_num)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    min_arg = "-m %d" % (m) if m else ""
    cmd1 = "{raider} -q -c {f} {min_arg} {seed} {input_file} {output_dir}".format(raider = Locations['raider'], f = f, min_arg = min_arg, seed = seed, input_file = input_file, output_dir = output_dir)

    out_file = raider_dir + "/" + input_base + ".s" + str(seed_num) + ".raider_consensus.txt"
    lib_file = raider_dir + "/" + input_base + ".s" + str(seed_num) + ".raider_consensus.fa"

    cmd2 = "{python} consensus_seq.py -s {seq_file} -e {elements_dir}/elements {output_file} {fa_file}".format(python = Locations['python'], seq_file = input_file, elements_dir = output_dir, output_file = out_file, fa_file = lib_file)

    if show_progress:
        sys.stderr.write("\nLaunching raider:\n%s\n%s\n" % (cmd1, cmd2))
        sys.stderr.flush()
    progress_fp.write("\nLaunching raider:\n%s\n%s\n" % (cmd1, cmd2))
    progress_fp.flush()

    batch_name = raider_dir + "/" + input_base + ".raider.batch"
    job_name = "raider.%d" % get_job_index("raider")
    #progress_fp.write("Sim batch: " + batch_name + "\n")
    p = pbsJobHandler(batch_file = batch_name, executable = cmd1 + "; " + cmd2, job_name = job_name,
                      stdout_file = input_base + ".raider.stdout", stderr_file = input_base + ".raider.stderr",
                      output_location = output_dir, walltime = time_limit)

    p.submit()
    p.tool_resources = [0]*4

    p.description = "raider"
    p.tools_resources = [0]*4
    p.seed = seed
    p.seed_num = seed_num
    p.seq_file = input_file
    p.lib_file = lib_file

    return p

def run_composites_finder(elements_file, seq_file, compositesFinderDir):
    input_base = file_base(elements_file)
    output_dir = compositesFinderDir + "/" + input_base.upper()
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    compositesDiscover = compositesFinderDir + "/" + "CompositesDiscover"
    slimComFinder = compositesFinderDir + "/" + "SlimComFinder.py"
    cmd1 = "{compositesFinder} {input_file}".format(compositesFinder = compositesDiscover, input_file = elements_file)
    cmd2 = "{python} {slim_composites_finder} {elements} {sequence_file} {output_file}".format(python = "python", slim_composites_finder = slimComFinder,
                        elements = elements_file, sequence_file = seq_file, output_file = output_dir + "/" + "ConsensusSequences")

    if show_progress:
        sys.stderr.write("\nLaunching composites finder:\n%s\n%s\n" % (cmd1, cmd2))
        sys.stderr.flush()
    progress_fp.write("\nLaunching composites finder:\n%s\n%s\n" % (cmd1, cmd2))
    progress_fp.flush()

    batch_name =  compositesFinderDir + "/" + input_base + ".composites finder.batch"
    job_name = "composites finder.%d" % get_job_index("composites finder")
    p = pbsJobHandler(batch_file = batch_name, executable = cmd1 + "; " + cmd2, job_name = job_name,
                      stdout_file = input_base + ".comFinder.stdout", stderr_file = input_base + ".comFinder.stderr",
                      output_location = output_dir, walltime = time_limit)

    p.submit()
    p.description = "composites.finder"
    p.elementsFile = elements_file
    p.seqFile = seq_file

    return p

def run_araider(seed, seed_num, f, m, input_file, araider_dir):
    """Given raider parameters and an input file, run RAIDER and put the output into
    the directory specified in output_dir (creating a random name is none is
    specified."""

    input_base = file_base(input_file).rstrip(".fa")
    output_dir = araider_dir + "/" + input_base.upper() + ".s" + str(seed_num)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    min_arg = "-m %d" % (m) if m else ""
    cmd1 = "{araider} -q -c {f} {min_arg} {seed} {input_file} {output_dir}".format(araider = Locations['araider'], f = f, min_arg = min_arg, seed = seed, input_file = input_file, output_dir = output_dir)

    out_file = araider_dir + "/" + input_base + ".s" + str(seed_num) + ".araider_consensus.txt"
    lib_file = araider_dir + "/" + input_base + ".s" + str(seed_num) + ".araider_consensus.fa"

    cmd2 = "{python} consensus_seq.py -s {seq_file} -e {elements_dir}/elements {output_file} {fa_file}".format(python = Locations['python'], seq_file = input_file, elements_dir = output_dir, output_file = out_file, fa_file = lib_file)

    if show_progress:
        sys.stderr.write("\nLaunching araider:\n%s\n%s\n" % (cmd1, cmd2))
        sys.stderr.flush()
    progress_fp.write("\nLaunching araider:\n%s\n%s\n" % (cmd1, cmd2))
    progress_fp.flush()

    batch_name = araider_dir + "/" + input_base + ".araider.batch"
    job_name = "araider.%d" % get_job_index("araider")
    #progress_fp.write("Sim batch: " + batch_name + "\n")
    p = pbsJobHandler(batch_file = batch_name, executable = cmd1 + "; " + cmd2, job_name = job_name,
                      stdout_file = input_base + ".araider.stdout", stderr_file = input_base + ".araider.stderr",
                      output_location = output_dir)

    p.submit()
    p.tool_resources = [0]*4

    p.description = "araider"
    p.tools_resources = [0]*4
    p.seed = seed
    p.seed_num = seed_num
    p.seq_file = input_file
    p.lib_file = lib_file

    return p

def run_bigfoot(input_file, bigfoot_dir, L, C, I, T):
    """Runs BIGFOOT and returns a submitted pbs object with specific attributes used to run RepeatMasker.
    * input_file: The name of the .fa sequence file being searched.
    * bigfoot_dir: The name of the directory that will contain all files from this run.
    """
    input_base = file_base(input_file).rstrip(".fa")     # The name of the inputfile -- which I've been using as a basis for all file names
    output_dir = bigfoot_dir + "/" + input_base.upper() # If bigfoot creates its own directory for information, use this as the name of that directory.
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    cmd1 = "{bigfoot} -l {L} -c {C} --I {I} --T {T} {input_file} {output_dir}".format(bigfoot = Locations['bigfoot'], L = L, C = C, I = I, T = T, output_dir = output_dir, input_file = input_file)   # Put the command-line executable for for bigfoot here.  Use input_file for the input file name, and put any output into bigfoot_dir
    cmd2 = "cp {output_dir}/seeds {bigfoot_dir}/{input_base}.seeds".format(output_dir=output_dir, bigfoot_dir=bigfoot_dir, input_base=input_base)
    cmd = cmd1 + "; " + cmd2;

    if show_progress:
        sys.stderr.write("\nLaunching bigfoot:\n%s\n" % (cmd))
        sys.stderr.flush()
    progress_fp.write("\nLaunching bigfoot:\n%s\n" % (cmd))
    progress_fp.flush()

    
    lib_file = bigfoot_dir + "/" + input_base + ".seeds"
    batch_name = bigfoot_dir + "/" + input_base + ".bigfoot.batch"    # This is the batch fils for the qsub command.
    job_name = "bigfoot.%d" % get_job_index("bigfoot")                # This is the redhawk jobname.  get_job_index just assigned the next unused number (for then running multiple jobs)
    stdout_file = input_base + ".bigfoot.stdout"                      # Anything bigfoot prints to stdout will be redirected here
    stderr_file = input_base + ".bigfoot.stderr"                      # Anything bigfoot prints to stderr will be redirected here
    p = pbsJobHandler(batch_file = batch_name, executable = cmd, job_name = job_name,
                      stdout_file = stdout_file, stderr_file = stderr_file,
                      output_location = output_dir, walltime = time_limit)    


    p.submit()
    p.description = "bigfoot"
    p.tool_resources = [0]*4
    p.seq_file = input_file     # Required by run_repeat_masker -- uses this as the source sequence.
    p.lib_file = lib_file     # This should be set to the file name that will be the library for the repeatmasker run

    return p 

def run_piler(input_file, piler_dir):
    """Runs Piler and returns a submitted pbs object with specific attributes used to run RepeatMasker."""
    input_base = file_base(input_file).rstrip(".fa")
    lib_file = input_base + ".lib"
    if not os.path.exists(piler_dir):
        os.makedirs(piler_dir)


    cmd = "{python} run_piler.py {input_file} {piler_dir} {output_file}".format(python = Locations['python'], input_file = input_file, piler_dir = piler_dir, output_file = lib_file);
    
    if show_progress:
        sys.stderr.write("\nLaunching Piler:\n%s\n" % (cmd))
        sys.stderr.flush()
    progress_fp.write("\nLaunching Piler:\n%s\n" % (cmd))
    progress_fp.flush()

    batch_name = piler_dir + "/" + input_base + ".piler.batch";
    job_name = "piler%d" % get_job_index("piler")
    stdout_file = input_base + ".piler.stdout";
    stderr_file = input_base + ".piler.stderr";
    p = pbsJobHandler(batch_file = batch_name, executable = cmd, job_name = job_name,
                      stdout_file = stdout_file, stderr_file = stderr_file,
                      output_location = piler_dir, walltime = time_limit)

    p.submit()
    p.description = "piler"
    p.tool_resources = [0]*4
    p.seq_file = input_file
    p.lib_file = piler_dir + "/" + lib_file

    return p

# def create_raider_consensus(p, output):
#     """Given the pbs object (from redhawk.py) used to start a RAIDER job, this
#     waits until the job is done, then invokes consensus_seq.py on the ouput and
#     writes the results to the output directory."""
#     p.wait();
#     cmd = "python3.3 consensus_seq.py -s %s -e %s/elements %s" % (p.seq_file, p.raider_output, output)
#     #cmd = "python3.3 consensus_seq.py -s %s -e %s/elements %s" % (p.file, p.raider_output, output)
#     if show_progress:
#     p2.submit()
#     p2.seq_file = p.seq_file
#     p2.output = output
#     return p2


def run_scout(input_file, output_dir, min_freq, length, use_first_filter):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    input_name= file_base(input_file)
    
    # First: run build_lmer_table
    lmer_output = output_dir + "/" + input_name.rstrip(".fa") + ".freq.fa"
    cmd1 = "{build_lmer_table_exe} -min {min} -sequence {sequence} -freq {freq}".format(build_lmer_table_exe=Locations['build_lmer_table'], min=min_freq,  
                                                                                               sequence = input_file, freq = lmer_output)

    # Next: Run RepeatScout
    rptscout_output = output_dir + "/" + input_name.rstrip(".fa") + ".repscout.fa"
    cmd2 = "{RptScout_exe} -sequence {sequence} -freq {freq} -output {output}".format(RptScout_exe = Locations['RptScout'], sequence = input_file, freq = lmer_output, output = rptscout_output)

    # Next: Run filter-stage-1
    if use_first_filter:
        filter_stage_output = output_dir + "/" + input_name.rstrip(".fa") + ".repscout.filtered.fa"
        cmd3 = "cat {input} | perl {filter} > {filter_output}".format(input=rptscout_output, filter = Locations['filter_stage-1'], filter_output = filter_stage_output)
    else:
        cmd3 = ""

    if show_progress:
        sys.stderr.write("\nRepeatScout:\n%s\n%s\n%s\n" % (cmd1, cmd2, cmd3))
        sys.stderr.flush()
    progress_fp.write("\nRepeatScout:\n%s\n%s\n%s\n" % (cmd1, cmd2, cmd3))
    progress_fp.flush()
        
    batch_name = output_dir + "/" + file_base(input_file) + ".repscout1.batch"
    job_name = "rptscout.%d" % (get_job_index("repscout"))
    #progress_fp.write("Sim batch: " + batch_name + "\n")
    p = pbsJobHandler(batch_file = batch_name, executable = cmd1 + "; " + cmd2 + "; " + cmd3, job_name = job_name,
                      stdout_file = file_base(rptscout_output) + ".stdout", stderr_file = file_base(rptscout_output) + ".stderr",
                      output_location = output_dir, walltime = time_limit)

    p.submit()
    p.description = "rep_scout"
    p.tool_resources = [0,0,0,0]

    p.seq_file = input_file
    p.lib_file = filter_stage_output if use_first_filter else rptscout_output

    return p

    # output_file = re.sub(".fa$", ".scout.fa",os.path.basename(input_file))
    
    # output_path = "%s/%s" %(curr_dir, output_file)
    # cmd = "./RunRepeatScout.sh {sequence} {length} {output} {min}".format(sequence=input_file, length=m, output=output_path, min=f)
    # print(cmd)
    # p = pbsJobHandler(batch_file = "%s.batch" % output_file, executable = cmd)
    # p.submit()
    # p.seq_file = input_file
    # p.scout_output = output_path
    # p.curr_dir = curr_dir
    # return p

def scout_second_filter(p, min_freq):
    """NOT CURRENTLY WORKING!!! Does not run correctly, and does not properly adjust time"""
    p.wait()

    filter2_stage_output = p.seq_file.rstrip(".fa") + ".repscout.filtered2.fa"
    cmd = "cat {output} | perl {filter} --cat={cat} --thresh={thresh} > {final}".format(output = p.lib_file, filter = Locations['filter_stage-2'], cat = p.rm_output, thresh = min_freq, final = filter2_stage_output)
    
    if show_progress:
        sys.stderr.write("\nRepeatScout Filter2:\n%s\n" % cmd)
        sys.stderr.flush()
    progress_fp.write("\nRepeatScout Filter2:\n%s\n" % cmd)
    progress_fp.flush()

    batch_name = file_dir(p.rm_output) + "/" + file_base(p.seq_file).rstrip(".fa") + ".repscout2.fa"
    job_name = "filter2%d" % get_job_index("filter2")

    #print("Sim batch: " + batch_name + "\n")
    p2 = pbsJobHandler(batch_file = batch_name, executable = cmd, job_name = job_name,
                       stdout_file = file_base(p.seq_file) + ".repscout2.stdout", stderr_file = file_base(p.seq_file) + ".repscout2.stderr",
                       output_location = file_dir(p.seq_file), walltime = time_limit)

    p2.submit()
    p2.description = "rep_scout"
    p2.time_resources = p.time_resources + p.getResources()
    p2.lib_file = p.lib_file
    p2.seq_file = p.seq_file
    p2.lib_file = filter2_stage_output

    return p2



def run_repeat_masker(p, num_processors):
    """Given the pbs object used to start a consensus sequence job as well as
    repeatmasker arguments, wait until the job is done and then call repeatmasker 
    on the output and put results in masker_dir (current dir if unspecified)"""
    p.wait()
    p.loadResources()

    input_base = file_base(p.seq_file)  # Base name of the file used for input

    output_dir = file_dir(p.lib_file)  # Output will fo in the same directory as the lib file
    cmd = "RepeatMasker -nolow -lib {library} -pa {pa} -dir {dir} {seq_file}".format(library = p.lib_file, pa = num_processors, dir = output_dir, seq_file = p.seq_file)

    if show_progress:
        sys.stderr.write("\nLaunch repeatmasker (%s):\n%s\n" % (p.description, cmd))
        sys.stderr.flush()
    progress_fp.write("\nLaunch repeatmasker (%s):\n%s\n" % (p.description, cmd))
    progress_fp.flush()

    batch_name = p.lib_file.rstrip(".fa") + ".rm.batch"
    job_name = "repmask.%d" % get_job_index("repmask")
    #print("Sim batch: " + batch_name + "\n")
    p2 = pbsJobHandler(batch_file = batch_name, executable = cmd, ppn = num_processors, RHmodules = ["RepeatMasker", "python-3.3.3"],
                       job_name = job_name, stdout_file = input_base + ".repmask.stdout", stderr_file = input_base + ".repmask.stderr",
                       output_location = output_dir, walltime = time_limit);
    p2.submit(preserve=True)

    p2.description = "RptMasker"
    p2.seed = p.seed if hasattr(p, "seed") else "NA"
    p2.seed_num = p.seed_num if hasattr(p, "seed_num") else "NA"
    p2.dir = output_dir
    p2.lib_file = p.lib_file
    p2.seq_file = p.seq_file
    p2.rm_output = output_dir + "/" + file_base(p.seq_file) + ".out"
    p2.tool_resources = p.resources
    p2.tool_description = p.description
    return p2





def performance_sum(stats_jobs, stats_dir, curr_dir, test):
    """Given a list of all of the statistics jobs, uses the statistics output files to
    generate a summary file indicative of overall performance. Put results in stats_dir
    (Current dir if unspecified)"""
    tps = 0
    tns = 0
    fps = 0
    fns = 0
    for p in stats_jobs:
        sf = open(p.stats_output, "r")
        tps += int(re.split("\s+", sf.readline())[1])
        fps += int(re.split("\s+", sf.readline())[1])
        tns += int(re.split("\s+", sf.readline())[1])
        fns += int(re.split("\s+", sf.readline())[1])
        sf.close()
    stats_file = "summary.%s.stats" % (test)
    smry_path = "%s/%s" % (stats_dir, stats_file) if stats_dir else "%s/%s" %(curr_dir, stats_file) if curr_dir else stats_file
    
    smry = open(smry_path, 'w')
    smry.write("Evaluation completed.\n")
    smry.write("True Positives (TP): \t %d \n" % (tps))
    smry.write("False Positives (FP): \t %d \n" % (fps))
    smry.write("True Negatives (TN): \t %d \n" % (tns))
    smry.write("False Negatives (FN): \t %d \n" % (fns))
    smry.write("\nPerformance of Repeat Classification Tool\n")
    smry.write("Sensitivity (TPR): \t\t %f %%\n" % (tps/(tps + fns)))
    smry.write("Specificity (TNR): \t\t %f %%\n" % (tns/(fps + tns)))
    smry.write("Precision (PPV): \t\t %f %%\n" % (tps/(tps + fps)))
    smry.write("Neg. Pred. Val. (NPV): \t %f %%\n" % (tns/(tns + fns)))
    smry.write("Fall-Out (FPR): \t\t %f %%\n" % (fps/(fps + tns)))
    smry.write("False Disc. Rate (FDR): \t  %f %%\n" % (fps/(tps + fps)))
    smry.close()



test_tools = ["raider", "bigfoot", "piler", "rep_scout", "araider"]  # List of implemented tools 
############################################################################################
if __name__ == "__main__":
    args = parse_params(sys.argv[1:])

    ####
    # Currently: We check for the RepeatScout executable at the location on my Mac; if found, we
    # assume we are running on the Mac.  If not, we check for in the Redhawk location, and if found
    # assume we are running on redhawk.  Otherwise we print and error and quit.
    if os.path.exists(MacLocations["RptScout"]):
        Locations = MacLocations
    elif os.path.exists(RedhawkLocations['RptScout']):
        Locations = RedhawkLocations
    else:
        assert False, "Could not determine host."
    ###

    if args.nuke and os.path.exists(args.results_dir):
        subprocess.call("rm -r %s" % args.results_dir, shell = True)
    if not os.path.exists(args.results_dir):
        os.makedirs(args.results_dir)


    ### Generate simulated file(s) and run to completion

    ### Set up the debugging log file (if needed)
    progress_fp = open(args.results_dir + "/debug.txt", "w")

 
    data_dir = args.results_dir + "/" + args.data_dir
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)



    # First: we put the chromosomes (simulated or real) into data_dir
    if args.subparser_name == "chrom_sim":
        # Launch the jobs

        f = lambda i: simulate_chromosome(chromosome = args.chromosome, repeat = args.chromosome + ".out",
                                          rng_seed = args.rng_seed, length = args.length, 
                                          neg_strand = args.negative_strand, fam_file = args.family_file, 
                                          data_dir = args.results_dir + "/" + args.data_dir, output_file = args.output, file_index = i, 
                                          k = args.k, mc_file = args.mc_file, mi = args.max_interval,
                                          retain_n = args.retain_n, num_repeats = args.num_repeats, low_complexity = args.low_complexity)
        J = [f(i) for i in range(args.num_sims)]

        # Run jobs to completion
        [j.wait() for j in J]    # Let all the simulations finish

        # Get the list of simulated file names
        file_list = [j.sim_output for j in J]

    else:
        file_list = []
        for file in args.seq_files:
            file_list.append(data_dir + "/" + file_base(file))
            shutil.copy(file, file_list[-1])
            shutil.copy(file + ".out", file_list[-1] + ".out")

    ### Start running each tool.  Each tool should run, creating the repeat masker library (putting the file name
    ### in the pbs lib_file attribute), then run repeat masker (putting the output file name in the pbs
    ### rm_output job.

    ############## First: Launch tools
    jobs = []
    if args.run_raider:
        seed_list = [seed for line in open(args.seed_file) for seed in re.split("\s+", line.rstrip()) if seed] if args.seed_file else [args.seed]
        jobs += [run_raider(seed = seed, seed_num = i, f = args.f, m = args.raider_min, input_file = file, 
                            raider_dir = args.results_dir + "/" + args.raider_dir) for i,seed in enumerate(seed_list)
                 for file in file_list]

    if args.run_araider:
        seed_list = [seed for line in open(args.seed_file) for seed in re.split("\s+", line.rstrip()) if seed] if args.seed_file else [args.seed]
        jobs += [run_araider(seed = seed, seed_num = i, f = args.f, m = args.raider_min, input_file = file, 
                            araider_dir = args.results_dir + "/" + args.araider_dir) for i,seed in enumerate(seed_list)
                 for file in file_list]
    if args.run_repscout:
        jobs += [run_scout(input_file = file, output_dir = args.results_dir + '/' + args.rptscout_dir, min_freq = args.f, length = args.repscout_min, use_first_filter = args.use_first_filter) for file in file_list]

    if args.run_bigfoot:
        bigfoot_dir = args.results_dir + "/" + args.bigfoot_dir    # Name of the directory all bigfoot files will go into
        if not os.path.exists(bigfoot_dir):
           os.makedirs(bigfoot_dir)
        jobs += [run_bigfoot(input_file = file, bigfoot_dir = bigfoot_dir, L = args.bigfoot_L, C = args.bigfoot_min, I = args.bigfoot_I, T = args.bigfoot_T) for file in file_list]

    if args.run_piler:
        piler_dir = args.results_dir + "/" + args.piler_dir    # Name of the directory all piler files will go into
        if not os.path.exists(piler_dir):
           os.makedirs(piler_dir)
        jobs +=[run_piler(input_file = file, piler_dir = piler_dir) for file in file_list]


    ############## Second: Launch repeatmasker jobs
    job_set = {j for j in jobs}
    RM_jobs = set()
    while job_set:
        finished_jobs = set()
        for j in job_set:
            if not j.isJobRunning():
                new_job = run_repeat_masker(j, args.pa)
                RM_jobs.add(new_job)
                finished_jobs.add(j)
        job_set = job_set - finished_jobs

    # Let them run to completion
    job_dic = {tool:[] for tool in test_tools}
               
    for j in RM_jobs:
        j.wait();
        job_dic[j.tool_description].append(j)
    job_dic['raider'].sort(key = lambda x: x.seed_num)
    job_dic['araider'].sort(key = lambda x: x.seed_num)
        

    # Print output files log
    with open(args.results_dir + "/file_log.txt", "w") as fp:
        for i in range(len(file_list)):
            fp.write("%d simulation_file %s\n" % (i, file_list[i]))
        for k in test_tools:
            fp.write(k + "\n")
            for j in job_dic[k]:
                fp.write(j.rm_output + "\n")

    ######
    # Create copy of seed file (if RAIDER is being used)
    if job_dic['raider']:
        with open(args.results_dir + "/seed_file.txt", "w") as fp:
            fp.write("\n".join(["{index:<5}{seed}".format(index=i,seed=s) for i,s in enumerate(seed_list)]) + "\n")
            
    if job_dic['araider']:
        with open(args.results_dir + "/seed_file.txt", "w") as fp:
            fp.write("\n".join(["{index:<5}{seed}".format(index=i,seed=s) for i,s in enumerate(seed_list)]) + "\n")

    ######
    # Calculate statistics (not bothering with parallelization yet)
    print_str = "{:<12}" + "{:<5}" + "".join("{:<14}"*4) + "".join("{:<14}"*6) + "".join("{:<14}"*8) + "\n"
    with open(args.results_dir + "/" + args.stats_file, "w") as fp:
        fp.write(print_str.format("#tool", "seed", "tp", "fp", "fn", "tn", "tpr", "tnr", "ppv", "npv", "fpr", "fdr","ToolCpuTime", "ToolWallTime", "ToolMem", "ToolVMem", "RMCpuTime", "RMWallTime", "RMMem", "RMVMem"))

        for key in test_tools:
            for p in job_dic[key]:
                Counts, Stats, Sets = perform_stats.perform_stats(p.seq_file + ".out", p.rm_output, None)
                Stats = [round(x,5) for x in Stats]
                fp.write(print_str.format(*([key, p.seed_num] + list(Counts) + list(Stats) + list(p.tool_resources) + list(p.getResources()))))

        # for i in range(len(J)):
        #     if RAIDER_JOBS:
        #         for j in range(len(RAIDER_JOBS[i])):
        #             p = RAIDER_JOBS[i][j]
        #             Counts, Stats, Sets = perform_stats.perform_stats(J[i].sim_output + ".out", p.rm_output, None)
        #             Stats = [round(x,5) for x in Stats]
        #             fp.write(print_str.format(*(["raider", p.seed_num] + list(Counts) + list(Stats) + list(p.tool_resources) + list(p.getResources()))))
        #     if SCOUT_JOBS:
        #         p = SCOUT_JOBS[i]
        #         CountSJ, StatsSJ, SetsSJ = perform_stats.perform_stats(J[i].sim_output + ".out", p.rm_output, None)
        #         StatsSJ = [round(x,5) for x in StatsSJ]
        #         fp.write(print_str.format(*(["repscout", "NA"] + list(CountSJ) + list(StatsSJ) + list(p.tool_resources) + list(p.getResources()))))
        #     if BIGFOOT_JOBS:
        #         p = BIGFOOT_JOBS[i]
        #         CountBF, StatsBF, SetsBF = perform_stats.perform_stats(J[i].sim_output + ".out", p.rm_output, None)
        #         StatsBF = [round(x,5) for x in StatsBF]
        #         fp.write(print_str.format(*(["bigfoot", "NA"] + list(CountBF) + list(StatsBF) + list(p.tool_resources) + list(p.getResources()))))
        #     if PILER_JOBS:
        #         p = PILER_JOBS[i]
        #         CountP, StatsP, SetsP = perform_stats.perform_stats(J[i].sim_output + ".out", p.rm_output, None)
        #         StatsP = [round(x,5) for x in StatsP]
        #         fp.write(print_str.format(*(["piler", "NA"] + list(CountP) + list(StatsP) + list(p.tool_resources) + list(p.getResources()))))
                            

            
