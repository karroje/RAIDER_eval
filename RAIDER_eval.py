import sys
import subprocess
import os
import os.path
import shutil
import argparse
import tempfile
import re
import perform_stats
import time
import pickle
from parse_pra_output import parse_pra_output
from checkpoint_jobs import *

try:
    fp = open("locations.txt")
except: 
    sys.stderr.println("locations.txt file not present.")
    exit(1);

location = fp.readline().strip()
if not (location.lower() in {'redhawk', 'oakley', 'osx'}):
    sys.stderr.println("location.txt: Bad content")
    exit(1)


wait_time = 100        # Amount of time to spin before checking job progress in a wait() call.
                       # Need to be higher on Oakley.
sleep_pause = 60
    
#################################################################
# The following global variables are related to debugging issues.
show_progress = False
stats_only = False
job_index = {}
default_time_limit = "4:00:00"
#default_time_limit = "00:20:00"
rm_time_limit = "25:00:00"
#rm_time_limit = "2:00:00"

time_limit = default_time_limit

timing = False
timing_jobs = False
start_time = None
quit_time = None
prog_walltime = None
safety_margin = None
continue_prev = False
check_fname = 'reval.dat'
log_fname = 'reval.log'

flist_start = "START_FILE_LIST"
flist_end = "END_FILE_LIST"
csjobs_start = "START_CHROM_SIM_JOBS"
csjobs_end = "END_CHROM_SIM_JOBS"
tjobs_start = "START_TOOL_JOBS"
tjobs_end = "END_TOOL_JOBS"
rmjobs_start = "START_REPMASK_JOBS"
rmjobs_end = "END_REPMASK_JOBS"
prajobs_start = "START_PRA_JOBS"
prajobs_end = "END_PRA_JOBS"
jobdic_start = "START_JOB_DICT"
jobdic_end = "END_JOB_DICT"
blast_db_start = "START_BLAST_DB"
blast_db_end = "END_BLAST_DB"
stats_start = "START_STATS_JOBS"
stats_end = "END_STATS_JOBS"


#################################################################
def print_time():
    return time.strftime("%x %X", time.localtime())
                         


#################################################################
# These global variables have to do with executable locations.
MacLocations = {'build_lmer_table':'/usr/local/RepeatScout/build_lmer_table',
                'RptScout':'/usr/local/RepeatScout/RepeatScout',
                'filter_stage-1':'/usr/local/RepeatScout/filter-stage-1.prl',
                'filter_stage-2':'/usr/local/RepeatScout/filter-stage-2.prl',
                'raider':'./raider',
                'raider_pre':'./raider_pre',
                'bigfoot':'./bigfoot',
                'python':'python3.4',
                'araider':'./araider',
                'raider2': './phRAIDER',
                'rm_modules': None,
                'RepeatMasker' : 'RepeatMasker',
                'proc_per_node' : 1,
                'basic_arch_type' : None,
                'high_mem_arch' : None}
RedhawkLocations = {'build_lmer_table':'./build_lmer_table',
                    'RptScout':'./RepeatScout',
                    'filter_stage-1':'./filter-stage-1.prl',
                    'filter_stage-2':'./filter-stage-2.prl',
                    'raider':'./raider',
                    'raider_pre':'./raider_pre',
                    'bigfoot':'./bigfoot',
                    'python':'python3.3',
                    'araider':'./araider',
                    'raider2': './phRAIDER',
                    'rm_modules' : ['RepeatMasker', 'python-3.3.3'],
                    'RepeatMasker' : 'RepeatMasker',
                    'proc_per_node' : 4,
                    'basic_arch_type' : ["n09","bigmem"],
                    'high_mem_arch' : 'redhawk'}
OakleyLocations = {'build_lmer_table':'./build_lmer_table',
                   'RptScout':'./RepeatScout',
                   'filter_stage-1':'./filter-stage-1.prl',
                   'filter_stage-2':'./filter-stage-2.prl',
                   'raider':'./raider',
                   'raider_pre':'./raider_pre',
                   'bigfoot':'./bigfoot',
                   'python':'python',
                   'araider':'./araider',
                   'raider2': './phRAIDER',
                   'rm_modules' : None,
                   'RepeatMasker' : 'RepeatMasker',
                   'proc_per_node' : 12,
                   'basic_arch_type' : None,
                   'high_mem_arch' : 'oakley'}
Locations = None;    # This will be set to one of the above two, and references to find exectuable locations.


#########
# Utility functions
def sum_resources(T1, T2):
    if T1[0] == -1 or T2[0] == -1:
        return [-1]*4
    return [T1[0] + T2[0], T1[1] + T2[1], max(T1[2], T2[2]), max(T1[3], T2[3])] 

def get_job_index(s):
    global job_index
    if s not in job_index:
        job_index[s] = 0
    v = job_index[s]
    job_index[s] += 1
    return v

def file_base(file):
    return os.path.basename(file)

def file_dir(file):
    return file.rstrip(file_base(file)).rstrip("/")

def parse_redhawk_time(time_str):
    """Parse time limit string for redhawk (format HH:MM:SS) into seconds amount"""
    secs = sum(int(x) * 60 ** i for i,x in enumerate(reversed(time_str.split(":"))))
    #print(time_str, '/t', secs)
    return secs

def convert_seed(seed):
    """Convert an abriviated seed to a full seed (e.g. "1{2}0{3}1{2}" => "1100011" """
    i = 0
    while (i < len(seed)-1):
        if seed[i+1] == '^':
            j = i+2
            assert seed[j] == "{"
            k = j+1
            while seed[k] != '}':
                k += 1
            n = int(seed[j+1:k])
            seed = seed[:i] + seed[i]*n + seed[k+1:]
        i += 1
    return seed

def parse_params(args):
    """Parse command line arguments using the argparse library"""
    parser = argparse.ArgumentParser(description = "Evaluate RAIDER against RepeatScout")
    
    # GENERAL ARGUMENTS
    #parser2 = parser.add_mutually_exclusive_group()
    #parser2.add_argument('--organize', action = "store_true", help = "Create directory for all Raider Eval output", default = False)
    #parser2.add_argument('--no', '--named_organize', dest = "named_organize", help = "Organize under a named directory", default = None)


    # TOOL SELECTION
    parser_tools = parser.add_argument_group("tool selection (all on by default)")
    parser_tools.add_argument('-R', '--raider_on', dest = 'run_raider', action = 'store_true', help = 'Turn RAIDER on', default = False)
    parser_tools.add_argument('--R2', '--raider2_on', dest = 'run_raider2', action = 'store_true', help = 'Turn RAIDERV2 on', default = False)
    parser_tools.add_argument('--AR', '--araider_on', dest = 'run_araider', action = 'store_true', help = 'Turn ARAIDER on', default = False)
    parser_tools.add_argument('--RS', '--repscout_on', dest = 'run_repscout', action = 'store_true', help = 'Turn RAIDER on', default = False)
    parser_tools.add_argument('-B', '--bigfoot_on', dest = 'run_bigfoot', action = 'store_true', help = 'Turn BIGFOOT on', default = False)
    parser_tools.add_argument('-P', '--piler_on', dest = 'run_piler', action = 'store_true', help = 'Turn PILER on', default = False)
    parser_tools.add_argument('-A', '--all_tools', dest = 'all_tools', action = 'store_true', help = 'Turn all tools on (overide all other tool arguments)', default = False)
    parser_tools.add_argument('--A2', '--all_tools2', dest = 'all_tools2', action = 'store_true', help = 'Turn all tools on except araider (overide all other tool arguments)', default = False)
    parser_tools.add_argument('--tl', '--time_limit', dest = 'time_limit', help = 'Redhawk time limit (max: 400:00:00 default: 4:00:00)', default = default_time_limit)
    parser_tools.add_argument("--mn", '--max_nodes', dest = "max_nodes", action="store_true", help="Reserve all nodes of a processor for each tool (disabled by default).", default=False)

    # Will later add: RepeatModeler, RECON, PILER (other?)


    # I/O ARGUMENTs
    parser_io = parser.add_argument_group("i/o arguments")
    parser_io.add_argument('-r', '--results_dir', dest = "results_dir", help = "Directory containing all results", default = "EVAL")
    parser_io.add_argument('--nuke', dest ='nuke', action = "store_true", help = "Nuke the results directory", default = False)
    parser_io.add_argument('--rd', '--raider_dir', dest = "raider_dir", help = "Subdirectory containing raider results", default = "RAIDER")
    parser_io.add_argument('--ard', '--araider_dir', dest = "araider_dir", help = "Subdirectory containing araider results", default = "ARAIDER")
    parser_io.add_argument('--r2d', '--raider2_dir', dest = "raider2_dir", help = "Subdirectory containing araider results", default = "RAIDERV2")
    parser_io.add_argument('--rsd', '--rptscout_dir', dest = 'rptscout_dir', help = "Subdirectory containing rpt scout results", default = "RPT_SCT")
    parser_io.add_argument('--bfd', '--bigfoot_dir', dest = 'bigfoot_dir', help = "Subdirectory containing bigfoot results", default = "BIGFOOT")
    parser_io.add_argument('--pd', '--pilder_dir', dest = 'piler_dir', help = "Subdirectory containing piler results", default = "PILER")
    parser_io.add_argument('--dd', '--data_dir', dest = 'data_dir', help = "Directory containing the resulting simulated chromosome", default = "SOURCE_DATA")
    parser_tools.add_argument('--hj', '--hooke_jeeves', dest = 'hooke_jeeves', action = 'store_true', help = 'Simply print the tp+tn statistics counts', default = False)
    
    # RAIDER ARGUMENTS
    raider_argument = parser.add_argument_group("RAIDER parameters")
    raider_argument.add_argument('-f', type = int, help = "E.R. occurrence threshold", default = 5)
    raider_argument.add_argument('-d', '--output_dir', help = "Raider output directory", default = None)
    raider_argument.add_argument('-e', '--output_ext', help = "Output Extension", default = None)
    raider_argument.add_argument('-C', '--cleanup_off', dest = "cleanup", action = "store_false", help = "Turn off file cleanup", default = True)
    raider_argument.add_argument('--raider_min', '--raider_min', type = int, help = "Minimum repeat length. Defaults to pattern length.", default = None)
    raider_argument.add_argument('--pre', '--pre_scan', action = 'store_true', help = "Use pre-scan version of raider", default = False)
    raider_argument.add_argument('--mem', action = 'store_true', help = "Use large memory-nodes", default = False);
    seed_group = raider_argument.add_mutually_exclusive_group(required = False)     
    seed_group.add_argument('-s', '--seed', dest = "seed", help = "Spaced seed string", default = "111111111111111111111111111111")    
    seed_group.add_argument('--sf', '--seed_file', dest = 'seed_file', help = 'File containing raider seeds', default = None)
   
    # RAIDER2 ARGUMENTS
    raider2_argument = parser.add_argument_group("RAIDER2 parameters")
    raider2_argument.add_argument('--age', type = int, help="Use older version of raider2", default=1) 
    raider2_argument.add_argument('--aa', '--all_ages', dest="all_ages", action="store_true", help="Run all ages of raider2", default=False) # type = int, help="Use older version of raider", default=0)
    #raider2_argument.add_argument('--multi', '--multi_seed', dest="multi_seed", action="store_true", help="Run all seeds in seed file concurrently",default=False)
    raider2_argument.add_argument('--na', '--no_family_array', dest="family_array", action="store_false", help="Disable family array in Raider2", default=True)
    raider2_argument.add_argument('--ex', '--excise', dest="excising", action="store_true", help="Enable excising in RAIDER2", default=False)
    raider2_argument.add_argument('--no', '--no_overlaps', dest="overlaps", action="store_false", help="Do not require overlaps in RAIDER2", default=True)
    raider2_argument.add_argument('--tu', '--tie_up', dest="tieup", action="store_true", help="Enable alternative tie ups", default=False)
    raider2_argument.add_argument('--ps', '--prosplit', dest="prosplit", action="store_true", help="Enable proactive splitting(disabled by default).", default=False)
    raider2_argument.add_argument("--pf", '--prevfam', dest="prevfam", action="store_true", help="Enable pointers to prev family (disabled by default).", default=False)

    # REPSCOUT ARGUMENTS
    repscout_argument = parser.add_argument_group("REPSCOUT parameters")
    repscout_argument.add_argument('--repscout_min', type = int, help = "Minimum repeat length for repscout.", default = 10)
    repscout_argument.add_argument('--rs_min_freq', type = int, help = "Minimum repeat length for repscout.", default = 3)
    repscout_argument.add_argument('--rs_filters', type = int, dest = "rs_filters", help = "Specify how many RS filters to use {0,1,2}. 3 specifies to run all versions", default = 0)
    #raider_argument.add_argument('--uff', '--use_first_filter', dest = "use_first_filter", action = "store_true", help = "Use the first RepScout filter", default = True)
    #raider_argument.add_argument('--usf', '--use_second_filter', dest = "use_second_filter", action = "store_true", help = "Use the second RepScout filter", default = True)

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
    stats_group.add_argument('--stats_dir', dest = 'stats_dir', help = "Statistics output directory", default = "STATS_OUTPUT")
    stats_group.add_argument('--stats_file', dest = 'stats_file', help = "Statistics output file", default = "stats.txt")
    stats_group.add_argument('--stats_only', dest = 'stats_only', action = 'store_true', help = 'Remove files not involved in stats analysis', default = False)
    #stats_group.add_argument('--print_reps', action = "store_true", help = "Print out repeats in statistics file", default = False)

    # DEBUGGING ARGUMENTS
    debug_group = parser.add_argument_group(title = "debugging")
    debug_group.add_argument('--sp', '--show_progress', dest = 'show_progress', action = 'store_true', help = "Print reports on program progress to stderr", default = False)
    debug_group.add_argument('--so', '--simulate_only', dest = 'simulate_only', action = 'store_true', help = "Quit after creating simulated file", default = False)

    # ANALYSIS 
    parser_analysis = parser.add_argument_group("Analysis options")
    parser_analysis.add_argument('--PRA', '--pre_rm_analysis_off', dest = 'pra', action = 'store_false', help = 'Turn off pre-RM stats. analysis', default = True) 
    parser_analysis.add_argument('--RA', '--rm_analysis_off', dest = 'repmask', action = 'store_false', help = 'Turn off RM stats. analysis', default = True) 
    parser_analysis.add_argument('--ce', '--class_exclude', dest = 'exclude', action = 'store', help = 'File of family classes to exclude from PRA analysis', default = 'exclude.txt')
    ### KARRO END


    # CONTINUE PREVIOUS RUN ARGUMENTS
    cont_group = parser.add_argument_group(title = "continuing previous")
    cont_group.add_argument('--timing_jobs', dest = 'timing_jobs', action = 'store_true', help = "Set up timed jobs", default = False)
    cont_group.add_argument('--pwt', '--prog_walltime', dest = 'prog_walltime', help = 'Redhawk time limit for program', default = None)
    cont_group.add_argument('--cp', '--continue_prev', dest = 'continue_prev', action = 'store_true', help = "Continue previously started job", default = False)
    cont_group.add_argument('--sm', '--safe_marg', dest = 'safety_margin', help = "Amount of time left on clock (secs) when start to save run state", default = None)

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
    parser_chrom.add_argument('--st', '--sim_type', dest = 'sim_type', type = int, help = "0 = use mdern sequence exactly (default); 1 = use ancestor fragment; 2 = preserve mutations, but not indels; 3 = preserve indels, but not mutations", default = 0)
    parser_chrom.add_argument('chromosome', help = "Template chromosome file")



    arg_return =  parser.parse_args(args)

    global time_limit
    time_limit = arg_return.time_limit

    global show_progress
    show_progress = arg_return.show_progress

    global stats_only
    stats_only = arg_return.stats_only


 
    ###
    # Update global vars related to continuing previous jobs
    global prog_walltime
    prog_walltime = parse_redhawk_time(arg_return.prog_walltime) if arg_return.prog_walltime else None

    global timing
    timing = True if prog_walltime else False

    global timing_jobs
    timing_jobs = True if timing else arg_return.timing_jobs

    global safety_margin
    safety_margin = arg_return.safety_margin if arg_return.safety_margin else prog_walltime/10.0 if timing else None

    global continue_prev
    continue_prev = arg_return.continue_prev if timing else False

    global check_fname
    check_fname = arg_return.results_dir + "/" + check_fname
    
    global log_fname
    log_fname = arg_return.results_dir + "/" + log_fname

    if arg_return.all_tools or arg_return.all_tools2:
        arg_return.run_raider = True
        arg_return.run_repscout = True
        arg_return.run_piler = True

    if arg_return.all_tools:
        arg_return.run_araider = True
        arg_return.run_raider2 = True
        arg_return.run_bigfoot = True


    #### The following is to set the global debugging variables 
    if arg_return.simulate_only:    # Set to supress all tools
        arg_return.run_raider = False
        arg_return.run_araider = False
        arg_return.run_raider2 = False
        arg_return.run_repscout = False
        arg_return.run_bigfoot = False
        arg_return.run_piler = False


    return arg_return


############################################################
# Main functions 
def simulate_chromosome(chromosome_file, rng_seed, length, neg_strand, fam_file, data_dir, output_file, file_index, k, mc_file, mi, retain_n, num_repeats, low_complexity, sim_type):
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
    seq_arg = chromosome_file
    repeat_arg = chromosome_file + ".out"

    output_file = (output_file if output_file else re.sub(".fa$", ".sim.%d.fa" % (file_index), file_base(chromosome_file)))
    output_path = "%s/%s" % (data_dir, output_file)

    mc = "--mc %s" % mc_file if mc_file else ""

    if k == 0:  # Really only for debugging
        cmd = "{python} chromosome_simulator3.py {mi} {length} {mc} {k} {seed} {neg} {fam} {retain_n} {num_repeats} {lc} {seq} {repeat} {output}".format(python = Locations['python'], mi=mi, mc=mc, length=length_arg, k=k_arg, seed=seed_arg, neg=neg_arg, fam=fam_arg, retain_n=retain_n, num_repeats=num_repeats, lc=low_complexity, seq = seq_arg, repeat=repeat_arg, output=output_path)
    elif sim_type == 0 and os.path.isfile(repeat_arg):
        cmd = "{python} chromosome_simulator.py {mi} {length} {mc} {k} {seed} {neg} {fam} {retain_n} {num_repeats} {lc} {seq} {repeat} {output}".format(python = Locations['python'], mi=mi, mc=mc, length=length_arg, k=k_arg, seed=seed_arg, neg=neg_arg, fam=fam_arg, retain_n=retain_n, num_repeats=num_repeats, lc=low_complexity, seq = seq_arg, repeat=repeat_arg, output=output_path)
    else:
        sim_type = "--st %d" % (sim_type)
        cmd = "{python} chromosome_simulator2.py {sim_type} {mi} {length} {mc} {k} {seed} {neg} {fam} {retain_n} {num_repeats} {lc} {seq} {output}".format(python = Locations['python'], sim_type = sim_type, mi=mi, mc=mc, length=length_arg, k=k_arg, seed=seed_arg, neg=neg_arg, fam=fam_arg, retain_n=retain_n, num_repeats=num_repeats, lc=low_complexity, seq=seq_arg, output=output_path)


    if show_progress:
        sys.stderr.write("Creating simulation:\n%s\n" % (cmd))
        sys.stderr.flush()
    ##progress_fp.write(print_time() + "\n")
    progress_fp.write("Creating simulation:\n%s\n" % (cmd))
    progress_fp.flush()

    batch_name = data_dir + "/" + output_file + ".sim.batch"
    job_name = "simulation.%d" % (get_job_index("simulation"))
    
    p = pbsJobHandler(batch_file = batch_name, executable = cmd, job_name = job_name,
                      stdout_file = output_file + ".stdout", stderr_file = output_file + ".stderr", 
                      output_location = data_dir, walltime = time_limit, arch_type = Locations['basic_arch_type'])
    if not timing_jobs:
        p.submit(preserve=True, delay = wait_time)
    else:
        p.submit_timed_job(preserve=True, delay = wait_time)

    p.output_file = output_file
    p.seq_file = file_base(output_file)
    p.sim_output = output_path
    p.index = file_index
    return p


def run_raider(seed, seed_num, f, m, input_file, raider_dir, mem, max_nodes):
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

    ##progress_fp.write(print_time() + "\n")
    progress_fp.write("\nLaunching raider:\n%s\n%s\n" % (cmd1, cmd2))
    progress_fp.flush()

    batch_name = raider_dir + "/" + input_base + ".raider.batch"
    job_name = "R.{input}.{seed}.{num}".format( num = get_job_index("raider") , input=re.sub("hg18.","",input_base), seed=seed_num)
    p = pbsJobHandler(batch_file = batch_name, executable = cmd1 + "; " + cmd2, job_name = job_name,
                      stdout_file = input_base + ".raider.stdout", stderr_file = input_base + ".raider.stderr",
                      output_location = output_dir, walltime = time_limit, mem = Locations['high_mem_arch'] if mem else False, ppn = Locations['proc_per_node'] if max_nodes else 1,
                      arch_type = Locations['basic_arch_type'] if not mem else False)
    if not timing_jobs:
        p.submit(preserve=True, delay = wait_time)
    else:
        p.submit_timed_job(preserve=True, delay = wait_time)

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

    ##progress_fp.write(print_time() + "\n")
    progress_fp.write("\nLaunching composites finder:\n%s\n%s\n" % (cmd1, cmd2))
    progress_fp.flush()

    batch_name =  compositesFinderDir + "/" + input_base + ".composites finder.batch"
    job_name = "composites finder.%d" % get_job_index("composites finder")
    p = pbsJobHandler(batch_file = batch_name, executable = cmd1 + "; " + cmd2, job_name = job_name,
                      stdout_file = input_base + ".comFinder.stdout", stderr_file = input_base + ".comFinder.stderr",
                      output_location = output_dir, walltime = time_limit)

    if not timing_jobs:
        p.submit(preserve=True, delay = wait_time)
    else:
        p.submit_timed_job(preserve=True, delay = wait_time)

    p.description = "composites.finder"
    p.elementsFile = elements_file
    p.seqFile = seq_file

    return p

def run_raider2(seed, seed_num, f, m, input_file, raider2_dir, family_array, excise, overlaps, tieup, prosplit, prevfam, age, age_only, max_nodes, mem):
    """Given raider parameters and an input file, run RAIDER and put the output into
    the directory specified in output_dir (creating a random name is none is
    specified."""

    input_base = file_base(input_file).rstrip(".fa")
    #raider2_dir += "NO_FA." if not family_array else "FA."
    #raider2_dir += "EXC." if excise else "NO_EXC."
    #raider2_dir += "NO_OV." if not overlaps else "OV."
    #raider2_dir += "TU." if tieup else "NO_TU."
    #raider2_dir += "PS" if prosplit else "NO_PS."
    #raider2_dir += "PF" if prevfam else "NO_PF"
    output_dir = raider2_dir + "/" + input_base.upper() + ".s" + str(seed_num)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    min_arg = "-m %d" % (m) if m else ""
    
    if type(seed) is list:
        seed_string = "-s " + " -s ".join(seed)
    else:
        seed_string = "-s {seed}".format(seed = seed)

    opt_str = ""
    if not age_only:
        opt_str += "--na " if not family_array else ""
        opt_str += "--e " if excise else ""
        opt_str += "--no " if not overlaps else ""
        opt_str += "--t " if tieup else ""
        opt_str += "--ps " if prosplit else ""
        opt_str += "--pf " if prevfam else ""
    else:
        opt_str += "--age " + str(age) 

    cmd1 = "{raider2} -q -c {f} {version} {min_arg} {seed} {input_file} {output_dir}".format(raider2 = Locations['raider2'], f = f, version = opt_str, min_arg = min_arg, seed = seed_string, input_file = input_file, output_dir = output_dir)

    out_file = raider2_dir + "/" + input_base + ".s" + str(seed_num) + ".raider2_consensus.txt"
    lib_file = raider2_dir + "/" + input_base + ".s" + str(seed_num) + ".raider2_consensus.fa"

    cmd2 = "{python} consensus_seq.py -s {seq_file} -e {elements_dir}/elements {output_file} {fa_file}".format(python = Locations['python'], seq_file = input_file, elements_dir = output_dir, output_file = out_file, fa_file = lib_file)

    element_file = output_dir + "/elements"
    family_file = output_dir + "/families"

    cmd3 = "rm {elements}; rm {family}".format(elements = element_file, family = family_file ) if stats_only else ""
    
    if show_progress:
        sys.stderr.write("\nLaunching raider2:\n%s\n%s\n\n" % (cmd1, cmd2))
        sys.stderr.flush()

    #progress_fp.write(print_time() + "\n")
    progress_fp.write("Launching raider2:\n%s\n%s\n\n" % (cmd1, cmd2))
    progress_fp.flush()

    batch_name = raider2_dir + "/" + input_base + ".s" + str(seed_num) +  ".raider2.batch"
    job_name = "R2.{input}.{seed}.{num}".format( num = get_job_index("raider2") , input=re.sub("hg18.","",input_base), seed=seed_num)
    p = pbsJobHandler(batch_file = batch_name, executable = cmd1 + "; " + cmd2 + "; " + cmd3, job_name = job_name,
                      stdout_file = input_base + ".raider2.stdout", stderr_file = input_base + ".raider2.stderr",
                      output_location = output_dir, walltime= time_limit, ppn = Locations['proc_per_node'] if max_nodes else 1, mem = Locations['high_mem_arch'] if mem else False)

    if not timing_jobs:
        p.submit(preserve=True, delay = wait_time)
    else:
        p.submit_timed_job(preserve=True, delay = wait_time)

    p.tool_resources = [0]*4

    p.description = "raider2"
    p.tools_resources = [0]*4
    p.seed = seed
    p.seed_num = seed_num
    p.seq_file = input_file
    p.lib_file = lib_file

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
    
    element_file = output_dir + "/elements"
    family_file = output_dir + "/families"

    cmd3 = "rm {elements}; rm {family}".format(elements = element_file, family = family_file ) if stats_only else ""

    if show_progress:
        sys.stderr.write("\nLaunching araider:\n%s\n%s\n" % (cmd1, cmd2))
        sys.stderr.flush()

    #progress_fp.write(print_time() + "\n")
    progress_fp.write("\nLaunching araider:\n%s\n%s\n" % (cmd1, cmd2))
    progress_fp.flush()

    batch_name = araider_dir + "/" + input_base + ".araider.batch"
    job_name = "araider.%d" % get_job_index("araider")
    p = pbsJobHandler(batch_file = batch_name, executable = cmd1 + "; " + cmd2 + "; " + cmd3, job_name = job_name,
                      stdout_file = input_base + ".araider.stdout", stderr_file = input_base + ".araider.stderr",
                      output_location = output_dir, walltime = time_limit)

    if not timing_jobs:
        p.submit(preserve=True, delay = wait_time)
    else:
        p.submit_timed_job(preserve=True, delay = wait_time)

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

    #progress_fp.write(print_time() + "\n")
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


    if not timing_jobs:
        p.submit(preserve=True, delay = wait_time)
    else:
        p.submit_timed_job(preserve=True, delay = wait_time)

    p.description = "bigfoot"
    p.tool_resources = [0]*4
    p.seq_file = input_file     # Required by run_repeat_masker -- uses this as the source sequence.
    p.lib_file = lib_file     # This should be set to the file name that will be the library for the repeatmasker run

    return p 

def run_piler(input_file, piler_dir, max_nodes):
    """Runs Piler and returns a submitted pbs object with specific attributes used to run RepeatMasker."""
    input_base = file_base(input_file).rstrip(".fa")
    lib_file = input_base + ".lib"
    if not os.path.exists(piler_dir):
        os.makedirs(piler_dir)


    cmd = "{python} run_piler.py {input_file} {piler_dir} {output_file}".format(python = Locations['python'], input_file = input_file, piler_dir = piler_dir, output_file = lib_file);
    
    if show_progress:
        sys.stderr.write("\nLaunching Piler:\n%s\n" % (cmd))
        sys.stderr.flush()

    #progress_fp.write(print_time() + "\n")
    progress_fp.write("\nLaunching Piler:\n%s\n" % (cmd))
    progress_fp.flush()

    batch_name = piler_dir + "/" + input_base + ".piler.batch";
    job_name = "piler%d" % get_job_index("piler")
    stdout_file = input_base + ".piler.stdout";
    stderr_file = input_base + ".piler.stderr";
    p = pbsJobHandler(batch_file = batch_name, executable = cmd, job_name = job_name,
                      stdout_file = stdout_file, stderr_file = stderr_file,
                      output_location = piler_dir, walltime = time_limit, ppn = Locations['proc_per_node'] if max_nodes else 1)

    if not timing_jobs:
        p.submit(preserve=True, delay = wait_time)
    else:
        p.submit_timed_job(preserve=True, delay = wait_time)

    p.description = "piler"
    p.tool_resources = [0]*4
    p.seq_file = input_file
    p.lib_file = piler_dir + "/" + lib_file

    return p


def run_scout(input_file, output_dir, min_freq, length, use_first_filter, use_second_filter, threshold, max_nodes, mem):
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
        cmd3 = "{filter} {input} > {filter_output}".format(input=rptscout_output, filter = Locations['filter_stage-1'], filter_output = filter_stage_output)
    else:
        cmd3 = ""

    if show_progress:
        sys.stderr.write("\nRepeatScout:\n%s\n%s\n%s\n" % (cmd1, cmd2, cmd3))
        sys.stderr.flush()

    #progress_fp.write(print_time() + "\n")
    progress_fp.write("\nRepeatScout:\n%s\n%s\n%s\n" % (cmd1, cmd2, cmd3))
    progress_fp.flush()
        
    batch_name = output_dir + "/" + file_base(input_file) + ".repscout1.batch"
    job_name = "rptscout.{input}.{num}".format( num = get_job_index("repscout") , input=file_base(input_file))
    p = pbsJobHandler(batch_file = batch_name, executable = cmd1 + "; " + cmd2 + "; " + cmd3, job_name = job_name, RHmodules = Locations['rm_modules'],
                      stdout_file = file_base(rptscout_output) + ".stdout", stderr_file = file_base(rptscout_output) + ".stderr",
                      output_location = output_dir, walltime = time_limit, arch_type = Locations['basic_arch_type'] if not mem else [], ppn = Locations['proc_per_node'] if max_nodes else 1,
                      mem = Locations['high_mem_arch'] if mem else False)

    if not timing_jobs:
        p.submit(preserve=True, delay = wait_time)
    else:
        p.submit_timed_job(preserve=True, delay = wait_time)
    p.description = "rep_scout" if not use_first_filter else "rep_scout1" if not use_second_filter else "rep_scout12"
    p.tool_resources = [0,0,0,0]

    p.seq_file = input_file
    p.should_filter_stage2 = use_second_filter
    p.input_name= input_name
    p.threshold = threshold
    p.stage = "1"
    p.lib_file = filter_stage_output if use_first_filter else rptscout_output
    return p


def run_scout_second_filter_RM(p, num_processors):
    if p.should_filter_stage2:
        p2 = run_repeat_masker(p, num_processors)
        p2.should_filter_stage2 = p.should_filter_stage2
        p2.threshold = p.threshold
        p2.input_name = p.input_name
        p2.stage = "RM"
        p2.description = p.description#"rep_scout"
        #print("F1 resources : " + str(p2.tool_resources))
        return p2
        
    else:
        return p

def run_scout_second_filter(p):
    if p.should_filter_stage2:
        filter_stage2_output = p.dir + "/" + p.input_name.rstrip(".fa") + ".repscout.filtered2.fa"
        cmd = "cat {filtered}  | {filter} --cat={rm_output} --thresh={thresh} > {filter_output}".format(filtered=p.lib_file, filter=Locations['filter_stage-2'], 
            filter_output= filter_stage2_output, thresh=p.threshold, rm_output = p.rm_output) #RM_output_dir + "/" + file_base(input_file) + ".out")
        if show_progress:
            sys.stderr.write("\nRepeatScout Filter2:\n%s\n" % cmd)
            sys.stderr.flush()

        #progress_fp.write(print_time() + "\n")
        progress_fp.write("\nRepeatScout Filter2:\n%s\n" % cmd)
        progress_fp.flush()

        batch_name = file_dir(p.rm_output) + "/" + p.input_name.rstrip(".fa") + ".repscout2.fa"
        job_name = "filter2.%d" % get_job_index("filter2")

        p2 = pbsJobHandler(batch_file = batch_name, executable = cmd, job_name = job_name,
                       stdout_file = file_base(p.seq_file) + ".repscout2.stdout", stderr_file = file_base(p.seq_file) + ".repscout2.stderr",
                       output_location = file_dir(p.seq_file), walltime = time_limit)
        if not timing_jobs:
            p2.submit(preserve=True, delay = wait_time)
        else:
            p2.submit_timed_job(preserve=True, delay = wait_time)
        p2.description = p.description#"rep_scout"
        p2.stage = "2"
        #print("RM resources : " + str(p.getResources(cleanup=False)))

        p2.tool_resources = [x + y for x, y in zip(p.tool_resources, p.getResources(cleanup=False))]
        #print("F1 + RM resources : " + str(p2.tool_resources))
        p2.seq_file = p.seq_file
        p2.lib_file = filter_stage2_output
        p2.should_filter_stage2 = p.should_filter_stage2
        return p2
    else:
        return p

def scout_second_filter(p, min_freq):
    """NOT CURRENTLY WORKING!!! Does not run correctly, and does not properly adjust time"""
    p.wait(wait_time)

    filter2_stage_output = p.seq_file.rstrip(".fa") + ".repscout.filtered2.fa"
    cmd = "cat {output} | perl {filter} --cat={cat} --thresh={thresh} > {final}".format(output = p.lib_file, filter = Locations['filter_stage-2'], cat = p.rm_output, thresh = min_freq, final = filter2_stage_output)
    
    if show_progress:
        sys.stderr.write("\nRepeatScout Filter2:\n%s\n" % cmd)
        sys.stderr.flush()

    #progress_fp.write(print_time() + "\n")
    progress_fp.write("\nRepeatScout Filter2:\n%s\n" % cmd)
    progress_fp.flush()

    batch_name = file_dir(p.rm_output) + "/" + file_base(p.seq_file).rstrip(".fa") + ".repscout2.fa"
    job_name = "filter2%d" % get_job_index("filter2")
    p2 = pbsJobHandler(batch_file = batch_name, executable = cmd, job_name = job_name,
                       stdout_file = file_base(p.seq_file) + ".repscout2.stdout", stderr_file = file_base(p.seq_file) + ".repscout2.stderr",
                       output_location = file_dir(p.seq_file), walltime = time_limit)
    if not timing_jobs:
        p2.submit(preserve=True, delay = wait_time)
    else:
        p2.submit_timed_job(preserve=True, delay = wait_time)
    p2.description = "rep_scout"
    p2.time_resources = p.time_resources + p.getResources(cleanup=False)
    p2.lib_file = p.lib_file
    p2.seq_file = p.seq_file
    p2.lib_file = filter2_stage_output

    return p2



def run_repeat_masker(p, num_processors):
    """Given the pbs object used to start a consensus sequence job as well as
    repeatmasker arguments, wait until the job is done and then call repeatmasker 
    on the output and put results in masker_dir (current dir if unspecified)"""
    p.wait(wait_time)
    p.loadResources()

    input_base = file_base(p.seq_file)  # Base name of the file used for input

    output_dir = file_dir(p.lib_file) + "/" + file_base(p.lib_file).upper() + ".RM"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    cmd = "{RepeatMasker} -nolow -lib {library} -pa {pa} -dir {dir} {seq_file}".format(RepeatMasker = Locations['RepeatMasker'], library = p.lib_file, pa = num_processors, dir = output_dir, seq_file = p.seq_file)

    if show_progress:
        sys.stderr.write("\nLaunch repeatmasker (%s):\n%s\n" % (p.description, cmd))
        sys.stderr.flush()

    #progress_fp.write(print_time() + "\n")
    progress_fp.write("\nLaunch repeatmasker (%s):\n%s\n" % (p.description, cmd))
    progress_fp.flush()

    batch_name = p.lib_file.rstrip(".fa") + ".rm.batch"

    tool_name = p.description if not "raider" in p.description else "R" if p.description == "raider" else "R2" 
    job_name = "RM.{input}.{tool}.{num}".format( num = get_job_index("repmask") , input=re.sub("hg18.","",input_base), tool=tool_name)
    ppn_arg = 4*num_processors if num_processors == 1 else num_processors

    p2 = pbsJobHandler(batch_file = batch_name, executable = cmd, nodes = 1, ppn = ppn_arg, RHmodules = ["RepeatMasker", "python-3.3.3"],
                       job_name = job_name, stdout_file = input_base + ".repmask.stdout", stderr_file = input_base + ".repmask.stderr",
                       output_location = output_dir, walltime = rm_time_limit, always_outputs=False);
    if not timing_jobs:
        p2.submit(preserve=True, delay = wait_time)
    else:
        p2.submit_timed_job(preserve=True, delay = wait_time)

    p2.description = "RptMasker"
    p2.seed = p.seed if hasattr(p, "seed") else "NA"
    p2.seed_num = p.seed_num if hasattr(p, "seed_num") else "NA"
    p2.dir = output_dir
    p2.lib_file = p.lib_file
    p2.seq_file = p.seq_file
    p2.rm_output = output_dir + "/" + file_base(p.seq_file) + ".out"
    p2.tool_resources = [x + y for x, y in zip(p.tool_resources, p.getResources(cleanup=False))]
    p2.tool_description = p.description
    return p2


def run_perform_stats(p, exclusion_file = None):
    """Given the pbs object used to start a consensus sequence job as well as
    repeatmasker arguments, wait until the job is done and then call repeatmasker 
    on the output and put results in masker_dir (current dir if unspecified)"""
    p.wait(wait_time)
    p.loadResources()

    input_base = file_base(p.seq_file)  # Base name of the file used for input

    known_repeats = p.seq_file + ".out"
    found_repeats = p.rm_output
    output_dir = p.dir
    output_path = output_dir + "/" + file_base(p.seq_file) + ".stats"
    exclusion_part = "-e {exclude}".format(exclude = exclusion_file) if exclusion_file else ""
    cmd = "python perform_stats.py {exclusions} {known} {found} {output}".format(exclusions=exclusion_part, known = known_repeats, found = found_repeats, output = output_path)

    if show_progress:
        sys.stderr.write("\nLaunch perform_stats: %s\n" % (cmd))
        sys.stderr.flush()

    #progress_fp.write(print_time() + "\n")
    progress_fp.write("\nLaunch perform_stats: %s\n" % (cmd))
    progress_fp.flush()

    batch_name = p.lib_file.rstrip(".fa") + ".stats.batch"
    job_name = "stats.%d" % get_job_index("stats")
    p2 = pbsJobHandler(batch_file = batch_name, executable = cmd,
                       job_name = job_name, stdout_file = input_base + ".stats.stdout", stderr_file = input_base + ".stats.stderr",
                       output_location = output_dir, walltime = time_limit, arch_type = Locations['basic_arch_type'])
    if not timing_jobs:
        p2.submit(preserve=True, delay = wait_time)
    else:
        p2.submit_timed_job(preserve=True, delay = wait_time)

    p2.description = "Stats"
    p2.seed = p.seed if hasattr(p, "seed") else "NA"
    p2.seed_num = p.seed_num if hasattr(p, "seed_num") else "NA"
    p2.dir = p.dir
    p2.lib_file = p.lib_file
    p2.seq_file = p.seq_file
    p2.rm_output = p.rm_output #output_dir + "/" + file_base(p.seq_file) + ".out"
    p2.tool_resources = p.tool_resources
    p2.tool_description = p.tool_description
    return p2
    

#def performance_sum(job_dic, PRA_jobs):
#    """Given a list of all of the statistics jobs, uses the statistics output files to
#    generate a summary file indicative of overall performance. Put results in stats_dir
#    (Current dir if unspecified)"""
#    ######
#    # Calculate statistics (not bothering with parallelization yet)
#    print_str = "{:<12}" + "{:<5}" + "".join("{:<14}"*4) + "".join("{:<14}"*6) + "".join("{:<14}"*8) + "{:<14}" + "\n"
#    stats_jobs = set()
#    for key in test_tools:
#        for p in job_dic[key]:
#            stats_jobs.add(run_perform_stats(p))
#
#
#    with open(args.results_dir + "/" + args.stats_file, "w") as fp:
#        fp.write(print_str.format("#tool", "seed", "tp", "fp", "fn", "tn", "tpr", "tnr", "ppv", "npv", "fpr", "fdr","ToolCpuTime", "ToolWallTime", "ToolMem", "ToolVMem", "RMCpuTime", "RMWallTime", "RMMem", "RMVMem", "coverage"))
#        
#        for key in test_tools:
#            for p in job_dic[key]:
#                try:
#                
#                    s = run_perform_stats(p)
#
#                except Exception as E:
#                    progress_fp.write("performance Exception: " + str(E) + "\n");
#                    fp.write("\t".join([str(key), str(p.seed_num) if hasattr(p, "seed_num") else "NA", "INCOMPLETE\n"]))
#
#    ### KARRO
#    # Finally: we should not terminate until all the pra jobs are done.  (If pra is off, this list will be empty.)
#    for p in PRA_jobs:
#        p.timed_wait()        # KARRO: Is this the correct method to use to ensure resubmission of needed
#    ### KARRO END
#    regex = re.compile("(?<=\# Average consensus coverage: )\d+.\d+")
#    >>> m = regex.findall(text)
#    >>> m
#    ['0.0063']
#
#    tps = 0
#    tns = 0
#    fps = 0
#    fns = 0
#    for p in stats_jobs:
#        sf = open(p.stats_output, "r")
#        tps += int(re.split("\s+", sf.readline().rstrip())[1])
#        fps += int(re.split("\s+", sf.readline().rstrip())[1])
#        tns += int(re.split("\s+", sf.readline().rstrip())[1])
#        fns += int(re.split("\s+", sf.readline().rstrip())[1])
#        sf.close()
#    stats_file = "summary.%s.stats" % (test)
#    smry_path = "%s/%s" % (stats_dir, stats_file) if stats_dir else "%s/%s" %(curr_dir, stats_file) if curr_dir else stats_file
#    
#    smry = open(smry_path, 'w')
#    smry.write("Evaluation completed.\n")
#    smry.write("True Positives (TP): \t %d \n" % (tps))
#    smry.write("False Positives (FP): \t %d \n" % (fps))
#    smry.write("True Negatives (TN): \t %d \n" % (tns))
#    smry.write("False Negatives (FN): \t %d \n" % (fns))
#    smry.write("\nPerformance of Repeat Classification Tool\n")
#    smry.write("Sensitivity (TPR): \t\t %f %%\n" % (tps/(tps + fns)))
#    smry.write("Specificity (TNR): \t\t %f %%\n" % (tns/(fps + tns)))
#    smry.write("Precision (PPV): \t\t %f %%\n" % (tps/(tps + fps)))
#    smry.write("Neg. Pred. Val. (NPV): \t %f %%\n" % (tns/(tns + fns)))
#    smry.write("Fall-Out (FPR): \t\t %f %%\n" % (fps/(fps + tns)))
#    smry.write("False Disc. Rate (FDR): \t  %f %%\n" % (fps/(tps + fps)))
#    smry.close()

# def run_pra_analysis(jobs, BLAST_DATABASE):
#     """This takes a list of the jobs, and a list of the BLAST_DATABASE jobs.  For each 
#     it makes sure the BLAST_DATABASE job is done for the corresponding datafile, then launches 
#     the analysis tool.  Returns a list of the analysis tool jobs.  RepeatMasker is NOT dependent
#     on these jobs -- it can be launched immediately."""

#     submitted_jobs =[]

#     cmd = "{python} blast_consensus.py {consensus_file} {rm_fa_file} {database_file} {output_file}"
#     for j in jobs:
#         j.wait()
#         o = BLAST_DATABASE[j.seq_file]  # This is the create_database job that was launched on this sequence file.
#         o.wait()    # Wait until the DATABASE file has been created.  
#                     # (This should be parallilzed, but probably not worth the effort)

#         analysis_cmd = cmd.format(python = Locations["python"], consensus_file = j.lib_file, 
#                                   rm_fa_file = o.rm_seq_file, database_file = o.rm_seq_file, 
#                                   output_file = j.lib_file.rstrip(".fa") + ".pra.txt")
        
#         if show_progress:
#             sys.stderr.write("\nLaunching pre-rm analysis:\n%s\n" % (analysis_cmd))
#             sys.stderr.flush()
#         progress_fp.write("pre-rm analysis:\n%s\n" % (analysis_cmd))
#         progress_fp.flush()

#         job_name = "pra.%d" % get_job_index("pra")

#         base_name = file_base(j.lib_file)[:-3] + ".pra"
#         batch_name = base_name + ".batch"
#         stdout_file = base_name + ".stdout"
#         stderr_file = base_name + ".stderr"
#         location = file_dir(j.lib_file)
#         print("Y: " + base_name + " " + batch_name + " " + job_name + " " + stdout_file + " " + stderr_file + " " + analysis_cmd)
#         p = pbsJobHandler(batch_file = batch_name, executable = analysis_cmd, job_name = job_name,
#                           stdout_file = stdout_file, stderr_file = stderr_file,
#                           output_location = location, walltime = time_limit, RHmodules = ["blast+"]);
#         p.submit_timed_job()    # KARRO: What parameters should be used here for resubmission?
#         submitted_jobs.append(p)

#     return submitted_jobs

def create_blast_db(file_list):
    # Launch the create_database.py tool on each chromosome.  (Needed for pre-RM result analysis; not needed for any tool.)
    # BLAST_DATABASE will be a dictionary mapping the (simulated) chromosome to the pbs object
    BLAST_DATABASE = {}
    blast_database_command = "{python} ./create_database.py {seq_file} {rm_file} {out_file}"
    for i,file_name in enumerate(file_list):
        file = file_name.rstrip(".fa")
        seq_file = file_name
        rm_file = file + ".fa.out"
        rm_seq_file = file + ".rptseq.fa"
        cmd = blast_database_command.format(python = Locations['python'], seq_file = seq_file, rm_file = rm_file, out_file = rm_seq_file)

        if show_progress:
            sys.stderr.write("\nLaunching create_database: %s\n\n" % (cmd))
            sys.stderr.flush()

        #progress_fp.write(print_time() + "\n")
        progress_fp.write("\nLaunching create_database: %s\n\n" % (cmd))
        progress_fp.flush()

        batch_name = file + ".blast_db.batch"
        job_name = "create_db.%d" % (i)
        stdout_file = file + ".blast_db.stdout"
        stderr_file = file + ".blast_db.stderr"
        
        o =  pbsJobHandler(batch_file = batch_name, executable = cmd, job_name = job_name,
                           ppn = 2,   walltime = "00:20:00", RHmodules = ["blast+"],
                           stdout_file = stdout_file, stderr_file = stderr_file)

        o.seq_file = seq_file
        o.rm_file = rm_file
        o.rm_seq_file = rm_seq_file

        BLAST_DATABASE[file_name] = o

    for o in BLAST_DATABASE.values():
        #o.submit_timed_job()   # KARRO: Highly unlikely this will ever exceed 20 minutes (or even 5 minutes) -- so I just took the default parameters.
        if not timing_jobs:
            o.submit(preserve=True, delay = wait_time)
        else:
            o.submit_timed_job(preserve=True, delay = wait_time)

    return BLAST_DATABASE


### KARRO
def run_pra_analysis(tool_job, database_job):
    """This launchs a pra_analysis jobs and returns the job object.
    * tool_job: the pbsJob for one of the de novo search tool.  This will use tool_job.lib file as the query sequence set.
    * database_job: the database job for the sequence on which tool_job was run.
    This function will wait on completion of both jobs."""

    cmd = "./pra_analysis {consensus_file} {rm_fa_file} {database_file} {output_file}"
    database_job.wait(100)


    analysis_cmd = cmd.format(python = Locations["python"], consensus_file = tool_job.lib_file, 
                              rm_fa_file = database_job.rm_seq_file, database_file = database_job.rm_seq_file, 
                              output_file = tool_job.lib_file.rstrip(".fa") + ".pra.txt", walltime = time_limit)
        
    if show_progress:
        sys.stderr.write("\nLaunching pre-rm analysis:\n%s\n" % (analysis_cmd))
        sys.stderr.flush()

    #progress_fp.write(print_time() + "\n")
    progress_fp.write("\npre-rm analysis:\n%s\n" % (analysis_cmd))
    progress_fp.flush()

    job_name = "pra.%d" % get_job_index("pra")

    location = file_dir(tool_job.lib_file)
    base_name = file_base(tool_job.lib_file)[:-3] + ".pra"
    batch_name = location + "/" + base_name + ".batch"
    stdout_file = base_name + ".stdout"
    stderr_file = base_name + ".stderr"
    p = pbsJobHandler(batch_file = batch_name, executable = analysis_cmd, job_name = job_name,
                      stdout_file = stdout_file, stderr_file = stderr_file,
                      output_location = location, walltime = time_limit, RHmodules = ["blast+"]);
    #p.submit_timed_job()    # KARRO: What parameters should be used here for resubmission?
    if not timing_jobs:
        p.submit(preserve=True, delay = wait_time)
    else:
        p.submit_timed_job(preserve=True, delay = wait_time)

    p.description = "PraAnalysis"
    p.seed = tool_job.seed if hasattr(tool_job, "seed") else "NA"
    p.seed_num = tool_job.seed_num if hasattr(tool_job, "seed_num") else "NA"
    p.lib_file = tool_job.lib_file
    p.seq_file = tool_job.seq_file
    p.pra_output = tool_job.lib_file.rstrip(".fa") + ".pra.txt"
    #p.tool_resources = tool_job.getResources(cleanup=False)
    p.tool_resources = [x + y for x, y in zip(tool_job.tool_resources, tool_job.getResources(cleanup=False))] #j.tool_resources + j.getResources(cleanup = False)
    #p2.tool_resources = p.resources
    p.tool_description = tool_job.description
    return p
### KARRO END


def exit_now():
    """If we submitted the evaluation as a PBS job with a set walltime, use this to write to the logging
    file an indicator that the evalation program has not yet completed, then exit the program"""
    logging_fp.write("CONTINUE\n")
    sys.exit(0)

def have_time_for_another_run(last_run_time):
    """If we submitted the evaluation as a PBS job with a set walltime, check the amount
    of time left and see if we can fit in another time interval of length at least 'last_run_time'"""
    time_left = quit_time - time.time() - last_run_time
    if quit_time - time.time() - last_run_time >= 0:
        return True
    else:
        logging_fp.write("Running out of time. Only have {t} left. Dumping data to new checkpoint file\n".format(t=time_left))
        return False



def run_timed_chrom_sim_jobs(jobs, flist=[]):
    """Given a set of chromosome simulation jobs and a list of paths to finished simulation files. 
    Three cases: 
        (1) If we submitted the evaluation as a PBS job with a set walltime, keep checking to
            see if we have reached point to save work and exit.
        2) If we are 'timing_jobs', keep calling timed_wait on each job until all have finished.
            This ensures that if a job is running out of time it will resubmit with a longer
            walltime (all work for this is done in redhawk.py)
        3) Otherwise, call wait() as usual. 
    if/when all jobs complete, returns list of paths to resulting simulation files.""" 
    chrom_job_set = {j for j in jobs}
    finished_jobs = set()
    time_est = None
    for j in chrom_job_set: 
        t1 = time.time()
        if timing and not time_est:
            time_est = parse_redhawk_time(j.walltime)
            if not have_time_for_another_run(time_est):
                save_timed_chrom_sim_jobs(chrom_job_set, finished_jobs, flist)            
                exit_now()
        if not timing_jobs:
            j.wait(100)
        else:
            j.timed_wait(100)
        t2 = time.time()
        finished_jobs.add(j)
        time_est = t2 - t1
        # Make sure won't run out of time if continue to next iteration
        if timing and not have_time_for_another_run(time_est):
            save_timed_chrom_sim_jobs(chrom_job_set, finished_jobs, flist)            
            exit_now()
    return [j.sim_output for j in finished_jobs]

def run_timed_tool_jobs(jobs, run_rm, pa, run_pra, blast_db, RM_jobs=None, PRA_jobs=None):
    """Given a set of repeat finding tool jobs and repmask jobs (with pa info), keep track of
    what tool jobs have completed and submit corresponding repmask job upon tool job completion.
    We call isJobRunning on each tool job -- if we are 'timing_jobs', this information is saved
    in the job object and redhawk.py will handle whether jobs need to be resubmitted with more time.
    If/when all tool jobs complete, returns list of repmask jobs (some of which are still running). 
    Note: If we submitted the evaluation as a PBS job with a set walltime, keep checking to
    see if we have reached point to save work and exit."""
    
    job_set = {j for j in jobs}
    if not RM_jobs:
        RM_jobs = set()
    if not PRA_jobs:
        PRA_jobs = set()
    time_est = None
    while job_set:
        finished_jobs = set()
        added_jobs = set()
        t1 = time.time()
        for j in job_set:
            if not j.isJobRunning():
                finished_jobs.add(j)
                if "rep_scout" in j.description and j.should_filter_stage2 and j.stage != "2":
                    if j.stage == "1":
                        added_jobs.add(run_scout_second_filter_RM(j,pa))
                    elif j.stage == "RM":
                        added_jobs.add(run_scout_second_filter(j))
                else:
                    if "rep_scout" in j.description and j.stage == "2":
                        #print("F2 resources : " + str(list(j.getResources(cleanup=False))))
                        #j.tool_resources = j.tool_resources + j.getResources(cleanup=False) #should we include the RM inside RS in timing? #p.getResources(cleanup=False)
                        j.tool_resources = [x + y for x, y in zip(j.tool_resources, j.getResources(cleanup=False))] #j.tool_resources + j.getResources(cleanup = False)
                        #print("F1 + RM + F2 resources : " + str(j.tool_resources))
                    rm_job = None
                    pra_job = None
                    if run_rm:
                        rm_job = run_repeat_masker(j,pa)
                        RM_jobs.add(rm_job)
                    if run_pra:
                        pra_job = run_pra_analysis(j, blast_db[j.seq_file])
                        PRA_jobs.add(pra_job)
                    if rm_job:
                        rm_job.pra_job = pra_job if pra_job else None
                    if pra_job:
                        pra_job.rm_job = rm_job if rm_job else None
        job_set = job_set - finished_jobs
        job_set = job_set | added_jobs
        time_est = time.time() - t1
        if timing and not have_time_for_another_run(time_est):
            save_timed_tool_jobs(job_set, RM_jobs, PRA_jobs, blast_db)
            exit_now()

        time.sleep(sleep_pause)  # We have checked all the jobs; lets sleep for ten minutes before checking again.
                                 # Important for keeping time down on the head node -- required on Oakley.
    return RM_jobs, PRA_jobs



test_tools = ["raider", "bigfoot", "piler", "rep_scout", "rep_scout1", "rep_scout12", "araider", "raider2", "raider2.0", "raider2.1", "raider2.2"]  # List of implemented tools 
def run_timed_analysis_jobs(run_rm, run_pra, RM_jobs, PRA_jobs, results_dir , stats_jobs=None, job_dic=None):
    """Given a set of repmask jobs and a working job dictionary, keep track of what repmask jobs 
    have completed and add completed jobs to job dictionary under appropriate tool name.
    We call isJobRunning on each repmask job -- if we are 'timing_jobs', this information is saved
    in the job object and redhawk.py will handle whether jobs need to be resubmitted with more time.
    If/when all repmask jobs complete, returns job dict to compute statistics on results. 
    Note: If we submitted the evaluation as a PBS job with a set walltime, keep checking to
    see if we have reached point to save work and exit."""
    job_dic = job_dic if job_dic else {tool:[] for tool in test_tools}
    if not stats_jobs:
        stats_jobs = set()
    for tool in test_tools:
        if tool not in job_dic.keys():
            job_dic[tool] = []


    pra_job_set = {}
    rm_job_set = {}
    if RM_jobs or PRA_jobs:
        pra_job_set = {j for j in PRA_jobs}
        rm_job_set = {j for j in RM_jobs}
        #finished_jobs = set()
        while rm_job_set or pra_job_set:
            finished_rm_jobs = set()
            finished_pra_jobs = set()
            time_est = None
            t1 = time.time()
            for j in rm_job_set:
                if not j.isJobRunning():
                    job_dic[j.tool_description].append(j)
                    finished_rm_jobs.add(j)
                #else:
            for j in pra_job_set:
                if not j.isJobRunning():
                    if run_pra and not run_rm:
                        j.pra_resources = list(j.getResources(cleanup=False))
                        job_dic[j.tool_description].append(j)
                    finished_pra_jobs.add(j)
            t2 = time.time() 
            time_est = t2 - t1
            if timing and not have_time_for_another_run(time_est):
                save_timed_PRA_jobs(pra_job_set - finished_pra_jobs)
                save_timed_RM_jobs(rm_job_set - finished_rm_jobs, stats_jobs, results_dir, job_dic)
                exit_now()
            pra_job_set = pra_job_set - finished_pra_jobs
            rm_job_set = rm_job_set-finished_rm_jobs

    return job_dic, stats_jobs, pra_job_set




############################################################################################
if __name__ == "__main__":
    start = time.time()
    args = parse_params(sys.argv[1:])

    start_time = start
    quit_time = prog_walltime + start_time - safety_margin if prog_walltime else None 

    ####
    # Currently: We check for the RepeatScout executable at the location on my Mac; if 
    # found, we assume we are running on the Mac.  If not, we check for in the Redhawk 
    # location, and if found assume we are running on redhawk.  Otherwise we print and 
    # error and quit.

    if location.lower() == "oakley":
        Locations = OakleyLocations;
        assert 1 <= args.pa <= 2, "Make sure you set the --pa parameter to a value between 1 and 4 on oakley (%d)" % (args.pa)
    elif location.lower() == "osx":
        Locations = MacLocations
    elif location.lower() == "redhawk":
        Locations = RedhawkLocations
        assert 1 <= args.pa <= 2, "Make sure you set the --pa parameter to a value between 1 and 4 on redhawk (%d)" % (args.pa)
    else:
        sys.stderr.println("locations.txt file: bad content")
        exit(1);
    ###

    
    data_dir = args.results_dir + "/" + args.data_dir    
      
    if not continue_prev:
        if args.nuke:
            if os.path.exists(args.results_dir):
                subprocess.call("rm -r %s" % args.results_dir, shell = True)
        else:
            if os.path.exists(args.results_dir):
                sys.stderr.write("%s exists; need to use --nuke option" % args.results_dir)
        if not os.path.exists(args.results_dir):
            os.makedirs(args.results_dir)
        if timing:
            checkpoint_fp = open(check_fname, "w")
            logging_fp = open(log_fname, "w")
            logging_fp.write("Starting new run from scratch\n")
        
        ### Generate simulated file(s) and run to completion
        ### Set up the debugging log file (if needed)
        progress_fp = open(args.results_dir + "/debug.txt", "w")
        #progress_fp.write(print_time() + "\n")        
        progress_fp.write(" ".join(sys.argv) + "\n\n");

        if not os.path.exists(data_dir):
            os.makedirs(data_dir)

        # First: we put the chromosomes (simulated or real) into data_dir
        if args.subparser_name == "chrom_sim":
            # Launch the jobs
            f = lambda i: simulate_chromosome(chromosome_file = args.chromosome, 
                                              rng_seed = args.rng_seed, length = args.length, 
                                              neg_strand = args.negative_strand, fam_file = args.family_file, 
                                              data_dir = args.results_dir + "/" + args.data_dir, output_file = args.output, file_index = i, 
                                              k = args.k, mc_file = args.mc_file, mi = args.max_interval,
                                              retain_n = args.retain_n, num_repeats = args.num_repeats, low_complexity = args.low_complexity, 
                                              sim_type = args.sim_type)
            J = [f(i) for i in range(args.num_sims)]
            if args.family_file:
                family_list = [fam for line in open(args.family_file) for fam in re.split("\s+", line.rstrip()) if fam]
                with open(args.results_dir + "/family_file.txt", "w") as fp:
                    fp.write("\n".join(["{fam}".format(fam=f) for f in family_list]) + "\n")
            # Get the list of simulated file names
            file_list = run_timed_chrom_sim_jobs(J) #[j.sim_output for j in J]

        else:
            file_list = []
            for file in args.seq_files:
                file_list.append(data_dir + "/" + file_base(file))
                shutil.copy(file, file_list[-1])
                shutil.copy(file + ".out", file_list[-1] + ".out")
        if timing:
            write_flist_to_checkpoint(file_list)
    

    else:
        if timing:
            if os.path.exists(check_fname):
                os.rename(check_fname, check_fname + ".old")
                old_checkpoint_fp = open(check_fname + ".old", "r")
            if os.path.exists(log_fname):
                os.rename(log_fname, log_fname + ".old")
            checkpoint_fp = open(check_fname, "w")
            logging_fp = open(log_fname, "w")
            logging_fp.write("Continuing previous run\n")
            flush_files()

        progress_fp = open(args.results_dir + "/debug.txt", "a")
        #progress_fp.write(print_time() + "\n")
        progress_fp.write(" ".join(sys.argv) + "\n\n");
        
        # Get file list from old checkpoint file
        file_list = []
        next_step = old_checkpoint_fp.readline().rstrip() 
        if next_step == flist_start:
            file_list, next_step = recover_file_list()

        if next_step == csjobs_start:
            chrom_job_set, next_step = recover_sim_jobs()
            run_timed_chrom_sim_jobs(chrom_job_set, file_list)
        flush_files()


    if not continue_prev or next_step == '':
        BLAST_DATABASE = create_blast_db(file_list) if args.pra else {} #CARLY: Moved this into a method to make main method (slightly) easier to follow
        if timing:
            write_blast_db_to_checkpoint(BLAST_DATABASE, args.results_dir)
    
    else:
        if next_step == blast_db_start:
            BLAST_DATABASE, next_step = recover_blast_db()
        if timing:
            write_blast_db_to_checkpoint(BLAST_DATABASE, args.results_dir)


    if not continue_prev or next_step == '':
        ### Start running each tool.  Each tool should run, creating the repeat masker library (putting the file name
        ### in the pbs lib_file attribute), then run repeat masker (putting the output file name in the pbs
        ### rm_output job.

        ############## Second: Launch tools
        ############## Need to initially launch all of the tool jobs
        jobs = []
        if args.pre:
            Locations['raider'] = Locations['raider_pre']
        if args.run_raider:
            seed_list = [seed for line in open(args.seed_file) for seed in re.split("\s+", line.rstrip()) if seed] if args.seed_file else [args.seed]
            jobs += [run_raider(seed = convert_seed(seed), seed_num = i, f = args.f, m = args.raider_min, input_file = file, 
                                raider_dir = args.results_dir + "/" + args.raider_dir, mem = args.mem, max_nodes = args.max_nodes) for i,seed in enumerate(seed_list)
                     for file in file_list]

        if args.run_araider:
            seed_list = [seed for line in open(args.seed_file) for seed in re.split("\s+", line.rstrip()) if seed] if args.seed_file else [args.seed]
            jobs += [run_araider(seed = convert_seed(seed), seed_num = i, f = args.f, m = args.raider_min, input_file = file, 
                                araider_dir = args.results_dir + "/" + args.araider_dir, max_nodes = args.max_nodes) for i,seed in enumerate(seed_list)
                     for file in file_list]
        
        if args.run_raider2:
            #raider2_ages = [0,1,2]
            
            seed_list = [seed for line in open(args.seed_file) for seed in re.split("\s+", line.rstrip()) if seed] if args.seed_file else [args.seed]
            
            
            jobs += [run_raider2(seed = convert_seed(seed), seed_num = i, f = args.f, m = args.raider_min, input_file = file,
                                 raider2_dir = args.results_dir + "/" + args.raider2_dir, family_array = args.family_array, excise = args.excising, 
                                 overlaps = args.overlaps, tieup = args.tieup, prosplit=args.prosplit, prevfam=args.prevfam,
                                 age=args.age, age_only=False, max_nodes=args.max_nodes, mem=args.mem) for i,seed in enumerate(seed_list) for file in file_list]
            #if args.all_ages:
            #    if not args.multi_seed:
            #        jobs += [run_raider2(seed = convert_seed(seed), seed_num = i, f = args.f, m = args.raider_min, input_file = file, 
            #                    raider2_dir = args.results_dir + "/" + args.raider2_dir + "." + str(curr_age), age=curr_age) for i,seed in enumerate(seed_list) for curr_age in raider2_ages
            #                    for file in file_list]
            #    else:
            #        jobs += [run_raider2(seed = seed_list, seed_num = "all", f = args.f, m = args.raider_min, input_file = file,
            #                    raider2_dir = args.results_dir + "/" + args.raider2_dir + "." + str(curr_age), age=curr_age) for curr_age in raider2_ages
            #                    for file in file_list]
            #else:
            #    if not args.multi_seed:
            #        jobs += [run_raider2(seed = convert_seed(seed), seed_num = i, f = args.f, m = args.raider_min, input_file = file, 
            #                    raider2_dir = args.results_dir + "/" + args.raider2_dir + "." + str(args.age), age=args.age) for i,seed in enumerate(seed_list)
            #                    for file in file_list]
            #    else:
            #        jobs += [run_raider2(seed = seed_list, seed_num = "all", f = args.f, m = args.raider_min, input_file = file, 
            #                    raider2_dir = args.results_dir + "/" + args.raider2_dir + "." + str(args.age), age=args.age) for file in file_list]

        if args.run_repscout:
            if args.rs_filters == 3:
                jobs += [run_scout(input_file = file, output_dir = args.results_dir + '/' + args.rptscout_dir, min_freq = args.rs_min_freq, length = len(args.seed) if args.seed else args.repscout_min,
                                        use_first_filter = False, use_second_filter = False, threshold = args.f, max_nodes = args.max_nodes, mem = args.mem) for file in file_list]
                jobs += [run_scout(input_file = file, output_dir = args.results_dir + '/' + args.rptscout_dir + "1", min_freq = args.rs_min_freq, length = len(args.seed) if args.seed else args.repscout_min,
                                        use_first_filter = True, use_second_filter = False, threshold = args.f, max_nodes = args.max_nodes, mem = args.mem) for file in file_list]
                jobs += [run_scout(input_file = file, output_dir = args.results_dir + '/' + args.rptscout_dir + "12", min_freq = args.rs_min_freq, length = len(args.seed) if args.seed else args.repscout_min,
                                        use_first_filter = True, use_second_filter = True, threshold = args.f, max_nodes = args.max_nodes, mem=args.mem) for file in file_list]
            else:
                use_first_filter = (args.rs_filters >= 1)
                use_second_filter = (args.rs_filters >= 2)
                dir_addon = ".F1F2" if use_second_filter else ".F1" if use_first_filter else ""
                jobs += [run_scout(input_file = file, output_dir = args.results_dir + '/' + args.rptscout_dir + dir_addon, min_freq = args.rs_min_freq, length = len(args.seed) if args.seed else args.repscout_min, 
                                   use_first_filter = use_first_filter, use_second_filter = use_second_filter, threshold = args.f, max_nodes = args.max_nodes, mem = args.mem) for file in file_list]

        if args.run_bigfoot:
            bigfoot_dir = args.results_dir + "/" + args.bigfoot_dir    # Name of the directory all bigfoot files will go into
            if not os.path.exists(bigfoot_dir):
                os.makedirs(bigfoot_dir)
            jobs += [run_bigfoot(input_file = file, bigfoot_dir = bigfoot_dir, L = args.bigfoot_L, C = args.bigfoot_min, I = args.bigfoot_I, T = args.bigfoot_T) for file in file_list]

        if args.run_piler:
            piler_dir = args.results_dir + "/" + args.piler_dir    # Name of the directory all piler files will go into
            if not os.path.exists(piler_dir):
                os.makedirs(piler_dir)
            jobs +=[run_piler(input_file = file, piler_dir = piler_dir, max_nodes = args.max_nodes) for file in file_list]

            
        ############## Third: Launch repeatmasker jobs
        job_set = {j for j in jobs}
        RM_jobs, PRA_jobs = run_timed_tool_jobs(jobs, args.repmask, args.pa, args.pra, BLAST_DATABASE)

    else:
        ############# Didn't finish processing all of the tool jobs
        jobs = []
        RM_jobs = set()
        PRA_jobs = set()
        if next_step == tjobs_start: 
            jobs, next_step = recover_tool_jobs()
        if next_step == prajobs_start:
            PRA_jobs, next_step = recover_pra_jobs()
        if next_step == rmjobs_start:
            RM_jobs, next_step = recover_rm_jobs()
        flush_files()
        RM_jobs, PRA_jobs = run_timed_tool_jobs(jobs, args.repmask, args.pa, args.pra, BLAST_DATABASE)

    if not continue_prev or next_step == '':
        ########## Need to run analysis jobs

        job_dic, stats_jobs, PRA_jobs = run_timed_analysis_jobs(args.repmask, args.pra, RM_jobs, PRA_jobs, args.results_dir)
    else:
        PRA_jobs = set()
        RM_jobs = set()
        stats_jobs = set()
        job_dic = None
        if next_step == prajobs_start:
            PRA_jobs, next_step = recover_pra_jobs()
        if next_step == rmjobs_start:
            RM_jobs, next_step = recover_rm_jobs()
        if next_step == stats_start:
            stats_jobs, next_step = recover_stats_jobs()
        if next_step == jobdic_start:
            old_job_dic, next_step = recover_job_dic()
            if not RM_jobs:
                job_dic = old_job_dic
            else:
                job_dic, stats_jobs, PRA_jobs = run_timed_analysis_jobs(args.repmask, args.pra, RM_jobs, PRA_jobs, args.results_dir, stats_jobs, old_job_dic)
        else:
            flush_files()
            job_dic, stats_jobs, PRA_jobs = run_timed_analysis_jobs(args.repmask, args.pra, RM_jobs, PRA_jobs, args.results_dir, stats_jobs)

    
    job_dic['raider'].sort(key = lambda x: x.seed_num)
    job_dic['raider2'].sort(key = lambda x: x.seed_num)
    job_dic['araider'].sort(key = lambda x: x.seed_num)
    
    # Print output files log
    with open(args.results_dir + "/file_log.txt", "w") as fp:
        for i in range(len(file_list)):
            fp.write("%d simulation_file %s\n" % (i, file_list[i]))
        for k in test_tools:
            fp.write(k + "\n")
            for j in job_dic[k]:
                if args.repmask:
                    fp.write(j.rm_output + "\n")
    
    ######
    # Create copy of seed file (if RAIDER is being used)
    if job_dic['raider'] or job_dic['araider'] or job_dic['raider2']:
        with open(args.results_dir + "/seed_file.txt", "w") as fp:
            fp.write("\n".join(["{index:<5}{seed}".format(index=i,seed=s) for i,s in enumerate(seed_list)]) + "\n")
            
    ### KARRO
    # Finally: we should not terminate until all the pra jobs are done.  (If pra is off, this list will be empty.)
    #for p in PRA_jobs:
    #    p.timed_wait()        # KARRO: Is this the correct method to use to ensure resubmission of needed
    ### KARRO END

    regex = re.compile("(?<=\# Average consensus coverage: )\d+.\d+")
    
    ######
    # Calculate statistics (not bothering with parallelization yet)

    #print_str: tool       seed      tp/fp/fn/tn           tpr - fdr              ToolCpuTime - RVMem  Con/QuCoverage
    print_str = "{:<12}" + "{:<5}" + "".join("{:<14}"*4) + "".join("{:<14}"*6) + "".join("{:<14}"*8) + "{:<14}"*2 + "\n"
    with open(args.results_dir + "/" + args.stats_file, "w") as fp:
        fp.write(print_str.format("#tool", "seed", "tp", "fp", "fn", "tn", "tpr", "tnr", "ppv", "npv", "fpr", "fdr","ToolCpuTime", "ToolWallTime", "ToolMem", "ToolVMem", "RMCpuTime", "RMWallTime", "RMMem", "RMVMem", "ConCoverage", "QuCoverage"))
        
        for key in test_tools:
            for p in job_dic[key]:
                Counts = [0,0,0,0]
                Stats = [0,0,0,0,0,0]
                RMResources = [0,0,0,0]
                Coverage = 0
                CoverageResources = [0,0,0,0]
                try:
                    if args.repmask:
                        #progress_fp.write(print_time() + "\n")                        
                        progress_fp.write("Calling: perform_stats.py(%s, %s, None)" % (p.seq_file + ".out", p.rm_output))
                        try:
                            Counts, Stats, Sets = perform_stats.perform_stats(p.seq_file + ".out", p.rm_output, None)
                            Stats = [round(x,5) for x in Stats]
                            RMResources = list(p.getResources(cleanup=False))
                            if args.hooke_jeeves:
                                print(Counts[1]+Counts[2])
                        except Exception as E:
                            #progress_fp.write(print_time() + "\n")
                            progress_fp.write("performance Exception: " + str(E) + "\n");
                            fp.write("\t".join([str(key), str(p.seed_num) if hasattr(p, "seed_num") else "NA", "INCOMPLETE\t"]))
                            fp.write("Resources: \t" + "\t".join(p.getResources(cleanup=False)) + "\n");
                            continue
                        
                        if p.pra_job:
                            try:
                                progress_fp.write("parse_pra_outpt: %s %s\n" % (p.pra_job.pra_output, args.exclude))
                                consensus_coverage, query_coverage, Used = parse_pra_output(p.pra_job.pra_output, args.exclude)
                                #matches = regex.findall(open(p.pra_job.pra_output, "r").read())
                                #if len(matches) > 0:
                                #    Coverage = matches[0]
                            except Exception as E:
                                #progress_fp.write(print_time() + "\n")
                                progress_fp.write("PRA Parsing Exception: " + str(E) + "\n");
                                progress_fp.flush()
                    else:
                        try:
                            progress_fp.write("parse_pra_outpt: %s %s\n" % (p.pra_output, args.exclude))
                            consensus_coverage, query_coverage, Used = parse_pra_output(p.pra_output, args.exclude)
                            #matches = regex.findall(open(p.pra_output, "r").read())
                            #if len(matches) > 0:
                            #    Coverage = matches[0]
                            #CoverageResources = list(p.getResources(cleanup=False))
                        except Exception as E:
                            #progress_fp.write(print_time() + "\n")
                            progress_fp.write("PRA Parsing Exception: " + str(E) + "\n");
                            progress_fp.flush()

                    fp.write(print_str.format(*([key, p.seed_num] + list(Counts) + list(Stats) + list(p.tool_resources) + list(RMResources) + [consensus_coverage] + [query_coverage])))

                except Exception as E:
                    ##progress_fp.write(print_time() + "\n")
                    progress_fp.write("performance Exception: " + str(E) + "\n");
                    fp.write("\t".join([str(key), str(p.seed_num) if hasattr(p, "seed_num") else "NA", "INCOMPLETE\n"]))

                            

            
