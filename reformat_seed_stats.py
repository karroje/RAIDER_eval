import seed_fun
import sys
import re

def main(set, print_header = False):
    fp = open("seeds%d.stats.txt" % set)
    wp = open("seeds%s.stats2.txt" % set, "w")
    header = fp.readline()
    

    H = re.split("\s+", header.strip())
    D = {x:i for i,x in enumerate(H)}

    wp.write("{:<15}".format("seed_set") + "\t" + "\t".join(["{:<15}".format(x) for x in H]) + "\n")
    if print_header:
        wp2.write("{:<15}".format("seed_set") + "\t" + "\t".join(["{:<15}".format(x) for x in H]) + "\n")
    for line in fp:
        A = {H[i]:v for i,v in enumerate(re.split("\s+", line.rstrip()))}
        if A['seed'] != 'NA':
            A['seed'] = seed_fun.compact_seed(A['seed'])
        A['ToolMem'] = 'NA' if A['ToolMem'] == 'NA' else str( int(A['ToolMem'])/1000000.0 ) 
        A['ToolVMem'] = 'NA' if A['ToolVMem'] == 'NA' else str( int(A['ToolVMem'])/1000000.0 ) 
        wp.write("{:<15}".format(set) + "\t" + "\t".join("{:<15}".format(A[h]) for h in H) + "\n")
        if A['tool'] == 'phRAIDER':
            wp2.write("{:<15}".format(set) + "\t" + "\t".join("{:<15}".format(A[h]) for h in H) + "\n")

wp2 = open("all.stats2.txt", "w")
for i in range(12,14):
    main(i, i==1)
