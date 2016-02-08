import numpy as np
import sys
import glob
import re

marker = 1
stat = sys.argv[marker]
marker+=1
result_list = sys.argv[marker:]

def read_seed_file(results_dir):
    file = "{results_dir}/seed.txt".format(results_dir=results_dir)
    return [line.rstrip() for line in open(file)]

def create_seed_dic(results_list):
    return {results_dir:read_seed_file(results_dir) for results_dir in results_list}

def collect_stats(results_dir, seed_list):
    results = []
    chr_set = set()
    for file in glob.glob("{results_dir}/chr*/stats.txt".format(results_dir=results_dir)):
        print(file)
        chr  = re.search("(chr\d+)", file).group(1)
        chr_set.add(chr)
        with open(file) as fp:
            fp.readline()
            for line in fp:
                A = re.split("\s+", line.rstrip())
                if len(A) != 3 and A[1] != 'NA':
                    # tool, seed_index, seed_list, "tpr, tnr, ppv, npv, fpr, fdr, cpu, mem, rm_cpu, rm_mem
                    results.append({'chr':chr, 'tool':A[0], 'seed_index':int(A[1]), 'seed':seed_list[int(A[1])], 
                                    'tpr':float(A[6]), 'tnr':float(A[7]), 'ppv':float(A[8]), 'npv':float(A[9]), 
                                    'fpr':float(A[10]), 'fdr':float(A[11]), 
                                    'cpu':int(A[12]), 'mem':int(A[14]), 
                                    'rm_cpu':int(A[16]), 'rm_mem':int(A[18])})

    return results

def rank_chr(results, chr_set = None, stat = 'tpr'):
    if not chr_set:
        chr_set = {r['chr'] for r in results}
    for chr in chr_set:
        L = sorted([r for r in results if r['chr'] == chr], key = lambda r: r[stat])
        j = 0
        while L[j][stat] == 0:
            j += 1
        for i in range(len(L)):
            L[i]['improvement'] = L[i][stat] / float(L[j][stat]) if i >= j else 0
            L[i]['improvement2'] = L[i][stat] / float(L[i-1][stat]) if i >= j+1 else 0
        for l in L:
            print("{chr:<10} {seedIndex:5} {improve:<0.4f} {last:<0.4f} {absolute:0.4f} {seed}".format(chr=l['chr'], seedIndex=l['seed_index'], improve=l['improvement'], last=l['improvement2'], absolute=l[stat], seed=l['seed']))
        print("-------")


D = create_seed_dic(result_list)
R = [collect_stats(dir, D[dir]) for dir in result_list]
rank_chr(R[0],stat=stat)
