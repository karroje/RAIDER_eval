import subprocess

chr_list = [19,18]
cmd = "python3.3 RAIDER_eval.py --nuke -r TEST1/chr{chr}.sim1.all -R --sf seed.TEST1.txt --mem chrom_sim --st 0 data/chr{chr}.fa > nohup.{chr}.out"

P = []
for i in chr_list:
    c = cmd.format(chr=i)
    print(c)
    P.append(subprocess.Popen(c, shell = True))

for p in P:
    p.wait()

