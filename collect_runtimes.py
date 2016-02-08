import glob
import re

ToolsDic = {'RS':'RepeatScout', 'phRA':'phRAIDER', 'RA':'RAIDER', 'prephRA':'pre-phRAIDER', 'na':'naieve'}
SizeDic = {org + "." + chr : size for line in open("sizes.txt") for org,chr,size in [re.split("\s+", line.rstrip())]}

def compact_seed(s):
    while True:
        r = re.search("(111+)", s)
        if r:
            w = r.group(1)
            s = re.sub("1"*len(w), "1^{" + str(len(w)) + "}", s)
        else:
            break

    while True:
        r = re.search("(000+)", s)
        if r:
            w = r.group(1)
            s = re.sub("0"*len(w), "0^{" + str(len(w)) + "}", s)
        else:
            break
    return s

def expand_seed(seed):
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


def create_seed_dic():
   D = {}
   for line in open("seed_file.txt"):
      seed_num, seed_comp = re.split("\s+", line.strip())
      seed = expand_seed(seed_comp)
      w = sum([1 for i in seed if i=="1"])
      l = len(seed)
      D[seed_num] = (seed,w,l,round(w/l,3))
   return D;

SeedMap = create_seed_dic()

def readFile(file):
   for line in open(file):
      r = re.search("User time \(seconds\):\s+(\d+\.\d+)", line)
      if r:
         time= float(r.group(1))
         continue
      r = re.search("Maximum resident set size \(kbytes\):\s+(\d+)", line)
      if r:
         mem = float(r.group(1))/4
         return (time,mem)
   return ('NA', 'NA')

def readFileOakley(file):
   for line in open(file):
      r = re.search("cput=(\d\d?):(\d\d?):(\d\d?)", line)
      if r:
         time = 60*60*int(r.group(1))  + 60*int(r.group(2)) + int(r.group(3))
         continue
      r = re.search("^mem=(\d+(?:\.\d+)?)\s+GB", line)
      if r:
         mem = float(r.group(1))
         return (time,mem)
   return ('NA', 'NA')


R = {}
for file in glob.glob("*/*/*.stderr"):
   r = re.search("\w+/(.*)\.(?:CHR)?(\w+)\.S(\d+).F(\d+)/(\w+)", file)
   if not r:
      continue
   org, chr, seed, f, tool = r.group(1,2,3,4,5);
   if org.startswith("ARAB"):
      org = "arab"
   else:
      org = org.lower()
   chr = "chr" + chr

   if org + "." + chr not in SizeDic:
      continue
   key = (tool,org,chr,seed,f)

   time,mem = readFile(file)

   if time == 'NA' or mem == 'NA':
      continue

   if key not in R:
      R[key] = (time,mem)
   else:
      t,m = R[key]
      R[key] = (t+time, max(m,mem))


print("\t".join(["tool", "label", "seed", "f", "size", "time", "mem", "len", "weight", "d", "seed"]))
for key in sorted(R.keys()):
   seed = key[3]
   print("\t".join([ToolsDic[key[0]] if key[0] in ToolsDic else key[0], key[1] + "." + key[2], str(key[3]), str(key[4]), SizeDic[key[1] + "." + key[2]], \
                    str(R[key][0]), str(R[key][1]), str(SeedMap[seed][1]), str(SeedMap[seed][2]), str(SeedMap[seed][3]), SeedMap[seed][0]]))
