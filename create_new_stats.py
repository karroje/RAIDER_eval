import sys
import re

def seed2num(s):
    t = 0
    for i in s:
        t = 2*t + int(i)
    return t

def num2seed(v):
    s = ""
    while v > 0:
        s = str(v %2) + s
        v = v // 2
    return s

if __name__ == "__main__":
    fp = open(sys.argv[1])
    wp = open(sys.argv[2], "w")

    for line in fp:
        index, seed = re.split("\s+", line.rstrip())
        wp.write("{index:<15}{new_index:<30}{seed}\n".format(index=index, new_index=seed2num(seed), seed=seed))
        
