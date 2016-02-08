import random
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

def weight(key):
    return sum([1 for s in key if s == "1"])

def check_weights(file):
    i = 0
    for line in open(file):
        seed = expand_seed(line.rstrip())
        w = weight(seed)
        print(len(seed), w, float(w)/len(seed))
        i += 1
    return i

def gen_seed(L):
    s = ""
    t = 1
    for i in L:
        s += str(t) + "^{" + str(i) + "}"
        t = 1 - t

    return s + "\n"

def random_seed(w,p):
    """generate a random seed of weight w, with each character
    being a 0 with probability p."""
    c = 1
    s = "1"
    while c < w-1:
        if random.uniform(0,1) < p:
            s = s + "0"
        else:
            s = s + "1"
            c = c + 1
    return s + "1"


def random_seed2(w,l):
    """Generate a random seed of length l and weight w"""
    L = ["1"]*(w-2) + ["0"]*(l-w)
    random.shuffle(L)
    return "1" + "".join(L) + "1"

def random_seed3(w,l):
    """Generate a random palindomic seed o length l and weight w"""
    L = ["1"]*(w//2 - 1) + ["0"]*((l-w)//2)
    random.shuffle(L)
    return "1" + "".join(L + L[::-1]) + "1"

def generate_list(w, p_lower, p_upper, p_step, n):
    debug = 0
    M = {}

    p = p_lower
    while p <= p_upper:
        l = int(round(float(w)/p))
        for k in range(n):
            num_tries = 0
            while True:
                s = random_seed2(w,l)
                if s not in M:
                    M[s] = True
                    break
                num_tries += 1
                if (num_tries >= 100):
                    break
            debug += 1
        p += p_step

    return sorted(M.keys(), key=len)

def create1(): 
    fp = open("seeds40.1.txt", "w")
    L = generate_list(40, 0.40, 0.9751, 0.025, 20)
    fp.write("\n".join([compact_seed(s) for s in L]))
    fp.close()

def create2(): 
    fp = open("seeds40.2.txt", "w")
    L = generate_list(40, 0.7, 0.9, 0.025, 30)
    fp.write("\n".join([compact_seed(s) for s in L]))
    fp.close()

def create3(): 
    fp = open("seeds40.3.txt", "w")
    L = generate_list(40, 0.65, 0.851, 0.025, 100)
    fp.write("\n".join([compact_seed(s) for s in L]))
    fp.close()
    
def main():
    create3()

if __name__ == "__main__":
    main()
