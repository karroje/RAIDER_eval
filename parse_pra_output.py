import re
import sys

r = re.compile("^# Average consensus coverage: (\d+(?:\.\d+)?)\s*$")
r2 = re.compile("^([^\|]*)\|([^\|]*)\|(\S*)\s+(\d+)\s+(\d+)\s+(\d+(?:\.\d+)?)\s*$")
def parse_pra_output(file, CLASS_Exclude = None):
    try:
        fp = open(file)
    except:
        return -1, -1, {}
    line = fp.readline()
    while (line and not r.match(line)):
        line = fp.readline()

    if not line:
        return -1, -1, {}

    consensus_cover = float(r.match(line).group(1))

    found = 0;
    total = 0;

    fp.readline()
    line = fp.readline()

    CLASSES_Used = set()  # For debugging
    while (line[0] != '#'):
        rfamily, rclass, rloc, b1, b2, p = r2.match(line).group(1,2,3,4,5,6)
        if rclass not in CLASS_Exclude:
            b1 = int(b1)
            b2 = int(b2)
            CLASSES_Used.add(rclass)
        
            found += int(b1)
            total += int(b2)
        line = fp.readline()
        if not line:
            return -1, -1, CLASSES_Used

    return consensus_cover, (1.0*found) / total if total > 0 else -1, CLASSES_Used

if __name__ == "__main__":
    Exclude = {line.rstrip() for line in open(sys.argv[2] if len(sys.argv) > 2 else "exclude.txt")}
    c1, c2, S = parse_pra_output(sys.argv[1], Exclude)
    print("Used: " + " ".join([x for x in S]))
    print("Consensus: " + str(c1))
    print("Query: " + str(c2))



