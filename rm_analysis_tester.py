import rm_analysis
import subprocess
import re

header = "   SW  perc perc perc  query      position in query           matching       repeat              position in  repeat\nscore  div. del. ins.  sequence    begin     end    (left)    repeat         class/family         begin  end (left)   ID\n\n"

line = "  197  13.5  5.4  0.0  chr22     {} {} ({}) +  {}           SINE/MIR                86  124  (144)      1\n"

def generate_test(tool, real, chrLen):
    tool_name = "tool.fa.out"
    real_name = "real.fa.out"
    real_fa   = "real.fa"

    open(real_fa, "w").write(">bogus\n" + "A"*chrLen);
    fp = open(tool_name, "w")
    fp.write(header)
    for start,stop in tool:
        fp.write(line.format(start, stop, chrLen - stop, "XXX"))
    fp.close()
    
    fp = open(real_name, "w")
    fp.write(header)
    for start,stop in real:
        fp.write(line.format(start, stop, chrLen - stop, "ALU"))
    fp.close()

    return tool_name, real_name, real_fa

cmd = "python rm_analysis.py {} {} {}"
parse_re1 = "Overall sensitivity: (\d+)\s+(\d+)\s+(\d+(?:\.\d+)?)"
parse_re2 = "Total uncovered bases\:\s*(\d+)"
parse_re3 = "Total false positives:\s*(\d+)"
parse_re4 = "Total FPs[^\:]*\: (\d+)"

def run_test(tool, real, chrLen, fp_dist, total, tp, fp, fp_d, neg, debug = False):
    """
       Problem parameters:
         Tool: Intervals produced by tool
         Real: Real intervals
         chrLen: Total length
         fp_dist: Limit on distance from a real for the fp_d count
       Correct answers:
         total: Total bases in family
         tp: Total number of true positives
         fp: Total number of false positives
         fp_d: Total number of fp that are within fp_dist of a real interval
         neg: Total number of bases not in family
    """
    tool_name, real_name, real_fa = generate_test(tool, real, chrLen)

    c = cmd.format(real_name,tool_name,fp_dist)
    if debug:
        print("cmd: " + c)
    p = subprocess.Popen(c, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    o, e = p.communicate()

    e = e.decode()
    if len(e) > 0:
        if debug:
            print("ERROR: " + e.rstrip())
        return False
    o = o.decode()
    if debug:
        print("OUTPUT:\n" + o.rstrip())

    r1 = re.search(parse_re1, o)    
    r2 = re.search(parse_re2, o)
    r3 = re.search(parse_re3, o)
    r4 = re.search(parse_re4, o)

    if (debug):
        print("**********")
        print("TP: " + str(int(r1.group(1)) == tp) + " " + str(r1.group(1)) + " " + str(tp))
        print("TOTAL: " + str(int(r1.group(2)) == total) + " " + str(r1.group(2)) + " " + str(total))
        print("NEG: " + str(int(r2.group(1)) == neg) + " " + r2.group(1) + " " + str(neg))
        print("FP: " + str(int(r3.group(1)) == fp) + " " + r3.group(1) + " "  + str(fp))
        print("FP_D: " + str(int(r4.group(1)) == fp_d) + " " + r4.group(1) + " " + str(fp_d))

    return int(r1.group(1)) == tp and int(r1.group(2)) == total and int(r2.group(1)) == neg and int(r3.group(1)) == fp and int(r4.group(1)) == fp_d


def test1(debug = False):
    """Single intervals, tool on left"""
    tool = [[11,20]]   # Tool reported intevals
    real = [[31,40]]   # Real intervals
    chrLen = 50        # Length of chromosome
    fp_dist = 15       # Designated fp_d limit
    total = 10         # Total number of  TE bases (sum of length of real intervals)
    tp = 0             # Actual true positives
    fp = 10;           # Actual false positives
    fp_d = 5           # Actual FPs within fp_dist bases of an actual positive
    neg = 40           # Actual number of uncovered bases
    result = run_test(tool, real, chrLen, fp_dist, total, tp, fp, fp_d, neg, debug)
    return result

def test16(debug = False):
    """Single intervals, tool on left"""
    tool = [[11,20]]
    real = [[31,40]]
    chrLen = 50
    fp_dist = 25
    total = 10
    tp = 0
    fp = 10;
    fp_d = 10
    neg = 40
    result = run_test(tool, real, chrLen, fp_dist, total, tp, fp, fp_d, neg, debug)
    return result

def test2(debug = False):
    """Single intervals, tool on right"""
    tool = [[31, 40]]
    real = [[11, 20]]
    chrLen = 50
    fp_dist = 15
    total = 10
    tp = 0
    fp = 10
    fp_d = 5
    neg = 40
    result = run_test(tool, real, chrLen, fp_dist, total, tp, fp, fp_d, neg, debug)
    return result

def test17(debug = False):
    """Single intervals, tool on right"""
    tool = [[31, 40]]
    real = [[11, 20]]
    chrLen = 50
    fp_dist = 25
    total = 10
    tp = 0
    fp = 10
    fp_d = 10
    neg = 40
    result = run_test(tool, real, chrLen, fp_dist, total, tp, fp, fp_d, neg, debug)
    return result

def test3(debug = False):
    """Single intervals, tool overlap on left"""
    tool = [[11, 30]]
    real = [[21, 40]]
    chrLen = 100
    fp_dist = 5
    total = 20
    tp = 10
    fp = 10
    fp_d = 5
    neg = 80
    result = run_test(tool, real, chrLen, fp_dist, total, tp, fp, fp_d, neg, debug)
    return result

def test4(debug = False):
    """Single intervals, tool overlap on right"""
    tool = [[21, 40]]
    real = [[11, 30]]
    chrLen = 50
    fp_dist = 5
    total = 20
    tp = 10
    fp = 10
    fp_d = 5
    neg = 30
    result = run_test(tool, real, chrLen, fp_dist, total, tp, fp, fp_d, neg, debug)
    return result

def test5(debug = False):
    """Single intervals, tool contains real"""
    tool = [[21, 50]]
    real = [[31, 40]]
    chrLen = 70
    total = 10
    fp_dist = 5
    tp = 10
    fp = 20
    fp_d = 10
    neg = 60
    result = run_test(tool, real, chrLen, fp_dist, total, tp, fp, fp_d, neg, debug)
    return result

def test6(debug = False):
    """Single intervals, tool contained in real"""
    tool = [[31, 40]]
    real = [[21, 50]]
    chrLen = 70
    total = 30
    fp_dist = 5
    tp = 10
    fp = 0
    fp_d = 0
    neg = 40
    result = run_test(tool, real, chrLen, fp_dist, total, tp, fp, fp_d, neg, debug)
    return result

def test7(debug = False):
    """Intervals the same"""
    tool = [[31, 40]]
    real = [[31, 40]]
    chrLen = 50
    fp_dist = 5
    total = 10
    tp = 10
    fp = 0
    fp_d = 0
    neg = 40
    result = run_test(tool, real, chrLen, fp_dist, total, tp, fp, fp_d, neg, debug)
    return result

def test8(debug = False):
    """Single tool containing two real"""
    tool = [[21,70]]
    real = [[31, 40], [51,60]]
    chrLen = 70
    fp_dist = 3
    total = 20
    tp = 20
    fp = 30
    fp_d = 12
    neg = 50
    result = run_test(tool, real, chrLen, fp_dist, total, tp, fp, fp_d, neg, debug)
    return result

def test9(debug = False):
    """Single tool overlapping two real"""
    tool = [[41,70]]
    real = [[31, 50], [61,80]]
    chrLen = 100
    fp_dist = 3
    total = 40
    tp = 20
    fp = 10
    fp_d = 6
    neg = 60
    result = run_test(tool, real, chrLen, fp_dist, total, tp, fp, fp_d, neg, debug)
    return result

def test10(debug = False):
    """Single real containing two tool"""
    real = [[21,70]]
    tool = [[31, 40], [51,60]]
    chrLen = 100
    fp_dist = 5
    total = 50
    tp = 20
    fp = 0
    fp_d = 0
    neg = 50
    result = run_test(tool, real, chrLen, fp_dist, total, tp, fp, fp_d, neg, debug)
    return result

def test11(debug = False):
    """Single real overlapping two tools"""
    real = [[41,70]]
    tool = [[31, 50], [61,80]]
    chrLen = 100
    fp_dist = 5
    total = 30
    tp = 20
    fp = 20
    fp_d = 10
    neg = 70
    result = run_test(tool, real, chrLen, fp_dist, total, tp, fp, fp_d, neg, debug)
    return result

def test12(debug = False):
    """Single real overlapping, two overlapping tools"""
    real = [[41,60]]
    tool = [[31, 55], [45,80]]
    chrLen = 100
    fp_dist = 5
    total = 20
    tp = 20
    fp = 30
    fp_d = 10
    neg = 80
    result = run_test(tool, real, chrLen, fp_dist, total, tp, fp, fp_d, neg, debug)
    return result

def test13(debug = False):
    """Left-over reals"""
    tool = [[11,20],[51,60]]
    real = [[21,30],[71,80],[91,100]]
    chrLen = 200
    total = 30
    fp_dist = 5
    tp = 0
    fp = 20
    fp_d = 5
    neg = 170
    result = run_test(tool, real, chrLen, fp_dist, total, tp, fp, fp_d, neg, debug)
    return result

def test14(debug = False):
    """Three-way-overlap"""
    tool = [[11,30], [21,40], [31,50]]
    real = [[11,50]]
    chrLen = 100
    total = 40
    fp_dist = 5
    tp = 40
    fp = 0
    fp_d = 0
    neg = 60
    result = run_test(tool, real, chrLen, fp_dist, total, tp, fp, fp_d, neg, debug)
    return result

def test15(debug=False):
    """Overlap between multiple real"""
    real = [[21,40], [61,80]]
    tool = [[11,50], [31,70], [51, 90]]
    chrLen = 100
    fp_dist = 5
    total = 40
    tp = 40
    fp = 40
    fp_d = 20
    neg = 60
    result = run_test(tool, real, chrLen, fp_dist, total, tp, fp, fp_d, neg, debug)
    return result

def test18(debug=False):
    """Testing overlap in adjacent regions when reals are adjacent"""
    real = [[21,30],[31,40]]
    tool = [[11,50]]
    chrLen = 60
    fp_dist = 5
    total = 20
    tp = 20
    fp = 20
    fp_d = 10
    neg = 40
    result = run_test(tool, real, chrLen, fp_dist, total, tp, fp, fp_d, neg, debug)
    return result    

def test19(debug=False):
    """Testing overlap in adjacent regions"""
    real = [[51,100],[151,200]]
    tool = [[101,150]]
    chrLen = 250
    fp_dist = 100
    total = 100
    tp = 0
    fp = 50
    fp_d = 50
    neg = 150
    result = run_test(tool, real, chrLen, fp_dist, total, tp, fp, fp_d, neg, debug)
    return result

def test20(debug=False):
    """Testing overlap in adjacent regions"""
    real = [[51,100],[131,150]]
    tool = [[81,120]]
    chrLen = 200
    fp_dist = 200
    total = 70
    tp = 20
    fp = 20
    fp_d = 20
    neg = 130
    result = run_test(tool, real, chrLen, fp_dist, total, tp, fp, fp_d, neg, debug)
    return result

def test21(debug=False):
    """Testing overlap in adjacent regions"""
    real = [[51,100],[131,150]]
    tool = [[101,140]]
    chrLen = 200
    fp_dist = 200
    total = 70
    tp = 10
    fp = 30
    fp_d = 30
    neg = 130
    result = run_test(tool, real, chrLen, fp_dist, total, tp, fp, fp_d, neg, debug)
    return result


def test22(debug=False):
    """First tool encompases second tool, both within real"""
    real = [[101,200]]
    tool = [[121,180], [141,150]]
    chrLen = 300
    fp_dist = 10
    total = 100
    tp = 60
    fp = 0
    fp_d = 0
    neg = 200
    result = run_test(tool, real, chrLen, fp_dist, total, tp, fp, fp_d, neg, debug)
    return result    

def test23(debug=False):
    """First tool encompases second tool, overlapping real"""
    real = [[101,200]]
    tool = [[121,220], [141,150]]
    chrLen = 300
    fp_dist = 10
    total = 100
    tp = 80
    fp = 20
    fp_d = 10
    neg = 200
    result = run_test(tool, real, chrLen, fp_dist, total, tp, fp, fp_d, neg, debug)
    return result    

def test24(debug=False):
    """First tool encompases second tool, overlapping real"""
    real = [[101,200]]
    tool = [[81,180], [141,150]]
    chrLen = 300
    fp_dist = 10
    total = 100
    tp = 80
    fp = 20
    fp_d = 10
    neg = 200
    result = run_test(tool, real, chrLen, fp_dist, total, tp, fp, fp_d, neg, debug)
    return result

def test25(debug=False):
    """First tool encompases second tool, overlapping real"""
    real = [[51,250]]
    tool = [[101,200], [121,160], [141,180], [161,220]]
    chrLen = 300
    fp_dist = 10
    total = 200
    tp = 120
    fp = 0
    fp_d = 0
    neg = 100
    result = run_test(tool, real, chrLen, fp_dist, total, tp, fp, fp_d, neg, debug)
    return result    

def test26(debug=False):
    """Check fp_d under second case"""
    real = [[11,20]]
    tool = [[21,30]]
    chrLen = 40
    fp_dist = 5
    total = 10
    tp = 0
    fp = 10
    fp_d = 5
    neg = 30
    result = run_test(tool, real, chrLen, fp_dist, total, tp, fp, fp_d, neg, debug)
    return result    


def run_all(debug = False):
    try:
        i = 1
        while True:
            r = eval("test%d" % (i))(debug)
            print("Test " + str(i) + ": " + str(r))
            i += 1
    except Exception as E:
        if debug:
            raise E
        pass


if __name__ == "__main__":
    run_all()
