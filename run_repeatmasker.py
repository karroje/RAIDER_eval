from redhawk import *
import sys
import subprocess

if relaunch(sys.argv, walltime = "400:00:00"):
    exit(0)

cmd = "RepeatMasker chr%d.fa -species mammals -a -ps 8"
for i in range(2,3):
    o = pbsJobHandler("chr%d.batch" % (i), cmd % (i), nodes = 2, ppn = 4, walltime = "400:00:00", RHmodules = ["RepeatMasker"])
    o.submit()
