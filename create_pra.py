import sys
import redhawk
import glob
import os
import re
import redhawk

reg = re.compile("^(.*\/)([^\/]*)\.s(\d+)\.f(\d+).blast.6.txt(.bz2)?$")

def main(DIR):
    for file in glob.glob("{DIR}/*/*/*.blast.6.txt*".format(DIR=DIR)):
        print(file)
        r = reg.search(file)
        if not r:
            sys.stderr.write("Bad file: " + file + "\n")
            continue

        D, org, seed_num, f = r.group(1,2,3,4)

        pra_output =  D + "/" + "{org}.s{seed}.f{f}.pra.txt".format(org=org, seed=seed_num, f=f)

        if os.path.exists(pra_output):
            continue

        if file.endswith(".bz2"):
            cmd = "bzcat {blast_output} | ./pra_analysis2 {output}".format(blast_output=file, output=pra_output)
        else:
            cmd = "cat {blast_output} | ./pra_analysis2 {output}".format(blast_output=file, output=pra_output)

        p = redhawk.pbsJobHandler(batch_file = D + "pra_job", executable = cmd, job_name = "pra.{org}.s{seed}.f{f}".format(org=org, seed=seed_num, f=f),
                                  walltime = "3:00:00", output_location = D)

        p.submit(preserve=True)
        
if __name__ == "__main__":
    main(sys.argv[1])
