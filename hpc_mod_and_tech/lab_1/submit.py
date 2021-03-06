import os
import sys

class Submission:
    def __init__(self, typeStr1, typeStr2):
        self.typeStr1 = typeStr1
        self.typeStr2 = typeStr2
        self.command = r"mpisubmit.pl = -g --stdout=" + self.typeStr1 + "." 
        pass

## mpisubmit.pl -g --stdout=$(1).RMAT.out --stderr=$(1).RMAT.err ./color -- --graph-scale=$(1) --graph-edgefactor=32 --test-run=0 --graph-type=rmat
if __name__ == "__main__":
    submit = len(sys.argv) > 2 and "submit" == sys.argv[2]
    sub = list()
    with open("config.txt") as file:
        n, a, b, t = [int(x) for x in file.readline().strip().split()]
        for i in range(t):
            typeStr1, typeStr2 = file.readline().strip().split()
            sub.append(Submission(typeStr1, typeStr2))

    for s in sub:
        for scale in range(a, b + 1):
            for rep in range(n):
                
                scale_str = "{:02d}".format(scale) 
                rep_str = "{:02d}".format(rep)
                id_str = s.typeStr1 + "." + scale_str + "." + rep_str
                outfile = "outputs/" + id_str
                errfile = "errors/" + id_str
                graphfile = "inputs/" + s.typeStr1 + "." + scale_str

                cmd0 = "bsub -gpu \"num=1:mode=exclusive_process\" -oo " + outfile + " -eo " + errfile + " -q normal"
                if s.typeStr1 == "SSCA2":
                    cmd1 = "./color --graph-type=market" + " --graph-file=" + graphfile + " --vertex-start-from-zero=0 --undirected --validation=none"
                else:
                    cmd1 = "./color --graph-edgefactor=32 --graph-type=" + s.typeStr2 + " --graph-scale=" + str(scale) + " --undirected --validation=none"
                cmd = cmd0 + " " + cmd1
                print(cmd)
                if submit:
                    os.system(cmd)
         
