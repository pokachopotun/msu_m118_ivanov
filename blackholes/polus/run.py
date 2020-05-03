import os
import subprocess

class Submission:
    def __init__(self, typeStr1, typeStr2):
        self.typeStr1 = typeStr1
        self.typeStr2 = typeStr2
        self.command = r"mpisubmit.pl = -g --stdout=" + self.typeStr1 + "."
        pass

if __name__ == "__main__":
    sub = list()
    with open("config.txt") as file:
        n, a, b, t = [int(x) for x in file.readline().strip().split()]
        for i in range(t):
            typeStr1, typeStr2 = file.readline().strip().split()
            sub.append(Submission(typeStr1, typeStr2))

    for scale in range(a, b + 1, 2):
        print("echo scale = " + str(scale))
        for s in sub:
            scale_str = "{:02d}".format(scale)
            id_str = s.typeStr1 + "." + scale_str
            outfile = "outputs/" + id_str
            errfile = "errors/" + id_str
            inputfile = "graphs/" + id_str + ".bin"

            chinese = "../chinese_cli " + inputfile + " 0 0 0"
            topsort = "../topsort_cli " + inputfile + " 0 0 0"

            cmd1 = "bsub -n 1 -W 30 -o " + outfile + ".chinese " + chinese
            cmd2 = "bsub -n 1 -W 30 -o " + outfile + ".topsort " + topsort

            print(cmd1)
            print(cmd2)

#define graph
#	./generator -s $(1) -directed -weighted -file inputs/$(1).$(2).graph.bin -type $(2)
#	./graph_to_file inputs/$(1).$(2).graph.bin inputs/$(1).$(2).graph
#endef
