#!/usr/bin/python3

import sys
import os
import subprocess

class Submission:
    def __init__(self, typeStr1, typeStr2):
        self.typeStr1 = typeStr1
        self.typeStr2 = typeStr2
        self.command = r"mpisubmit.pl = -g --stdout=" + self.typeStr1 + "."
        pass

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("use python3 run.py config.txt")
        exit(0)

    configFileName = sys.argv[1]

    # run params
    binaryInput = 1
    printDebugInfo = 0
    useDivideAndConquer = 1
    useCondensation = 0
    skipModeList =[0,]
    useOpenMP = []
    # run params

    sub = list()
    with open(configFileName, 'r') as file:
        runName = file.readline().strip()
        utils = [x.strip() for x in file.readline().strip().split()]
        n, a, b, step, t = [int(x) for x in file.readline().strip().split()]
        for i in range(t):
            typeStr1, typeStr2 = file.readline().strip().split()
            sub.append(Submission(typeStr1, typeStr2))

        useOpenMP = [int(x) for x in file.readline().strip().split()[1:]]
        useCondensation = [int(x) for x in file.readline().strip().split()[1:]]
        maxBHSize = [int(x) for x in file.readline().strip().split()[1:]]
        skipModeList = [int(x) for x in file.readline().strip().split()[1:]]

    for scale in range(a, b + 1, step):
        print("echo scale = " + str(scale))
        for ompThreads in useOpenMP:
            for useCond in useCondensation:
                for bh_size in maxBHSize:
                    for skipMode in skipModeList:
                        for s in sub:
                            scale_str = "{:02d}".format(scale)
                            id_str = "{}.{}".format(s.typeStr1, scale_str)
                            file_str = id_str + ".{}.{}.{}".format(ompThreads, useCond, skipMode)
                            outfile = "outputs/" + file_str
                            errfile = "errors/" + file_str
                            inputfile = "../../graphs/" + id_str + ".bin"

                            for utility in utils:

                                if utility in ["chinese", "topsort_over_chinese"]:
                                    cmd1 = "../../{}_cli ".format(utility) + "{} {} {} {} {} {} {}".format(inputfile, binaryInput, bh_size, useCond, ompThreads, skipMode, printDebugInfo)
                                if utility == "topsort":
                                    cmd1 = "../../{}_cli ".format(utility) + "{} {} {} {} {}".format(inputfile, binaryInput, useDivideAndConquer, ompThreads, skipMode)

                                cmd = "bsub -n 1 -W 30 -oo {}.{} -eo {}".format(outfile, utility, errfile) + ".{} {}".format(utility, cmd1)

                                print(cmd)

#define graph
#	./generator -s $(1) -directed -weighted -file inputs/$(1).$(2).graph.bin -type $(2)
#	./graph_to_file inputs/$(1).$(2).graph.bin inputs/$(1).$(2).graph
#endef
