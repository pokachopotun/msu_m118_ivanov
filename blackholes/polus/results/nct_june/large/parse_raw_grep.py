#!/usr/bin/python3
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

import sys
import os

def PlotGraphFromScale(data, graph, figid, measure):
    fig = plt.figure(figid)
    plt.rcParams.update({'font.size': 14})
    fig.suptitle("{}".format(graph), fontsize=16)
    fig.tight_layout()
    top = fig.add_subplot(1, 1, 1)
    top.set_xlabel('Graph scale')
    top.set_ylabel(measure)
    for algo in ["topsort", "chinese", "topsort_over_chinese"]:
        if not algo in data:
            continue
        dg = data[algo][graph]
        if algo == "topsort_over_chinese":
            algo = "TopOver"
        if algo == "topsort":
            algo = "TopSort_BH"
        if algo == "chinese":
            algo = "iBlackholeDC"
        for threads in [1, 8]:
            for cond in ["cond", "nocond"]:
                if algo  == "iBlackholeDC":
                    skips = ["NoSkip",]
                else:
                    skips = ["SkipFast", "SkipPrecise"]
                for skip in skips:
                    if not threads in dg:
                        continue
                    dgc = dg[threads][cond][skip]
                    scales = dgc.keys()
                    times = dgc.values()
                    if len(scales) > 0 and len(scales) == len(times):
                        scales, times = zip(*sorted(zip(scales, times)))
                        top.plot(scales, times, label = "{}.{}.{}".format(algo, cond, skip), linewidth = 3)
                        top.grid(True)
                        fontP = FontProperties()
                        fontP.set_size(10)
                        top.legend(prop = fontP)
    top.set_yscale('symlog')

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("use python parse_raw_grep.py raw.txt")
        exit(0)

    inputFileName = sys.argv[1]
    measure = sys.argv[1]

    if inputFileName == "time.txt":
        measure = "Time (s)"
    if inputFileName == "bh_count.txt":
        measure = "Found blackholes"
    if inputFileName == "bh_filtered.txt":
        measure = "Filtered blackholes"

    # data -> graph -> threads -> scale -> time
    data = dict()

    with open(inputFileName, 'r') as file:
        for line in file:
            descr, seconds = line.strip().split()[:2]
            seconds = float(seconds)
            vals =  [x.strip() for x in descr.strip().split('.')]
            cond = "nocond"
            graph, size, threads, useCond, skip, algo = vals

            if int(useCond) == 1:
                cond = "cond"

            if int(skip) == 0:
                skip = "SkipPrecise"
            else:
                skip = "SkipFast"

            if algo == "chinese":
                skip = "NoSkip"

            size = int(size)
            threads = int(threads)

            if not algo in data:
                data[algo] = dict()
                data[algo]["RMAT"] = dict()
                data[algo]["SSCA2"] = dict()
                data[algo]["UR"] = dict()

            dg = data[algo][graph]
            if not threads in dg:
                dg[threads] = dict()
                dg[threads]["cond"] = dict()
                dg[threads]["nocond"] = dict()
                dg[threads]["cond"]["SkipFast"] = dict()
                dg[threads]["cond"]["SkipPrecise"] = dict()
                dg[threads]["cond"]["NoSkip"] = dict()
                dg[threads]["nocond"]["SkipFast"] = dict()
                dg[threads]["nocond"]["SkipPrecise"] = dict()
                dg[threads]["nocond"]["NoSkip"] = dict()

            dgc = dg[threads][cond][skip]

            dgc[size] = seconds

    plotConfig = dict()

    #plotConfig["RMAT"] = list([[4,5], [6,7], [8,8], [9,10]])
    #plotConfig["SSCA2"] = list([[4,5],[6,6]])
    #plotConfig["UR"] = list([[4,5], [6,7], [8,10], [11,12]])

    figid = 0
    for graph in ["RMAT", "SSCA2", "UR"]:
        figid += 1
        PlotGraphFromScale(data, graph, figid, measure)

    figid = 0
    for graph in ["RMAT", "SSCA2", "UR"]:
        figid += 1
        plt.figure(figid).savefig(inputFileName + "{}_{}.pdf".format(figid, graph))

    plt.show()
