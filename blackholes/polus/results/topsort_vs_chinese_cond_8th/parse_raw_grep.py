import matplotlib.pyplot as plt

import sys
import os

def PlotGraphFromScale(data, graph, figid):
    fig = plt.figure(figid)
    fig.suptitle("{}".format(graph), fontsize=16)
    fig.tight_layout()
    top = fig.add_subplot(1, 1, 1)
    top.set_xlabel('Graph scale')
    top.set_ylabel('Time (s)')
    for algo in ["topsort", "chinese"]:
        dg = data[algo][graph]
        if algo == "topsort":
            algo = "TopSort_BH"
        if algo == "chinese":
            algo = "iBlackholeDC"
        for threads in [1, 8]:
            for cond in ["cond", "nocond"]:
                if not threads in dg:
                    continue
                dgc = dg[threads][cond]
                scales = dgc.keys()
                times = dgc.values()
                if len(scales) > 0 and len(scales) == len(times):
                    top.plot(scales, times, label = "{}".format(algo))
                    top.grid(True)
                    top.legend()
#    top.set_yscale('symlog')

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("use python parse_raw_grep.py raw.txt")
        exit(0)

    inputFileName = sys.argv[1]

    # data -> graph -> threads -> scale -> time
    data = dict()

    with open(inputFileName, 'r') as file:
        for line in file:
            descr, seconds, _ = line.strip().split()
            seconds = float(seconds)
            vals =  [x.strip() for x in descr.strip().split('.')]
            cond = "nocond"
            graph, size, threads, useCond, algo = vals

            if int(useCond) == 1:
                cond = "cond"

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

            dgc = dg[threads][cond]

            dgc[size] = seconds

    plotConfig = dict()

    #plotConfig["RMAT"] = list([[4,5], [6,7], [8,8], [9,10]])
    #plotConfig["SSCA2"] = list([[4,5],[6,6]])
    #plotConfig["UR"] = list([[4,5], [6,7], [8,10], [11,12]])

    figid = 0
    for graph in ["RMAT", "SSCA2", "UR"]:
        figid += 1
        PlotGraphFromScale(data, graph, figid)

    figid = 0
    for graph in ["RMAT", "SSCA2", "UR"]:
        figid += 1
        plt.figure(figid).savefig(inputFileName + "{}_{}.pdf".format(figid, graph))

    plt.show()
