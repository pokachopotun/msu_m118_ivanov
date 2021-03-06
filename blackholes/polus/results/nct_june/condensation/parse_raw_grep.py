import matplotlib.pyplot as plt

import sys
import os

def PlotGraphFromThreads(data, graph, figid, size_a, size_b):
    fig = plt.figure(figid)
    fig.suptitle(graph, fontsize=16)
    fig.tight_layout()
    top = fig.add_subplot(1, 1, 1)
    top.set_xlabel('Threads')
    top.set_ylabel('Time (s)')
    dg = data[graph]
    for size in range(size_a, size_b + 1):
        for cond in ["cond", "nocond"]:
            dgc = dg[size][cond]
            time = dgc["time"]
            threads = dgc["threads"]
            top.plot(threads, time, label = "scale={},{}".format(size, cond))
            top.grid(True)
            top.legend()
    top.set_yscale('symlog')

def PlotGraphSimple(data, graph, figid):
    fig = plt.figure(figid)
    fig.suptitle("{}".format(graph), fontsize=16)
    fig.tight_layout()
    top = fig.add_subplot(1, 1, 1)
    top.set_xlabel('Graph scale')
    top.set_ylabel('Time (s)')
    dg = data[graph]
    scales = dg.keys()
    times = dg.values()
    if len(scales) > 0 and len(scales) == len(times):
        top.plot(scales, times, label = "Condensation")
        top.grid(True)
        top.legend()
    top.set_yscale('symlog')

def PlotGraphFromScale(data, graph, figid):
    fig = plt.figure(figid)
    fig.suptitle("{}".format(graph), fontsize=16)
    fig.tight_layout()
    top = fig.add_subplot(1, 1, 1)
    top.set_xlabel('Graph scale')
    top.set_ylabel('Time (s)')
    dg = data[graph]
    for threads in [1, 8]:
        for cond in ["cond", "nocond"]:
            dgc = dg[threads][cond]
            scales = dgc.keys()
            times = dgc.values()
            if len(scales) > 0 and len(scales) == len(times):
                top.plot(scales, times, label = "Condensation,threads={},{}".format(threads, cond))
                top.grid(True)
                top.legend()
    top.set_yscale('symlog')

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("use python parse_raw_grep.py raw.txt")
        exit(0)

    inputFileName = sys.argv[1]

    # data -> graph -> threads -> scale -> time
    data = dict()
    data["RMAT"] = dict()
    data["SSCA2"] = dict()
    data["UR"] = dict()

    with open(inputFileName, 'r') as file:
        for line in file:
            descr, seconds, _ = line.strip().split()
            seconds = float(seconds)
            vals =  [x.strip() for x in descr.strip().split('.')]
            graph, size, algo = vals

            size = int(size)

            data[graph][size] = seconds

    plotConfig = dict()

    #plotConfig["RMAT"] = list([[4,5], [6,7], [8,8], [9,10]])
    #plotConfig["SSCA2"] = list([[4,5],[6,6]])
    #plotConfig["UR"] = list([[4,5], [6,7], [8,10], [11,12]])
    plotConfig["RMAT"] = list([[4,10]])
    plotConfig["SSCA2"] = list([[4,6]])
    plotConfig["UR"] = list([[4,12]])

    figid = 0
    for graph in ["RMAT", "SSCA2", "UR"]:
        figid += 1
        PlotGraphSimple(data, graph, figid)

    figid = 0
    for graph in ["RMAT", "SSCA2", "UR"]:
        figid += 1
        plt.figure(figid).savefig(inputFileName + "{}_{}.pdf".format(figid, graph))

    plt.show()
