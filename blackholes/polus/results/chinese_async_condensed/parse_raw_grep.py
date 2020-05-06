import matplotlib.pyplot as plt

import sys
import os

def PlotGraph(data, graph, figid, size_a, size_b):
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
      #  top.set_yscale('symlog')

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("use python parse_raw_grep.py raw.txt")
        exit(0)

    inputFileName = sys.argv[1]

    # data -> graph -> size -> threads -> time
    data = dict()
    data["RMAT"] = dict()
    data["SSCA2"] = dict()
    data["UR"] = dict()

    with open(inputFileName, 'r') as file:
        for line in file:
            descr, seconds, _ = line.strip().split()
            seconds = float(seconds)
            vals =  [x.strip() for x in descr.strip().split('.')]
            cond = "nocond"
            graph, size, threads, algo = vals[:4]
            if len(vals) > 4:
                cond = vals[4]

            size = int(size)
            threads = int(threads)

            dg = data[graph]
            if not size in dg:
                dg[size] = dict()
                dg[size]["cond"] = dict()
                dg[size]["nocond"] = dict()
                dg[size]["cond"]["time"] = list()
                dg[size]["cond"]["threads"] = list()
                dg[size]["nocond"]["time"] = list()
                dg[size]["nocond"]["threads"] = list()

            dgc = dg[size][cond]

            dgc["time"].append(seconds)
            dgc["threads"].append(threads)

    plotConfig = dict()

    plotConfig["RMAT"] = list([[4,5], [6,7], [8,8], [9,10]])
    plotConfig["SSCA2"] = list([[4,5],[6,6]])
    plotConfig["UR"] = list([[4,5], [6,7], [8,10], [11,12]])

    figid = 0
    for graph in plotConfig:
        for a, b in plotConfig[graph]:
            figid += 1
            PlotGraph(data, graph, figid, a, b)

    figid = 0
    for graph in ["RMAT", "SSCA2", "UR"]:
        for a, b in plotConfig[graph]:
            figid += 1
            plt.figure(figid).savefig(inputFileName + "{}_{}.pdf".format(figid, graph))

    plt.show()
