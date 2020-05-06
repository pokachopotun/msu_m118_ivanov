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
        time = dg[size]["time"]
        threads = dg[size]["threads"]
        top.plot(threads, time, label = str(size))
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
            graph, size, threads, algo = [x.strip() for x in descr.strip().split('.')]
            size = int(size)
            threads = int(threads)

            dg = data[graph]
            if not size in dg:
                dg[size] = dict()
                dg[size]["time"] = list()
                dg[size]["threads"] = list()

            dg[size]["time"].append(seconds)
            dg[size]["threads"].append(threads)

    PlotGraph(data, "RMAT", 1, 8, 10)
    PlotGraph(data, "RMAT", 2, 4, 7)
    PlotGraph(data, "SSCA2", 3, 4, 5)
    PlotGraph(data, "SSCA2", 4, 6, 6)
    PlotGraph(data, "UR", 5, 4, 7)
    PlotGraph(data, "UR", 6, 8, 12)

    plt.figure(1).savefig(inputFileName + "{}_RMAT.pdf".format(1))
    plt.figure(2).savefig(inputFileName + "{}_RMAT.pdf".format(2))
    plt.figure(3).savefig(inputFileName + "{}_SSCA2.pdf".format(3))
    plt.figure(4).savefig(inputFileName + "{}_SSCA2.pdf".format(4))
    plt.figure(5).savefig(inputFileName + "{}_UR.pdf".format(5))
    plt.figure(6).savefig(inputFileName + "{}_UR.pdf".format(6))

    plt.show()
