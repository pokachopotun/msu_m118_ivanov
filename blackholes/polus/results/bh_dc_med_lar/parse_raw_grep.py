import matplotlib.pyplot as plt

import operator

import sys
import os

def PlotGraphFromScale(data, graph, figid):
    fig = plt.figure(figid)
    fig.suptitle("{}".format(graph), fontsize=16)
    fig.tight_layout()
    top = fig.add_subplot(1, 1, 1)
    top.set_xlabel('Graph scale')
    top.set_ylabel('Blackholes found')
    dg = data[graph]
    for threads in [1, 8]:
        if not threads in dg:
            continue
        for cond in ["cond", "nocond"]:
            dgc = dg[threads][cond]
            scales = dgc.keys()
            values = dgc.values()
            if len(scales) > 0 and len(scales) == len(values):
                scales, values = zip(*sorted(zip(scales, values)))
                scales = scales[:4]
                values = values[:4]
                print(threads, cond, scales, values)
                top.plot(scales, values, label = "iBlackholeDC,threads={},{}".format(threads, cond))
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
            descr, blackholes = line.strip().split()
            blackholes = int(blackholes)
            vals =  [x.strip() for x in descr.strip().split('.')]
            print(vals)
            cond = "nocond"
            graph, size, threads, useCond, algo = vals

            if int(useCond) == 1:
                cond = "cond"

            size = int(size)
            threads = int(threads)

            dg = data[graph]
            if not threads in dg:
                dg[threads] = dict()
                dg[threads]["cond"] = dict()
                dg[threads]["nocond"] = dict()

            dgc = dg[threads][cond]
            dgc[size] = blackholes

    figid = 0
    for graph in ["RMAT"]:
        figid += 1
        PlotGraphFromScale(data, graph, figid)

    figid = 0
    for graph in ["RMAT"]:
        figid += 1
        plt.figure(figid).savefig(inputFileName + "{}_{}.pdf".format(figid, graph))

    plt.show()
