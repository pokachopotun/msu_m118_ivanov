#!/usr/bin/python3
import matplotlib.pyplot as plt

import operator

import sys

value = "count"

def PlotGraphFromScale(data, graph, figid):
    fig = plt.figure(figid)
    fig.suptitle("{}".format(graph), fontsize=16)
    fig.tight_layout()
    top = fig.add_subplot(1, 1, 1)
    top.set_xlabel('Threads')
    top.set_ylabel('Blackholes {}'.format(value))
    for algo in ["chinese"]:
        dga = data[algo][graph]
        if algo == "chinese":
            algo = "iBlackholeDC"
        if algo == "topsort":
            algo = "TopSort_BH"
        for cond in ["cond",]:
            for binding in [0, 1]:
                dgc = dga[18][cond][binding]
                if binding == 0:
                    binding = "no binding"
                if binding == 1:
                    binding = "core binding"
                scales = dgc.keys()
                values = dgc.values()
                if len(scales) > 0 and len(scales) == len(values):
                    scales, values = zip(*sorted(zip(scales, values)))
                    scales = scales
                    values = values
                    print(cond, scales, values)
                    top.plot(scales, values, label = "{}.{}".format(algo, binding))
                    top.grid(True)
                    top.legend()
#    top.set_xscale('symlog')

def PlotSpeedup(data, graph, figid):
    fig = plt.figure(figid)
    fig.suptitle("{}".format(graph), fontsize=16)
    fig.tight_layout()
    top = fig.add_subplot(1, 1, 1)
    top.set_xlabel('Threads')
    top.set_ylabel('Speedup'.format(value))
    for algo in ["chinese"]:
        dga = data[algo][graph]
        if algo == "chinese":
            algo = "iBlackholeDC"
        if algo == "topsort":
            algo = "TopSort_BH"
        for cond in ["cond",]:
            for binding in [0, 1]:
                dgc = dga[18][cond][binding]
                if binding == 0:
                    binding = "no binding"
                if binding == 1:
                    binding = "core binding"
                scales = dgc.keys()
                values = dgc.values()
                if len(scales) > 0 and len(scales) == len(values):
                    scales, values = zip(*sorted(zip(scales, values)))
                    single = values[1] / 2
                    values = [float(x) / single for x in values]
                    print(cond, scales, values)
                    top.plot(scales, values, label = "{}.{}".format(algo, binding))
#    top.set_xscale('symlog')
    scales = values = [x for x in range(1, 21)]
    print(cond, scales, values)
    top.plot(scales, values, label = "Linear speedup")
    top.grid(True)
    top.legend()

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("use python parse_raw_grep.py raw.txt value")
        exit(0)

    inputFileName = sys.argv[1]
    value = sys.argv[2]

    # data -> graph -> threads -> scale -> time
    data = dict()
    for algo in ["chinese","topsort"]:
        data[algo] = dict()
        data[algo]["RMAT"] = dict()
        data[algo]["SSCA2"] = dict()
        data[algo]["UR"] = dict()

    with open(inputFileName, 'r') as file:
        for line in file:
            descr, blackholes = line.strip().split()
            blackholes = int(blackholes)
            vals =  [x.strip() for x in descr.strip().split('.')]
            print(vals)
            cond = "nocond"
            graph, size, threads, useCond, binding, algo = vals

            binding = int(binding)
            if int(useCond) == 1:
                cond = "cond"

            size = int(size)
            threads = int(threads)

            dg = data[algo][graph]
            if not size in dg:
                dg[size] = dict()
                dg[size]["cond"] = dict()
                dg[size]["nocond"] = dict()
                dg[size]["cond"][0] = dict()
                dg[size]["cond"][1] = dict()
                dg[size]["nocond"][0] = dict()
                dg[size]["nocond"][1] = dict()

            dgc = dg[size][cond][binding]
            dgc[threads] = blackholes

    figid = 0
    for graph in ["RMAT"]:
        figid += 1
        PlotGraphFromScale(data, graph, figid)
        figid += 1
        PlotSpeedup(data, graph, figid)

    figid = 0
    for graph in ["RMAT"]:
        figid += 1
        plt.figure(figid).savefig(inputFileName + "{}_{}.pdf".format(figid, graph))
        figid += 1
        plt.figure(figid).savefig(inputFileName + "{}_{}_speedup.pdf".format(figid, graph))

    plt.show()
