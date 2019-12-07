import matplotlib.pyplot as plt

class Res:
    def __init__(self, sz, bh, t):
        self.sz = sz
        self.bh = bh
        self.t = t

def PlotOne(figid, graph, name, results, form):
    time = list()
    bh_found = list()
    size = list()
    for res in results:
        sz = res.sz
        bh = res.bh
        t = res.t

        size.append(sz)
        bh_found.append(bh)
        time.append(t)

    fig = plt.figure(figid)
    fig.tight_layout()
    top = fig.add_subplot(2, 1, 1)
    bot = fig.add_subplot(2, 1, 2)

    top.set_yscale('symlog')
    top.plot(size, time, form, label = name)
    top.grid(True)
    top.legend()
    top.set_xlabel('Graph scale')
    top.set_ylabel('Time (s)')

    if graph != 'UR':
        bot.set_yscale('symlog')
    bot.plot(size, bh_found, form, label = name)
    bot.grid(True)
    bot.legend()
    bot.set_xlabel('Graph scale')
    bot.set_ylabel('Number of blackholes')


def PlotGraph(graph, res, figid):
    f = ['b', 'r']
    keys = list(res.keys());
    for i in range(len(keys)):
        alg = keys[i]
        PlotOne(figid, graph, alg, res[alg], f[i])

if __name__ == "__main__":
    

    rmat = dict()

    rmat["iBlackhole"] = list()
    rmat["TopSort"] = list()

    ssca = dict()

    ssca["iBlackhole"] = list()
    ssca["TopSort"] = list()

    ur = dict()

    ur["iBlackhole"] = list()
    ur["TopSort"] = list()

    with open("results.txt", 'r') as file:
        for line in file:
            graph, sz, alg, bh, t = [x.strip() for x in line.strip().split()]
            sz = int(sz)
            t = float(t)
            bh = int(bh)
    
            if graph == "RMAT":
                rmat[alg].append(Res(sz, bh, t))
            if graph == "SSCA2":
                ssca[alg].append(Res(sz, bh, t))
            if graph == "UR":
                ur[alg].append(Res(sz, bh, t)) 

    #plot rmat

    

    PlotGraph("RMAT", rmat, 1)
    PlotGraph("SSCA2", ssca, 2)
    PlotGraph("UR", ur, 3)

    plt.figure(1).savefig("rmat.pdf")
    plt.figure(2).savefig("SSCA2.pdf")
    plt.figure(3).savefig("UR.pdf")

    plt.show()
