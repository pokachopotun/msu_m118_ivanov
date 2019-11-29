import os

class RunData:
    def __init__(self):
        self.time = list()
        self.teps = list()
        self.edges = list()

if __name__ == "__main__":

    os.system("grep \"non-zero\|total time\|Degree Histogram\" outputs/* | tr \,\:\(\) \"\\ \" | tr \"abcdefghijklmnopqrstuvwxyz/HD\" \"\\ \" > rep.txt")
    with open("rep.txt") as file:
        lines = [x.strip().split() for x in file.readlines()]
    data = dict()
    for i in range(len(lines)):
        line = lines[i]
        
        name = line[0]
        name = name[:-3]

        if name not in data:
            data[name] = RunData()

        if len(line) == 2:
            data[name].time.append(float(line[1]))
        if len(line) == 3:
            _, vert, edges = line
            data[name].edges.append(int(edges))
    arr = list()
    for name in data:
        teps_avg = 0
        time_avg = 0
        edges_avg = 0
        for t, e in zip(data[name].time, data[name].edges):
            teps = e / t * 1000
            data[name].teps.append(teps)
            teps_avg += teps
            time_avg += t
            edges_avg += e
        if len(data[name].time) == 0:
            name1, scale = name.split('.')
            arr.append((name, scale, 0, 0, 0))
            continue
        teps_avg /= len(data[name].time)
        time_avg /= len(data[name].time)
        edges_avg /= len(data[name].time)
        name1, scale = name.split('.')
        scale = int(scale)
        arr.append((name, scale, teps_avg, time_avg, edges_avg))
    arr.sort()
    with open("report.txt", 'w') as file:
        file.write("{} {} {} {} {}\n".format("name", "scale", "teps_avg", "time_avg_ms", "edges_avg"))
        for x in arr:
            file.write("{} {} {} {} {}\n".format(x[0], x[1], x[2], x[3], x[4]))
