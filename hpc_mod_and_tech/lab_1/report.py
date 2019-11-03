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
            data[name].edges.append(int(line[2]))
    arr = list()
    for name in data:
        teps_avg = 0
        for t, e in zip(data[name].time, data[name].edges):
            teps = e / t * 1000
            data[name].teps.append(teps)
            teps_avg += teps
        teps_avg /= len(data[name].time)
        name1, scale = name.split('.')
        scale = int(scale)
        #print(name1, scale, teps_avg)
        arr.append((name, scale, teps_avg))
    arr.sort()
    for x in arr:
        print(x)
#   print(lines)
