import os
import sys

def sc(x):
    return str(x) #'{:.2e}'.format(x)

if __name__ == "__main__":
    
    resFile = sys.argv[1]

    base = dict()
    with open(resFile, 'r') as file:
        for line in file:
            a = line.strip().split()

            pc = a[0]
            nodes = int(a[1])
            grid = int(a[2])
            omp = int(a[4])
            time = float(a[6])
            delta = float(a[14]) 

            if nodes == 1 and omp == 1:
                base[grid] = time


    with open(resFile, 'r') as file:
        delim = " & "
        s = "HPC" + delim 
        s += "Nodes" + delim
        s += "Grid" + delim
        s += "OMP Threads" + delim
        s += "Wall time" + delim
        s += "Delta" + delim
        s += "Speedup" + delim
        s += "Efficiency" + r" \\"
        print(s)
        for line in file:
            a = line.strip().split()
            pc = a[0]
            nodes = int(a[1])
            grid = int(a[2])
            omp = int(a[4])
            time = float(a[6])
            delta = float(a[14]) 
            speedup = base[grid] / time
            expected = nodes * omp
            efficiency = speedup  / expected

            lst = [pc, nodes, grid, omp, time, delta, speedup, efficiency]
            delim = " & "
            s = ""
            s += pc + delim 
            s += str(nodes) + delim
            s += str(grid) + delim
            s += str(omp) + delim
            s += sc(time) + delim
            s += sc(delta) + delim
            s += sc(speedup) + delim
            s += sc(efficiency) + r" \\"

            print(s)
