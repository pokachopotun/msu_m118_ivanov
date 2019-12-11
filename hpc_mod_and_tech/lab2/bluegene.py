import sys 
import os

if __name__ == "__main__":
    
    submit = False
    if len(sys.argv) == 2:
        submit = sys.argv[1] == "submit"

    proc = [(8, 16), (16, 16), (16, 32)]
    grid = [8192, 16384, 32768]
    for ompNumProc in range(1, 4): 
        for p in proc:
            totalProc = p[0] * p[1]
            for g in grid:
                filename = "bluegene." + str(totalProc) + "." + str(g) + "." + str(ompNumProc)
                errfile = "errors/" + filename
                outfile = "outputs/" + filename
                params = str(p[0]) + " " + str(p[1]) + " " + str(g) + " " + str(g) + " " + str(20) + " " + str(1.0) + " " + str(1.0) + " " + str(ompNumProc)
                if totalProc == 128:
                    mins = 15
                if totalProc == 256:
                    mins = 10
                if totalProc == 512:
                    mins = 5 
                cmd = "mpisubmit.bg -n " + str(totalProc) + " -m smp -w " + str(mins) + ":00 --stdout " + outfile + " --stderr " + errfile + " ./bluegene -- " + params
                print(cmd)
                if submit:
                    os.system(cmd)
