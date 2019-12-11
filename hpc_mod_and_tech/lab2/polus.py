import sys
import os

if __name__ == "__main__":
    run = False
    if len(sys.argv) == 2:
        run  = sys.argv[1] == "submit"

    proc = [(1, 1), (2, 5), (4, 5), (5, 8)]
    grid = [2048, 4096, 8192]
    for p in proc:
        totalProc = p[0] * p[1]
        for g in grid:
            outfile = "outputs/polus." + str(totalProc) + "." + str(g)
            errfile = "errors/polus." + str(totalProc) + "." + str(g)
            params = "{0} {1} {2} {3} {4} {5} {6} {7}".format(p[0], p[1], g, g, 20, 1.0, 1.0, 1)
            cmd = "mpisubmit.pl -p {0} -w 00:30 --stdout {1} --stderr {2} ./polus -- {3}".format(totalProc, outfile, errfile, params)
            print(cmd)
            if run:
                os.system(cmd)    
