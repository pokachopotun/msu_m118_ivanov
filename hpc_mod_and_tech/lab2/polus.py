

if __name__ == "__main__":
    proc = [(1, 1), (2, 5), (4, 5), (5, 8)]
    grid = [2048, 4096, 8192]
    for p in proc:
        totalProc = p[0] * p[1]
        for g in grid:
            outfile = "outputs/polus." + str(totalProc) + "." + str(g)
            errfile = "errors/polus." + str(totalProc) + "." + str(g)
            params = "{0} {1} {2} {3} {4} {5} {6} {7}".format(g, g, p[0], p[1], 20, 1.0, 1.0, 1)
            cmd = "mpisubmit.pl -p {0} -w 60:00 --stderr {1} --stdout {2} ./polus -- {3}".format(totalProc, outfile, errfile, params)
            print(cmd)
