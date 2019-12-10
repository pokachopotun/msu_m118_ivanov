

if __name__ == "__main__":
    proc = [(8, 16), (16, 16)]
    grid = [8192, 16384, 32768]
    for p in proc:
        totalProc = p[0] * p[1]
        for g in grid:
            outfile = "outputs/bluegene." + str(totalProc) + "." + str(g)
            errfile = "errors/bluegene." + str(totalProc) + "." + str(g)
            params = "{0} {1} {2} {3} {4} {5} {6} {7}".format(g, g, p[0], p[1], 20, 1.0, 1.0, 1)
            if totalProc == 128:
                mins = 15
            else:
                mins = 10
            cmd = "mpisubmit.pl -n {0} -m smp -w {4}:00 --stderr {1} --stdout {2} ./bluegene -- {3}".format(totalProc, outfile, errfile, params, mins)
            print(cmd)
