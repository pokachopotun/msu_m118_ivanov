class Submission:
    def __init__(self, typeStr1, typeStr2):
        self.typeStr1 = typeStr1
        self.typeStr2 = typeStr2
        self.command = r"mpisubmit.pl = -g --stdout=" + self.typeStr1 + "." 
        pass

## mpisubmit.pl -g --stdout=$(1).RMAT.out --stderr=$(1).RMAT.err ./color -- --graph-scale=$(1) --graph-edgefactor=32 --test-run=0 --graph-type=rmat
if __name__ == "__main__":
    sub = list()
    with open("config.txt") as file:
        n, a, b, t = [int(x) for x in file.readline().strip().split()]
        for i in range(t):
            typeStr1, typeStr2 = file.readline().strip().split()
            sub.append(Submission(typeStr1, typeStr2))

    for s in sub:
        for scale in range(a, b + 1):
            for rep in range(n):
                
                scale_str = "{:02d}".format(scale) 
                rep_str = "{:02d}".format(rep)
                id_str = s.typeStr1 + "." + scale_str + "." + rep_str
                outfile = "outputs/" + id_str
                errfile = "errors/" + id_str
                graphfile = "inputs/" + s.typeStr1 + "." + scale_str

                cmd0 = "mpisubmmit.pl -g --stdout=" + outfile + " --stderr=" + errfile
                if s.typeStr2 == "matrix":
                    cmd1 = "./color -- --test-run=0 --graph-type=" + s.typeStr2 + " --graph-file=" + graphfile
                else:
                    cmd1 = "./color -- --graph-edgefactor=32 --test-run=0" + " --graph-scale=" + scale_str + " --graph-type=" + s.typeStr2  
                cmd = cmd0 + " " + cmd1
                print(cmd)
                os.system(cmd)
         
