import os

class Submission:
    def __init__(self, typeStr1, typeStr2):
        self.typeStr1 = typeStr1
        self.typeStr2 = typeStr2
        self.command = r"mpisubmit.pl = -g --stdout=" + self.typeStr1 + "." 
        pass

if __name__ == "__main__":
    sub = list()
    with open("config.txt") as file:
        n, a, b, t = [int(x) for x in file.readline().strip().split()]
        for i in range(t):
            typeStr1, typeStr2 = file.readline().strip().split()
            sub.append(Submission(typeStr1, typeStr2))
    
    for s in sub:
        for scale in range(a, b + 1):
            scale_str = "{:02d}".format(scale) 
            id_str = s.typeStr1 + "." + scale_str
            outfile = "outputs/" + id_str
            errfile = "errors/" + id_str
            graphfile = "graphs/" + id_str + ".bin"
            inputfile = "inputs/" + id_str

            cmd1 = "./generator/generator -s " + scale_str + " -unweighted -file " + graphfile + " -type " + s.typeStr1
            #cmd2 = "./graph_to_file/graph_to_file " + graphfile + " " + inputfile
            print(cmd1)
            #print(cmd2)

            os.system(cmd1)
            #os.system(cmd2)


#define graph
#	./generator -s $(1) -directed -weighted -file inputs/$(1).$(2).graph.bin -type $(2)
#	./graph_to_file inputs/$(1).$(2).graph.bin inputs/$(1).$(2).graph
#endef
