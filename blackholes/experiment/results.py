import os


if __name__ == "__main__":
    
    for filename in os.listdir("outputs"):
        path = os.path.join("outputs", filename)
        vals = [x.strip() for x in filename.split('.')]
        graph, size, algo = vals
        res = 0
        time = str(60 * 20)
        with open(path, 'r') as file:
            for line in file:
                line = line.strip()
                vals = line.split()
                if len(vals) == 2:
                    s, cnt = vals
                    if s == "bhcount":
                        res = max(res, int(cnt))
                else:
                    s = vals[0]
                    if s == "It":
                        time = vals[3]
        out = graph + "\t" + str(size) + "\t" + algo + "\t" + str(res) + "\t" + str(time)
        print(out) 
