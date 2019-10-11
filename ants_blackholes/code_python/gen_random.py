import os
import sys
import random
if __name__ == "__main__":
	if len(sys.argv) < 3:
		print("Use python gen_random.py nNodes density")
		exit()
	nNodes = int(sys.argv[1])
	density = float(sys.argv[2])
	outputFileName = "random_" + str(nNodes) + ".graph"
	with open(outputFileName, 'w') as file:
		for i in range(nNodes):
			s = str(i)
			for j in range(nNodes):
				if i == j:
					continue
				if random.random() < density:
					s += " " + str(j)
			file.write( s + os.linesep )
