import os
import sys
if __name__ == "__main__":
	if len(sys.argv) < 2:
		print("Use python gen_circle.py nNodes")
		exit()
	nNodes = int(sys.argv[1])
	outputFileName = "circle_" + str(nNodes) + ".graph"
	with open(outputFileName, 'w') as file:
		for i in range(nNodes):
			file.write(str(i) + " " + str( (i + 1) % nNodes ) + os.linesep )
