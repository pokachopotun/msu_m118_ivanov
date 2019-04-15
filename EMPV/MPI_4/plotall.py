import os
import sys
import matplotlib.pyplot as plt

if __name__ == "__main__":
	curDir = os.getcwd()
	targetDir = os.path.join( curDir, sys.argv[1] )
	files = os.listdir( targetDir )
	for fileName in files:
		filePath = os.path.join( targetDir, fileName )
		with open( filePath, 'r' ) as file:
			lines = [ x.strip().split(' ') for x in file.readlines() ]
		data = list();
		nProc, L1, *rest = [ x.strip() for x in fileName.strip().split('.') ]
		nProc = int(nProc)
		L1 = int(L1)
		for i in range( nProc ):
			data.append( list( [ 0 for i in range(20000) ] ) )
		for line in lines:
			if line[0] != 'final:':
				dummy, dummy, it, dummy2, rank, dummy3, load = line
			
			it = int(it)
			rank = int(rank)
			load = int(load)
			data[ rank ][ it ] = load
		for i in range( nProc ):
			plt.plot( data[i])
		plt.xlabel('Iteration')
		plt.ylabel('Load')
		plt.grid(True)
		title = "Load(Iter) nProc " + str(nProc) + ", System L" + str(L1)
		plt.title(title)
		plt.savefig( title  + '.png')
		plt.show();

