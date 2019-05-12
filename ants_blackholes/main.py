import networkx as nx
import matplotlib.pyplot as plt
import itertools
import numpy
thermoSpeed = 0.25
cooling = 0.95
antsHeating = 0.005

class Node:
	def __init__(self, G, idx):
		self.id = idx
		self.heat = 0
		self.newHeat = self.heat
		self.antsCnt = 0
		self.G = G
		self.recvAntsCnt = 0
		
	def GetLabel(self):
		return str(self.id) + " / " + str( self.antsCnt ) #+ " / " + str(self.heat)

	def CalcNewHeat(self):
		sumHeat = self.heat
		cnt = 1
		
		for neighb in self.G.successors(self):
			if neighb.heat < self.heat:
				sumHeat += neighb.heat
				cnt += 1
		
		for neighb in self.G.predecessors(self):
			if neighb.heat > self.heat:
				sumHeat += neighb.heat
				cnt += 1

		#print("Cnt " + str(cnt) + " sumHeat " + str(sumHeat) )	
		
		targetHeat = sumHeat / cnt
		self.newHeat = self.heat + ( targetHeat - self.heat ) * thermoSpeed
		
#		self.newHeat = min( self.newHeat, 1.0 )
#		self.newHeat = max( self.newHeat, 0.0 )

	def Cool(self):
		self.heat *= cooling

	def HeatAnts(self):
		self.heat += self.antsCnt * antsHeating

	def UpdateHeat(self):
		self.heat = self.newHeat;

	def Spawn(self, cnt):
		self.antsCnt += cnt

	def GetColor(self):
		return ( min( 1.0, self.heat ), 0.0, max(0.0, 1.0 - self.heat))	
	
	def SendAnts(self):
		totalHeat = 0
		maxHeat = 0
		cnt = 0
		ns = list(self.G.successors(self))
		if len(ns) == 0:
			return
		for neighb in ns:
			totalHeat += neighb.heat
			maxHeat = max( maxHeat, neighb.heat )
		p = list()
		psum = 0.0
		for neighb in ns:
			p.append( maxHeat - neighb.heat )
			psum += p[-1]
		if psum == 0:
			p = [ 1.0 / len(p) for x in p ]
		else:
			p = [ x / psum for x in p ]
		for _ in range(self.antsCnt):
			i = numpy.random.choice( range(len(ns)), p=p)
			ns[i].recvAntsCnt += 1
		self.antsCnt = 0
			
	def RecvAnts(self):
		self.antsCnt += self.recvAntsCnt
		self.recvAntsCnt = 0
class Solution:
	def __init__(self):
		self.Iter = 0
		self.G = nx.DiGraph()
		
		self.nodes = [ Node(self.G, x) for x in range(10) ]
		self.G.add_nodes_from( self.nodes )

		edges = [ (0, 1), (1, 2), (2,3), (3,0), (0, 4), (5, 0), (5,6), (6,7), (7,5), (7,1), (7,8), (8,6) ]
		self.G.add_edges_from( [ (self.nodes[x[0]], self.nodes[x[1]]) for x in edges] )
		self.colorMap = [ (0, 0, 0) for node in nx.nodes(self.G) ]
		self.UpdateColors()
		self.labels = dict()
		self.UpdateLabels()
		self.CheatInit()

	def CheatInit(self):
		#self.nodes[0].heat = 1.0
		#self.nodes[0].antsCnt = 10
		for node in nx.nodes(self.G):
			if node.id >= 5:
				node.Spawn(10)
	
	def UpdateLabels(self):
		for node in nx.nodes(self.G):
			self.labels[node] = node.GetLabel()

	def UpdateColors(self):
		for node in nx.nodes(self.G):
			self.colorMap[node.id] = node.GetColor()

	def Show(self):
		self.UpdateColors()
		self.UpdateLabels()
		#print(self.colorMap)
		nx.draw_shell(self.G, with_labels = True, labels = self.labels, node_color = self.colorMap, font_color = 'k' )
	#	plt.show()

	def Diffuse(self):
		for node in nx.nodes(self.G):
			node.CalcNewHeat()
		for node in nx.nodes(self.G):
			node.UpdateHeat()		

	def Cool(self):
		for node in nx.nodes(self.G):
			node.Cool()		
	
	def MoveAnts(self):
		for node in nx.nodes(self.G):
			node.SendAnts()
		for node in nx.nodes(self.G):
			node.RecvAnts()	
	def Heat(self):
		for node in nx.nodes(self.G):
			node.HeatAnts()

	def RunStep(self):
		self.Iter += 1
		self.MoveAnts()
		self.Heat()
		self.Diffuse()
		self.Cool()

	def SaveFigure(self):
		plt.savefig("vis_{:04d}.png".format( self.Iter ) )
		plt.clf()
	def Run(self):
		while True:
			self.RunStep()
			self.Show()
			self.SaveFigure()
		
if __name__ == "__main__":
	Solution().Run()
