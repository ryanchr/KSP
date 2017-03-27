#!/usr/bin/python
#Name: Ren Chen
#Date: 03/22/2017

import getopt
import sys
import copy
import random
import collections

def genVFile(numVertex):
	fileName = "./test_in/vertex_"+str(numVertex)+".csv"

	with open(fileName, "w") as vFile:
		vFile.write("NodeID, NodeName\n")
		for i in range(numVertex):
			vFile.write(str(i)+", node_"+str(i)+"\n")


def genEFile(numVertex, numEdge):
	fileName = "./test_in/edge_"+str(numVertex)+".csv"
	maxBW, maxCost, maxPhyD = 10, 100, 100
	edgeSet = set()
	conMap = collections.defaultdict(set)

	with open(fileName, "w") as eFile:
		eFile.write("LinkID, SrcID, DstID, MaxBW, UpResidue, DownResidue, Cost, PhyDistance\n")

		while not allConnect(conMap, numVertex):
			for i in range(len(edgeSet), len(edgeSet)+numEdge):
				link_id = i
				while True:
					src_id, dst_id = random.sample(range(numVertex), 2) 
					if str(src_id)+"_"+str(dst_id) not in edgeSet:
						edgeSet.add(str(src_id)+"_"+str(dst_id))
						break
				conMap[src_id] = conMap[src_id] | conMap[dst_id] | set([dst_id])
				max_bw = 10
				up_residue, down_residue = random.randint(1, maxBW), random.randint(1, maxBW)
				cost, phy_distance = random.randint(1, maxCost), random.randint(1, maxPhyD)
				printOut = [link_id, src_id, dst_id, max_bw, up_residue, down_residue, cost, phy_distance]
	
				eFile.write(','.join([str(x) for x in printOut])+"\n")


def allConnect(conMap, numVertex):
	res = True
	for x in conMap:
		if len(conMap[x]) < numVertex-1:
			res = False
	return res and len(conMap.keys()) == numVertex

	
##def genEFile(numVertex, numEdge):
##	fileName = "./test_in/edge_"+str(numVertex)+".csv"
##	maxBW, maxCost, maxPhyD = 10, 100, 100
##	edgeSet = set()
##
##	with open(fileName, "w") as vFile:
##		vFile.write("LinkID, SrcID, DstID, MaxBW, UpResidue, DownResidue, Cost, PhyDistance\n")
##		for i in range(numVertex):
##			src_id = i
##			for j in range(numVertex):
##				if j != i:
##					dst_id = j
##					link_id = i*numVertex+j
##					max_bw = 10
##					up_residue, down_residue = random.randint(1, maxBW), random.randint(1, maxBW)
##					cost, phy_distance = random.randint(1, maxCost), random.randint(1, maxPhyD)
##					printOut = [link_id, src_id, dst_id, max_bw, up_residue, down_residue, cost, phy_distance]
##					vFile.write(','.join([str(x) for x in printOut])+"\n")

######Run 
numVertex, numEdge = 400, 400
genVFile(numVertex)
genEFile(numVertex, numEdge) 


		