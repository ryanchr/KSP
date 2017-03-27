#!/usr/bin/python
#name: Ren Chen
#date: Feb/27/2017

import getopt
import sys
import copy 
import matplotlib.pyplot as plt
import numpy as np

def readFiles():
    args = sys.argv[1:]
    linkFileName = args[0]
    pathFileName = args[1]
    
    with open(linkFileName, 'r') as inputFile:
    	inputFile.readline()
        in_link = inputFile.read().splitlines()
    
    with open(pathFileName, 'r') as inputFile:
        in_agg_path = inputFile.read().splitlines()

    return in_link, in_agg_path
        
def getEdgesCost(in_link):
	edge_cost = {}
	for x in in_link:
		x_list = [y.strip() for y in x.split(',')]
		edge_cost[(x_list[1],x_list[2])] = int(x_list[3])
	return edge_cost

def getPathCost(edge_cost, path):
	l, res = len(path), 0
	for i in range(l-1):
		res += edge_cost[(path[i],path[i+1])] if (path[i],path[i+1]) in edge_cost  \
											  else edge_cost[(path[i+1],path[i])]
	return res

def getAggPairsCost(edge_cost, in_agg_path):
	agg_pairs_cost = {}
	for x in in_agg_path:
		x_list = [y.strip() for y in x.split(',')]
		if "AggDst" not in x_list[1]:
			path = x_list[2].split('|'); path.pop()
			pair = (x_list[0],x_list[1])
			agg_pairs_cost[pair] = agg_pairs_cost.get(pair,[]) + [getPathCost(edge_cost, path)]
	return agg_pairs_cost

def getCostHist(agg_pairs_cost):
	cost_hist = {}; 
	for x in agg_pairs_cost:
		for y in agg_pairs_cost[x]:
			cost_hist[y] = cost_hist.get(y,0) + 1
	plt.bar(range(len(cost_hist)), cost_hist.values(), align='center')
	plt.xticks(range(len(cost_hist)), cost_hist.keys())
	plt.show()

def getPairCostHist(pari, pair_cost):
	num_bins = 30
	plt.hist(pair_cost, num_bins, color='green')
	plt.show()


def test():
	in_link, in_agg_path = readFiles()
	edge_cost = getEdgesCost(in_link)
	agg_pairs_cost = getAggPairsCost(edge_cost, in_agg_path)
	#getCostHist(agg_cost)

test()


