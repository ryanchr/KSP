#!/usr/bin/python 
#name: Ren Chen 
#date: Feb/27/2017 

import getopt 
import sys 
import copy  
import matplotlib
matplotlib.use('Agg')

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

def getPathCost(edge_cost, path, links): 
	l, res = len(path), 0 
	miss_link = False
	for i in range(l-1): 
		if (path[i],path[i+1]) in edge_cost:
			res += edge_cost[(path[i],path[i+1])]  
		elif (path[i+1],path[i]) in edge_cost:
			res += edge_cost[(path[i+1],path[i])]  
		else:
			miss_link = True
			links.add((path[i],path[i+1]))
		
	return res, miss_link

def getAggPairsCost(edge_cost, in_agg_path): 
	agg_pairs_cost = {}; agg_pairs_miss_link = {}; missed_links = set(); #pairs_miss_link = {}
	for x in in_agg_path: 
		x_list = [y.strip() for y in x.split(',')] 
		if "AggDst" not in x_list[1]: 
			path = x_list[2].split('|'); path.pop() 
			pair = (x_list[0],x_list[1]) 
			cost, miss_link = getPathCost(edge_cost, path, missed_links)
			agg_pairs_cost[pair] = agg_pairs_cost.get(pair,[]) + [cost] 
			agg_pairs_miss_link[pair] = agg_pairs_miss_link.get(pair,[]) + [miss_link]
			##if pair not in pairs_miss_link
			#pairs_miss_link[pair] = True if miss_link else False
	return agg_pairs_cost, agg_pairs_miss_link, missed_links

def getCostHist(agg_pairs_cost): 
	cost_hist = {};  
	for x in agg_pairs_cost: 
		for y in agg_pairs_cost[x]: 
			cost_hist[y] = cost_hist.get(y,0) + 1 
	##plt.bar(range(len(cost_hist)), cost_hist.values(), align='center') 
	##plt.xticks(range(len(cost_hist)), cost_hist.keys()) 
	##plt.show() 
	return cost_hist

def getPairCostHist(pair, pair_cost): 
	num_bins = 30 
	plt.hist(pair_cost, num_bins, range=[0,600000], color='green')
	plt.xlabel("Cost")
	plt.ylabel("Frequency")
	
	min_val = min(pair_cost)
	max_val = max(pair_cost)
	stdev_val =  format(np.std(pair_cost),'.3f')
	avg_val = format(np.mean(pair_cost),'.3f')
	to_print = "Min: "+str(min_val)+", Max: "+str(max_val)+" ,Ave: "+str(avg_val)+ " ,Std: " + str(stdev_val)
	
	plt.title("Cost distribution for AGG pair ("+pair[0]+","+pair[1]+") (Pathsize: "+str(len(pair_cost))+")\n "+to_print)
	plt.savefig("agg_"+pair[0]+"_"+pair[1]+'.png')
	plt.clf()
	#plt.show() 

def getPathsMissingLinks(agg_pairs_miss_link):
	pair_names = []; pair_miss = []; xvals = [i for i in range(1,len(agg_pairs_miss_link)+1)]
	for x in agg_pairs_miss_link:
		pair_names.append(x)
		pair_miss.append(agg_pairs_miss_link[x].count(True)/len(agg_pairs_miss_link[x]))
	plt.plot(xvals, pair_miss)
	plt.xticks(xvals, pair_names, rotation='vertical')
	plt.margins(0.7)
	plt.savefig("missing_paths"+'.png')
	plt.clf()
	
def test(): 
	in_link, in_agg_path = readFiles() 
	edge_cost = getEdgesCost(in_link) 
	agg_pairs_cost, agg_pairs_miss_link, missed_links = getAggPairsCost(edge_cost, in_agg_path) 
	
	#print len(missed_links)
	##
	#getPathsMissingLinks(agg_pairs_miss_link)
	##calculate number of missing links
	#pairs = [('0','2'), ('0','334'), ('2','334'), ('2','394'), ('1','2'), ('0','3')]
	#pairs = []
	#for x in agg_pairs_cost:
	#	if all(y==False for y in agg_pairs_miss_link[x]):
	#		pairs.append(x)
	#	if len(pairs) > 5:
	#		break
	#print pairs
	
	#######Get all pair_cost distribution
	#cost_hist = getCostHist(agg_pairs_cost)
	all_costs = []
	for x in agg_pairs_cost:
	#	#if 
		all_costs.extend(agg_pairs_cost[x])
	getPairCostHist(("for","all"), all_costs)
	
	#getPairCostHist(cost_hist[('0','1')])
	#getPairCostHist("qs_0_2", [x[0] for x in agg_pairs_cost[('0','2')]])
	
	##Draw cost distribution for all pairs
	for x in agg_pairs_cost:
		getPairCostHist(x, agg_pairs_cost[x])
	
	##Print cost distribution value for all pairs
	#for x in agg_pairs_cost:
		#pair_cost = agg_pairs_cost[x]
		#min_val = min(pair_cost)
		#max_val = max(pair_cost)
		#stdev_val =  format(np.std(pair_cost),'.3f')
		#avg_val = format(np.mean(pair_cost),'.3f')
		#print x
	
	#getPairCostHist([x[0] for x in agg_pairs_cost[('0','1')]])
	#getCostHist(agg_cost) 
 
test() 
