#import python packages
import os as os
import pickle
import argparse as argparse
import numpy as np
import csv
import math
from random import randint
import scipy
import pandas as pd
import requests
import networkx as nx
import operator
# import other pieces of our software
import networkConstructor as nc
import utils as utils

class modelHolder:
	def __init__(self,valueList):
		[ size, nodeList, individualParse, andNodeList, andNodeInvertList, andLenList,nodeList, nodeDict, initValueList]=valueList
		self.size=size
		self.individualParse=individualParse
		self.andNodeList=andNodeList
		self.andNodeInvertList=andNodeInvertList
		self.andLenList=andLenList
		self.nodeList=nodeList
		self.nodeDict=nodeDict
		self.initValueList=initValueList

def outputGraphs(pathway, RAval, comparator, pathImportances):
	original=pathway[3].copy()

	
	nx.set_node_attributes(original,'Display Name',{k: k for k in original.nodes()})
	nx.set_node_attributes(original,'andNode',{k: 0 for k in original.nodes()})
	nx.set_node_attributes(original,'RA',{k: RAval[k] for k in original.nodes()})
	nx.set_node_attributes(original,'IS',{k: float(pathImportances[k]) for k in original.nodes()})

	ruleGraph=utils.Get_expanded_network(pathway[2][0].split('\n'),equal_sign='*=')
	nx.set_node_attributes(ruleGraph,'RA',{k: RAval[k] if k in original.nodes() else 0. for k in ruleGraph.nodes()})
	nx.set_node_attributes(ruleGraph,'IS',{k: float(pathImportances[k]) if k in original.nodes() else 0. for k in ruleGraph.nodes()})
	nx.write_graphml(original,comparator+'/'+pathway[0]+'.graphml')
	nx.write_graphml(ruleGraph,comparator+'/'+pathway[0]+'_rules.graphml')


def findPersistentRules(ruleList):
	finalRules=[]
	ruleSets=[ruleSet.split('\n') for ruleSet in ruleList]
	for i in range(len(ruleSets[0])):
		ruleDict={}
		for j in range(len(ruleSets)):
			if ruleSets[j][i] in ruleDict:
				ruleDict[ruleSets[j][i]]=ruleDict[ruleSets[j][i]]+1
			else:
				ruleDict[ruleSets[j][i]]=1
		finalRules.append(max(ruleDict.iteritems(), key=operator.itemgetter(1))[0])
	return finalRules
def retrievePathKey():
	pathDict={}

	requester='http://rest.kegg.jp/list/pathway'
	r=requests.get(requester)
	lines=r.text
	for line in lines.split('\n'):
		pieces=line.split('\t')
		if(len(pieces[0])>8):
			pathDict[str(pieces[0][8:])]=str(pieces[1])
		else:
			print(pieces)
	return pathDict

def makeRA(data,comparison,groups):
	RAdict={}
	group1=comparison[0]
	group2=comparison[1]
	for element in data:
		mean1=np.mean([data[element][temp] for temp in groups[group1]])
		mean2=np.mean([data[element][temp] for temp in groups[group2]])
		if mean1<=0 or mean2<=0:
			print(element)
			print(mean1)
			print(mean2)
			print(data[element])
		else:
			differ=abs(math.log(mean1,2)-math.log(mean2,2))
			RAdict[element]=differ
	return RAdict

def findPathwayList():

	pathways=[]
	codes=[]
	for file in os.listdir("gpickles"):
		if file.endswith(".gpickle"):
			if os.path.isfile('RD_out/'+file[:-8]+'_1_output.pickle'):
				codes.append(file[:-8])
			else:
				print(file+' has no output')
	print(codes)
	for code in codes:
		pathVals=[]
		rules=[]
		for i in range(1,6):
			[bruteOut1,dev,storeModel, storeModel3, equivalents]=pickle.Unpickler(open( 'RD_out/'+code+'_'+str(i)+'_local1.pickle', "rb" )).load()
			model=modelHolder(storeModel3)
			pathVals.append(pickle.Unpickler(open( 'RD_out/'+code+'_'+str(i)+'_scores.pickle', "rb" )).load())
			rules.append(utils.writeModel(bruteOut1, model))
		graph = nx.read_gpickle("gpickles/"+code+".gpickle")
		ImportanceVals={}
		for node in range(len(storeModel[1])):
			ImportanceVals[storeModel[1][node]]=math.log(np.mean([pathVals[i][node] for i in range(5)]),2)
		pathways.append([code,ImportanceVals, rules, graph])
	return pathways

def readFpkm(dataName,delmited):
	data=[]
	with open(dataName) as csvfile:
		data={}
		reader = csv.reader(csvfile, delimiter=delmited)
		firstline=reader.next()
		for row in reader:
			data[row[0]]=[float(row[k]) for k in range(1,len(row))]
	firstline.pop(0)
	
	# identify positions of each sample in the data
	colNums={}
	for item in range(len(firstline)):
		colNums[firstline[item]]=item
	return data, colNums

def readDiffs(diffName,delmited):
	# open the set of differences to be considered
	with open(diffName) as csvfile:
		diffs=[]
		reader = csv.reader(csvfile, delimiter=delmited)
		for row in reader:
			diffs.append(row)
	return diffs

def readMatrix(matrixName,delmited, colNums):
	# open and analyze matrix of group memberships
	with open(matrixName) as csvfile:
		matrix=[]
		reader = csv.reader(csvfile, delimiter=delmited)
		groupNames=reader.next()
		groupNames.pop(0)
		groups={}
		for name in groupNames:
			groups[name]=[]
		for row in reader:
			groupname=row.pop(0)
			for i in range(len(row)):
				if int(row[i])==1:
					groups[groupNames[i]].append(colNums[groupname])
	return groups

def analyze_pathways(diffName, matrixName, dataName, delmited):

	pathList=findPathwayList()
	data, colNums= readFpkm(dataName,delmited)	# read in fpkm data
	diffs= readDiffs(diffName,delmited) # read in difference to be considered
	groups=readMatrix(matrixName,delmited, colNums) # read design matrix

	# create an index of relative activities for all comparisons
	csvmaker=[]
	RAvals=[]
	comparisonStrings=[]
	pathDict=retrievePathKey()
	for comparison in diffs:
		comparisonStrings.append(comparison[0]+'-'+comparison[1])
		RAvals.append(makeRA(data,comparison,groups))
		if not os.path.exists(comparison[0]+'-'+comparison[1]):
			os.makedirs(comparison[0]+'-'+comparison[1])
	# iterate over pathways and calculate scores for pathway
	for pathway in pathList:
		print(pathway[0])
		for RAval, comparator in zip(RAvals, comparisonStrings):
			outputGraphs(pathway, RAval, comparator, pathway[1])
		z_scores=[]
		# iterate over comparisons for that pathway
		for RAval in RAvals:
			z_scores.append(scorePathway(RAval,pathway[1]))
		pvals=scipy.stats.norm.sf(map(abs,z_scores))
		tempdict={}
		tempdict['pathway']=pathDict[pathway[0][3:]]
		tempdict['code']=pathway[0][3:]
		for i in range(len(comparisonStrings)):
			tempdict[comparisonStrings[i]]=-math.log(pvals[i],10)
		csvmaker.append(tempdict)
	df=pd.DataFrame(csvmaker,columns=['pathway','code'].extend(comparisonStrings))
	df.to_csv(path_or_buf='pvalues.csv')


def scorePathway(RAs,pathImportances):
	score=0
	allNodes=RAs.keys()
	for node in pathImportances:
		score+=RAs[node]*pathImportances[node]
	print(np.mean([RAs[value] for value in allNodes]))
	randomScores=[]
	for i in range(1000):
		tempscore=0
		for node in pathImportances:
			tempscore+=RAs[allNodes[randint(0,len(allNodes)-1)]]*pathImportances[node]
		randomScores.append(tempscore)
	meaner=np.mean(randomScores)
	stdev=np.std(randomScores)
	zscore=(score-meaner)/stdev
	return zscore

if __name__ == '__main__':
	import time
	start_time = time.time()
	parser = argparse.ArgumentParser(prog='BONITA') 
	parser.set_defaults(sep=',')
	parser.add_argument("-sep", "--sep", metavar="seperator", help="How are columns in datafile specified")	
	parser.add_argument("-t", action='store_const',const='\t', dest="sep",help="Tab delimited?")	
	parser.add_argument("data")

	parser.add_argument("matrix")
	parser.add_argument("diffName")
	results = parser.parse_args()
	matrixName= results.matrix

	analyze_pathways(results.diffName, matrixName, results.data, results.sep)

	print("--- %s seconds ---" % (time.time() - start_time))

