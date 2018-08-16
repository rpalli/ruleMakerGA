# import python packages needed in this file
import pickle
from deap import base, creator, gp, tools
from deap import algorithms as algo
from random import random, seed, shuffle, randint
import networkx as nx
import numpy as numpy
import copy as copy
import operator
import argparse as argparse
import gc as gc

# import other pieces of our software
#import networkConstructor as nc
import simulation as sim
import GA as ga
from utils import genInitValueList, synthesizeInputs
# make empty list representing no knockouts or knockins
def setupEmptyKOKI(samples):
	knockoutLists=[]
	knockinLists=[]
	for q in range(samples):
		knockoutLists.append([])
		knockinLists.append([])
	return knockoutLists, knockinLists

# make a list of dictionaries giving values at each node from list of values across samples and a dictionary structure with random numbers
def genSampleList(output, sampleDict, samples):
	newSampleList=[]
	for k in range(0,samples):
		newSample=copy.deepcopy(sampleDict[k])
		for j in range(0,len(model.nodeList)):
			newSample[model.nodeList[j]]=output[k][j]
		newSampleList.append(newSample)
	return newSampleList


# run an experiment comparing 
def runExperiment(graph, name, samples, noise, edgeNoise, individual, params):
	#creates a model, runs simulations, then tests reverse engineering capabilities of models in a single function
	#samples is the number of different initial conditions to provide per trial
	#graph specifies the network we are testing. 
	# does everything except select the initial random bitstring and setup parameters

	sampleList=synthesizeInputs(graph,samples) # get empty list of inputs
	model=sim.modelClass(graph,sampleList, True) # generate empty model

	initModel=[(model.size), list(model.nodeList), list(model.individualParse), list(model.andNodeList) , list(model.andNodeInvertList), list(model.andLenList),	list(model.nodeList), dict(model.nodeDict), list(model.initValueList)]
	knockoutLists, knockinLists= setupEmptyKOKI(samples)
		
	# generate some simulated samples
	output=ga.runProbabilityBooleanSims(individual[1], model, samples, params.cells, params, knockoutLists, knockinLists)

	# add noise in omics data


	# add noise in RPKN
	if edgeNoise > 0:
		newgraph=graph.copy()
		edgelist=newgraph.edges()
		nodelist=newgraph.nodes()
		for newer in range(edgeNoise): # add edgeNoise FP edges
			rand1=randint(0,len(nodelist)-1)
			rand2=randint(0,len(nodelist)-1)
			edgeCandidate=(nodelist[rand1],nodelist[rand2])
			while edgeCandidate in edgelist or edgeCandidate[0]==edgeCandidate[1]:
				rand1=randint(0,len(nodelist)-1)
				rand2=randint(0,len(nodelist)-1)
				edgeCandidate=(nodelist[rand1],nodelist[rand2])
			if random()<.5:
				activity1='a'
			else:
				activity1='i'
			print(edgeCandidate)
			newgraph.add_edge(nodelist[rand1],nodelist[rand2], signal=activity1)
			edgelist.append((nodelist[rand1],nodelist[rand2]))

		print(edgelist)
		print(newgraph.edges())
	
	# output the initial generated data
	pickle.dump( output, open( name+"_input.pickle", "wb" ) )

	# copy simulated data into right format
	newSampleList=genSampleDict(output, sampleList, samples)
	testModel=sim.modelClass(graph,newSampleList, False)

	# put initial values into correct format, add to model
	newInitValueList=genInitValueList(newSampleList,testModel)
	testModel.initValueList=newInitValueList
	
	#find rules
	testModel, dev, bruteOut =ga.GAsearchModel(testModel, newSampleList, params, knockoutLists, knockinLists, name) # run GA
	bruteOut, equivalents = ga.localSearch(testModel, bruteOut, newSampleList, params, knockoutLists, knockinLists) # run local search
	storeModel3=[(testModel.size), list(testModel.nodeList), list(testModel.individualParse), list(testModel.andNodeList) , list(testModel.andNodeInvertList), list(testModel.andLenList),	list(testModel.nodeList), dict(testModel.nodeDict), list(testModel.initValueList)]

	outputList=[individual[1],bruteOut,storeModel, storeModel3, equivalents]
	pickle.dump( outputList, open( name+"_local1.pickle", "wb" ) )

def sampleTester(graph, name, samples):
	individual=ga.genBits(model) #generate random set of logic rules to start with
	params=sim.paramClass() # load in parameters
	runExperiment(graph, name, samples, 0., 0,individual, params)

def omicsNoiseTester(graph, name, noise):
	individual=ga.genBits(model) #generate random set of logic rules to start with
	params=sim.paramClass() # load in parameters
	runExperiment(graph, name, params.samples, noise,0, individual, params)

def RPKNnoiseTester(graph, name, noiseEdges):
	# runs experiment using graph and rule from Liu et al. 2016 along with additional false positive edges
	# loop over number of times we want to generate fake data and perform sequence of events
	params=sim.paramClass() # load in parameters
	# set up rule from Liu network
	individual=ga.genBits(model)
	individual[1]=	[0,1,1,1,1,1,0,0,1,0,1,1,1]
	print(individual[1])
	print(graph.edges())
	runExperiment(graph, name, params.samples, 0. , noiseEdges, individual, params)

def transformTest(graph,name,filename):
	# can't fit a rule to only one node
	if len(graph.nodes())<2:
		print('not enough overlap')
		return

	# load data, params, make empty knockout and knockin lists (no KO or KI in transform tests)
	sampleDict = constructBinInput(fileName)
	params=sim.paramClass()
	print(graph.nodes())
	knockoutLists, knockinLists= setupEmptyKOKI(samples)


	# generate turn sample dict into sample list (list of dicts instead of dict of lists)
	keyList=sampleDict.keys()
	sampleList=[{} for i in range(len(sampleDict[keyList[0]]))]
	for i in range(len(sampleList)):
		for key in keyList:
			if key in graph.nodes():
				sampleList[i][key]=sampleDict[key][i]
	
	# generate model encompassing graph and samples to do rule inference on
	model=sim.modelClass(graph,sampleList, False)

	# cpy data into correct order for simulation 
	newInitValueList=genInitValueList(sampleList,model)
	model.initValueList=newInitValueList
	print('setup successful')

	# find the rules
	model, dev1, bruteOut =ga.GAsearchModel(model, sampleList, params, knockoutLists, knockinLists, name)
	bruteOut, equivalents, dev2 = ga.localSearch(model, bruteOut, sampleList, params, knockoutLists, knockinLists)
	pickle.dump( [[dev1],[dev2],[bruteOut],[model]], open( name+"_output.pickle", "wb" ) )

def findPathways(geneDict):
	returnerlist=[]
	aliasDict={}
	dict1={}
	nc.parseKEGGdicthsa('inputData/hsa00001.keg',aliasDict,dict1)
	dict2={}
	nc.parseKEGGdict('inputData/ko00001.keg',aliasDict,dict2)
	namelist=pa.find_overlaps('filtered.c2.cp.kegg.v3.0.symbols.gmt',geneDict)
	print('num of overlap nodes')
	print(len(namelist))
	for name in namelist:
		returnerlist.append(pa.retrieveGraph(name,aliasDict,dict1,dict2, geneDict))
	return returnerlist

def constructBinInput(filename):
	binput = open(filename, 'r')
        ssDict = {}
	lines = binput.readlines()
	for line in lines:
		geneCount = line.split('\t')
		i = len(geneCount) 
		ssDict[geneCount[0]] = [float(q) for q in geneCount[1:i]]
	return(ssDict)

if __name__ == '__main__':

	# use the command line to input the file that you want
	#Is this right?
	#python experiments.py graphname(to be given) binData\kMeans.bin iterNum(to be given)

	parser = argparse.ArgumentParser()
	parser.add_argument("graph")
	parser.add_argument("datafile")
	parser.add_argument("iterNum")

	results = parser.parse_args()
	graphName=results.graph
	fileName=results.datafile
	iterNum=int(results.iterNum)
	
	name=graphName[:-8]+fileName+'_'+results.iterNum ##

	# insert your code to get a ssDict here.
	
	# use fileName to read in the discretized file... you want to standardize the format
	#ssDict= SOMETHING
	
	# edit the below name to include the pertinent piece of the filename of discretized data
	#outname=graphName[:-8]+'_'+results.iterNum
	outname=graphName[:-8]+fileName+'_'+results.iterNum
	print(name) #
	graph = nx.read_gpickle(graphName)
	print(len(graph.nodes()))
	transformTest(graph,outname,fileName)