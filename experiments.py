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
import utils as utils
import simulation as sim
# import murphyReader as mr
import motif_cutter as mc
# import pathway_analysis as pa
import GA as ga

def setupEmptyKOKI(samples):
	knockoutLists=[]
	knockinLists=[]
	for q in range(samples):
		knockoutLists.append([])
		knockinLists.append([])
	return knockoutLists, knockinLists

def genSampleList(output, sampleDict, samples):
	newSampleList=[]
	for k in range(0,samples):
		newSample=copy.deepcopy(sampleDict[k])
		for j in range(0,len(model.nodeList)):
			newSample[model.nodeList[j]]=output[k][j]
		newSampleList.append(newSample)
	return newSampleList

def genInitValueList(newSampleList,model)
	newInitValueList=[]
	for j in range(0,len(newSampleList)):
		newInitValueList.append([])
	for j in range(0,len(model.nodeList)):
		for k in range(0,len(newSampleList)):
			ss=newSampleList[k]
			if  model.nodeList[j] in newSampleList[0]:
				newInitValueList[k].append(ss[model.nodeList[j]])
			else:
				newInitValueList[k].append(0.5)
	return newInitValueList

def paramExperiment(graph, name, samples, noise):
	#creates a model, runs simulations, then tests reverse engineering capabilities of models in a single function
	#samples is the number of different initial conditions to provide per trial
	#graph specifies the network we are testing. 

	params=sim.paramClass() # load in parameters
	sampleList=utils.synthesizeInputs(graph,samples) # get empty list of inputs
	model=sim.modelClass(graph,sampleList, True) # generate empty model

	iteratorDict={} # dummy variable passed around
	storeModel=[(model.size), list(model.nodeList), list(model.individualParse), list(model.andNodeList) , list(model.andNodeInvertList), list(model.andLenList),	list(model.nodeList), dict(model.nodeDict), list(model.initValueList)]
	knockoutLists, knockinLists= setupEmptyKOKI(samples)
	
	#generate random set of logic rules to start with
	individual=ga.genBits(model)
	truthList.append(list(individual[1]))
	
	# generate some simulated samples
	output=ga.runProbabilityBooleanSims(individual[1], model, samples, params.cells, params, knockoutLists, knockinLists, iteratorDict)

	# output the initial generated data
	pickle.dump( output, open( name+"_input.pickle", "wb" ) )

	# copy simulated data into right format
	newSampleList=genSampleDict(output, sampleList, samples)
	model=[]
	model=sim.modelClass(graph,newSSS, False)


	newInitValueList=
	model.initValueList=newInitValueList
	#perform GA run
	params.adaptive=True
	model, dev, bruteOut =ga.GAsearchModel(model, newSSS, params, knockoutLists, knockinLists, iteratorDict, name)
	bruteOut, equivalents = ga.localSearch(model, bruteOut, newSSS, params, knockoutLists, knockinLists)
	
	storeModel3=[(model.size), list(model.nodeList), list(model.individualParse), list(model.andNodeList) , list(model.andNodeInvertList), list(model.andLenList),	list(model.nodeList), dict(model.nodeDict), list(model.initValueList)]

	outputList=[individual[1],bruteOut,storeModel, storeModel3, equivalents]
	pickle.dump( outputList, open( name+"_local1.pickle", "wb" ) )

def sampleTester(graph, name, samples):
	#creates a model, runs simulations, then tests reverse engineering capabilities of models in a single function
	#trials is the number of different individuals to try
	#samples is the number of different initial conditions to provide per trial
	#graph specifies the network we are testing. 

	# set up params, sss, model, samples

	params=sim.paramClass()
	sss=utils.synthesizeInputs(graph,samples)
	model=sim.modelClass(graph,sss, True)

	# set up necessary lists for output
	truthList= [] 
	testListAdapt =[]
	devListAdapt=[]
	iteratorDict={}
	storeModel=[(model.size), list(model.nodeList), list(model.individualParse), list(model.andNodeList) , list(model.andNodeInvertList), list(model.andLenList),	list(model.nodeList), dict(model.nodeDict), list(model.initValueList)]
	knockoutLists, knockinLists= setupEmptyKOKI(samples)
	#generate random set of logic rules to start with
	individual=ga.genBits(model)
	truthList.append(list(individual[1]))
	# generate Boolean model for this trial
	output=ga.runProbabilityBooleanSims(individual[1], model, samples, params.cells, params, knockoutLists, knockinLists, iteratorDict)

	truthList.append(list(individual[1]))
	# copy new output into newSSS and initial values
	pickle.dump( output, open( name+"_input.pickle", "wb" ) )

	newSSS=[]
	for k in range(0,samples):
		newSS=copy.deepcopy(sss[k])
		for j in range(0,len(model.nodeList)):
			newSS[model.nodeList[j]]=output[k][j]
		newSSS.append(newSS)
	# print(newSSS)
	model=[]
	model=sim.modelClass(graph,newSSS, False)

	newInitValueList=[]
	for j in range(0,len(sss)):
		newInitValueList.append([])
	for j in range(0,len(model.nodeList)):
		for k in range(0,len(sss)):
			ss=newSSS[k]
			if  model.nodeList[j] in sss[0]:
				newInitValueList[k].append(ss[model.nodeList[j]])
			else:
				newInitValueList[k].append(0.5)
	model.initValueList=newInitValueList
	#perform GA run
	params.adaptive=True
	model, dev, bruteOut =ga.GAsearchModel(model, newSSS, params, knockoutLists, knockinLists, iteratorDict, name)
	bruteOut, equivalents = ga.localSearch(model, bruteOut, newSSS, params, knockoutLists, knockinLists)
	
	storeModel3=[(model.size), list(model.nodeList), list(model.individualParse), list(model.andNodeList) , list(model.andNodeInvertList), list(model.andLenList),	list(model.nodeList), dict(model.nodeDict), list(model.initValueList)]

	outputList=[[individual[1]],[bruteOut],storeModel, storeModel3, equivalents]
	pickle.dump( outputList, open( name+"_local1.pickle", "wb" ) )

	model, dev, bruteOut = ga.GASearchModel2(model, newSSS, params, knockoutLists, knockinLists, bruteOut, name)
	bruteOut, equivs2 = ga.localSearch(model, bruteOut, newSSS, params, knockoutLists, knockinLists)
	testListAdapt.append(list(bruteOut))
	devListAdapt.append(dev)
	iteratorDict={}
	# params.adaptive=False
	gc.collect()

	# # set up output and save as a pickle
	storeModel2=[(model.size), list(model.nodeList), list(model.individualParse), list(model.andNodeList) , list(model.andNodeInvertList), list(model.andLenList),	list(model.nodeList), dict(model.nodeDict), list(model.initValueList)]

	outputList=[truthList,testListAdapt,storeModel, storeModel2, equivs2]
	pickle.dump( outputList, open( name+"_output.pickle", "wb" ) )

def liuTester(graph, name, noiseEdges):
	#creates a model, runs simulations, then tests reverse engineering capabilities of models in a single function
	#trials is the number of different individuals to try
	#samples is the number of different initial conditions to provide per trial
	#graph specifies the network we are testing. 

	# set up params, sss, model, samples

	params=sim.paramClass()
	# rewire graph if necessary
	#if params.rewire:
	#graph=nc.rewireNetwork(graph)

	sss=utils.synthesizeInputs(graph,params.samples)
	model=sim.modelClass(graph,sss, True)
	samples=params.samples

	# set up necessary lists for output
	truthList= [] 
	testListGA =[]
	devListGA=[]
	testListAdapt =[]
	devListAdapt=[]
	iteratorDict={}
	knockoutLists=[]
	knockinLists=[]
	testListAdaptglobal=[]
	devListAdaptglobal=[]
	testListGAglobal=[]
	devListGAglobal=[]
	testListLocal1=[]
	testListLocal2=[]

	storeModel=[(model.size), list(model.nodeList), list(model.individualParse), list(model.andNodeList) , list(model.andNodeInvertList), list(model.andLenList),	list(model.nodeList), dict(model.nodeDict), list(model.initValueList)]
	for q in range(samples):
		temRand=randint(0,len(model.nodeList))
		knockoutLists.append([])
		knockinLists.append([])

	# loop over number of times we want to generate fake data and perform sequence of events
	for i in range(0,params.trials):
		iteratorDict={}
		#generate random set of logic rules to start with
		individual=ga.genBits(model)
		individual[1]=	[0,1,1,1,1,1,0,0,1,0,1,1,1]
		print(individual[0].nodeList)
		print(individual[0].individualParse)
		print(individual[0].andNodeList)
		print(individual[1])
		print(graph.edges())
		newgraph=graph.copy()
		edgelist=newgraph.edges()
		nodelist=newgraph.nodes()
		for newer in range(noiseEdges):
			edger=newgraph.edges()
			edger=edger[0]
			while edger in edgelist or edger[0]==edger[1]:
				rand1=randint(0,len(nodelist)-1)
				rand2=randint(0,len(nodelist)-1)
				edger=(nodelist[rand1],nodelist[rand2])
			if random()<.5:
				activity1='a'
			else:
				activity1='i'
			print(edger)

			newgraph.add_edge(nodelist[rand1],nodelist[rand2], signal=activity1)
			edgelist.append((nodelist[rand1],nodelist[rand2]))

		print(edgelist)
		print(newgraph.edges())


		truthList.append(list(individual[1]))
		# generate Boolean model for this trial
		output=ga.runProbabilityBooleanSims(individual[1], model, samples, params.cells, params, knockoutLists, knockinLists, iteratorDict)
		# copy new output into newSSS and initial values
		pickle.dump( output, open( name+"_input.pickle", "wb" ) )

		newSSS=[]
		for k in range(0,samples):
			newSS=copy.deepcopy(sss[k])
			for j in range(0,len(model.nodeList)):
				newSS[model.nodeList[j]]=output[k][j]
			newSSS.append(newSS)
		# print(newSSS)
		model=[]
		model=sim.modelClass(newgraph,newSSS, False)

		newInitValueList=[]
		for j in range(0,len(sss)):
			newInitValueList.append([])
		for j in range(0,len(model.nodeList)):
			for k in range(0,len(sss)):
				ss=newSSS[k]
				if  model.nodeList[j] in sss[0]:
					newInitValueList[k].append(ss[model.nodeList[j]])
				else:
					newInitValueList[k].append(0.5)
		model.initValueList=newInitValueList
		#perform GA run
		params.adaptive=True
		model, dev, bruteOut =ga.GAsearchModel(model, newSSS, params, knockoutLists, knockinLists, iteratorDict, name)
		bruteOut, equivalents = ga.localSearch(model, bruteOut, newSSS, params, knockoutLists, knockinLists)
		
		storeModel3=[(model.size), list(model.nodeList), list(model.individualParse), list(model.andNodeList) , list(model.andNodeInvertList), list(model.andLenList),	list(model.nodeList), dict(model.nodeDict), list(model.initValueList)]

		outputList=[[individual[1]],[bruteOut],storeModel, storeModel3, equivalents]
		pickle.dump( outputList, open( name+"_local1.pickle", "wb" ) )

		# model, dev, bruteOut =ga.GAsearchModel2(model, newSSS, params, knockoutLists, knockinLists, iteratorDict, name)
		# dev=ga.evaluateByNode(bruteOut, 1000, model,  newSSS, params, knockoutLists, knockinLists,iteratorDict)
		# for badnode in range(len(dev)):
		# 	if (dev[badnode]>.2) and (len(model.possibilityList[badnode])>3):
		# 		truth, model=ga.checkUpstreamPossibilities(badnode, bruteOut, newSSS, 1000, model,params, knockoutLists, knockinLists,iteratorDict )
		model, dev, bruteOut = ga.GASearchModel2(model, newSSS, params, knockoutLists, knockinLists, bruteOut, name)
		bruteOut, equivs2 = ga.localSearch(model, bruteOut, newSSS, params, knockoutLists, knockinLists)
		testListAdapt.append(list(bruteOut))
		devListAdapt.append(dev)
		iteratorDict={}
		# params.adaptive=False
		gc.collect()

	# # set up output and save as a pickle
	storeModel2=[(model.size), list(model.nodeList), list(model.individualParse), list(model.andNodeList) , list(model.andNodeInvertList), list(model.andLenList),	list(model.nodeList), dict(model.nodeDict), list(model.initValueList)]

	outputList=[truthList,testListAdapt,storeModel, storeModel2, equivs2]
	pickle.dump( outputList, open( name+"_output.pickle", "wb" ) )

def partitionTester(graph,name,ssDict):
	# first trim the graph to just what is in the dataset and the dataset to just what is in the graph

	if len(graph.nodes())>1:	
		partitionTest(graph, ssDict, name) # divide by max
	else:
		print('not enough overlap')

def partitionTest(graph,  ssDict, name):

	# ssDict, valueLister, ones, zeros= mr.constructOmicsInput(geneDict, maxRem, divider)
	params=sim.paramClass()
	
	#graph=nc.simplifyNetwork(graph, ssDict) # simplify the graph
	#graph=mc.cutMotifs(graph)
	print(graph.nodes())
	keyList=ssDict.keys()
	# generate SSS which is a list of SSs where SS is a dict pointing to some values
	sss=[{} for i in range(len(ssDict[keyList[0]]))]
	newInitValueList=[[] for i in range(len(ssDict[keyList[0]]))]
	for i in range(len(sss)):
		for key in keyList:
			if key in graph.nodes():
				sss[i][key]=ssDict[key][i]
	
	model=sim.modelClass(graph,sss, False)
	# generate empty knockout and knockin lists for support
	knockoutLists=[]
	knockinLists=[]
	for q in range(len(sss)):
		knockoutLists.append([])
		knockinLists.append([])

	# generate a set of initial values which are just the values from the SS in order. 
	newInitValueList=[]
	for j in range(0,len(sss)):
		newInitValueList.append([])
	for j in range(0,len(model.nodeList)):
		for k in range(0,len(sss)):
			ss=sss[k]
			if  model.nodeList[j] in sss[0]:
				newInitValueList[k].append(ss[model.nodeList[j]])
			else:
				newInitValueList[k].append(0.5)
	model.initValueList=newInitValueList
	print('setup successful')
	#perform GA run
	# dev, bruteOut= ga.GAsearchModel(model, sss, params, knockoutLists, knockinLists)
	iteratorDict={}
	model, dev1, bruteOut =ga.GAsearchModel(model, sss, params, knockoutLists, knockinLists, iteratorDict, name)
	bruteOut, equivalents, dev2 = ga.localSearch(model, bruteOut, sss, params, knockoutLists, knockinLists)
	# storeModel3=[(model.size), list(model.nodeList), list(model.individualParse), list(model.andNodeList) , list(model.andNodeInvertList), list(model.andLenList),	list(model.nodeList), dict(model.nodeDict), list(model.initValueList)]
	# outputList=[[individual[1]],[bruteOut],storeModel, storeModel3, equivalents]
	# pickle.dump( outputList, open( name+"_local1.pickle", "wb" ) )
	#model, dev, bruteOut = ga.GASearchModel2(model, sss, params, knockoutLists, knockinLists, bruteOut, name)
	#bruteOut, equivs2 = ga.localSearch(model, bruteOut, sss, params, knockoutLists, knockinLists)

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
		# ssDict[geneCount[0]] = [ float(q) for q in geneCount[1:-1]]
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
	ssDict = constructBinInput(fileName)
	
	# edit the below name to include the pertinent piece of the filename of discretized data
	#outname=graphName[:-8]+'_'+results.iterNum
	outname=graphName[:-8]+fileName+'_'+results.iterNum
	print(name) #
	graph = nx.read_gpickle(graphName)
	print(len(graph.nodes()))
	partitionTester(graph,outname,ssDict)