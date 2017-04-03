import pickle
from deap import base
from deap import creator
from deap import gp
from deap import tools
from deap import algorithms as algo
from random import random
from random import seed
import utils as utils
import simulation as sim
import networkx as nx
import networkConstructor as nc
import numpy as numpy
import copy as copy
import operator
import matplotlib.pyplot as plt
import liu_networks as liu
from random import shuffle
import math as math
from sets import Set
from joblib import Parallel, delayed


class probInitSeqClass:
	def __init__(self):
		self.startNodes=[]
		self.startProbs=[]
		self.nodeOrder=[]

def buildToolbox( individualLength, bitFlipProb):
	# # #setup toolbox
	toolbox = base.Toolbox()
	pset = gp.PrimitiveSet("MAIN", arity=1)
	pset.addPrimitive(operator.add, 2)
	pset.addPrimitive(operator.sub, 2)
	pset.addPrimitive(operator.mul, 2)
	weightTup=(-1.0,)
	creator.create("FitnessMin", base.Fitness, weights=weightTup) # make a fitness minimization function
	creator.create("Individual", list, fitness=creator.FitnessMin)	# create a class of individuals that are lists

	#register our bitsring generator and how to create an individual, population
	toolbox.register("genRandomBitString", utils.genRandBits, individualLength=individualLength)
	toolbox.register("individual", tools.initIterate, creator.Individual, toolbox.genRandomBitString)
	toolbox.register("population", tools.initRepeat, list , toolbox.individual)
	
	#create statistics toolbox and give it functions
	stats = tools.Statistics(key=lambda ind: ind.fitness.values)
	stats.register("avg", numpy.mean)
	stats.register("std", numpy.std)
	stats.register("min", numpy.min)
	stats.register("max", numpy.max)
	
	# finish registering the toolbox functions
	toolbox.register("mate", tools.cxTwoPoint)
	toolbox.register("mutate", tools.mutFlipBit, indpb=bitFlipProb)
	toolbox.register("select", tools.selNSGA2)
	toolbox.register("similar", numpy.array_equal)
	
	return toolbox, stats
			
def evaluate(individual, params, model, simulator, sss):
	SSEs=[]
	for j in range(0,len(sss)):
		ss=sss[j]
		initValues=model.initValueList[j]
		SSE=0
		boolValues, addnodeNums=sim.runModel(individual, model,simulator, model.initValueList[j])	
		for i in range(0, len(model.evaluateNodes)):
			SSE+=(boolValues[model.evaluateNodes[i]]-ss[model.nodeList[model.evaluateNodes[i]]])**2
		SSEs.append(SSE)
	summer=0
	for i in range(0,len(SSEs)):
		summer+=SSEs[i]
	summer=summer/len(SSEs)
	likelihood=1-summer/len(model.andLenList)
	for i in range(len(model.nodeList)):
		if i==len(model.nodeList)-1:
			end= model.size
		else:
			end=model.individualParse[i+1]	 
		if sum(individual[model.individualParse[i]:model.individualParse[i+1]])==0:
			likelihood=.001	
	if params.IC==1:
		return addnodeNums-math.log(likelihood),
	elif params.IC==2:
		return addnodeNums*math.log(len(sss))-2*math.log(likelihood),
	else:
		return summer,
	

# generates a random set of samples made up of cells by using parameteris from probInit seq
# to set up then iterating using strict Boolean modeling. 
def runProbabilityBooleanSims(individual, model, sampleNum, cells):
	samples=[]
	seeds=[]
	for i in range(0,sampleNum):
		seeds.append(random)
	samples=Parallel(n_jobs=min(6,sampleNum))(delayed(sampler)(individual, model, sampleNum, seeds[i]) for i in range(sampleNum))
	# counter=0
	# for sample in range(len(samples)):
	# 	if sum(samples(sample))==0:
	# 		temp=samples.pop(sample) 
	# newSamples=Parallel(n_jobs=min(8,sampleNum))(delayed(sampler)(individual, model, sampleNum, seeds[i]) for i in range(sampleNum))
	# print(samples)
	return samples
	

def sampler(individual, model, cells, seeder):
	seed(seeder)
	cellArray=[]
	sampleProbs=[]
	simulator=sim.simulatorClass('bool')
	params=sim.paramClass()
	for j in range(0,len(model.nodeList)):
		sampleProbs.append(random())
	#samples= Parallel(n_jobs=8)(delayed(sampler)(individual, model, cells, randos[i]) for i in range(sampleNum))
	for j in range(0,cells):
		# shuffle nodes to be initially called.... 
		#simulations that are truly random starting states should have all nodes added to this list
		#get initial values for all nodes
		initValues=genPBNInitValues(individual, model,sampleProbs)
		# run Boolean simulation with initial values and append
		vals, nums=sim.runModel(individual, model, simulator, initValues)
		cellArray.append(vals)
	return [sum(col) / float(cells) for col in zip(*cellArray)]

def genPBNInitValues(individual, model,sampleProbs):
	initValues=[False for x in range(0,len(model.nodeList))]
	for node in range(0,len(sampleProbs)):
		randNum=random()
		if randNum<sampleProbs[node]:
			initValues[node]=True
	return initValues

def updateInitSeq(individual, model, probInitSeq):
	#looks through list of nodes and adds to list of initial updated nodes if there are no specified inputs
	initSeq=probInitSeqClass()
	initSeq.startNodes=list(probInitSeq.startNodes)
	initSeq.startProbs=list(probInitSeq.startProbs)
	initSeq.nodeOrder=list(probInitSeq.nodeOrder)
	for i in range(0,len(model.nodeList)):
		parser=model.individualParse[i]
		andList=model.andNodeList[i] # find the list of possible input combinations for the node we are on 
		
		if (model.andLenList[i]==0 or (len(andList)<2 and individual[0]==0)) and i not in initSeq.startNodes:
			# if there are no inputs or the rule excludes a single input then add to initial node list if not already there
			initSeq.startNodes.append(i)
			initSeq.nodeOrder.pop(i)
			initSeq.startProbs.append(1)		
	return initSeq

def GAautoSolver(model, sss, propSimulator):
	params=sim.paramClass()
	toolbox, stats=buildToolbox(model.size,params.bitFlipProb)
	# reset simSteps to 1 so we just see first step in Sim...
	propSimulator.simSteps=1
	toolbox.register("evaluate", evaluate, params=params,model=model,simulator=propSimulator,sss=sss)
	population=toolbox.population(n=params.popSize)
	hof = tools.HallOfFame(params.hofSize, similar=numpy.array_equal)
	output=algo.eaMuCommaLambda(population, toolbox, mu=params.mu, lambda_=params.lambd, stats=stats, cxpb=params.crossoverProb, mutpb=params.mutationProb, ngen=params.generations, verbose=False, halloffame=hof)
	# stats = tools.Statistics(key=lambda ind: ind.fitness.values)
	# stats.register("avg", numpy.mean)
	# stats.register("std", numpy.std)
	# stats.register("min", numpy.min)
	# stats.register("max", numpy.max)
	# toolbox.register("evaluate", evaluate, params=params,model=model,simulator=propSimulator,sss=newSSS)
	return output, hof	

def evalNode(individual, currentNode, params, model, simulator, sss):
	SSEs=[]
	for j in range(0,len(sss)):
		ss=sss[j]
		initValue=model.initValueList[j]
		SSE=0
		value= sim.updateNode(currentNode,initValue,individual,  model,simulator)
		SSE+=(value-ss[model.nodeList[model.evaluateNodes[currentNode]]])**2
		SSEs.append(SSE)
	summer=0
	for i in range(0,len(SSEs)):
		summer+=SSEs[i]
	summer=summer/len(SSEs)
	likelihood=1-summer
	for i in range(len(model.nodeList)):
		if i==len(model.nodeList)-1:
			end= model.size
		else:
			end=model.individualParse[i+1]	 
		if sum(individual[model.individualParse[i]:model.individualParse[i+1]])==0:
			likelihood=.001	
	totalNodes=0
	for bit in range(len(individual)):
		if individual[bit]==1:
			totalNodes=totalNodes+len(model.andNodeInvertList[currentNode][bit])
	if params.IC==1:
		return totalNodes-math.log(likelihood),
	elif params.IC==2:
		return totalNodes*math.log(len(sss))-2*math.log(likelihood),
	elif params.IC==3:
		return totalNodes-4*math.log(likelihood),
	else:
		return summer,

def piecewiseGASolver(model, sss, propSimulator):
	params=sim.paramClass()
	# reset simSteps to 1 so we just see first step in Sim...
	propSimulator.simSteps=1
	bestList=[]
	for i in range(len(model.nodeList)):
		if model.andLenList[i]>1:
			toolbox, stats=buildToolbox(model.andLenList[i],params.bitFlipProb)
			toolbox.register("evaluate", evalNode,currentNode=i, params=params,model=model,simulator=propSimulator,sss=sss)
			population=toolbox.population(n=params.popSize)
			hof = tools.HallOfFame(params.hofSize, similar=numpy.array_equal)
			output=algo.eaMuCommaLambda(population, toolbox, mu=params.mu, lambda_=params.lambd, stats=stats, cxpb=params.crossoverProb, mutpb=params.mutationProb, ngen=params.generations, verbose=False, halloffame=hof)
			bestList.append(hof[0])
		elif  model.andLenList[i]==1:
			bestList.append([1])
		else:
			bestList.append([])
	return [item for sublist in bestList for item in sublist]

#exhaustively search boolean networks for best option going node by node for rules
def bruteForceSearchModel(model, sss1, simulator):
	initValueList=[]
	for j in range(0,len(sss1)):
		initValueList.append([])
	for i in range(0,len(model.nodeList)):
		for j in range(0,len(sss1)):
			ss=sss1[j]
			if  model.nodeList[i] in sss1[0].keys():
				initValueList[j].append(ss[model.nodeList[i]])
			else:
				initValueList[j].append(0.5)
	#print(model.initValueList)
	model.initValueList=initValueList
	bestlist=Parallel(n_jobs=6)(delayed(singleNodeBF)(model, sss1, simulator, i) for i in range(len(model.nodeList)))
	return [item for sublist in bestlist for item in sublist]

def singleNodeBF(model, sss1, simulator, i):
	best=[]
	if model.andLenList[i]>0:
		currentDev=10000*len(sss1)
		best=utils.bitList(0,model.andLenList[i])
		if model.andLenList[i]<11:
			checkRange=2**(model.andLenList[i])
		else:
			checkRange=2**(11)
		for j in range(1,checkRange):
			bits=[]
			bits=utils.bitList(j,model.andLenList[i] )
			deviation=0
			for steadyStateNum in range(0,len(model.initValueList)):
				derivedVal=sim.updateNode(i,model.initValueList[steadyStateNum],bits, model, simulator)
				deviation=deviation+(derivedVal-sss1[steadyStateNum][model.nodeList[i]])**2
			if(deviation<currentDev):
				best=bits
				currentDev=deviation
	return best


def compareIndividualsNodeWise(truthList, testList, model):
	nodesensitivity=[]
	nodespecificity=[]
	netOnes=[]
	netZeros=[]
	netNegOnes=[]
	for node in range(len(model.nodeList)):
		ones=[]
		zeros=[]
		negones=[]
		# get indices to extract individual for this node
		if node==len(model.nodeList)-1:
			end=len(model.nodeList)
		else:
			end=model.individualParse[node+1]
		start=model.individualParse[node]
		sumindividual=[]
		#loop over individuals provided and calculate relevant values
		for i in range(len(truthList)):
			truth= truthList[i][start:end]
			test= testList[i][start:end]
			sumindividual.append(sum(truth))
			newindividual=[a_i - b_i for a_i, b_i in zip(truth, test)]
			ones.append(newindividual.count(1))
			zeros.append(newindividual.count(0))
			negones.append(newindividual.count(-1))
		# append node-wise breakdowns to list of breakdowns for the model as a whole
		netOnes.append(sum(ones))
		netZeros.append(sum(zeros))
		netNegOnes.append(sum(negones))
		
		# calculate sensitivity and specificity for the node
		temp=[100 if zeros[i]==0 else (1.*zeros[i])/(zeros[i]+ones[i]) for i in range(0,len(ones))]
		temp=filter(lambda a: a != 100, temp)
		if len(temp)==0:
			sensitivity=100
		else:
			sensitivity=(sum(temp)/len(temp))
		temp=[100 if (len(newindividual)-sumindividual[i])==0 else (1.*len(newindividual)-sumindividual[i]-negones[i])/(len(newindividual)-sumindividual[i]) for i in range(0,len(ones))]
		temp=filter(lambda a: a != 100, temp)
		if len(temp)==0:
			specificity=100
		else:
			specificity=(sum(temp)/len(temp))


		# add to list of sensitivity and specificity by node
		nodesensitivity.append(sensitivity)
		nodespecificity.append(specificity)
	#calculate sensitivity and specificity on the network as a whole
	ones=[]
	zeros=[]
	negones=[]
	for i in range(len(truthList)):
		truth= truthList[i]
		test= testList[i]
		sumindividual.append(sum(truth))
		newindividual=[a_i - b_i for a_i, b_i in zip(truth, test)]
		ones.append(newindividual.count(1))
		zeros.append(newindividual.count(0))
		negones.append(newindividual.count(-1))
	temp=[100 if zeros[i]==0 else (1.*zeros[i])/(zeros[i]+ones[i]) for i in range(0,len(ones))]
	temp=filter(lambda a: a != 100, temp)
	if len(temp)==0:
		sensitivity=100
	else:
		sensitivity=(sum(temp)/len(temp))
	temp=[100 if (len(newindividual)-sumindividual[i])==0 else (1.*len(newindividual)-sumindividual[i]-negones[i])/(len(newindividual)-sumindividual[i]) for i in range(0,len(ones))]
	temp=filter(lambda a: a != 100, temp)
	if len(temp)==0:
		specificity=100
	else:
		specificity=(sum(temp)/len(temp))
	return sensitivity, specificity, nodesensitivity, nodespecificity



def genBits(model):
	startInd=utils.genRandBits(model.size)
	for node in range(0,len(model.nodeList)):
		if node==(len(model.nodeList)-1):
			end=len(startInd)-1
		else:
			end=model.individualParse[node+1]
		start=model.individualParse[node]
		truth=startInd[start:end]
		if len(truth)>1:
			for i in range(len(truth)):
				if random()<(1./len(truth)):
					truth[i]=1
				else:
					truth[i]=0
			counter=0
			while sum(truth)>5 and counter < 100000:
				indices = [i for i in range(len(truth)) if truth[i] == 1]
				chosen=math.floor(random()*len(indices))
				truth[indices[int(chosen)]]=0
				counter+=1
			startInd[start:end]=truth
		elif len(truth)==1:
			truth[0]=1
			startInd[start:end]=truth
	return startInd
		
def simTester(model, sss, simClass):
	#creates a model, runs simulations, then tests reverse engineering capabilities of models in a single function
	#trials is the number of different individuals to try
	#samples is the number of different initial conditions to provide per trial
	#graph specifies the network we are testing. 

	# set up empty lists and call parameters
	params=sim.paramClass()
	samples=params.samples
	trials=params.trials
	cells=params.cells
	truthIndividuals=[]
	truthAvgs=[]
	hofScores=[]
	newSSS=[]
	trueNodeList=[]
	zeros=[]
	ones=[]
	negones=[]
	truthCounter=0
	sumindividual=[]
	truthList= [] 
	testList =[]
	# loop over number of times we want to generate fake data and perform sequence of events
	for i in range(0,trials):

		#generate random set of logic rules to start with
		individual=genBits(model)
		for node in range(0,len(model.nodeList)):
			if model.andLenList[node]>0:
				if node==len(model.nodeList)-1:
					end=len(model.nodeList)
				else:
					end=model.individualParse[node+1]
				if sum(individual[model.individualParse[node]:end])==0:
					individual[model.individualParse[node]]=1

		# generate Boolean model for this trial
		output=runProbabilityBooleanSims(individual, model, samples, cells)
		# copy new output into newSSS and initial values
		newSSS=[]
		for k in range(0,samples):
			newSS=copy.deepcopy(sss[k])
			for j in range(0,len(model.nodeList)):
				newSS[model.nodeList[j]]=output[k][j]
			newSSS.append(newSS)
		newInitValueList=[]
		for j in range(0,len(sss)):
			newInitValueList.append([])
		# print(len(newSSS))
		# print(samples)
		# print(len(sss))
		for j in range(0,len(model.nodeList)):
			for k in range(0,len(sss)):
				ss=newSSS[k]
				if  model.nodeList[j] in sss[0].keys():
					newInitValueList[k].append(ss[model.nodeList[j]])
				else:
					newInitValueList[k].append(0.5)
		model.initValueList=newInitValueList


		# set up PBN-based simulator
		propSimulator=sim.simulatorClass(simClass)
		propSimulator.trainingData=model.initValueList
		propSimulator.train(model)
		#perform brute force search
		bruteOut=bruteForceSearchModel(model, newSSS,propSimulator)
		truthList.append(individual)
		testList.append(bruteOut)
	tuple1=compareIndividualsNodeWise(truthList, testList, model)
	sensitivities=[]
	specificities=[]
	inEdgeNums=[]
	# run correction that simplifies truth rules
	TPsum=0
	TNsum=0
	FNsum=0
	FPsum=0
	newtruths=[]
	for k in range(len(truthList)):
		newtruths.append([])
	for node in range(0,len(model.nodeList)):
		FP=0
		TP=0
		TN=0
		FN=0
		if node==(len(model.nodeList)-1):
			end=len(truthList[0])-1
		else:
			end=model.individualParse[node+1]
		start=model.individualParse[node]
		andNodeList=model.andNodeList[node]
		inEdges=[]
		for lister in andNodeList:
			inEdges.append(set(lister))
		for k in range(len(truthList)):
			truth=truthList[k][start:end]
			test= testList[k][start:end]
			for i in range(len(truth)):
				if truth[i]==1:
					for j in range(len(truth)):
						if truth[j]==1:
							if inEdges[i].issubset(inEdges[j]):
								truth[j]=0
			newtruths[k].extend(truth)
		for k in range(len(truthList)):
			truth=newtruths[k][start:end]
			test= testList[k][start:end]
			truthSet=Set([])
			testSet=Set([])
			baseSet=Set([])
			for i in range(0,len(truth)):
				if truth[i]==1:
					for nodeToAdd in andNodeList[i]:
						truthSet.add(nodeToAdd)
			for i in range(0,len(test)):
				if test[i]==1:
					for nodeToAdd in andNodeList[i]:
						testSet.add(nodeToAdd)
				for nodeToAdd in andNodeList[i]:
					baseSet.add(nodeToAdd)
			FP+=1.*len(testSet.difference(truthSet))
			TP+=1.*len(testSet.intersection(truthSet))
			TN+=1.*(len(baseSet)-len(truthSet))
			FN+=1.*len(truthSet.difference(testSet))
		if len(truthSet)>0:
			sensitivity=TP/(len(truthSet))
		else:
			sensitivity=100
		if TN+FP>0:
			specificity=TN/(TN+FP)
		else:
			specificity=100
		sensitivities.append(sensitivity)
		specificities.append(specificity)
		TPsum+=TP
		TNsum+=TN
		FNsum+=FN
		FPsum+=FP
		inEdgeNums.append(len(baseSet))
	tuple2=compareIndividualsNodeWise(newtruths, testList, model)
	if (TPsum+FNsum)>0:
		sensitivity=TPsum/(TPsum+FNsum)
	else:
		sensitivity=100
	if (FPsum+TNsum)>0:
		specificity= TNsum/(FPsum+TNsum)
	else:
		specificity=100
	tuple3= (sensitivity, specificity, sensitivities, specificities)


	return tuple1, tuple2, tuple3, inEdgeNums,  [testList, truthList,  newtruths]

def ifngStimTest(bioReplicates):
	aliasDict={}
	dict1={}
	nc.parseKEGGdicthsa('hsa00001.keg',aliasDict,dict1)
	dict2={}
	nc.parseKEGGdict('ko00001.keg',aliasDict,dict2)
	

	# read in list of codes then load them into network
	inputfile = open('ID_filtered.c2.cp.kegg.v3.0.symbols.txt', 'r')
	lines = inputfile.readlines()
	data=dict(utils.loadFpkms('Hela-C-1.count'))
	sensitivities=[]
	specificities=[]
	lines.pop(0)
	nodesensitivities=[[],[],[]]
	nodespecificities=[[],[],[]]
	truthholder=[]
	edgeDegree=[]
	lines.pop(0)
	for code in lines:
		graph=nx.DiGraph()
		coder=str('ko'+code[:-1])
		nc.uploadKEGGcodes([coder], graph, dict2)
		coder=str('hsa'+code[:-1])
		nc.uploadKEGGcodes_hsa([coder], graph,dict1, dict2)
		if(len(graph.edges())>1):
			graph=nc.simplifyNetwork(graph, data)
		
		if(len(graph.edges())>1):
			print(coder)
			print(len(graph.edges()))
			
			nx.write_graphml(graph,coder+'.graphml')
			tempsensitivities=[[],[],[]]
			tempspecificities=[[],[],[]]
			truthlists=[[],[],[]]
			for i in range(bioReplicates):
				tuple1, tuple2, tuple3, inEdgeNums, truthlisttemp=rewireSimTest(graph)
				sensitivity1, specificity1, nodesensitivity1, nodespecificity1 = tuple1
				sensitivity2, specificity2, nodesensitivity2, nodespecificity2 = tuple2
				sensitivity3, specificity3, nodesensitivity3, nodespecificity3 = tuple3
				tempsensitivities[0].append(sensitivity1)
				tempsensitivities[1].append(sensitivity2)
				tempsensitivities[2].append(sensitivity3)
				tempspecificities[0].append(specificity1)
				tempspecificities[1].append(specificity2)
				tempspecificities[2].append(specificity3)
				nodesensitivities[0].extend(nodesensitivity1)
				nodesensitivities[1].extend(nodesensitivity2)
				nodesensitivities[2].extend(nodesensitivity3)
				nodespecificities[0].extend(nodespecificity1)
				nodespecificities[1].extend(nodespecificity2)
				nodespecificities[2].extend(nodespecificity3)
				edgeDegree.extend(inEdgeNums)
				truthlists[0].extend(truthlisttemp[0])
				truthlists[1].extend(truthlisttemp[1])
				truthlists[2].extend(truthlisttemp[2])
			for i in range(len(tempsensitivities)):
				tempsensitivities[i]=filter(lambda a: a != 100, tempsensitivities[i])
				if len(tempsensitivities[i])==0:
					tempsensitivities[i].append(0.)
			for tempholder in tempspecificities:
				tempspecificities[i]=filter(lambda a: a != 100, tempspecificities[i])
				if len(tempspecificities[i])==0:
					tempspecificities[i].append(0.)
			sensitivity=[sum(tempsensitivities[0])/len(tempsensitivities[0]),sum(tempsensitivities[1])/len(tempsensitivities[1]),sum(tempsensitivities[2])/len(tempsensitivities[2])]
			specificity=[sum(tempspecificities[0])/len(tempspecificities[0]),sum(tempspecificities[1])/len(tempspecificities[1]),sum(tempspecificities[2])/len(tempspecificities[2])]
			truthholder.append(truthlists)
			specificities.append(specificity)
			sensitivities.append(sensitivity)
			print(sensitivity)
			print(specificity)
			print(nodesensitivities)
			print(nodespecificities)
	nodeLookup={}
	for number in edgeDegree:
		nodeLookup[number]=[[],[],[],[],[],[]]
	for i in range(len(edgeDegree)):
		nodeLookup[edgeDegree[i]][0].append(nodesensitivities[0][i])
		nodeLookup[edgeDegree[i]][1].append(nodesensitivities[1][i])
		nodeLookup[edgeDegree[i]][2].append(nodesensitivities[2][i])
		nodeLookup[edgeDegree[i]][3].append(nodespecificities[0][i])
		nodeLookup[edgeDegree[i]][4].append(nodespecificities[1][i])
		nodeLookup[edgeDegree[i]][5].append(nodespecificities[2][i])
	finalNodeData=[]
	for key in nodeLookup.keys():
		templisting=[]
		for lister in nodeLookup[key]:
			newlist=filter(lambda a: a != 100, lister)
			if len(newlist)==0:
				templisting.append(0.)
			else:
				templisting.append(sum(newlist)/len(newlist))
		finalNodeData.append(templisting)
	print(nodeLookup.keys())
	print(finalNodeData)
	pickle.dump( finalNodeData, open( "node_by_node_data.pickle", "wb" ) )
	pickle.dump( [sensitivities,specificities], open( "network_by_network_data.pickle", "wb" ) )
	pickle.dump( truthholder, open( "expt_true_corrected_bits.pickle", "wb" ) )
	print(sensitivities)
	print(specificities)
	print(finalNodeData)
def rewireSimTest(graph):
	graph2=nc.rewireNetwork(graph)
	params=sim.paramClass()
	sss=utils.synthesizeInputs(graph,params.samples)
	model=sim.modelClass(graph,sss)
	#print(model.nodeList)
	#print(model.size)
	#print(model.individualParse)
	return simTester(model, sss,'prop')

if __name__ == '__main__':
	import time
	start_time = time.time()
	ifngStimTest(10)
	print("--- %s seconds ---" % (time.time() - start_time))
