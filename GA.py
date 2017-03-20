
from deap import base
from deap import creator
from deap import gp
from deap import tools
from deap import algorithms as algo
from random import random
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
	simulator=sim.simulatorClass('fuzzy')
	params=sim.paramClass()
	for i in range(0,sampleNum):
		cellArray=[]
		sampleProbs=[]
		for j in range(0,len(model.nodeList)):
			sampleProbs.append(random())
		for j in range(0,cells):
			# shuffle nodes to be initially called.... 
			#simulations that are truly random starting states should have all nodes added to this list
			#get initial values for all nodes
			initValues=genPBNInitValues(individual, model,sampleProbs)
			# run Boolean simulation with initial values and append
			vals, nums=sim.runModel(individual, model, simulator, initValues)
			cellArray.append(vals)
		samples.append([sum(col) / float(cells) for col in zip(*cellArray)])
	return samples

def genPBNInitValues(individual, model,sampleProbs):
	simulator=sim.simulatorClass('fuzzy')
	initValues=[0 for x in range(0,len(model.nodeList))]
	for node in range(0,len(sampleProbs)):
		if random()<sampleProbs[node]:
			initValues[node]=1 
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
	params=sim.paramClass()
	bestList=[]
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
	for i in range(0,len(model.nodeList)):
		# print(model.nodeList[i])
		# print(model.individualParse[i])
		currentDev=10000*len(sss1)
		best=[]
		if model.andLenList[i]>0:
			if model.andLenList[i]<10:
				checkRange=2**(model.andLenList[i])
			else:
				checkRange=3*(model.andLenList[i])
			for j in range(1,checkRange):
				bits=[]
				bits=utils.bitList(j,model.andLenList[i] )
				deviation=0
				for steadyStateNum in range(0,len(model.initValueList)):
					derivedVal=sim.updateNode(i,model.initValueList[steadyStateNum],bits, model, simulator)
					deviation=deviation+(derivedVal-sss1[steadyStateNum][model.nodeList[i]])**2
				# print(utils.writeBruteNode(i,bits,model))
				# print(bits)
				# print(deviation)
				if(deviation<currentDev):
					# print("best")
					best=bits
					currentDev=deviation	
		bestList.append(best)
		# print(model.nodeList[i])
		# print(currentDev)

	return [item for sublist in bestList for item in sublist]

def simTester(model, sss, simClass):
	#creates a model, runs simulations, then tests reverse engineering capabilities of models in a single function
	#trials is the number of different individuals to try
	#samples is the number of different initial conditions to provide per trial
	#graph specifies the network we are testing. 
	#probInitSeq gives the order in which the network should be built... 
	# probInitSeq should contain all nodes for a totally random setup. 

	# set up empty lists and call parameters
	params=sim.paramClass()
	samples=params.samples
	trials=params.trials
	cells=params.cells
	avgSSEs=[]
	hofs=[]
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
	
	# loop over number of times we want to generate fake data and perform sequence of events
	for i in range(0,trials):

		#generate random set of logic rules to start with
		individual=utils.genRandBits(model.size)
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
		# print(individual)
		# print(bruteOut)
		sumindividual.append(sum(individual))
		newindividual=[a_i - b_i for a_i, b_i in zip(individual, bruteOut)]
		ones.append(newindividual.count(1))
		zeros.append(newindividual.count(0))
		negones.append(newindividual.count(-1))
		truth=utils.writeModel(individual, model)
		BF=utils.writeModel(bruteOut,model)
		if truth==BF:
			truthCounter=truthCounter+1
		else:
			truthlines=truth.split('\n')
			newlines=BF.split('\n')
			trueNodes=0
			for k in range(0,len(truthlines)):
				if truthlines[k]==newlines[k]:
					trueNodes=trueNodes+1
				#else:
					# print("incorrect pair: true then test")
					# print(truthlines[k])
					# print(newlines[k])
			trueNodeList.append(trueNodes)
		# print(i)
		# print(trueNodeList)

	# print('# true of # trials')
	# print(truthCounter)
	# print(trials)
	# print("for the incorrects, by node")
	# print(trueNodeList)
	# print(len(truthlines))
	# print(zeros)
	# print(ones)
	# print(negones)

	temp=[(zeros[i])/(zeros[i]+ones[i]) for i in range(0,len(ones))]
	sensitivity=(sum(temp)/len(temp))
	temp=[(1.*len(individual)-sumindividual[i]-negones[i])/(len(individual)-sumindividual[i]) for i in range(0,len(ones))]
	specificity=(sum(temp)/len(temp))
	return sensitivity,specificity

def uploadIFNGstimData(filename):
	inputfile = open(filename, 'r')
	lines = inputfile.readlines()
	filename.split()


def ifngStimTest(bioReplicates):
	aliasDict={}
	dict1={}
	nc.parseKEGGdicthsa('hsa00001.keg',aliasDict,dict1)
	dict2={}
	nc.parseKEGGdict('ko00001.keg',aliasDict,dict2)
	graph=nx.DiGraph()

	# read in list of codes then load them into network
	inputfile = open('IFNG_hsa_pathways.txt', 'r')
	lines = inputfile.readlines()
	uploadList=[line[:-1] for line in lines]
	nc.uploadKEGGcodes_hsa(uploadList,graph,dict1,dict2)
	nx.write_graphml(graph,'IFNG.graphml')
	
	data=dict(utils.loadFpkms('Hela-C-1.count'))
	graph=nc.simplifyNetwork(graph, data)
	nx.write_graphml(graph,'IFNG_simplified.graphml')
	sensitivities=[]
	specificities=[]
	for i in range(bioReplicates):
		tempsens,tempspec=rewireSimTest(graph)
		sensitivities.append(tempsens)
		specificities.append(tempspec)
	sensitivity=sum(sensitivities)/(1.*len(sensitivities))
	specificity=sum(specificities)/(1.*len(specificities))
	print(sensitivity)
	print(specificity)

def rewireSimTest(graph):
	graph2=nc.rewireNetwork(graph)
	params=sim.paramClass()
	sss=utils.synthesizeInputs(graph,params.samples)
	model=sim.modelClass(graph,sss)
	print(model.nodeList)
	print(model.size)
	print(model.individualParse)
	return simTester(model, sss,'prop')

if __name__ == '__main__':
	ifngStimTest(10)
