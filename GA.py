
from deap import base
from deap import creator
from deap import gp
from deap import tools
from deap import algorithms as algo
from random import random
import utils as utils
import liu_networks as liu
import simulation as sim
import networkx as nx
import numpy as numpy
import copy as copy
import operator
import matplotlib.pyplot as plt
import liu_networks as liu
from random import shuffle
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
		if params.async:
			boolValues=sim.iterateModel(individual, params, model, simulator, model.initValueList[j], False)	
		else:
			boolValues=sim.runModel(individual, params, model,simulator, model.initValueList[j], False)	
		for i in range(0, len(model.evaluateNodes)):
			SSE=SSE+(boolValues[model.evaluateNodes[i]]-ss[model.nodeList[model.evaluateNodes[i]]])**2
		SSEs.append(SSE)
	edgeDegree=0
	if params.complexPenalty:
		for i in range(0,len(model.nodeList)):
			if model.possibilityNumList[i]>0:
				edgeDegree=edgeDegree+len(model.inputOrderList[i][bit2int(individual[model.individualParse[i][0]:model.individualParse[i][1]])%model.possibilityNumList[i]])
	SSEs.append(edgeDegree/1.)
	summer=0.
	for i in range(0,len(SSEs)):
		summer+=SSEs[i]
	return summer,
	

# generates a random set of samples made up of cells by using parameteris from probInit seq
# to set up then iterating using strict Boolean modeling. 
def runProbabilityBooleanSims(individual, model, probInitSeq, sampleNum, cells):
	samples=[]
	simulator=sim.simulatorClass('fuzzy')
	params=sim.paramClass()
	for i in range(0,sampleNum):
		cellArray=[]
		for j in range(0,len(probInitSeq.startProbs)):
			probInitSeq.startProbs[j]=random()
		for j in range(0,cells):
			# shuffle nodes to be initially called.... 
			#simulations that are truly random starting states should have all nodes added to this list
			startnodes= probInitSeq.startNodes
			shuffle(startnodes)
			#get initial values for all nodes
			initValues=genPBNInitValues(individual, model, probInitSeq)
			# run Boolean simulation with initial values and append
			cellArray.append(sim.runModel(individual, model, simulator, initValues))
		samples.append([sum(col) / float(cells) for col in zip(*cellArray)])
	return samples

def genPBNInitValues(individual, model, probInitSeq):
	simulator=sim.simulatorClass('fuzzy')
	initValues=[0 for x in range(0,len(model.nodeList))]
	for node in range(0,len(probInitSeq.startNodes)):
		if random()<probInitSeq.startProbs[node]:
			initValues[probInitSeq.startNodes[node]]=1 
	for node in range(0,len(probInitSeq.nodeOrder)):
		if(sim.updateNet(probInitSeq.nodeOrder[node],initValues,individual, model.individualParse[probInitSeq.nodeOrder[node]],model, simulator)==1):
			initValues[probInitSeq.nodeOrder[node]]=1
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
	toolbox, stats=buildToolbox( )
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
		currentDev=10*len(sss1)
		best=[]
		if model.andLenList[i]>0:
			for j in range(1,2**(model.andLenList[i])):
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


def simTester(model, probInitSeq, simClass):
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
	truthCounter=0
	
	# loop over number of times we want to generate fake data and perform sequence of events
	for i in range(0,trials):

		#generate random set of logic rules to start with
		individual=utils.genRandBits(model.size)
		# doublecheck that all nodes that have no upstream contributors are part of probInitSeq
		initSeq=updateInitSeq(individual, model, probInitSeq)
		# generate Boolean model for this trial
		output=runProbabilityBooleanSims(individual, model, probInitSeq, samples, cells)
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
		print(len(newSSS))
		print(samples)
		print(len(sss))
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
		
		#perform brute force search
		bruteOut=bruteForceSearchModel(model, newSSS,propSimulator)
		


		print(individual)
		print(bruteOut)
		truth=utils.writeModel(individual, model)
		BF=utils.writeModel(bruteOut,model)
		# print(BF)





		if truth==BF:
			truthCounter=truthCounter+1
		else:
			truthlines=truth.split('\n')
			newlines=BF.split('\n')
			trueNodes=0
			for k in range(0,len(truthlines)):
				if truthlines[k]==newlines[k]:
					trueNodes=trueNodes+1
				else:
					print("incorrect pair: true then test")
					print(truthlines[k])
					print(newlines[k])
			trueNodeList.append(trueNodes)
		print(i)
	print(trueNodeList)
	print(truthCounter)
		#print(avgs)

	# 	output, hof=GAsolver(model,params, sss,simClass, toolbox, propSimulator, stats)
		
	# 	for j in range(0,10):
	# 		bestRun=(utils.writeModel(hof[j], model))
	# 		if truth==bestRun:
	# 			truthCounter[j]+=1
	# 			break
	# 		elif j==0:
	# 			truthlines=truth.split('\n')
	# 			newlines=bestRun.split('\n')
	# 			trueNodes=0
	# 			for k in range(0,len(truthlines)):
	# 				if truthlines[k]==newlines[k]:
	# 					trueNodes=trueNodes+1
	# 				else:
	# 					print("incorrect pair: true then test")
	# 					print(truthlines[k])
	# 					print(newlines[k])
	# 			trueNodeList.append(trueNodes)
	# 	print(trueNodeList)
	# 	avgs=[output[1][k]['min'] for k in range(0,len(output[1]))]
	# 	#print(avgs)
	# 	hofs.append(hof)
	# 	temp=[]
	# 	for hofind in range(0,len(hof)):
	# 		tempVal=0
	# 		for bit in range(0,len(hof[hofind])):
	# 			if not hof[hofind][bit]==individual[bit]:
	# 				tempVal=tempVal+1
	# 		temp.append(tempVal)
	# 	hofScores.append(temp)
	# 	truthIndividuals.append(individual)
	# 	truthAvgs.append(boolValues)
	# f = open('hof_differences_'+str(nodeNoise)+str(networkNoise)+'.txt', 'w')
	# g = open('hof_individuals'+str(nodeNoise)+str(networkNoise)+'.txt',  'w')
	# h = open('truth_individuals'+str(nodeNoise)+str(networkNoise)+'.txt',  'w')
	# f.write(str(hofScores))
	# f.write('\n')
	# for hofy in range(0,len(hofs)):
	# 	g.write(str(hofs[hofy]))
	# 	g.write('\n')
	# h.write(str(truthIndividuals))
	# h.write('\n')
	# f.close()
	# g.close()
	# h.close()
	# print('# true of # trials')
	# for k in range(0,len(truthCounter)):
	# 	print(truthCounter[k])
	# print(trials)
	# print("for the incorrects, by node")
	# print(trueNodeList)
	# print(len(truthlines))
	




if __name__ == '__main__':
	params=sim.paramClass()
	graph=liu.LiuNetwork1Builder()
	graph.add_edge('a','f', signal='a')
	sss=utils.synthesizeInputs(graph,params.samples)
	model=sim.modelClass(graph,sss)
	probInitSeq=liu.liu1probInitSeqBuilder(model)
	simTester(model, probInitSeq,'prop')
	nx.write_graphml(graph,'Liu1.graphml')