
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
import scipy.stats as regress
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
	
#exhaustively search boolean networks for best option going node by node for rules
def bruteForceSearchModel(model, simulator,graph,sss1):
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
	print(model.initValueList)

	for i in range(0,len(model.nodeList)):
		# print(model.nodeList[i])
		# print(model.individualParse[i])
		currentDev=10*len(sss1)
		best=[]
		if model.possibilityNumList[i]>0:
			for j in range(0,int(2**int(model.individualParse[i][2]-model.individualParse[i][0]))):
				bits=[]
				bits=utils.bitList(j)
				while len(bits)<(model.individualParse[i][2]-model.individualParse[i][0]):
					bits.insert(0,0)
				deviation=0
				for steadyStateNum in range(0,len(model.initValueList)):
					derivedVal=sim.updateNode(i,model.initValueList[steadyStateNum],bits, model, simulator)
					deviation=deviation+(derivedVal-sss1[steadyStateNum][model.nodeList[i]])**2
				print(utils.writeBruteNode(i,bits,model))
				print(bits)
				print(deviation)
				if(deviation<currentDev):
					print("best")
					best=bits
					currentDev=deviation	
		bestList.append(best)
		print(model.nodeList[i])
		print(currentDev)

	return [item for sublist in bestList for item in sublist]


def calcDeviation(i,model,bits,simulator,sss1):
	oldTriple=model.individualParse[currentNode]
	triple=[]
	triple.append(0)
	triple.append(oldTriple[1]-oldTriple[0])
	triple.append(oldTriple[2]-oldTriple[0])

	inputOrder=model.inputOrderList[currentNode] # find the list of possible input combinations for the node we are on 
	inputOrderInvert=model.inputOrderInvertList[currentNode] #find list of lists of whether input nodes need to be inverted (corresponds to inputOrder)
	if model.possibilityNumList[currentNode]>0:
		logicOperatorFlags=list(individual[triple[1]:triple[2]]) # find demarcations on individual for bits we need
		inputOrder=list(inputOrder[utils.bit2int(individual[triple[0]:triple[1]])%model.possibilityNumList[currentNode]]) # determine what order of inputs is
		inputOrderInvert=inputOrderInvert[utils.bit2int(individual[triple[0]:triple[1]])%model.possibilityNumList[currentNode]] # lookup which inputs need to be inverted ('not')
		if len(inputOrder)==0:
			for oldValue in model.initValueList:
				value.append(oldValue[currentNode]) #if no inputs, maintain value
			return value
		elif len(inputOrder)==1:
			#if only one input, then can either affect or not affect the node. so either keep the value or update to the single input's value
			for oldValue in model.initValueList:
				value.append(Inv(oldValue[inputOrder[0]],inputOrderInvert[0]))
			return value
		else:
			#update nodes with more than one input
			# update and then or
			upstreamValueHolder=[]
			for oldValue in model.initValueList:
				upstreamVals=[]
				for upstream in range(0,len(inputOrder)):
					upstreamVals.append(Inv(oldValue[inputOrder[upstream]],inputOrderInvert[upstream]))
				upstreamValueHolder.append(upstreamVals)
			counter =0
			
			# print(upstreamVals)
			# print(logicOperatorFlags)
			while counter < len(logicOperatorFlags) and counter+1<len(inputOrder):
				if logicOperatorFlags[counter]==0:

					node1=[]
					node2=[]
					for upstreamVals in upstreamValueHolder:
						node1.append(upstreamVals[counter])
						node2.append(upstreamVals[counter+1])
					slope, intercept, r_value, p_value, std_err = regress.linregress(node2,node1)
					AB=slope+intercept
					slope, intercept, r_value, p_value, std_err = regress.linregress(node1,node2)
					BA=slope+intercept
					for upstreamVals in upstreamValueHolder:
						tempVal=simulator.And(upstreamVals[counter],upstreamVals[counter+1],AB,BA)
						newPartialValues.append(tempVal)
						upstreamVals.pop(counter)
						upstreamVals.pop(counter)
						upstreamVals.insert(counter,tempVal)
					inputOrder.pop(counter)
					inputOrder.pop(counter)
					logicOperatorFlags.pop(counter)

				else:
					counter=counter+1
				# print(upstreamVals)

			#first one uses the initial logic operator flag to decide and vs or then combines the first two inputs
			while len(upstreamVals)>1:
				tempVal=simulator.Or(upstreamVals.pop(0),upstreamVals.pop(0),inputOrder.pop(0),inputOrder.pop(0),simulator.corrMat)
				upstreamVals.insert(0,tempVal)
				# print(upstreamVals)
			return upstreamVals[0]
	else:
		#returns savme value if now inputs
		return oldValue[currentNode]	


	for steadyStateNum in range(0,len(model.initValueList)):
		derivedVal=sim.updateNode(i,model.initValueList[steadyStateNum],bits, model, simulator)
		deviation=deviation+(derivedVal-sss1[steadyStateNum][model.nodeList[i]])**2
#exhaustively searcj noolean networks for best option while simultaneously solving all the rules over all values to continuously update table
def bruteForceSolveSearch(model, simulator,graph,sss1):
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
	print(model.initValueList)

	for i in range(0,len(model.nodeList)):
		# print(model.nodeList[i])
		# print(model.individualParse[i])
		currentDev=10*len(sss1)
		best=[]
		if model.possibilityNumList[i]>0:
			for j in range(0,int(2**int(model.individualParse[i][2]-model.individualParse[i][0]))):
				bits=[]
				bits=utils.bitList(j)
				while len(bits)<(model.individualParse[i][2]-model.individualParse[i][0]):
					bits.insert(0,0)
				deviation=0
				deviation=calcDeviation(i,model,bits,simulator,sss1)
				print(utils.writeBruteNode(i,bits,model))
				print(bits)
				print(deviation)
				if(deviation<currentDev):
					print("best")
					best=bits
					currentDev=deviation	
		bestList.append(best)
		print(model.nodeList[i])
		print(currentDev)
	return [item for sublist in bestList for item in sublist]

def runProbabilityBooleanSims(individual, model, probInitSeq, sampleNum, cells):
	samples=[]
	simulator=sim.simulatorClass('fuzzy')
	params=sim.paramClass()
	for i in range(0,sampleNum):
		cellArray=[]
		for j in range(0,len(probInitSeq.startProbs)):
			probInitSeq.startProbs[j]=random()
		for j in range(0,cells):
			initValues=genProbabilisticInitValues(individual, model, probInitSeq)
			cellArray.append(sim.runModel(individual, params, model, simulator, initValues, False))
		samples.append([sum(col) / float(cells) for col in zip(*cellArray)])
	return samples

def genProbabilisticInitValues(individual, model, probInitSeq):
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
	initSeq=probInitSeqClass()
	initSeq.startNodes=list(probInitSeq.startNodes)
	initSeq.startProbs=list(probInitSeq.startProbs)
	initSeq.nodeOrder=list(probInitSeq.nodeOrder)
	for i in range(0,len(model.individualParse)):
		triple=model.individualParse[i]
		inputOrder=model.inputOrderList[i] # find the list of possible input combinations for the node we are on 
		inputOrderInvert=model.inputOrderInvertList[i] #find list of lists of whether input nodes need to be inverted (corresponds to inputOrder)
		# print(i)
		# print(model.possibilityNumList[i])
		if model.possibilityNumList[i]>0:
			logicOperatorFlags=list(individual[triple[1]:triple[2]]) # find demarcations on individual for bits we need
			inputOrder=list(inputOrder[utils.bit2int(individual[triple[0]:triple[1]])%model.possibilityNumList[i]]) # determine what order of inputs is
			inputOrderInvert=inputOrderInvert[utils.bit2int(individual[triple[0]:triple[1]])%model.possibilityNumList[i]] # lookup which inputs need to be inverted ('not')
			# print(inputOrder)
			# print(model.nodeList[inputOrder[0]])
			if len(inputOrder)<1 and i not in initSeq.startNodes:
				initSeq.startNodes.append(i)
				initSeq.nodeOrder.pop(i)
				initSeq.startProbs.append(1)		
		else:
			# print('no possibilities')
			if i not in initSeq.startNodes:
				# print('added')
				initSeq.startNodes.append(i)
				initSeq.nodeOrder.pop(i)
				initSeq.startProbs.append(1)

	return initSeq

def simTester(trials, model, graph, samples, probInitSeq, simClass):
	#creates a model, runs simulations, then tests reverse engineering capabilities of models in a single function
	#trials is the number of different individuals to try
	#samples is the number of different initial conditions to provide per trial
	#graph specifies the network we are testing. 
	#probInitSeq gives the order in which the network should be built

	params=sim.paramClass()

	avgSSEs=[]
	hofs=[]
	truthIndividuals=[]
	truthAvgs=[]
	hofScores=[]
	newSSS=[]
	truthCounter=[0 for number in xrange(params.hofSize)]
	trueNodeList=[]

	
	# toolbox, stats=buildToolbox( model.size, params.bitFlipProb)
	print(model.inputOrderList)
	print(model.individualParse)
	print(model.possibilityNumList)
	print(model.size)
	print(model.possibilityNumList[0:2])
	
	for i in range(0,trials):
		individual=utils.genRandBits(model.size)
		initSeq=updateInitSeq(individual, model, probInitSeq)
		output=testPropInference(model, sss, params.cells, samples, initSeq, individual)
		newSSS=[]
		for k in range(0,samples):
			newSS=copy.deepcopy(sss[k])
			for j in range(0,len(model.nodeList)):
				newSS[model.nodeList[j]]=output[k][j]
			newSSS.append(newSS)
		newInitValueList=[]
		for j in range(0,len(sss)):
			newInitValueList.append([])
		for j in range(0,len(model.nodeList)):
			for k in range(0,len(sss)):
				ss=newSSS[k]
				if  model.nodeList[j] in sss[0].keys():
					newInitValueList[k].append(ss[model.nodeList[j]])
				else:
					newInitValueList[k].append(0.5)
		model.initValueList=newInitValueList

		#set up prop simulator and corr list
		matt=[]
		for j in range(0,len(model.nodeList)):
			templist=[]
			jarray=[output[q][j] for q in range(0,samples)]
			for k in range(0,len(model.nodeList)):
				karray=[output[q][k] for q in range(0,samples)]
				slope, intercept, r_value, p_value, std_err = regress.linregress(jarray,karray)
				# print(jarray)
				# print(karray)
				# print(slope)
				# print(intercept)
				# print(slope-intercept)
				if((slope+intercept)>0):
					if (slope+intercept)>1:
						templist.append(1)
					else:
						templist.append(slope+intercept)
				else:
					templist.append(0)
			matt.append(templist)
		propSimulator=sim.simulatorClass(simClass)
		propSimulator.corrMat=matt
		
		boolValues=evaluate(individual, params, model, propSimulator, newSSS)

		print('\n'.join([''.join(['{:10}'.format(item) for item in row]) for row in matt]))
		outs=bruteForceSolveSearch(model, propSimulator,graph,newSSS)
		truth=utils.writeModel(individual, model)
		print("truth")
		print(truth)
		BF=utils.writeModel(outs,model)
		print(BF)
		# stats = tools.Statistics(key=lambda ind: ind.fitness.values)
		# stats.register("avg", numpy.mean)
		# stats.register("std", numpy.std)
		# stats.register("min", numpy.min)
		# stats.register("max", numpy.max)
		# toolbox.register("evaluate", evaluate, params=params,model=model,simulator=propSimulator,sss=newSSS)
		


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
	

def testPropInference(model, sss, cells, samples, probInitSeq, individual):
	
	print(utils.writeModel(individual,model))
	output=runProbabilityBooleanSims(individual, model, probInitSeq, samples, cells)
	#for node in range(0,len(model.nodeList)):
	#	initValues.append(sss[0][model.nodeList[node]])
	#return initValues
	return output

def liu1probInitSeqBuilder(model):
	probInitSeq=probInitSeqClass()
	probInitSeq.startNodes.append(model.nodeDict['j'])
	probInitSeq.startProbs.append(.5)
	probInitSeq.startNodes.append(model.nodeDict['a'])
	probInitSeq.startProbs.append(.5)
	probInitSeq.startNodes.append(model.nodeDict['b'])
	probInitSeq.startProbs.append(.5)
	probInitSeq.nodeOrder.append(model.nodeDict['c'])
	probInitSeq.nodeOrder.append(model.nodeDict['h'])
	probInitSeq.nodeOrder.append(model.nodeDict['d'])
	probInitSeq.nodeOrder.append(model.nodeDict['f'])
	probInitSeq.nodeOrder.append(model.nodeDict['g'])
	probInitSeq.nodeOrder.append(model.nodeDict['k'])
	return probInitSeq

def GAautoSolver(model,params, sss,simClass):
	propSimulator=simulatorClass(simClass)
	toolbox, stats=buildToolbox( )
	toolbox.register("evaluate", evaluate, params=params,model=model,simulator=propSimulator,sss=sss)
	return 	GAsolver(model,params, sss,simClass, toolbox, propSimulator, stats)

def GAsolver(model,params, sss,simClass, toolbox, propSimulator, stats):
	population=toolbox.population(n=params.popSize)
	hof = tools.HallOfFame(params.hofSize, similar=numpy.array_equal)
	output=algo.eaMuCommaLambda(population, toolbox, mu=params.mu, lambda_=params.lambd, stats=stats, cxpb=params.crossoverProb, mutpb=params.mutationProb, ngen=params.generations, verbose=False, halloffame=hof)
	return output, hof

if __name__ == '__main__':
	samples=15
	trials=1
	graph=liu.LiuNetwork1Builder()
	sss=utils.synthesizeInputs(graph,samples)
	model=sim.modelClass(graph,sss)
	probInitSeq=liu1probInitSeqBuilder(model)
	simTester(trials, model, graph, samples, probInitSeq,'propE2')
	nx.write_graphml(graph,'Liu1.graphml')