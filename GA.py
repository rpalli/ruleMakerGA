
# import python packages
import pickle
from deap import base, creator, gp, tools
from deap import algorithms as algo
from random import random, seed, shuffle
import random as rand
import networkx as nx
import numpy as numpy
import copy as copy
import operator
import math as math
from sets import Set
from joblib import Parallel, delayed

# import other pieces of our software
import networkConstructor as nc
import liu_networks as liu
import utils as utils
import simulation as sim

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
			# for i in range(len(truth)):
			# 	if random()<(1./len(truth)):
			# 		truth[i]=1
			# 	else:
			# 		truth[i]=0
			counter=0
			while sum(truth)>5 and counter < 100000:
				indices = [i for i in range(len(truth)) if truth[i] == 1]
				chosen=math.floor(random()*len(indices))
				truth[indices[int(chosen)]]=0
				counter+=1
			startInd[start:end]=truth
			if numpy.sum(truth)==0:
				chosen=math.floor(random()*len(truth))
				truth[int(chosen)]=1
		elif len(truth)==1:
			truth[0]=1
			startInd[start:end]=truth
	return startInd

def buildToolbox( individualLength, bitFlipProb, model):
	# # #setup toolbox
	toolbox = base.Toolbox()
	
	pset = gp.PrimitiveSet("MAIN", arity=1)
	pset.addPrimitive(operator.add, 2)
	pset.addPrimitive(operator.sub, 2)
	pset.addPrimitive(operator.mul, 2)
	weightTup=(-1.0,)
	params=sim.paramClass()
	if params.adaptive:
		for i in range(len(model.evaluateNodes)-1):
			weightTup+=(-1.0,)
	creator.create("FitnessMin", base.Fitness, weights=weightTup) # make a fitness minimization function
	creator.create("Individual", list, fitness=creator.FitnessMin)	# create a class of individuals that are lists

	#register our bitsring generator and how to create an individual, population
	toolbox.register("genRandomBitString", genBits, model=model)
	toolbox.register("Individual", tools.initIterate, creator.Individual, toolbox.genRandomBitString)
	toolbox.register("population", tools.initRepeat, list , toolbox.Individual)
	#create statistics toolbox and give it functions
	stats = tools.Statistics(key=lambda ind: ind.fitness.values)
	stats.register("avg", numpy.mean)
	stats.register("std", numpy.std)
	stats.register("min", numpy.min)
	stats.register("max", numpy.max)
	
	# finish registering the toolbox functions
	toolbox.register("mate", tools.cxTwoPoint)
	toolbox.register("mutate", tools.mutFlipBit, indpb=bitFlipProb)
	toolbox.register("select", tools.selBest)
	toolbox.register("similar", numpy.array_equal)
	return toolbox, stats
			
def evaluate(individual, cells, model,  sss, params):
	SSEs=[]
	for j in range(0,len(sss)):
		boolValues=iterateBooleanModel(individual, model, cells, model.initValueList[j], params)
		SSE=numpy.sum([(boolValues[model.evaluateNodes[i]]-sss[j][model.nodeList[model.evaluateNodes[i]]])**2 for i in range(0, len(model.evaluateNodes))])
		SSEs.append(SSE)
	summer=1.*numpy.sum(SSEs)/len(SSEs)
	for j in range(len(model.nodeList)):
		if model.andLenList[j]>0:
			if j==len(model.nodeList)-1:
				end= model.size
			else:
				end=model.individualParse[j+1]
			temptruther=individual[model.individualParse[j]:end]
			if len(temptruther)==1:
				individual[model.individualParse[j]:end]=[1]
			elif len(temptruther)>1:
				if numpy.sum(temptruther)==0:
					summer+=1000
	return summer,
	
def evaluateByNode(individual, cells, model,  sss, params):
	SSEs=[]
	boolValues=Parallel(n_jobs=min(6,len(sss)))(delayed(iterateBooleanModel)(list(individual), model, cells, model.initValueList[i], params) for i in range(len(sss)))
	for i in range(0, len(model.evaluateNodes)):
		SSE= numpy.sum([(boolValues[j][model.evaluateNodes[i]]-sss[j][model.nodeList[model.evaluateNodes[i]]])**2 for j in range(0,len(sss))])
		SSEs.append(SSE)
	return tuple(SSEs)
# generates a random set of samples made up of cells by using parameteris from probInit seq
# to set up then iterating using strict Boolean modeling. 
def runProbabilityBooleanSims(individual, model, sampleNum, cells, params):
	samples=[]
	seeds=[]
	for i in range(0,sampleNum):
		seeds.append(random())
	samples=Parallel(n_jobs=min(6,sampleNum))(delayed(sampler)(individual, model, sampleNum, seeds[i], params) for i in range(sampleNum))
	# counter=0
	# for sample in range(len(samples)):
	# 	if sum(samples(sample))==0:
	# 		temp=samples.pop(sample) 
	# newSamples=Parallel(n_jobs=min(8,sampleNum))(delayed(sampler)(individual, model, sampleNum, seeds[i]) for i in range(sampleNum))
	# print(samples)
	return samples
	
def sampler(individual, model, cells, seeder, params):
	seed(seeder)
	cellArray=[]
	sampleProbs=[]
	simulator=sim.simulatorClass('bool')
	# generate random proportions for each node to start
	for j in range(0,len(model.nodeList)):
		sampleProbs.append(random())
	return iterateBooleanModel(individual, model, cells, sampleProbs, params)

def iterateBooleanModel(individual, model, cells, initProbs, params):
	cellArray=[]
	simulator=sim.simulatorClass('bool')
	simulator.simSteps=3*len(model.nodeList)
	for j in range(0,cells):
		# shuffle nodes to be initially called.... 
		#simulations that are truly random starting states should have all nodes added to this list
		#get initial values for all nodes
		initValues=genPBNInitValues(individual, model,initProbs)
		# run Boolean simulation with initial values and append
		vals=sim.runModel(individual, model, simulator, initValues, params)
		cellArray.append(vals)
	return [numpy.sum(col) / float(cells) for col in zip(*cellArray)]

def genPBNInitValues(individual, model,sampleProbs):
	#return [True if (random()<sampleProbs[node]) else False for node in range(0,len(sampleProbs))]
	initValues=[False for x in range(0,len(model.nodeList))]
	for node in range(0,len(sampleProbs)):
		if random()<sampleProbs[node]:
			initValues[node]=True
	return initValues

def varOrAdaptive(population, toolbox, model, lambda_, cxpb, mutpb):
	# algorithm for generating a list of offspring... copied and pasted from DEAP with modification for adaptive mutation
	assert (cxpb + mutpb) <= 1.0, ("The sum of the crossover and mutation "
		"probabilities must be smaller or equal to 1.0.")
	offspring = []
	for _ in xrange(lambda_):
		op_choice = random()
		if op_choice < cxpb:            # Apply crossover
			ind1, ind2 = map(toolbox.clone, rand.sample(population, 2))
			ind1, ind2 = toolbox.mate(ind1, ind2)
			del ind1.fitness.values
			offspring.append(ind1)
		elif op_choice < cxpb + mutpb:  # Apply mutation
			ind = toolbox.clone(rand.choice(population))
			ind, = mutFlipBitAdapt(ind,  .5, model)
			del ind.fitness.values
			offspring.append(ind)
		else:                           # Apply reproduction
			offspring.append(rand.choice(population))
	return offspring

def mutFlipBitAdapt(individual, indpb, model):
	# get errors
	errors=list(individual.fitness.values)
	# if error is relatively low, do a totally random mutation
	if numpy.sum(errors)<.1:
		for i in xrange(len(individual)):
			if random() < indpb:
				individual[i] = type(individual[i])(not individual[i])
	else:
		# if errors are relatively high, focus on nodes that fit the worst
		for j in xrange(len(errors)):
			if model.andLenList[model.evaluateNodes[j]]<2:
				errors[model.evaluateNodes[j]]=0
		# normalize errors to get a probability that the node  is modified
		normerrors=[error/numpy.sum(errors) for error in errors]
		probs=numpy.cumsum(normerrors)
		# randomly select a node to mutate
		randy=random()
		focusNode=next(i for i in xrange(len(probs)) if probs[i]>randy)
		# mutate on average 1.5 shadow ands from the node. 
		if model.andLenList[model.evaluateNodes[focusNode]]>1:
			# find ends of the node of interest in the individual
			start=model.individualParse[model.evaluateNodes[focusNode]]
			if model.evaluateNodes[focusNode]==len(model.nodeList)-1:
				end= model.size
			else:
				end=model.individualParse[model.evaluateNodes[focusNode]+1]
			for i in range(start,end):
				if random()< .5:
					individual[i] = 1
				else:
					individual[i] = 0
			#ensure that there is at least one shadow and node turned on
			if sum(individual[start:end])==0:
				individual[start+1]=1
	print(individual)
	return individual,

def eaMuPlusLambdaAdaptive(population, toolbox, model, mu, lambda_, cxpb, mutpb, ngen, stats=None, halloffame=None, verbose=__debug__):
	logbook = tools.Logbook()
	logbook.header = ['gen', 'nevals'] + (stats.fields if stats else [])

	# Evaluate the individuals with an invalid fitness
	invalid_ind = [ind for ind in population if not ind.fitness.valid]
	fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
	for ind, fit in zip(invalid_ind, fitnesses):
		ind.fitness.values = fit

	if halloffame is not None:
		halloffame.update(population)

	record = stats.compile(population) if stats is not None else {}
	logbook.record(gen=0, nevals=len(invalid_ind), **record)
	if verbose:
		print logbook.stream

	# Begin the generational process
	for gen in range(1, ngen+1):
		# Vary the population
		errors=[ind.fitness.values for ind in population]
		offspring = varOrAdaptive(population, toolbox, model, lambda_, cxpb, mutpb)
		
		# Evaluate the individuals with an invalid fitness
		invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
		fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
		#fitnesses=Parallel(n_jobs=min(6,len(invalid_ind)))(delayed(toolbox.evaluate)(indy) for indy in invalid_ind)
		for ind, fit in zip(invalid_ind, fitnesses):
			ind.fitness.values = fit
		


		# Update the hall of fame with the generated individuals
		if halloffame is not None:
			halloffame.update(offspring)

		# Select the next generation population
		population[:] = toolbox.select(population + offspring, mu)

		# Update the statistics with the new population
		record = stats.compile(population) if stats is not None else {}
		logbook.record(gen=gen, nevals=len(invalid_ind), **record)
		if verbose:
			print logbook.stream
		if numpy.sum(halloffame[0].fitness.values)<.4:
			break
	return population, logbook

def GAautoSolver(model, sss, params):
	# set up toolbox and run GA with or without adaptive mutations turned on
	toolbox, stats=buildToolbox(model.size,params.bitFlipProb, model)

	if params.adaptive:
		toolbox.register("evaluate", evaluateByNode, cells=params.cells,model=model,sss=sss, params=params)
	else:
		toolbox.register("evaluate", evaluate, cells=params.cells,model=model,sss=sss)
	population=toolbox.population(n=params.popSize)
	hof = tools.HallOfFame(params.hofSize, similar=numpy.array_equal)
	
	if params.adaptive:
		output=eaMuPlusLambdaAdaptive(population, toolbox, model, mu=params.mu, lambda_=params.lambd, stats=stats, cxpb=params.crossoverProb, mutpb=params.mutationProb, ngen=params.generations, verbose=False, halloffame=hof)
	else:
		output=algo.eaMuCommaLambda(population, toolbox, mu=params.mu, lambda_=params.lambd, stats=stats, cxpb=params.crossoverProb, mutpb=params.mutationProb, ngen=params.generations, verbose=False, halloffame=hof)
	return output, hof

def GAsearchModel(model, newSSS,params):
	output, hof=GAautoSolver(model, newSSS, params)
	return hof[0].fitness.values, hof[0]

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
			sumindividual.append(numpy.sum(truth))
			newindividual=[a_i - b_i for a_i, b_i in zip(truth, test)]
			ones.append(newindividual.count(1))
			zeros.append(newindividual.count(0))
			negones.append(newindividual.count(-1))
		# append node-wise breakdowns to list of breakdowns for the model as a whole
		netOnes.append(numpy.sum(ones))
		netZeros.append(numpy.sum(zeros))
		netNegOnes.append(numpy.sum(negones))
		
		# calculate sensitivity and specificity for the node
		temp=[100 if (sumindividual[i]-negones[i]+ones[i])==0 else (sumindividual[i]-negones[i])/(sumindividual[i]-negones[i]+ones[i]) for i in range(0,len(ones))]
		temp=filter(lambda a: a != 100, temp)
		if len(temp)==0:
			sensitivity=100
		else:
			sensitivity=(sum(temp)/len(temp))
		print(sensitivity)
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
	sumindividual=[]

	for i in range(len(truthList)):
		truth= truthList[i]
		test= testList[i]
		sumindividual.append(1.*sum(truth))
		newindividual=[a_i - b_i for a_i, b_i in zip(truth, test)]
		ones.append(newindividual.count(1))
		zeros.append(newindividual.count(0))
		negones.append(newindividual.count(-1))
	temp=[100 if (sumindividual[i])==0 else (sumindividual[i]-ones[i])/(sumindividual[i]) for i in range(0,len(ones))]
	temp=filter(lambda a: a != 100, temp)
	if len(temp)==0:
		sensitivity=100
	else:
		sensitivity=(numpy.sum(temp)/len(temp))
	temp=[100 if (len(newindividual)-sumindividual[i])==0 else (1.*len(newindividual)-sumindividual[i]-negones[i])/(len(newindividual)-sumindividual[i]) for i in range(0,len(ones))]
	temp=filter(lambda a: a != 100, temp)
	if len(temp)==0:
		specificity=100
	else:
		specificity=(sum(temp)/len(temp))
	return sensitivity, specificity, nodesensitivity, nodespecificity

		
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
	devList=[]
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
				if numpy.sum(individual[model.individualParse[node]:end])==0:
					individual[model.individualParse[node]]=1

		# generate Boolean model for this trial
		output=runProbabilityBooleanSims(individual, model, samples, cells, params)
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


		# # set up PBN-based simulator
		# propSimulator=sim.simulatorClass(simClass)
		# propSimulator.trainingData=model.initValueList
		# propSimulator.train(model)
		# #perform brute force search
		# bruteOut=bruteForceSearchModel(model, newSSS,propSimulator)
		
		#perform GA run
		dev, bruteOut=GAsearchModel(model, newSSS, params)
		truthList.append(individual)
		testList.append(list(bruteOut))
		devList.append(dev)
	print('devList')
	print(devList)

	tuple1=compareIndividualsNodeWise(truthList, testList, model)
	sensitivities=[]
	specificities=[]
	inEdgeNums=[]
	overlaps=[]
	# run correction that simplifies truth rules
	TPsum=0
	TNsum=0
	FNsum=0
	FPsum=0
	newtruths=[]
	for k in range(len(truthList)):
		newtruths.append([])
	newtests=[]
	for k in range(len(testList)):
		newtests.append([])
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
						if truth[j]==1 and not i==j:
							if inEdges[i].issubset(inEdges[j]):
								truth[j]=0
			for i in range(len(test)):
				if test[i]==1:
					for j in range(len(test)):
						if test[j]==1 and not i==j:
							if inEdges[i].issubset(inEdges[j]):
								test[j]=0
			newtests[k].extend(test)
			newtruths[k].extend(truth)
		for k in range(len(truthList)):
			truth=newtruths[k][start:end]
			test= newtests[k][start:end]
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
			TN+=1.*(len((baseSet.difference(truthSet)).difference(testSet)))
			FN+=1.*len(truthSet.difference(testSet))
		if (TP+FN)>0:
			sensitivity=TP/(TP+FN)
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
		overlaper=-1
		if len(baseSet)==2:
			overlapSet=Set([])
			lappersets=[]
			inneredge=[]
			for upstream in baseSet:
				lapperset=Set([])
				if upstream==(len(model.nodeList)-1):
					end=len(truthList[0])-1
				else:
					end=model.individualParse[upstream+1]
				start=model.individualParse[upstream]
				andNodeList=model.andNodeList[upstream]
				for lister in andNodeList:
					lapperset.add(Set(lister))
				lappersets.append(lapperset)
			overlaper=len(lappersets[0].intersection(lappersets[1]))
			if overlaper==0:
				for upstream in baseSet:
					if upstream in lappersets[0] or upstream in lappersets[1]:
						print("feed-forward")
						overlaper=1
		overlaps.append(overlaper)
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

	tuple4=compareIndividualsNodeWise(newtruths, newtests, model)
	for i in range(len(newtruths)):
		print("truth, found")
		print(utils.writeModel(newtruths[i], model))
		print(utils.writeModel(newtests[i], model))

	return tuple1, tuple2, tuple3, tuple4, devList, inEdgeNums, overlaps, [testList, truthList,  newtruths, newtests]

def simpleNetBuild():
	graph = nx.DiGraph()
	graph.add_edge('a1','a', signal='a')	
	graph.add_edge('a2','a', signal='a')
	graph.add_edge('a3','a', signal='a')
	graph.add_edge('a4','a', signal='a')
	graph.add_edge('b1','b', signal='a')
	graph.add_edge('b2','b', signal='a')
	graph.add_edge('b3','b', signal='a')
	graph.add_edge('a','d', signal='a')	
	graph.add_edge('b','d', signal='a')
	graph.add_edge('d','e', signal='a')
	return graph

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
	specificityStds=[]
	sensitivityStds=[]
	lines.pop(0)
	nodesensitivities=[[],[],[],[]]
	nodespecificities=[[],[],[],[]]
	truthholder=[]
	edgeDegree=[]
	overlapNodes=[]
	lines.pop(0)
	lines=[0]
	lines=[lines[0]]
	for code in lines:
		# graph=nx.DiGraph()
		# coder=str('ko'+code[:-1])
		# nc.uploadKEGGcodes([coder], graph, dict2)
		# coder=str('hsa'+code[:-1])
		# nc.uploadKEGGcodes_hsa([coder], graph,dict1, dict2)

		# if(len(graph.edges())>1):
		# 	graph=nc.simplifyNetwork(graph, data)
		graph = simpleNetBuild()
		#graph=liu.LiuNetwork1Builder()
		coder='unrolled'
		if(len(graph.edges())>1):
			print(coder)
			print(len(graph.edges()))
			
			nx.write_graphml(graph,coder+'.graphml')
			tempsensitivities=[[],[],[],[]]
			tempspecificities=[[],[],[],[]]
			truthlists=[[],[],[],[]]
			devLists=[]
			for i in range(bioReplicates):
				tuple1, tuple2, tuple3, tuple4, devList, inEdgeNums, overlaps, truthlisttemp=rewireSimTest(graph)
				sensitivity1, specificity1, nodesensitivity1, nodespecificity1 = tuple1
				sensitivity2, specificity2, nodesensitivity2, nodespecificity2 = tuple2
				sensitivity3, specificity3, nodesensitivity3, nodespecificity3 = tuple3
				sensitivity4, specificity4, nodesensitivity4, nodespecificity4 = tuple4
				tempsensitivities[0].append(sensitivity1)
				tempsensitivities[1].append(sensitivity2)
				tempsensitivities[2].append(sensitivity3)
				tempsensitivities[3].append(sensitivity4)
				tempspecificities[0].append(specificity1)
				tempspecificities[1].append(specificity2)
				tempspecificities[2].append(specificity3)
				tempspecificities[3].append(specificity4)
				nodesensitivities[0].extend(nodesensitivity1)
				nodesensitivities[1].extend(nodesensitivity2)
				nodesensitivities[2].extend(nodesensitivity3)
				nodesensitivities[3].extend(nodesensitivity4)
				nodespecificities[0].extend(nodespecificity1)
				nodespecificities[1].extend(nodespecificity2)
				nodespecificities[2].extend(nodespecificity3)
				nodespecificities[3].extend(nodespecificity4)
				devLists.append(devList)
				overlapNodes.extend(overlaps)
				edgeDegree.extend(inEdgeNums)
				truthlists[0].extend(truthlisttemp[0])
				truthlists[1].extend(truthlisttemp[1])
				truthlists[2].extend(truthlisttemp[2])
				truthlists[3].extend(truthlisttemp[3])
			for i in range(len(tempsensitivities)):
				tempsensitivities[i]=filter(lambda a: a != 100, tempsensitivities[i])
				if len(tempsensitivities[i])==0:
					tempsensitivities[i].append(100)
			for tempholder in tempspecificities:
				tempspecificities[i]=filter(lambda a: a != 100, tempspecificities[i])
				if len(tempspecificities[i])==0:
					tempspecificities[i].append(0.)
			sensitivity=[numpy.sum(tempsensitivities[0])/len(tempsensitivities[0]),numpy.sum(tempsensitivities[1])/len(tempsensitivities[1]),numpy.sum(tempsensitivities[2])/len(tempsensitivities[2]),numpy.sum(tempsensitivities[3])/len(tempsensitivities[3])]
			sensitivityStd=[numpy.std(tempsensitivities[0])/len(tempsensitivities[0]),numpy.std(tempsensitivities[1])/len(tempsensitivities[1]),numpy.std(tempsensitivities[2])/len(tempsensitivities[2]),numpy.std(tempsensitivities[3])/len(tempsensitivities[3])]
			specificity=[numpy.sum(tempspecificities[0])/len(tempspecificities[0]),numpy.sum(tempspecificities[1])/len(tempspecificities[1]),numpy.sum(tempspecificities[2])/len(tempspecificities[2]),numpy.sum(tempspecificities[3])/len(tempspecificities[3])]
			specificityStd=[numpy.std(tempspecificities[0])/len(tempspecificities[0]),numpy.std(tempspecificities[1])/len(tempspecificities[1]),numpy.std(tempspecificities[2])/len(tempspecificities[2]),numpy.std(tempspecificities[3])/len(tempspecificities[3])]
			truthholder.append(truthlists)
			specificities.append(specificity)
			sensitivities.append(sensitivity)
			specificityStds.append(specificityStd)
			sensitivityStds.append(sensitivityStd)
			print(sensitivity)
			print(sensitivityStd)
			print(specificity)
			print(specificityStd)
			print(nodesensitivities)
			print(nodespecificities)
			print(devLists)
	nodeLookup={}
	for number in edgeDegree:
		nodeLookup[number]=[[],[],[],[],[],[],[],[]]
	for i in range(len(edgeDegree)):
		nodeLookup[edgeDegree[i]][0].append(nodesensitivities[0][i])
		nodeLookup[edgeDegree[i]][1].append(nodesensitivities[1][i])
		nodeLookup[edgeDegree[i]][2].append(nodesensitivities[2][i])
		nodeLookup[edgeDegree[i]][3].append(nodesensitivities[3][i])
		nodeLookup[edgeDegree[i]][4].append(nodespecificities[0][i])
		nodeLookup[edgeDegree[i]][5].append(nodespecificities[1][i])
		nodeLookup[edgeDegree[i]][6].append(nodespecificities[2][i])
		nodeLookup[edgeDegree[i]][7].append(nodespecificities[3][i])
	overlapLookup={}
	for overlap in overlapNodes:
		overlapLookup[overlap]=[[],[],[],[],[],[],[],[]]
	for i in range(len(overlapNodes)):
		overlapLookup[overlapNodes[i]][0].append(nodesensitivities[0][i])
		overlapLookup[overlapNodes[i]][1].append(nodesensitivities[1][i])
		overlapLookup[overlapNodes[i]][2].append(nodesensitivities[2][i])
		overlapLookup[overlapNodes[i]][3].append(nodesensitivities[3][i])
		overlapLookup[overlapNodes[i]][4].append(nodespecificities[0][i])
		overlapLookup[overlapNodes[i]][5].append(nodespecificities[1][i])
		overlapLookup[overlapNodes[i]][6].append(nodespecificities[2][i])
		overlapLookup[overlapNodes[i]][7].append(nodespecificities[3][i])
	finalNodeData=[]
	finalExtendData=[]
	for key in nodeLookup.keys():
		templisting=[]
		tempExtended=[]
		for lister in nodeLookup[key]:
			newlist=filter(lambda a: a != 100, lister)
			tempExtended.append(newlist)
			print("number of nodes by key")
			print(len(newlist))
			if len(newlist)==0:
				templisting.append(0.)
			else:
				templisting.append(sum(newlist)/len(newlist))
		finalNodeData.append(templisting)
		finalExtendData.append(tempExtended)
	finalOverlapExtendData=[]
	for key in overlapLookup.keys():
		templisting=[]
		tempOverlap=[]
		for lister in overlapLookup[key]:
			newlist=filter(lambda a: a != 100, lister)
			tempOverlap.append(newlist)
			print("number of overlaps by key")
			print(len(newlist))
		finalOverlapExtendData.append(tempOverlap)
	print(nodeLookup.keys())
	print(overlapLookup.keys())
	print(finalNodeData)
	print(finalExtendData)
	pickle.dump( finalNodeData, open( "node_by_node_data.pickle", "wb" ) )
	pickle.dump( finalExtendData, open( "extended_node_by_node_data.pickle", "wb" ) )
	pickle.dump( finalOverlapExtendData, open( "extended_overlap_data.pickle", "wb" ) )
	
	pickle.dump( devLists, open( "devLists.pickle", "wb" ) )

	pickle.dump( [sensitivities,sensitivityStds,specificities, specificityStds], open( "network_by_network_data.pickle", "wb" ) )
	print(sensitivities)
	print(specificities)
	print(finalNodeData)
	pickle.dump( truthholder, open( "expt_true_corrected_bits.pickle", "wb" ) )
def rewireSimTest(graph):
	params=sim.paramClass()
	if params.rewire:
		graph=nc.rewireNetwork(graph)

	sss=utils.synthesizeInputs(graph,params.samples)
	model=sim.modelClass(graph,sss)
	#print(model.nodeList)
	#print(model.size)
	#print(model.individualParse)
	return simTester(model, sss,'prop')

if __name__ == '__main__':
	import time
	start_time = time.time()
	ifngStimTest(1)
	print("--- %s seconds ---" % (time.time() - start_time))
