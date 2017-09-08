
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
import argparse as argparse
import scipy.stats as stat
from sets import Set
import matplotlib.pyplot as plt
import requests

import sklearn.mixture as mix
# import other pieces of our software
import networkConstructor as nc
import utils as utils
import simulation as sim
import murphyReader as mr
import motif_cutter as mc

def genBits(model):
	startInd=list(utils.genRandBits(model.size))
	counter=0
	while numpy.sum(startInd)==0 and counter < 10000:
		startInd=list(utils.genRandBits(model.size))
		counter+=1
	for node in range(0,len(model.nodeList)):
		if node==(len(model.nodeList)-1):
			end= model.size
		else:
			end=model.individualParse[node+1]
		start=model.individualParse[node]
		if (end-start)>1:
			counter=0
			while numpy.sum(startInd[start:end])>5 and counter < 10000:
				chosen=math.floor(random()*(end-start))
				startInd[start+int(chosen)]=0
				counter+=1
			if numpy.sum(startInd[start:end])==0:
				chosen=math.floor(random()*(end-start))
				startInd[start+int(chosen)]=1
		elif (end-start)==1:
			startInd[start]=1
	return startInd

def cxTwoPointNode(ind1, ind2, model):
	"""Executes a two-point crossover on the input :term:`sequence`
    individuals. The two individuals are modified in place and both keep
    their original length. 
    
    :param ind1: The first individual participating in the crossover.
    :param ind2: The second individual participating in the crossover.
    :returns: A tuple of two individuals.

    This function uses the :func:`~random.randint` function from the Python 
    base :mod:`random` module.
    """
	size = len(model.nodeList)
	cxpoint1 = rand.randint(1, size)
	cxpoint2 = rand.randint(1, size - 1)
	if cxpoint2 >= cxpoint1:
		cxpoint2 += 1
	else: # Swap the two cx points
		cxpoint1, cxpoint2 = cxpoint2, cxpoint1
	errors1=list(ind1.fitness.values)
	errors2=list(ind2.fitness.values)
	sum1= numpy.sum(errors1[:cxpoint1])+numpy.sum(errors2[cxpoint1:cxpoint2])+numpy.sum(errors1[cxpoint1:])
	sum2= numpy.sum(errors2[:cxpoint1])+numpy.sum(errors1[cxpoint1:cxpoint2])+numpy.sum(errors2[cxpoint1:])
	cxpoint1=model.individualParse[cxpoint1]
	cxpoint2=model.individualParse[cxpoint2]
	ind1[cxpoint1:cxpoint2], ind2[cxpoint1:cxpoint2] = ind2[cxpoint1:cxpoint2], ind1[cxpoint1:cxpoint2]
	if sum1 < sum2:
		return ind1, ind2
	else:
		return ind2, ind1
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
		for i in range(len(model.nodeList)-1):
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
	toolbox.register("select", tools.selNSGA2)
	toolbox.register("similar", numpy.array_equal)
	return toolbox, stats
			
def evaluate(individual, cells, model,  sss, params, KOlist, KIlist):
	SSEs=[]
	for j in range(0,len(sss)):
		boolValues=iterateBooleanModel(individual, model, cells, model.initValueList[j], params,KOlist[i], KIlist[i])
		SSE=numpy.sum([(boolValues[i]-sss[j][model.nodeList[i]])**2 for i in range(0, len(model.nodeList))])
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
	
def evaluateByNode(individual, cells, model,  sss, params, KOlist, KIlist):
	# SSEs=[]
	#boolValues=Parallel(n_jobs=min(6,len(sss)))(delayed(iterateBooleanModel)(list(individual), model, cells, model.initValueList[i], params) for i in range(len(sss)))
	boolValues=[iterateBooleanModel(list(individual), model, cells, model.initValueList[i], params, KOlist[i], KIlist[i]) for i in range(len(sss))]
	SSEs= [numpy.sum([(boolValues[j][i]-sss[j][model.nodeList[i]])**2 for j in range(0,len(sss))]) for i in range(0, len(model.nodeList))]
	return tuple(SSEs)
# generates a random set of samples made up of cells by using parameteris from probInit seq
# to set up then iterating using strict Boolean modeling. 
def runProbabilityBooleanSims(individual, model, sampleNum, cells, params, KOlist, KIlist):
	samples=[]
	seeds=[]
	for i in range(0,sampleNum):
		seeds.append(random())
	samples=Parallel(n_jobs=min(24,sampleNum))(delayed(sampler)(individual, model, sampleNum, seeds[i], params, KOlist[i], KIlist[i]) for i in range(sampleNum))
	# counter=0
	# for sample in range(len(samples)):
	# 	if sum(samples(sample))==0:
	# 		temp=samples.pop(sample) 
	# newSamples=Parallel(n_jobs=min(8,sampleNum))(delayed(sampler)(individual, model, sampleNum, seeds[i]) for i in range(sampleNum))
	# print(samples)
	return samples
	
def sampler(individual, model, cells, seeder, params, KOs, KIs):
	seed(seeder)
	cellArray=[]
	sampleProbs=[]
	simulator=sim.simulatorClass('bool')
	# generate random proportions for each node to start
	for j in range(0,len(model.nodeList)):
		sampleProbs.append(random())
	return iterateBooleanModel(individual, model, cells, sampleProbs, params, KOs, KIs)

def iterateBooleanModel(individual, model, cells, initProbs, params, KOs, KIs):
	cellArray=[]
	simulator=sim.simulatorClass('bool')
	simulator.simSteps=3*len(model.nodeList)
	for j in range(0,cells):
		# shuffle nodes to be initially called.... 
		#simulations that are truly random starting states should have all nodes added to this list
		#get initial values for all nodes
		initValues=genPBNInitValues(individual, model,initProbs)
		# run Boolean simulation with initial values and append
		vals=sim.runModel(individual, model, simulator, initValues, params, KOs, KIs)
		cellArray.append(vals)
	return [1.*numpy.sum(col) / float(cells) for col in zip(*cellArray)]

def genPBNInitValues(individual, model,sampleProbs):
	#return [True if (random()<sampleProbs[node]) else False for node in range(0,len(sampleProbs))]
	initValues=[False for x in range(0,len(model.nodeList))]
	for node in range(0,len(sampleProbs)):
		if random()<sampleProbs[node]:
			initValues[node]=True
	return initValues

def varOrAdaptive(population, toolbox, model, lambda_, cxpb, mutpb, genfrac):
	# algorithm for generating a list of offspring... copied and pasted from DEAP with modification for adaptive mutation
	assert (cxpb + mutpb) <= 1.0, ("The sum of the crossover and mutation "
		"probabilities must be smaller or equal to 1.0.")
	offspring = []
	for _ in xrange(lambda_):
		op_choice = random()
		if op_choice < cxpb:            # Apply crossover
			ind1, ind2 = map(toolbox.clone, rand.sample(population, 2))
			ind1, ind2 = cxTwoPointNode(ind1, ind2, model)
			del ind1.fitness.values
			offspring.append(ind1)
		elif op_choice < cxpb + mutpb:  # Apply mutation
			ind = toolbox.clone(rand.choice(population))
			ind, = mutFlipBitAdapt(ind,  .5, model, genfrac)
			del ind.fitness.values
			offspring.append(ind)
		else:                           # Apply reproduction
			offspring.append(rand.choice(population))
	return offspring

def mutFlipBitAdapt(individual, indpb, model, genfrac):
	# get errors
	errors=list(individual.fitness.values)

	# get rid of errors in nodes that can't be changed
	for j in xrange(len(errors)):
		if model.andLenList[j]<2:
			errors[j]=0
	# print(errors)
	# if error is relatively low, do a totally random mutation
	if numpy.sum(errors)<.05*len(model.nodeList):
		focusNode=int(math.floor(random()*len(model.andLenList)))
	else:
		# if errors are relatively high, focus on nodes that fit the worst
		normerrors=[1.*error/numpy.sum(errors) for error in errors]# normalize errors to get a probability that the node  is modified
		probs=numpy.cumsum(normerrors)
		# print(probs)
		randy=random()# randomly select a node to mutate
		focusNode=next(i for i in range(len(probs)) if probs[i]>randy)
	# print(focusNode)
	if model.andLenList[focusNode]>1:
		# find ends of the node of interest in the individual
		start=model.individualParse[focusNode]
		if focusNode==len(model.nodeList)-1:
			end= model.size
		else:
			end=model.individualParse[focusNode+1]
		for i in range(start,end):
			if random()< 2/(end-start+1):
				individual[i] = 1
			else:
				individual[i] = 0
		#ensure that there is at least one shadow and node turned on
		if numpy.sum(individual[start:end])==0:
			individual[start]=1
	return individual,

def selNSGA2(individuals, k):
	"""Apply NSGA-II selection operator on the *individuals*. Usually, the
	size of *individuals* will be larger than *k* because any individual
	present in *individuals* will appear in the returned list at most once.
	Having the size of *individuals* equals to *k* will have no effect other
	than sorting the population according to their front rank. The
	list returned contains references to the input *individuals*. For more
	details on the NSGA-II operator see [Deb2002]_.
	
	:param individuals: A list of individuals to select from.
	:param k: The number of individuals to select.
	:returns: A list of selected individuals.
	
	.. [Deb2002] Deb, Pratab, Agarwal, and Meyarivan, "A fast elitist
	   non-dominated sorting genetic algorithm for multi-objective
	   optimization: NSGA-II", 2002.
	"""
	pareto_fronts = sortNondominatedAdapt(individuals, k)
	for front in pareto_fronts:
		assignCrowdingDist(front)
	
	chosen = list(chain(*pareto_fronts[:-1]))
	k = k - len(chosen)
	if k > 0:
		sorted_front = sorted(pareto_fronts[-1], key=attrgetter("fitness.crowding_dist"), reverse=True)
		chosen.extend(sorted_front[:k])
		
	return chosen

def sortNondominatedAdapt(individuals, k, first_front_only=False):
	"""Sort the first *k* *individuals* into different nondomination levels 
	using the "Fast Nondominated Sorting Approach" proposed by Deb et al.,
	see [Deb2002]_. This algorithm has a time complexity of :math:`O(MN^2)`, 
	where :math:`M` is the number of objectives and :math:`N` the number of 
	individuals.
	
	:param individuals: A list of individuals to select from.
	:param k: The number of individuals to select.
	:param first_front_only: If :obj:`True` sort only the first front and
							 exit.
	:returns: A list of Pareto fronts (lists), the first list includes 
			  nondominated individuals.

	.. [Deb2002] Deb, Pratab, Agarwal, and Meyarivan, "A fast elitist
	   non-dominated sorting genetic algorithm for multi-objective
	   optimization: NSGA-II", 2002.
	"""
	if k == 0:
		return []

	map_fit_ind = defaultdict(list)
	for ind in individuals:
		map_fit_ind[ind.fitness].append(ind)
	fits = map_fit_ind.keys()
	
	current_front = []
	next_front = []
	dominating_fits = defaultdict(int)
	dominated_fits = defaultdict(list)
	
	# Rank first Pareto front
	for i, fit_i in enumerate(fits):
		for fit_j in fits[i+1:]:
			if dominated(fit_i, fit_j):
				dominating_fits[fit_j] += 1
				dominated_fits[fit_i].append(fit_j)
			elif dominated(fit_j, fit_i):
				dominating_fits[fit_i] += 1
				dominated_fits[fit_j].append(fit_i)
		if dominating_fits[fit_i] == 0:
			current_front.append(fit_i)
	
	fronts = [[]]
	for fit in current_front:
		fronts[-1].extend(map_fit_ind[fit])
	pareto_sorted = len(fronts[-1])

	# Rank the next front until all individuals are sorted or 
	# the given number of individual are sorted.
	if not first_front_only:
		N = min(len(individuals), k)
		while pareto_sorted < N:
			fronts.append([])
			for fit_p in current_front:
				for fit_d in dominated_fits[fit_p]:
					dominating_fits[fit_d] -= 1
					if dominating_fits[fit_d] == 0:
						next_front.append(fit_d)
						pareto_sorted += len(map_fit_ind[fit_d])
						fronts[-1].extend(map_fit_ind[fit_d])
			current_front = next_front
			next_front = []
	
	return fronts

def dominated(ind1, ind2):
	"""Return true if each objective of *self* is not strictly worse than 
		the corresponding objective of *other* and at least one objective is 
		strictly better.

		:param obj: Slice indicating on which objectives the domination is 
					tested. The default value is `slice(None)`, representing
					every objectives.
	"""
	not_equal = False
	for self_wvalue, other_wvalue in zip(ind1.wvalues[obj], ind2.wvalues[obj]):
		if self_wvalue > other_wvalue + .2:
			not_equal = True
		elif self_wvalue < other_wvalue:
			return False                
	return not_equal

def assignCrowdingDist(individuals):
	"""Assign a crowding distance to each individual's fitness. The 
	crowding distance can be retrieve via the :attr:`crowding_dist` 
	attribute of each individual's fitness.
	"""
	if len(individuals) == 0:
		return
	
	distances = [0.0] * len(individuals)
	crowd = [(ind.fitness.values, i) for i, ind in enumerate(individuals)]
	
	nobj = len(individuals[0].fitness.values)
	
	for i in xrange(nobj):
		crowd.sort(key=lambda element: element[0][i])
		distances[crowd[0][1]] = float("inf")
		distances[crowd[-1][1]] = float("inf")
		if crowd[-1][0][i] == crowd[0][0][i]:
			continue
		norm = nobj * float(crowd[-1][0][i] - crowd[0][0][i])
		for prev, cur, next in zip(crowd[:-2], crowd[1:-1], crowd[2:]):
			distances[cur[1]] += 1.*(next[0][i] - prev[0][i]) / norm

	for i, dist in enumerate(distances):
		individuals[i].fitness.crowding_dist = dist

def eaMuPlusLambdaAdaptive(population, toolbox, model, mu, lambda_, cxpb, mutpb, ngen, stats=None, halloffame=None, verbose=__debug__):
	logbook = tools.Logbook()
	logbook.header = ['gen', 'nevals'] + (stats.fields if stats else [])

	# Evaluate the individuals with an invalid fitness
	invalid_ind = [ind for ind in population if not ind.fitness.valid]
	# fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
	fitnesses=Parallel(n_jobs=min(24,len(invalid_ind)))(delayed(toolbox.evaluate)(list(indy)) for indy in invalid_ind)

	for ind, fit in zip(invalid_ind, fitnesses):
		ind.fitness.values = fit

	if halloffame is not None:
		halloffame.update(population)

	record = stats.compile(population) if stats is not None else {}
	logbook.record(gen=0, nevals=len(invalid_ind), **record)
	if verbose:
		print(logbook.stream)

	breaker=False
	for ind in population:
		if numpy.sum(ind.fitness.values)< .01*len(ind.fitness.values):
			breaker=True
	if breaker:
		return population, logbook
	# Begin the generational process
	for gen in range(1, ngen+1):
		# Vary the population
		offspring = varOrAdaptive(population, toolbox, model, lambda_, cxpb, mutpb, (1.*gen/ngen))
		
		# Evaluate the individuals with an invalid fitness
		invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
		#fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
		fitnesses=Parallel(n_jobs=min(24,len(invalid_ind)))(delayed(toolbox.evaluate)(list(indy)) for indy in invalid_ind)
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
			print(logbook.stream)
		breaker=False
		for ind in population:
			if numpy.sum(ind.fitness.values)< .01*len(ind.fitness.values):
				breaker=True
				saveInd=ind
		if breaker:
			errorTemp=saveInd.fitness.values
			for value in errorTemp:
				if value> .1:
					breaker=False
		if breaker:
			break
	return population, logbook

def GAautoSolver(model, sss, params, KOlist, KIlist):
	# set up toolbox and run GA with or without adaptive mutations turned on
	toolbox, stats=buildToolbox(model.size,params.bitFlipProb, model)

	if params.adaptive:
		toolbox.register("evaluate", evaluateByNode, cells=params.cells,model=model,sss=sss, params=params, KOlist=KOlist, KIlist=KIlist)
	else:
		toolbox.register("evaluate", evaluate, cells=params.cells,model=model,sss=sss, KOlist=KOlist, KIlist=KIlist)
	population=toolbox.population(n=params.popSize)
	hof = tools.HallOfFame(params.hofSize, similar=numpy.array_equal)
	
	if params.adaptive:
		output=eaMuPlusLambdaAdaptive(population, toolbox, model, mu=params.mu, lambda_=params.lambd, stats=stats, cxpb=params.crossoverProb, mutpb=params.mutationProb, ngen=params.generations, verbose=params.verbose, halloffame=hof)
	else:
		output=algo.eaMuCommaLambda(population, toolbox, mu=params.mu, lambda_=params.lambd, stats=stats, cxpb=params.crossoverProb, mutpb=params.mutationProb, ngen=params.generations, verbose=params.verbose, halloffame=hof)
	return output

def GAsearchModel(model, sss,params, KOlist, KIlist):

	newSSS=[]
	newSSS.append({})
	newSSS.append({})
	newSSS.append({})
	newSSS.append({})
	newSSS.append({})
	newSSS.append({})
	for key in sss[0].keys():
		newSSS[0][key]=numpy.mean([sss[i][key] for i in range(4)])
		newSSS[1][key]=numpy.mean([sss[4+i][key] for i in range(3)])
		newSSS[2][key]=numpy.mean([sss[7+i][key] for i in range(4)])
		newSSS[3][key]=numpy.mean([sss[11+i][key] for i in range(4)])
		newSSS[4][key]=numpy.mean([sss[15+i][key] for i in range(4)])
		newSSS[5][key]=numpy.mean([sss[19+i][key] for i in range(4)])
	newInitValueList=[]
	knockoutLists=[]
	knockinLists=[]
	for q in range(len(sss)):
		knockoutLists.append([])
		knockinLists.append([])
	for j in range(0,len(sss)):
		newInitValueList.append([])
	for j in range(0,len(model.nodeList)):
		for k in range(0,len(sss)):
			ss=sss[k]
			if  model.nodeList[j] in sss[0].keys():
				newInitValueList[k].append(ss[model.nodeList[j]])
			else:
				newInitValueList[k].append(0.5)
	
	newInitValueList=[]
	knockoutLists=[]
	knockinLists=[]
	for q in range(len(sss)):
		knockoutLists.append([])
		knockinLists.append([])
	for j in range(0,len(sss)):
		newInitValueList.append([])
	for j in range(0,len(model.nodeList)):
		for k in range(0,len(sss)):
			ss=sss[k]
			if  model.nodeList[j] in sss[0].keys():
				newInitValueList[k].append(ss[model.nodeList[j]])
			else:
				newInitValueList[k].append(0.5)
	model.initValueList=newInitValueList
	population, logbook=GAautoSolver(model, newSSS, params, knockoutLists, knockinLists)
	minny=1000000
	saveVal=-1
	for i in range(len(population)):
		if numpy.sum(population[i].fitness.values)< minny:
			minny=numpy.sum(population[i].fitness.values)
			saveVal=i
	ultimate=list(population[saveVal])
	minvals=population[saveVal].fitness.values
	# get rid of redundant edges
	newultimate=[]
	for node in range(0,len(model.nodeList)):
		#get start and end indices for node in individual
		if node==(len(model.nodeList)-1):
			end=len(ultimate)
		else:
			end=model.individualParse[node+1]
		start=model.individualParse[node]
		# get all the in edges for each and node
		andNodeList=model.andNodeList[node]
		inEdges=[]
		for lister in andNodeList:
			inEdges.append(set(lister))
		truth=ultimate[start:end]
		# check if any nodes are redundant
		for i in range(len(truth)):
			if truth[i]==1:
				for j in range(len(truth)):
					if truth[j]==1 and not i==j:
						if inEdges[i].issubset(inEdges[j]):
							truth[j]=0
		newultimate.extend(truth)
	ultimate=newultimate
	for node in range(0,len(model.nodeList)):
		#get start and end indices for node in individual
		if node==(len(model.nodeList)-1):
			end=len(ultimate)
		else:
			end=model.individualParse[node+1]
		start=model.individualParse[node]
		# get all the in edges for each and node
		andNodeList=model.andNodeList[node]
		inEdges=[]
		for lister in andNodeList:
			inEdges.append(set(lister))
		truth=ultimate[start:end]
		invalid_ind=[]
		# generate a list of possibilities
		for i in range(1,2**(end-start)):
			tempultimate=list(newultimate)
			tempInd=utils.bitList(i, len(truth))
			tempultimate[start:end]=tempInd
			invalid_ind.append(tempultimate)
		# see which rule is the best
		if end-start>1:
			fitnesses=Parallel(n_jobs=min(24,len(invalid_ind)))(delayed(evaluateByNode)(indy, params.cells, model,  newSSS, params, KOlist, KIlist) for indy in invalid_ind)
			minny=1000
			mini=100
			for i in range(len(fitnesses)):
				currentsum= numpy.sum(fitnesses[i])
				if currentsum< minny:
					minny=currentsum
					mini=i
			newultimate=invalid_ind[mini]
	ultimate=newultimate

	return minvals, ultimate

def deltaTester(ultimate, i, model,  newSSS, params, minny):
	print(i)
	modifyFlag=False
	copied=list(ultimate)
	copied[i]=1-copied[i]
	for node in range(0,len(model.nodeList)):
		better=True
		#get start and end indices for node in individual
		if node==(len(model.nodeList)-1):
			end=len(ultimate)-1
		else:
			end=model.individualParse[node+1]
		start=model.individualParse[node]
		if model.andLenList[node]>0:
			if sum(copied[start:end])>0:
				newtot=numpy.sum(evaluateByNode(copied, params.cells, model,  newSSS, params, KOlist, KIlist))
				if newtot<minny:
					modifyFlag=True
	return modifyFlag

def simTester(graph, name):
	#creates a model, runs simulations, then tests reverse engineering capabilities of models in a single function
	#trials is the number of different individuals to try
	#samples is the number of different initial conditions to provide per trial
	#graph specifies the network we are testing. 

	# set up params, sss, model, samples


	params=sim.paramClass()
	# rewire graph if necessary
	if params.rewire:
		graph=nc.rewireNetwork(graph)

	# # remove two node bunches
	# i=0
	# while i < range(len(graph.nodes())):
	# 	nodes=graph.nodes()
	# 	node=nodes[i]
	# 	if len(graph.predecessors(node))==0:
	# 		successorless=True
	# 		for successor in graph.successors(node):
	# 			if graph.successors(successor)>0 or graph.predecessors(successor)>1:
	# 				successorless=False
	# 		if successorless:
	# 			graph.remove_node(node)
	# 			i=i-1
	# 	i+=1

	# # remove orphan nodes
	# for i in range(len(graph.nodes())):
	# 	nodes=graph.nodes()
	# 	node=nodes[i]
	# 	if len(graph.predecessors(node))==0 and len(graph.successors(node))==0:
	# 		graph.remove_node(node)



	sss=utils.synthesizeInputs(graph,params.samples)
	model=sim.modelClass(graph,sss)
	samples=params.samples

	# set up necessary lists for output
	truthList= [] 
	testList =[]
	devList=[]


	knockoutLists=[]
	knockinLists=[]

	for q in range(samples):
		temRand=rand.randint(0,len(model.nodeList))
		knockoutLists.append([])
		knockinLists.append([])


	# loop over number of times we want to generate fake data and perform sequence of events
	for i in range(0,params.trials):
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
		output=runProbabilityBooleanSims(individual, model, samples, params.cells, params, knockoutLists, knockinLists)
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
		for j in range(0,len(model.nodeList)):
			for k in range(0,len(sss)):
				ss=newSSS[k]
				if  model.nodeList[j] in sss[0].keys():
					newInitValueList[k].append(ss[model.nodeList[j]])
				else:
					newInitValueList[k].append(0.5)
		model.initValueList=newInitValueList
		#perform GA run
		dev, bruteOut=GAsearchModel(model, newSSS, params, knockoutLists, knockinLists)
		# add output of this trial to lists
		truthList.append(individual)
		testList.append(list(bruteOut))
		devList.append(dev)
	# set up output and save as a pickle
	outputList=[truthList,testList, devList,model.size, model.nodeList, model.individualParse,model.andNodeList ,model.andNodeInvertList, model.andLenList,	model.nodeList, model.nodeDict, model.initValueList]
	pickle.dump( outputList, open( name+"_output.pickle", "wb" ) )

def partitionTester(graph,name):
	# first trim the graph to just what is in the dataset and the dataset to just what is in the graph

	if len(graph.nodes())>1:	
		# dev1, bruteOut1, model1= partitionTest(graph, False, False) # complicated method with max included
		# dev2, bruteOut2, model2= partitionTest(graph, True, False) # complicated method after throwing out max
		dev3, bruteOut3, model3= partitionTest(graph, False, True) # divide by max
		dev4, bruteOut4, model4= partitionTest(graph, True, True) # divide by second highest
		# set up output and save as a pickle
		# outputList=[[dev1,dev2,dev3,dev4],[bruteOut1,bruteOut2,bruteOut3, bruteOut4],[model1,model2,model3,model4]]
		outputList=[[dev3,dev4],[bruteOut3,bruteOut4],[model3,model4]]

		pickle.dump( outputList, open( name+"_output.pickle", "wb" ) )
	else:
		print('not enough overlap')

def indDistFind(name):
	ssDict=pickle.Unpickler(open( 'data_F_T.pickle', "rb" )).load()
	[[dev1,dev3],[bruteOut1,bruteOut3],[model1,model3]]=pickle.Unpickler(open( 'pickles/'+name+"_output.pickle", "rb" )).load()
	# generate SSS which is a list of SSs where SS is a dict pointing to some values
	keyList=ssDict.keys()
	sss=[{} for i in range(len(ssDict[keyList[0]]))]
	newInitValueList=[[] for i in range(len(ssDict[keyList[0]]))]
	for i in range(len(sss)):
		for key in keyList:
			if key in graph.nodes():
				sss[i][key]=ssDict[key][i]
	params=sim.paramClass()

	simulator=sim.simulatorClass('bool')
	importanceScores=sim.calcImportance(bruteOut1,params,model1,simulator, sss)

	pickle.dump( importanceScores, open( name+"_scores.pickle", "wb" ) )

def partitionTest(graph, maxRem, divider):
	
	# first construct the input -omics dataset
	if maxRem:
		if divider:
			fileName='data_T_T.pickle'
		else:
			fileName='data_T_F.pickle'
	else:
		if divider:
			fileName='data_F_T.pickle'
		else:
			fileName='data_F_F.pickle'
	ssDict=pickle.Unpickler(open( fileName, "rb" )).load()
	# ssDict, valueLister, ones, zeros= mr.constructOmicsInput(geneDict, maxRem, divider)
	params=sim.paramClass()
	
	graph=nc.simplifyNetwork(graph, ssDict) # simplify the graph
	graph=mc.cutMotifs(graph)
	print(graph.nodes())
	keyList=ssDict.keys()
	# generate SSS which is a list of SSs where SS is a dict pointing to some values
	sss=[{} for i in range(len(ssDict[keyList[0]]))]
	newInitValueList=[[] for i in range(len(ssDict[keyList[0]]))]
	for i in range(len(sss)):
		for key in keyList:
			if key in graph.nodes():
				sss[i][key]=ssDict[key][i]
	
	model=sim.modelClass(graph,sss)
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
			if  model.nodeList[j] in sss[0].keys():
				newInitValueList[k].append(ss[model.nodeList[j]])
			else:
				newInitValueList[k].append(0.5)
	model.initValueList=newInitValueList
	print('setup successful')
	#perform GA run
	dev, bruteOut= GAsearchModel(model, sss, params, knockoutLists, knockinLists)
	return dev, bruteOut, model

def read_gmt(filename):
	gmt_dict={}
	inputfile = open(filename, 'r')
	lines = inputfile.readlines()
	for line in lines:
		newline=line.split('\t')
		gmt_dict[newline[0]]=Set(newline[2:])
	return gmt_dict

def find_overlaps(filename,geneDict):
	overlapsets=[]
	genes=Set(geneDict.keys())
	keggDict=read_gmt(filename)
	for key in keggDict.keys():
		if len(genes.intersection(keggDict[key]))>4:
			overlapsets.append(key)
			print(key)
			print(len(genes.intersection(keggDict[key])))
	return overlapsets

def retrieveGraph(name,aliasDict,dict1,dict2, geneDict):
	print(name)
	namelist=name.split('_')
	namelist.pop(0)
	requester='http://rest.kegg.jp/find/pathway/'+namelist.pop(0)
	for item in namelist:
		requester=requester+'+'+item

	# print(requester)

	r=requests.get(requester)
	genes=Set(geneDict.keys())

	lines=r.text
	print(lines)
	if len(lines.split('\n')[0].split(':'))>1:

		code=lines.split('\n')[0].split(':')[1][3:8] # KEGG number of overlapped pathway
		graph=nx.DiGraph()
		coder=str('ko'+code) 
		nc.uploadKEGGcodes([coder], graph, dict2)
		coder=str('hsa'+code)
		nc.uploadKEGGcodes_hsa([coder], graph,dict1, dict2)
		if len(list(nx.connected_component_subgraphs(graph.to_undirected() )))>0:
			newgraph = max(nx.connected_component_subgraphs(graph.to_undirected()), key=len)
			newOverlap=genes.intersection(set(newgraph.nodes()))
			print(len(newOverlap))
			graphLen= nx.average_shortest_path_length(newgraph, weight=None)
			graph=nc.simplifyNetwork(graph, geneDict)
			print('nodes, edges')
			print(len(graph.nodes()))
			print(len(graph.edges()))
			if len(newOverlap)>4:
				nx.write_graphml(graph,coder+'.graphml')
				nx.write_gpickle(graph,coder+'.gpickle')
		else:
			graphLen=0
		print('graphlen:')
		print(graphLen)
		#check is average path length is sufficiently long


		
	else:
		print('not found:')
		print(requester)
		print(lines)


def findPathways():
	geneDict=mr.readData() # get -omics data
	
	aliasDict={}
	dict1={}
	nc.parseKEGGdicthsa('inputData/hsa00001.keg',aliasDict,dict1)
	dict2={}
	nc.parseKEGGdict('inputData/ko00001.keg',aliasDict,dict2)
	namelist=find_overlaps('filtered.c2.cp.kegg.v3.0.symbols.gmt',geneDict)
	print('num of overlap nodes')
	print(len(namelist))
	for name in namelist:
		retrieveGraph(name,aliasDict,dict1,dict2, geneDict)
	# # read in list of codes then load them into network
	# #inputfile = open('ID_filtered.c2.cp.kegg.v3.0.symbols.txt', 'r')
	# for code in lines:
		
	# 	# if(len(graph.edges())>1):
	# 	# 	graph=simplifyNetwork(graph, data)
		
	# 	#graph = utils.simpleNetBuild()
	# 	#coder='unrolled'

	# 	#graph=utils.LiuNetwork1Builder()
	# 	# coder='liu'

	# 	# if the code has an interesting logic rule to find, run the algorithm... 
	# 	checker=False
	# 	for x in graph.in_degree():
	# 		if x>1:
	# 			checker=True
	# 	if(checker):
	# 		print(coder)
	# 		codelist.append(coder)			
	# 		nx.write_graphml(graph,coder+'_unsimplified.graphml')
	# 		nx.write_gpickle(graph,coder+'_unsimplified.gpickle')


def findClustDiff(lister):
	mixer= mix.GaussianMixture(n_components=2, covariance_type='full').fit([[a] for a in lister])
	return abs(mixer.means_[0]-mixer.means_[1])


def testDiscretizationSetup():

	findPathways()

	# geneDict=mr.readData() # get -omics data

	# ssDict, valueLister, ones, zeros= mr.constructOmicsInput(geneDict, True, True)
	# pickle.dump( ssDict, open( "data_T_T.pickle", "wb" ) )


	# ssDict, valueLister, ones, zeros= mr.constructOmicsInput(geneDict, False, True)
	# pickle.dump( ssDict, open( "data_F_T.pickle", "wb" ) )

	# ssDict, valueLister, ones, zeros= mr.constructOmicsInput(geneDict, True, False)
	# pickle.dump( ssDict, open( "data_T_F.pickle", "wb" ) )


	# ssDict, valueLister, ones, zeros= mr.constructOmicsInput(geneDict, False, False)
	# pickle.dump( ssDict, open( "data_F_F.pickle", "wb" ) )


if __name__ == '__main__':
	import time
	start_time = time.time()
	parser = argparse.ArgumentParser()
	parser.add_argument("graph")
	results = parser.parse_args()
	graphName=results.graph
	name=graphName[:-8]+'_1'

	graph = nx.read_gpickle(graphName)
	indDistFind(name)


	# testDiscretizationSetup()


	# for numgraph in ['04066','04350','04380','04612','04630','04650','04657','04659','04660']:
	# valueList=[]
	# for key in geneDict.keys():
	# 	print(key)
	# 	valueList.append(findClustDiff(geneDict[key]))
	
	# pickle.dump( valueList, open( "distance_data.pickle", "wb" ) )
	# valueList=pickle.Unpickler(open( "distance_data.pickle", "rb" )).load()
	# print(len(valueList))
	# valueList=[valueList[i] for i in range(0,len(valueList),50)]
	# print(len(valueList))
	# bins = numpy.linspace(0,6000, 500)
	# n, bins, patches = plt.hist(valueList, 20, range=(0,20),  normed=0, facecolor='green')

	# plt.xlabel('distance between clusters')
	# plt.ylabel('Probability Density')
	# plt.title('Histogram of cluster distances')
	# # plt.axis([40, 160, 0, 0.03])
	# plt.grid(True)

	# plt.savefig('histogram_cluster_dist_1.png', bbox_inches='tight')




	# for numgraph in ['04657','04659','04660']:
	# 	graphName='hsa'+numgraph+'_unsimplified.gpickle'
	# 	name=graphName[:-8]+'_1'
	# 	graph = nx.read_gpickle(graphName)
	# 	partitionTester(graph, name, geneDict)
	# 	print("--- %s seconds ---" % (time.time() - start_time))
