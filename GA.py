
# import python packages
import pickle
import copy as copy
from deap import base, creator, gp, tools
from deap import algorithms as algo
from random import random, seed, shuffle, randint, sample, choice
import numpy as numpy
import operator
import math as math
from sets import Set
from joblib import Parallel, delayed
import gc as gc
from collections import defaultdict, deque
from itertools import chain
from operator import attrgetter, itemgetter
# import other pieces of our software
import utils as utils
from simulation import paramClass, EBNbool

# generates random bitstring with at least one value for each node
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
	return [copy.deepcopy(model),startInd]

# finds the lowest error individual in a population
def findPopBest(population):
	saveVal=-1
	minny=100000
	for i in range(len(population)):
		if numpy.sum(population[i].fitness.values)< minny:
			minny=numpy.sum(population[i].fitness.values)
			saveVal=i
	ultimate=population[saveVal]
	minvals=population[saveVal].fitness.values
	# get rid of redundant edges
	return minvals, ultimate[1], ultimate[0]

# executes two point crossover at node junctions
def cxTwoPointNode(ind1, ind2):
	"""Executes a two-point crossover on the input :term:`sequence`
    individuals. The two individuals are modified in place and both keep
    their original length. 
    
    :param ind1: The first individual participating in the crossover.
    :param ind2: The second individual participating in the crossover.
    :returns: A tuple of two individuals.

    This function uses the :func:`~random.randint` function from the Python 
    base :mod:`random` module.
    """
	size = len(ind1[0].nodeList)
	cxpointer1 = randint(1, size)
	cxpointer2 = randint(1, size - 1)
	if cxpointer2 >= cxpointer1:
		cxpointer2 += 1
	else: # Swap the two cx points
		cxpointer1, cxpointer2 = cxpointer2, cxpointer1
	cxpoint1=ind1[0].individualParse[cxpointer1]
	cxpoint2=ind1[0].individualParse[cxpointer2]
	ind1[1][cxpoint1:cxpoint2], ind2[1][cxpoint1:cxpoint2] = ind2[1][cxpoint1:cxpoint2], ind1[1][cxpoint1:cxpoint2]
	ind1[0].andNodeList[cxpointer1:cxpointer2], ind2[0].andNodeList[cxpointer1:cxpointer2] = ind2[0].andNodeList[cxpointer1:cxpointer2], ind1[0].andNodeList[cxpointer1:cxpointer2]
	ind1[0].andNodeInvertList[cxpointer1:cxpointer2], ind2[0].andNodeInvertList[cxpointer1:cxpointer2] = ind2[0].andNodeInvertList[cxpointer1:cxpointer2], ind1[0].andNodeInvertList[cxpointer1:cxpointer2]
	return ind1, ind2

# sets up GA toolbox for deap/adaptive
def buildToolbox( individualLength, bitFlipProb, model, params):
	# # #setup toolbox
	toolbox = base.Toolbox()
	
	pset = gp.PrimitiveSet("MAIN", arity=1)
	pset.addPrimitive(operator.add, 2)
	pset.addPrimitive(operator.sub, 2)
	pset.addPrimitive(operator.mul, 2)
	weightTup=(-1.0,)
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
	toolbox.register("select", selNSGA2)
	toolbox.register("similar", numpy.array_equal)
	return toolbox, stats

def evaluator(individual, cells, model,  sss, params, KOlist, KIlist,iteratorDict, j):
	boolValues=EBNbool(individual, model, cells, model.initValueList[j], params,KOlist[j], KIlist[j])
	return numpy.sum([(boolValues[i]-sss[j][model.nodeList[i]])**2 for i in range(0, len(model.nodeList))])	

# fitness calculation based on simulation			
def evaluate(individual, cells, model,  sss, params, KOlist, KIlist,iteratorDict):
	SSEs=[]
	# for j in range(0,len(sss)):
	# 	boolValues=EBNbool(individual, model, cells, model.initValueList[j], params,KOlist[j], KIlist[j], iteratorDict)
	# 	SSE=numpy.sum([(boolValues[i]-sss[j][model.nodeList[i]])**2 for i in range(0, len(model.nodeList))])
	# 	SSEs.append(SSE)
	SSEs=Parallel(n_jobs=min(24,len(sss)))(delayed(evaluator)(list(individual), cells, model,  sss, params, KOlist, KIlist,iteratorDict, j) for j in range(len(sss)))
	summer=1.*numpy.sum(SSEs)/len(SSEs)
	return summer,

# calculate node-wise fitness
def evaluateByNode(individual, cells, model,  sss, params, KOlist, KIlist,iteratorDict):
	boolValues=[EBNbool(list(individual), model, cells, model.initValueList[i], params, KOlist[i], KIlist[i]) for i in range(len(sss))]
	return tuple([numpy.sum([(boolValues[j][i]-sss[j][model.nodeList[i]])**2 for j in range(0,len(sss))]) for i in range(0, len(model.nodeList))])

# generates a random set of samples made up of cells by using parameteris from probInit seq
# to set up then iterating using strict Boolean modeling. 
def runProbabilityBooleanSims(individual, model, sampleNum, cells, params, KOlist, KIlist,iteratorDict):
	samples=[]
	seeds=[]
	for i in range(0,sampleNum):
		seeds.append(random())
	#samples=Parallel(n_jobs=min(24,sampleNum))(delayed(sampler)(individual, model, sampleNum, seeds[i], params, KOlist[i], KIlist[i],iteratorDict) for i in range(sampleNum))
	samples=[sampler(individual, model, sampleNum, seeds[i], params, KOlist[i], KIlist[i],iteratorDict) for i in range(sampleNum)]
	return samples

# generates random seed samples... i.e. generates random starting states then runs EBN
def sampler(individual, model, cells, seeder, params, KOs, KIs,iteratorDict):
	seed(seeder)
	cellArray=[]
	sampleProbs=[]
	# generate random proportions for each node to start
	for j in range(0,len(model.nodeList)):
		sampleProbs.append(random())
	# print(sampleProbs)
	return EBNbool(individual, model, cells, sampleProbs, params, KOs, KIs)

def varOrAdaptive(population, toolbox, model, lambda_, cxpb, mutpb, genfrac, mutModel):
	# algorithm for generating a list of offspring... copied and pasted from DEAP with modification for adaptive mutation
	assert (cxpb + mutpb) <= 1.0, ("The sum of the crossover and mutation "
		"probabilities must be smaller or equal to 1.0.")
	offspring = []
	for _ in xrange(lambda_):
		op_choice = random()
		if op_choice < cxpb:            # Apply crossover
			ind1, ind2 = map(toolbox.clone, sample(population, 2))
			ind1, ind2 = cxTwoPointNode(ind1, ind2)
			del ind1.fitness.values
			offspring.append(ind1)
		elif op_choice < cxpb + mutpb:  # Apply mutation
			ind = toolbox.clone(choice(population))
			ind, = mutFlipBitAdapt(ind, genfrac, mutModel)
			del ind.fitness.values
			offspring.append(ind)
		else:                           # Apply reproduction
			offspring.append(choice(population))
	return offspring

def mutFlipBitAdapt(indyIn, genfrac, mutModel):
	# get errors
	errors=list(indyIn.fitness.values)
	individual=indyIn[1]
	model=indyIn[0]
	for i in range(len(errors)):
		temper=model.successorNums[i]
		if temper==0:
			errors[i]=errors[i]*len(model.possibilityList[i])
		else:
			errors[i]=errors[i]*len(model.possibilityList[i])*temper

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
		if len(model.possibilityList[focusNode])>3 and random()<mutModel:
			temppermup=[]
			upstreamAdders=list(model.possibilityList[focusNode])
			rvals=list(model.rvalues[focusNode])
			while len(temppermup)<3:
				randy=random()# randomly select a node to mutate
				tempsum=sum(rvals)
				if tempsum==0:
					addNoder=int(math.floor(random()*len(upstreamAdders)))
				else:
					recalc=numpy.cumsum([1.*rval/tempsum for rval in rvals])
					# print(sum(recalc))
					# print(recalc)
					# print(rvals)
					addNoder=next(i for i in range(len(recalc)) if recalc[i]>randy)
				temppermup.append(upstreamAdders.pop(addNoder))
				rvals.pop(addNoder)
			model.update_upstream(focusNode,temppermup )
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
		indyIn[0]=model
		indyIn[1]=individual
	else:
		print('did not actually check')
	return indyIn,
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
	mean1=numpy.mean(ind1.wvalues)
	mean2=numpy.mean(ind2.wvalues)
	std1=numpy.std(ind1.wvalues)
	if mean1 > mean2 :
		not_equal = True
	elif mean1 < mean2:
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

def eaMuPlusLambdaAdaptive( toolbox, model, mu, lambda_, cxpb, mutpb, ngen, namer, newSSS,KOlist, KIlist,iteratorDict, params ,stats=None, verbose=__debug__):
	modelNodes=params.modelNodes
	# population=[[copy.deepcopy(model),genBits(model)]for i in range(params.popSize)]
	population=toolbox.population(n=params.popSize)
	mutModel=params.mutModel
	logbook = tools.Logbook()
	logbook.header = ['gen', 'nevals'] + (stats.fields if stats else [])
	lastcheck=[]
	modellist=[]
	fitnesslist=[]
	popList=[]
	# Evaluate the individuals with an invalid fitness
	invalid_ind = [ind for ind in population if not ind.fitness.valid]
	fitnesses=Parallel(n_jobs=min(24,len(invalid_ind)))(delayed(evaluateByNode)(indy[1], params.cells, indy[0],  newSSS, params, KOlist, KIlist,iteratorDict) for indy in invalid_ind)
	for ind, fit in zip(invalid_ind, fitnesses):
		ind.fitness.values = fit
	fitnesslist.append([list(ind.fitness.values) for ind in population])
	popList.append([list(inder[1]) for inder in population])
	modellist.append([[(modeler[0].size), list(modeler[0].nodeList), list(modeler[0].individualParse), list(modeler[0].andNodeList) , list(modeler[0].andNodeInvertList), list(modeler[0].andLenList),	list(modeler[0].nodeList), dict(modeler[0].nodeDict), list(modeler[0].initValueList)] for modeler in population])

	record = stats.compile(population) if stats is not None else {}
	logbook.record(gen=0, nevals=len(invalid_ind), **record)
	if verbose:
		print(logbook.stream)

	breaker=False
	for ind in population:
		if numpy.sum(ind.fitness.values)< .01*len(ind.fitness.values):
			breaker=True
	if breaker:
		# outputList=[fitnesslist, popList, modellist]
		# pickle.dump( outputList, open( namer+"_pops.pickle", "wb" ) )
		return population, logbook

	# Begin the generational process
	for gen in range(1, ngen+1):
		offspring = varOrAdaptive(population, toolbox, model, lambda_, .5+.5*(1.-1.*gen/ngen), (.5*gen/ngen), (1.*gen/ngen),mutModel)
		# Evaluate the individuals with an invalid fitness
		invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
		fitnesses=Parallel(n_jobs=min(24,len(invalid_ind)))(delayed(evaluateByNode)(indy[1], params.cells, indy[0],  newSSS, params, KOlist, KIlist,iteratorDict) for indy in invalid_ind)
		for ind, fit in zip(invalid_ind, fitnesses):
			ind.fitness.values = fit

		# Select the next generation population
		population[:] = toolbox.select(population + offspring, mu)
		fitnesslist.append([list(ind.fitness.values) for ind in population])
		popList.append([list(inder[1]) for inder in population])
		modellist.append([[(modeler[0].size), list(modeler[0].nodeList), list(modeler[0].individualParse), list(modeler[0].andNodeList) , list(modeler[0].andNodeInvertList), list(modeler[0].andLenList),	list(modeler[0].nodeList), dict(modeler[0].nodeDict), list(modeler[0].initValueList)] for modeler in population])

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
			# outputList=[fitnesslist, popList, modellist]
			# pickle.dump( outputList, open( namer+"_pops.pickle", "wb" ) )
			return population, logbook

	# outputList=[fitnesslist, popList, modellist]
	# pickle.dump( outputList, open( namer+"_pops.pickle", "wb" ) )
	return population, logbook

def eaMuPlusLambdaAdaptive2( toolbox, model, mu, lambda_, cxpb, mutpb, ngen, namer, newSSS,KOlist, KIlist,seed, params ,stats=None, verbose=__debug__):
	modelNodes=params.modelNodes
	# population=[[copy.deepcopy(model),genBits(model)]for i in range(params.popSize)]
	population=toolbox.population(n=params.popSize)
	for item in population:
		for i in range(len(item[1])):
			item[1][i]=seed[i]
	mutModel=params.mutModel
	mutpb=.8
	cxpb=.2
	logbook = tools.Logbook()
	logbook.header = ['gen', 'nevals'] + (stats.fields if stats else [])
	lastcheck=[]
	modellist=[]
	fitnesslist=[]
	popList=[]
	iteratorDict={}
	# Evaluate the individuals with an invalid fitness
	invalid_ind = [ind for ind in population if not ind.fitness.valid]
	fitnesses=Parallel(n_jobs=min(24,len(invalid_ind)))(delayed(evaluateByNode)(indy[1], params.cells, indy[0],  newSSS, params, KOlist, KIlist,iteratorDict) for indy in invalid_ind)
	for ind, fit in zip(invalid_ind, fitnesses):
		ind.fitness.values = fit
	fitnesslist.append([list(ind.fitness.values) for ind in population])
	popList.append([list(inder[1]) for inder in population])
	modellist.append([[(modeler[0].size), list(modeler[0].nodeList), list(modeler[0].individualParse), list(modeler[0].andNodeList) , list(modeler[0].andNodeInvertList), list(modeler[0].andLenList),	list(modeler[0].nodeList), dict(modeler[0].nodeDict), list(modeler[0].initValueList)] for modeler in population])

	record = stats.compile(population) if stats is not None else {}
	logbook.record(gen=0, nevals=len(invalid_ind), **record)
	if verbose:
		print(logbook.stream)

	breaker=False
	for ind in population:
		if numpy.sum(ind.fitness.values)< .01*len(ind.fitness.values):
			breaker=True
	if breaker:
		# outputList=[fitnesslist, popList, modellist]
		# pickle.dump( outputList, open( namer+"_pops2.pickle", "wb" ) )
		return population, logbook
	# Begin the generational process
	for gen in range(1, ngen+1):
		offspring = varOrAdaptive(population, toolbox, model, lambda_,  (.5*gen/ngen),.5+.5*(1.-1.*gen/ngen), (1.*gen/ngen),mutModel)
		# Evaluate the individuals with an invalid fitness
		invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
		fitnesses=Parallel(n_jobs=min(24,len(invalid_ind)))(delayed(evaluateByNode)(indy[1], params.cells, indy[0],  newSSS, params, KOlist, KIlist,iteratorDict) for indy in invalid_ind)
		for ind, fit in zip(invalid_ind, fitnesses):
			ind.fitness.values = fit

		# Select the next generation population
		population[:] = toolbox.select(population + offspring, mu)
		fitnesslist.append([list(ind.fitness.values) for ind in population])
		popList.append([list(inder[1]) for inder in population])
		modellist.append([[(modeler[0].size), list(modeler[0].nodeList), list(modeler[0].individualParse), list(modeler[0].andNodeList) , list(modeler[0].andNodeInvertList), list(modeler[0].andLenList),	list(modeler[0].nodeList), dict(modeler[0].nodeDict), list(modeler[0].initValueList)] for modeler in population])

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
			# outputList=[fitnesslist, popList, modellist]
			# pickle.dump( outputList, open( namer+"_pops2.pickle", "wb" ) )
			return population, logbook
	# outputList=[fitnesslist, popList, modellist]
	# pickle.dump( outputList, open( namer+"_pops2.pickle", "wb" ) )
	return population, logbook
def GAautoSolver(model, sss, params, KOlist, KIlist,iteratorDict, namer):
	# set up toolbox and run GA with or without adaptive mutations turned on
	toolbox, stats=buildToolbox(model.size,params.bitFlipProb, model, params)
	if not params.adaptive:
		toolbox.register("evaluate", evaluate, cells=params.cells,model=model,sss=sss, params=params, KOlist=KOlist, KIlist=KIlist, iteratorDict=iteratorDict)
		population=toolbox.population(n=params.popSize)
	if params.adaptive:
		output=eaMuPlusLambdaAdaptive(toolbox, model, mu=params.mu, lambda_=params.lambd, stats=stats, cxpb=params.crossoverProb, mutpb=params.mutationProb, ngen=params.generations, namer=namer, newSSS= sss,KOlist=KOlist, KIlist=KIlist, params=params, iteratorDict=iteratorDict,  verbose=params.verbose)
	else:
		output=algo.eaMuCommaLambda(population, toolbox, mu=params.mu, lambda_=params.lambd, stats=stats, cxpb=params.crossoverProb, mutpb=params.mutationProb, ngen=params.generations, verbose=params.verbose)
	return output

def GAsearchModel(model, sss,params, KOlist, KIlist, iteratorDict, namer):

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
			if  model.nodeList[j] in sss[0]:
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
			if  model.nodeList[j] in sss[0]:
				newInitValueList[k].append(ss[model.nodeList[j]])
			else:
				newInitValueList[k].append(0.5)
	model.initValueList=newInitValueList
	population, logbook=GAautoSolver(model, sss, params, knockoutLists, knockinLists,iteratorDict, namer)
	out1, out2, model  = findPopBest(population)
	return model,out1,out2

def GASearchModel2(model, sss,params, KOlist, KIlist, seed, namer):
	
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
			if  model.nodeList[j] in sss[0]:
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
			if  model.nodeList[j] in sss[0]:
				newInitValueList[k].append(ss[model.nodeList[j]])
			else:
				newInitValueList[k].append(0.5)
	model.initValueList=newInitValueList
	toolbox, stats=buildToolbox(model.size,params.bitFlipProb, model, params)
	population, logbook=eaMuPlusLambdaAdaptive2(toolbox, model, mu=params.mu, lambda_=params.lambd, stats=stats, cxpb=params.crossoverProb, mutpb=params.mutationProb, ngen=params.generations, namer=namer+'_2_', newSSS= sss,KOlist=KOlist, KIlist=KIlist, params=params, seed=seed,  verbose=params.verbose)
	out1, out2, model  = findPopBest(population)
	return model,out1,out2

def localSearch(model, indy, newSSS, params, KOlist, KIlist):
	iteratorDict={}
	outputs=Parallel(n_jobs=min(24,len(model.nodeList)))(delayed(checkNodePossibilities)(node, indy, newSSS, params.cells, model,params, KOlist, KIlist,iteratorDict ) for node in range(len(model.nodeList)))
	equivs=[]
	individual=[]
	devs=[]
	for output in outputs:
		individual.extend(output[0])
		equivs.append(output[1])
		devs.append(output[2])
	return individual, equivs, devs

def checkNodePossibilities(node, indy, newSSS, cellNum, model,params, KOlist, KIlist,iteratorDict ):
	tol=.01*len(newSSS)
	if node==(len(model.nodeList)-1):
		end=len(indy)
	else:
		end=model.individualParse[node+1]
	start=model.individualParse[node]
	truth=list(indy[start:end])
	equivs=[truth] 	#add if
	if (end-start)==0:
		return truth, equivs, 0.
	indOptions=[]
	indErrors=[]
	for i in range(1,2**(end-start)):
		tempultimate=list(indy)
		tempInd=utils.bitList(i, len(truth))
		tempultimate[start:end]=tempInd
		truth=list(tempInd)
		currentsumtemp=evaluateByNode(tempultimate, cellNum, model,  newSSS, params,KOlist, KIlist,iteratorDict )
		currentsum=currentsumtemp[node]
		indOptions.append(truth)
		indErrors.append(currentsum)
		iteratorDict={}
		# devListGA.append(dev)
		gc.collect()
	minny= min(indErrors)
	equivs=[]
	for i in range(len(indOptions)):
		if indErrors[i]< minny+tol:
			equivs.append(indOptions[i])
	return truth, equivs, minny
