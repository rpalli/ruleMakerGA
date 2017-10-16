
# import python packages
import pickle
from deap import base, creator, gp, tools
from deap import algorithms as algo
from random import random, seed, shuffle, randint, sample, choice
import numpy as numpy
import operator
import math as math
from sets import Set
from joblib import Parallel, delayed

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
	return startInd

# executes two point crossover at node junctions
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
	cxpoint1 = randint(1, size)
	cxpoint2 = randint(1, size - 1)
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

# sets up GA toolbox for deap/adaptive
def buildToolbox( individualLength, bitFlipProb, model):
	# # #setup toolbox
	toolbox = base.Toolbox()
	
	pset = gp.PrimitiveSet("MAIN", arity=1)
	pset.addPrimitive(operator.add, 2)
	pset.addPrimitive(operator.sub, 2)
	pset.addPrimitive(operator.mul, 2)
	weightTup=(-1.0,)
	params=paramClass()
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

# fitness calculation based on simulation			
def evaluate(individual, cells, model,  sss, params, KOlist, KIlist,iteratorDict):
	SSEs=[]
	for j in range(0,len(sss)):
		boolValues=EBNbool(individual, model, cells, model.initValueList[j], params,KOlist[i], KIlist[i], iteratorDict)
		SSE=numpy.sum([(boolValues[i]-sss[j][model.nodeList[i]])**2 for i in range(0, len(model.nodeList))])
		SSEs.append(SSE)
	summer=1.*numpy.sum(SSEs)/len(SSEs)
	return summer,

# calculate node-wise fitness
def evaluateByNode(individual, cells, model,  sss, params, KOlist, KIlist,iteratorDict):
	boolValues=[EBNbool(list(individual), model, cells, model.initValueList[i], params, KOlist[i], KIlist[i],iteratorDict) for i in range(len(sss))]
	return tuple([numpy.sum([(boolValues[j][i]-sss[j][model.nodeList[i]])**2 for j in range(0,len(sss))]) for i in range(0, len(model.nodeList))])

# generates a random set of samples made up of cells by using parameteris from probInit seq
# to set up then iterating using strict Boolean modeling. 
def runProbabilityBooleanSims(individual, model, sampleNum, cells, params, KOlist, KIlist,iteratorDict):
	samples=[]
	seeds=[]
	for i in range(0,sampleNum):
		seeds.append(random())
	samples=Parallel(n_jobs=min(24,sampleNum))(delayed(sampler)(individual, model, sampleNum, seeds[i], params, KOlist[i], KIlist[i],iteratorDict) for i in range(sampleNum))
	return samples

# generates random seed samples... i.e. generates random starting states then runs EBN
def sampler(individual, model, cells, seeder, params, KOs, KIs,iteratorDict):
	seed(seeder)
	cellArray=[]
	sampleProbs=[]
	# generate random proportions for each node to start
	for j in range(0,len(model.nodeList)):
		sampleProbs.append(random())
	return iterateBooleanModel(individual, model, cells, sampleProbs, params, KOs, KIs,iteratorDict)

def varOrAdaptive(population, toolbox, model, lambda_, cxpb, mutpb, genfrac):
	# algorithm for generating a list of offspring... copied and pasted from DEAP with modification for adaptive mutation
	assert (cxpb + mutpb) <= 1.0, ("The sum of the crossover and mutation "
		"probabilities must be smaller or equal to 1.0.")
	offspring = []
	for _ in xrange(lambda_):
		op_choice = random()
		if op_choice < cxpb:            # Apply crossover
			ind1, ind2 = map(toolbox.clone, sample(population, 2))
			ind1, ind2 = cxTwoPointNode(ind1, ind2, model)
			del ind1.fitness.values
			offspring.append(ind1)
		elif op_choice < cxpb + mutpb:  # Apply mutation
			ind = toolbox.clone(choice(population))
			ind, = mutFlipBitAdapt(ind,  .5, model, genfrac)
			del ind.fitness.values
			offspring.append(ind)
		else:                           # Apply reproduction
			offspring.append(choice(population))
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
		# fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
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

def GAautoSolver(model, sss, params, KOlist, KIlist,iteratorDict):
	# set up toolbox and run GA with or without adaptive mutations turned on
	toolbox, stats=buildToolbox(model.size,params.bitFlipProb, model)
	if params.adaptive:
		toolbox.register("evaluate", evaluateByNode, cells=params.cells,model=model,sss=sss, params=params, KOlist=KOlist, KIlist=KIlist, iteratorDict=iteratorDict)
	else:
		toolbox.register("evaluate", evaluate, cells=params.cells,model=model,sss=sss, KOlist=KOlist, KIlist=KIlist, iteratorDict=iteratorDict)
	population=toolbox.population(n=params.popSize)
	hof = tools.HallOfFame(params.hofSize, similar=numpy.array_equal)
	
	if params.adaptive:
		output=eaMuPlusLambdaAdaptive(population, toolbox, model, mu=params.mu, lambda_=params.lambd, stats=stats, cxpb=params.crossoverProb, mutpb=params.mutationProb, ngen=params.generations, verbose=params.verbose, halloffame=hof)
	else:
		output=algo.eaMuCommaLambda(population, toolbox, mu=params.mu, lambda_=params.lambd, stats=stats, cxpb=params.crossoverProb, mutpb=params.mutationProb, ngen=params.generations, verbose=params.verbose, halloffame=hof)
	return output

def GAsearchModel(model, sss,params, KOlist, KIlist):
	iteratorDict={}

	newSSS=[]
	newSSS.append({})
	newSSS.append({})
	newSSS.append({})
	newSSS.append({})
	newSSS.append({})
	newSSS.append({})
	for key in sss[0]:
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
	population, logbook=GAautoSolver(model, newSSS, params, knockoutLists, knockinLists,iteratorDict)
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
			fitnesses=Parallel(n_jobs=min(24,len(invalid_ind)))(delayed(evaluateByNode)(indy, params.cells, model,  newSSS, params, KOlist, KIlist,iteratorDict) for indy in invalid_ind)
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
