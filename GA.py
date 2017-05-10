
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

# import other pieces of our software
import networkConstructor as nc
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
			while numpy.sum(truth)>5 and counter < 100000:
				indices = [i for i in range(len(truth)) if truth[i] == 1]
				chosen=math.floor(random()*len(indices))
				truth[indices[int(chosen)]]=0
				counter+=1
			#startInd[start:end]=truth
			if numpy.sum(truth)==0:
				chosen=math.floor(random()*len(truth))
				truth[int(chosen)]=1
		elif len(truth)==1:
			truth[0]=1
		startInd[start:end]=truth
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
	cxpoint1 = random.randint(1, size)
	cxpoint2 = random.randint(1, size - 1)
	if cxpoint2 >= cxpoint1:
		cxpoint2 += 1
	else: # Swap the two cx points
		cxpoint1, cxpoint2 = cxpoint2, cxpoint1
	cxpoint1=model.individualParse[cxpoint1]
	cxpoint2=model.individualParse[cxpoint2]
	ind1[cxpoint1:cxpoint2], ind2[cxpoint1:cxpoint2] = ind2[cxpoint1:cxpoint2], ind1[cxpoint1:cxpoint2]
	return ind1, ind2

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
	toolbox.register("select", tools.selNSGA2)
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
	# SSEs=[]
	#boolValues=Parallel(n_jobs=min(6,len(sss)))(delayed(iterateBooleanModel)(list(individual), model, cells, model.initValueList[i], params) for i in range(len(sss)))
	boolValues=[iterateBooleanModel(list(individual), model, cells, model.initValueList[i], params) for i in range(len(sss))]
	# for i in range(0, len(model.evaluateNodes)):
	# 	SSEs= numpy.sum([(boolValues[j][model.evaluateNodes[i]]-sss[j][model.nodeList[model.evaluateNodes[i]]])**2 for j in range(0,len(sss))])
	# 	SSEs.append(SSE)
	SSEs= [numpy.sum([(boolValues[j][model.evaluateNodes[i]]-sss[j][model.nodeList[model.evaluateNodes[i]]])**2 for j in range(0,len(sss))]) for i in range(0, len(model.evaluateNodes))]
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
	return [1.*numpy.sum(col) / float(cells) for col in zip(*cellArray)]

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
			ind1, ind2 = cxTwoPointNode(ind1, ind2, model)
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

	# get rid of errors in nodes that can't be changed
	for j in xrange(len(errors)):
		if model.andLenList[model.evaluateNodes[j]]<2:
			errors[model.evaluateNodes[j]]=0

	# if error is relatively low, do a totally random mutation
	if numpy.sum(errors)<.05:
		for i in xrange(len(individual)):
			if random() < indpb:
				individual[i] = type(individual[i])(not individual[i])
	else:
		# if errors are relatively high, focus on nodes that fit the worst

		# normalize errors to get a probability that the node  is modified
		normerrors=[1.*error/numpy.sum(errors) for error in errors]
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
	for self_wvalue, other_wvalue in zip(self.wvalues[obj], other.wvalues[obj]):
		if self_wvalue > other_wvalue + 1:
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
	fitnesses=Parallel(n_jobs=min(6,len(invalid_ind)))(delayed(toolbox.evaluate)(list(indy)) for indy in invalid_ind)

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
		offspring = varOrAdaptive(population, toolbox, model, lambda_, cxpb, mutpb)
		
		# Evaluate the individuals with an invalid fitness
		invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
		#fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
		fitnesses=Parallel(n_jobs=min(6,len(invalid_ind)))(delayed(toolbox.evaluate)(list(indy)) for indy in invalid_ind)
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
		if breaker:
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
		output=eaMuPlusLambdaAdaptive(population, toolbox, model, mu=params.mu, lambda_=params.lambd, stats=stats, cxpb=params.crossoverProb, mutpb=params.mutationProb, ngen=params.generations, verbose=params.verbose, halloffame=hof)
	else:
		output=algo.eaMuCommaLambda(population, toolbox, mu=params.mu, lambda_=params.lambd, stats=stats, cxpb=params.crossoverProb, mutpb=params.mutationProb, ngen=params.generations, verbose=params.verbose, halloffame=hof)
	return output

def GAsearchModel(model, newSSS,params):
	population, logbook=GAautoSolver(model, newSSS, params)
	minny=1000000
	saveVal=-1
	for i in range(len(population)):
		if numpy.sum(population[i].fitness.values)< minny:
			minny=numpy.sum(population[i].fitness.values)
			saveVal=i
	ultimate=list(population[saveVal])
	minvals=population[saveVal].fitness.values
	newultimate=[]
	if minny>.01*len(population[saveVal].fitness.values):
		#iterate over nodes
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
		print(len(ultimate))
		flags=Parallel(n_jobs=min(24,len(ultimate)))(delayed(deltaTester)(list(ultimate), i, model,  newSSS, params, minny) for i in range(len(ultimate)))
		for flag, i in zip(flags,range(len(ultimate))):
			if flag:
				ultimate[i]=1-ultimate[i]
		minny=numpy.sum(evaluateByNode(ultimate, params.cells, model,  newSSS, params))
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
				newtot=numpy.sum(evaluateByNode(copied, params.cells, model,  newSSS, params))
				if newtot<minny:
					modifyFlag=True
	return modifyFlag

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
		temp=[100 if sumindividual[i]==0 else 1.*(sumindividual[i]-ones[i])/(sumindividual[i]) for i in range(0,len(ones))]
		temp=filter(lambda a: a != 100, temp)
		if len(temp)==0:
			sensitivity=100
		else:
			sensitivity=(1.*numpy.sum(temp)/len(temp))
		temp=[100 if (len(newindividual)-sumindividual[i])==0 else (1.*len(newindividual)-sumindividual[i]-negones[i])/(len(newindividual)-sumindividual[i]) for i in range(0,len(ones))]
		temp=filter(lambda a: a != 100, temp)
		if len(temp)==0:
			specificity=100
		else:
			specificity=(1.*numpy.sum(temp)/len(temp))
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
		sumindividual.append(1.*numpy.sum(truth))
		newindividual=[a_i - b_i for a_i, b_i in zip(truth, test)]
		ones.append(newindividual.count(1))
		zeros.append(newindividual.count(0))
		negones.append(newindividual.count(-1))
	temp=[100 if (sumindividual[i])==0 else 1.*(sumindividual[i]-ones[i])/(sumindividual[i]) for i in range(0,len(ones))]
	temp=filter(lambda a: a != 100, temp)
	if len(temp)==0:
		sensitivity=100
	else:
		sensitivity=(1.*numpy.sum(temp)/len(temp))
	temp=[100 if (len(newindividual)-sumindividual[i])==0 else (1.*len(newindividual)-sumindividual[i]-negones[i])/(len(newindividual)-sumindividual[i]) for i in range(0,len(ones))]
	temp=filter(lambda a: a != 100, temp)
	if len(temp)==0:
		specificity=100
	else:
		specificity=(1.*numpy.sum(temp)/len(temp))
	return sensitivity, specificity, nodesensitivity, nodespecificity

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

	sss=utils.synthesizeInputs(graph,params.samples)
	model=sim.modelClass(graph,sss)
	samples=params.samples

	# set up necessary lists for output
	truthList= [] 
	testList =[]
	devList=[]

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
		output=runProbabilityBooleanSims(individual, model, samples, params.cells, params)
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
		dev, bruteOut=GAsearchModel(model, newSSS, params)
		# add output of this trial to lists
		truthList.append(individual)
		testList.append(list(bruteOut))
		devList.append(dev)
	# set up output and save as a pickle
	outputList=[truthList,testList, devList,model.size, model.evaluateNodes, model.individualParse, model.andNodeInvertList, model.andLenList,model.earlyEvalNodes,	model.nodeList, model.nodeDict, model.initValueList]
	pickle.dump( outputList, open( name+"_output.pickle", "wb" ) )

if __name__ == '__main__':
	import time
	start_time = time.time()
	parser = argparse.ArgumentParser()
	parser.add_argument("graph")
	parser.add_argument("iterNum")
	results = parser.parse_args()
	graphName=results.graph
	iterNum=int(results.iterNum)
	name=graphName[:-8]+'_'+results.iterNum
	graph = nx.read_gpickle(graphName)
	simTester(graph, name)
	print("--- %s seconds ---" % (time.time() - start_time))