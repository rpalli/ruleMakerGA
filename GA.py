
from deap import base
from deap import creator
from deap import gp
from deap import tools
from deap import algorithms as algo

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
	toolbox.register("genRandomBitString", genRandBits, individualLength=individualLength)
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
		if async:
			boolValues=simulator.iterateModel(individual, params, model, initValueList[j])	
		else:
			boolValues=simulator.runModel(individual, params, model,initValueList[j])	
		for i in range(0, len(model.evaluateNodes)):
			SSE=SSE+(boolValues[model.evaluateNodes[i]]-ss[model.nodeList[model.evaluateNodes[i]]])**2
		SSEs.append(SSE)
	edgeDegree=0
	if complexPenalty:
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
	params=paramClass()
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
	for i in range(0,len(model.nodeList)):
		print(model.nodeList[i])
		print(model.individualParse[i])
		currentDev=10*len(sss1)
		best=[]
		if model.possibilityNumList[i]>0:
			for j in range(0,int(2**int(model.individualParse[i][2]-model.individualParse[i][0]+1))):
				bits=[]
				bits=bitList(j)
				while len(bits)<(model.individualParse[i][2]-model.individualParse[i][0]+1):
					bits.insert(0,0)
				deviation=0
				for steadyStateNum in range(0,len(sss1)):
					derivedVal=updateNode(i,model.initValueList[steadyStateNum],bits, model, simulator)
					deviation=deviation+(derivedVal-sss1[steadyStateNum][model.nodeList[i]])**2
				print(bits)
				print(deviation)
				if(deviation<currentDev):
					print("best")
					best=bits
					currentDev=deviation	
		bestList.append(best)
	return [item for sublist in bestList for item in sublist]


def autoSimulator(async, iters, toolbox, stats, graph, trials, sss,individualLength,evaluateNodes, individualParse,nodeList,inputOrderList,inputOrderInvertList,possibilityNumList,initValueList, steps ,nodeNoise,networkNoise,popSize,crossoverProb,mutationProb,generations,earlyEvalNodes, mu,lambd, complexPenalty, genSteps):
	hofSize=10
	avgSSEs=[]
	hofs=[]
	truthIndividuals=[]
	truthAvgs=[]
	hofScores=[]
	newSSS=[]
	truthCounter=[0 for number in xrange(hofSize)]
	trueNodeList=[]
	for i in range(0,trials):
		individual=genRandBits(individualLength)
		newSSS=[]
		for k in range(0,len(sss)):
			if async:
				boolValues=iterateFuzzyModel(async, iters,individual, individualParse, nodeList, inputOrderList, inputOrderInvertList, possibilityNumList, initValueList[k], genSteps,nodeNoise,networkNoise,earlyEvalNodes)
			else:
				boolValues=runFuzzyModel(individual, async,individualParse, nodeList, inputOrderList, inputOrderInvertList, possibilityNumList, initValueList[k], genSteps ,nodeNoise,networkNoise,earlyEvalNodes)
			newSS=copy.deepcopy(sss[k])
			for j in range(0,len(evaluateNodes)):
				newSS[nodeList[evaluateNodes[j]]]=boolValues[evaluateNodes[j]]
			newSSS.append(newSS)
		#update eval nodes
		print(evaluate(individual, async, iters, newSSS, evaluateNodes, individualParse, nodeList, inputOrderList, inputOrderInvertList, possibilityNumList, initValueList, steps,earlyEvalNodes,complexPenalty))
		newInitValueList=[]
		for j in range(0,len(sss)):
			newInitValueList.append([])
		for j in range(0,len(nodeList)):
			for k in range(0,len(sss)):
				ss=newSSS[k]
				if  nodeList[j] in sss[0].keys():
					newInitValueList[k].append(ss[nodeList[j]])
				else:
					newInitValueList[k].append(0.5)
		#print(len(sss))
		toolbox.register("evaluate", evaluate, iters=iters, async=async, sss=newSSS, evaluateNodes=evaluateNodes,individualParse= individualParse, nodeList=nodeList, inputOrderList=inputOrderList, inputOrderInvertList=inputOrderInvertList, possibilityNumList=possibilityNumList, initValueList=newInitValueList, steps=steps, earlyEvalNodes=earlyEvalNodes, complexPenalty=complexPenalty)
		#toolbox.register("evaluate", evaluate, ss=ss, evaluateNodes=evaluateNodes,individualParse= individualParse, nodeList=nodeList, inputOrderList=inputOrderList, inputOrderInvertList=inputOrderInvertList, possibilityNumList=possibilityNumList, initValues=boolValues, steps=steps, h=h,p=p,hillOn=hillOn)
		population=toolbox.population(n=popSize)
		hof = tools.HallOfFame(hofSize, similar=numpy.array_equal)
		output = tools.Logbook()
		output=algo.eaMuCommaLambda(population, toolbox, mu=mu, lambda_=lambd, stats=stats, cxpb=crossoverProb, mutpb=mutationProb, ngen=generations, verbose=True, halloffame=hof)
		truth=(writeFuzzyModel(individual, individualParse, nodeList,inputOrderList, inputOrderInvertList, possibilityNumList))
		
		for j in range(0,10):
			bestRun=(writeFuzzyModel(hof[j], individualParse, nodeList,inputOrderList, inputOrderInvertList, possibilityNumList))
			if truth==bestRun:
				truthCounter[j]+=1
				break
			
			elif j==0:
				truthlines=truth.split('\n')
				newlines=bestRun.split('\n')
				trueNodes=0
				for k in range(0,len(truthlines)):
					if truthlines[k]==newlines[k]:
						trueNodes=trueNodes+1
					else:
						print("incorrect pair: true then test")
						print(truthlines[k])
						print(newlines[k])
				trueNodeList.append(trueNodes)
		# else:
			# print('WRONG')
			# print(truth)
			# print(bestRun)
		avgs=[output[1][k]['min'] for k in range(0,len(output[1]))]
		#print(avgs)
		hofs.append(hof)
		temp=[]
		for hofind in range(0,len(hof)):
			tempVal=0
			for bit in range(0,len(hof[hofind])):
				if not hof[hofind][bit]==individual[bit]:
					tempVal=tempVal+1
			temp.append(tempVal)
		hofScores.append(temp)
		truthIndividuals.append(individual)
		truthAvgs.append(boolValues)
		gc.collect()
	# prefix='fatimaTest_noise_'+str(nodeNoise)+'_networkNoise_'+str(networkNoise)
	# outfile = open(prefix+'_hof.pkl', 'wb')
	# pickle.dump(hof, outfile)
	# outfile.close()
	f = open('hof_differences_'+str(nodeNoise)+str(networkNoise)+'.txt', 'w')
	g = open('hof_individuals'+str(nodeNoise)+str(networkNoise)+'.txt',  'w')
	h = open('truth_individuals'+str(nodeNoise)+str(networkNoise)+'.txt',  'w')
	f.write(str(hofScores))
	f.write('\n')
	for hofy in range(0,len(hofs)):
		g.write(str(hofs[hofy]))
		g.write('\n')
	h.write(str(truthIndividuals))
	h.write('\n')
	f.close()
	g.close()
	h.close()
	print('# true of # trials')
	for k in range(0,len(truthCounter)):
		print(truthCounter[k])
	print(trials)
	print("for the incorrects, by node")
	print(trueNodeList)
	print(len(truthlines))
	
if __name__ == '__main__':
	
	graph = old.LiuNetwork3Builder()