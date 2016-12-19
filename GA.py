def buildToolbox( individualLength, bitFlipProb, samples):
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
			
def setupGAparams(graph, sss):
	repeat=True
	while(repeat):
		repeat=False
		for node in graph.nodes():
			if node in graph.successors(node):
				graph.remove_edge(node,node)
				repeat=True
	
	evaluateNodes=[] #list of nodes that need to be compared to the steady state values (sss)
	individualParse=[] #list of triples that contain the start and end of the noder order value and operator flag list for each node
	inputOrderList=[] #list of lists of possible node orders for each node
	inputOrderInvertList=[] # keeps track of which incoming nodes for each node need to be inverted
	possibilityNumList=[] # keeps track of length of above inputOrderList for each node
	earlyEvalNodes=[] # nodes that don't have initial value and need to be re-evaluated early on in the simulation
	
	nodeList=graph.nodes()#define the node list simply as the nodes in the graph. 
	nodeDict={} #identifies names of nodes with their index in the node list
	for i in range(0,len(nodeList)):
		nodeDict[nodeList[i]]=i #constructs the node dict so we can easily look up nodes
	
	counter=int(0) #keeps track of where we are in the generic individual
	initValueList=[] #starting states for nodes
	for j in range(0,len(sss)): #construct empty lists to stick values in later for intiial value list
			initValueList.append([])
	

	for i in range(0,len(nodeList)):
		preds=graph.predecessors(nodeList[i]) # get predecessors of node. 
		if len(preds)>15: #handle case where there are too many predecessors by truncation
			preds=preds[1:15]
		for j in range(0,len(preds)):
			preds[j]=nodeDict[preds[j]]
		
		# the followign section constructs a list of possible node orders
		withNones = zip(preds, itertool.repeat('empty'))
		possibilities=list(itertool.product(*withNones))
		for j in range(0,len(possibilities)):
			possibilities[j]=list(possibilities[j])
			while 'empty' in possibilities[j]:
				possibilities[j].remove('empty')
			while [] in possibilities[j]:
				possibilities[j].remove([])
		while [] in possibilities:
			possibilities.remove([])
		
		#store activities of nodes in correspondence with node order. 
		activities=[] #list to store activities of nodes (a vs i)
		for sequence in possibilities:
			activity=[]
			for node in sequence:
				if graph.edge[nodeList[node]][nodeList[i]]['signal']=='a':
					activity.append(False)
				else:
					activity.append(True)
			activities.append(activity)
		inputOrderList.append(possibilities)
		inputOrderInvertList.append(activities)
		possibilityNumList.append(len(possibilities))
		
		#deal with nodes with no incoming or only one incoming first
		if len(possibilities)==0:
			individualParse.append([counter,counter,counter])
		elif len(possibilities)==1:
			individualParse.append([counter,counter,counter])
			counter=int(counter+1)
		else:
			logNum=ceil(log(len(possibilities))/log(2))#determine number of bits necessary to store information of which node order is the one for that individual
			individualParse.append([int(counter),int(counter+logNum-1),int(counter+logNum+len(preds)-2)]) #set up limits: lognum-1 for list then len(preds)-1 for the activity flags
			counter=int(counter+logNum+len(preds)-1) #uptick the counter so we know where to start on the next node
		#set up lists of initial values and values that need to be evaluated early
		for j in range(0,len(sss)):
			ss=sss[j]
			if  nodeList[i] in sss[0].keys():
				initValueList[j].append(ss[nodeList[i]])
				evaluateNodes.append(i)
			else:
				initValueList[j].append(0.5)
				earlyEvalNodes.append(i)
	return counter-1, individualParse, nodeList, inputOrderList, initValueList, evaluateNodes, possibilityNumList, inputOrderInvertList, earlyEvalNodes
	

def evaluate(individual, async, iters, sss, evaluateNodes, individualParse, nodeList, inputOrderList, inputOrderInvertList, possibilityNumList, initValueList, steps,earlyEvalNodes,complexPenalty):
	SSEs=[]
	for j in range(0,len(sss)):
		ss=sss[j]
		initValues=initValueList[j]
		SSE=0
		if async:
			boolValues=iterateFuzzyModel(async, iters, individual, individualParse, nodeList, inputOrderList, inputOrderInvertList, possibilityNumList, initValueList[j], steps, 0,0,earlyEvalNodes)	
		else:
			boolValues=runFuzzyModel(individual, async, individualParse, nodeList, inputOrderList, inputOrderInvertList, possibilityNumList,initValueList[j], steps,0,0,earlyEvalNodes)	
		for i in range(0, len(evaluateNodes)):
			SSE=SSE+(boolValues[evaluateNodes[i]]-ss[nodeList[evaluateNodes[i]]])**2
		SSEs.append(SSE)
	edgeDegree=0
	if complexPenalty:
		for i in range(0,len(nodeList)):
			if possibilityNumList[i]>0:
				edgeDegree=edgeDegree+len(inputOrderList[i][bit2int(individual[individualParse[i][0]:individualParse[i][1]])%possibilityNumList[i]])
	SSEs.append(edgeDegree/1.)
	gc.collect()
	summer=0.
	for i in range(0,len(SSEs)):
		summer+=SSEs[i]
	return summer,

	
def writeFuzzyModel(individual, individualParse, nodeList,inputOrderList, inputOrderInvertList, possibilityNumList):
	#iterate over nodes to generate a BooleanNet representation for the entire model
	addString=''
	for i in range(0,len(nodeList)):
		addString=addString+writeFuzzyNode(i,individual, individualParse, nodeList,inputOrderList, inputOrderInvertList, possibilityNumList)
		addString=addString+'\n'
	return addString
