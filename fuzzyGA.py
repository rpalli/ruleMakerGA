#import stuff
from random import random
from random import shuffle
from random import gauss
import numpy as numpy
from math import floor, ceil, log
import operator
from deap import base
from deap import creator
from deap import gp
from deap import tools
from deap import algorithms as algo
import networkx as nx
import fuzzyGA_oldTests as tester
import fuzzyNetworkConstructor as constructor
import csv
import itertools as itertool
import copy as copy
import gc as gc

def returnSimParams():
	simSteps=100 # number of steps each individual is run when evaluating
	generations=50 # generations to run
	popSize=100 #size of population
	mu= 100#individuals selected
	lambd= 200#children produced
	bitFlipProb=.1 # prob of flipping bits inside mutation
	crossoverProb=.2 # prob of crossing over a particular parent
	mutationProb=.2 # prob of mutating a particular parent
	async=False # run in asynchronous mode
	iters=100 #number of simulations to try in asynchronous mode
	complexPenalty=False #penalize models with increased complexity
	genSteps=10000 # steps to find steady state with fake data
	return async, iters, simSteps, generations, popSize, bitFlipProb, crossoverProb,mutationProb, mu,lambd, complexPenalty, genSteps

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
	
def genRandBits(individualLength): #makes a random bitstring
	arr = numpy.random.randint(2, size=(int(individualLength),))
	return list(arr) 
	
def loadFatimaData(filename,tempDict): #reads data from differential expression file
	inputfile = open(filename, 'r')
	lines = inputfile.readlines()
	for line in lines:
		if line[0]!='#':
			kline=line.split(' ')
			tempDict[str.upper(kline[0][3:])]=logistic.cdf(float(kline[1]))

			
def loadFpkms(filename): #loads data from fpkms tab delimited csv file
	with open(filename) as csvfile:
		data=[]
		reader = csv.reader(csvfile, delimiter='\t')
		for row in reader:
			data.append(row)
		return data
def sortFpkms(data): #puts fpkms data into appropriate lists of steady state dictionaries following normalization to largest value (as denominator)
	sss=[]
	print(data[1])
	for j in range(1,len(data[1])):
		sss.append({})
	for i in range(1,len(data)):
		maxdata=0
		for j in range(1,len(data[i])):
			if float(data[i][j])>maxdata:
				maxdata=float(data[i][j])
		if maxdata==0:
			maxdata=1
		for j in range(0,len(data[i])-1):
			sss[j][str.upper(data[i][0])]=float(data[i][1])/maxdata
	return sss
def Inv(x, inverter): #inverts if necessary then applies hill fun
	if inverter:
		return (1-x)
	else:
		return (x)
def fuzzyAnd(num1, num2):
	return num1*num2
	#return min(num1,num2)
			
def fuzzyOr(num1, num2):
	return (num1+num2-(num1*num2))
	#return max(num1, num2)

def bit2int(bitlist): #converts bitstring to integer
	out = 0
	for bit in bitlist:
		out = (out << 1) | bit
	return out





def propUpdateNet(currentNode,oldValue,individual, triple, nodeList,inputOrderList, inputOrderInvertList, possibilityNumList):
	inputOrder=inputOrderList[currentNode] # find the list of possible input combinations for the node we are on 
	inputOrderInvert=inputOrderInvertList[currentNode] #find list of lists of whether input nodes need to be inverted (corresponds to inputOrder)
	if possibilityNumList[currentNode]>0:
		logicOperatorFlags=list(individual[triple[1]:triple[2]]) # find demarcations on individual for bits we need
		inputOrder=list(inputOrder[bit2int(individual[triple[0]:triple[1]])%possibilityNumList[currentNode]]) # determine what order of inputs is
		inputOrderInvert=inputOrderInvert[bit2int(individual[triple[0]:triple[1]])%possibilityNumList[currentNode]] # lookup which inputs need to be inverted ('not')
		if len(inputOrder)==0:
			value=oldValue[currentNode] #if no inputs, maintain value
			return value
		elif len(inputOrder)==1:
			#if only one input, then can either affect or not affect the node. so either keep the value or update to the single input's value
			if individual[triple[0]]==1:
				value=Inv(oldValue[inputOrder[0]],inputOrderInvert[0])
			else:
				value=oldValue[currentNode]
			return value
		else:
			#update nodes with more than one input
			# update and then or
			upstreamVals=[]
			for upstream in range(0,len(inputOrder)):
				upstreamVals.append(Inv(oldValue[inputOrder[upstream]],inputOrderInvert[upstream]))
			counter =0
			# print(upstreamVals)
			# print(logicOperatorFlags)
			while counter < len(logicOperatorFlags) and counter+1<len(inputOrder):
				if logicOperatorFlags[counter]==0:
					tempVal=fuzzyAnd(upstreamVals[counter],upstreamVals[counter+1])
					inputOrder.pop(counter)
					inputOrder.pop(counter)
					logicOperatorFlags.pop(counter)
					upstreamVals.pop(counter)
					upstreamVals.pop(counter)
					upstreamVals.insert(counter,tempVal)
				else:
					counter=counter+1
				# print(upstreamVals)

			#first one uses the initial logic operator flag to decide and vs or then combines the first two inputs
			while len(upstreamVals)>1:
				tempVal=fuzzyOr(upstreamVals.pop(0),upstreamVals.pop(0))
				upstreamVals.insert(0,tempVal)
				# print(upstreamVals)
			return upstreamVals[0]
	else:
		#returns savme value if now inputs
		return oldValue[currentNode]						

def propUpdateNode(currentNode,oldValue,individual, individualParse, nodeList,inputOrderList, inputOrderInvertList, possibilityNumList):
	
	oldTriple=individualParse[currentNode]
	triple=[]
	triple.append(0)
	triple.append(oldTriple[1]-oldTriple[0])
	triple.append(oldTriple[2]-oldTriple[0])

	retinue = propUpdateNet(currentNode,oldValue,individual, triple, nodeList,inputOrderList, inputOrderInvertList, possibilityNumList)
	return retinue
	# inputOrder=inputOrderList[currentNode]
	# inputOrderInvert=inputOrderInvertList[currentNode]
	# if possibilityNumList[currentNode]>0:
	# 	logicOperatorFlags=individual[triple[1]:triple[2]]
	# 	inputOrder=inputOrder[bit2int(individual[triple[0]:triple[1]])%possibilityNumList[currentNode]]
	# 	inputOrderInvert=inputOrderInvert[bit2int(individual[triple[0]:triple[1]])%possibilityNumList[currentNode]]
	# 	if len(inputOrder)==0:
	# 		value=oldValue[currentNode]
	# 	elif len(inputOrder)==1:
	# 		if individual[0]==1:
	# 			value=Inv(oldValue[inputOrder[0]],inputOrderInvert[0])
	# 		else:
	# 			value=oldValue[currentNode]
	# 	else:
	# 		counter =0
	# 		if logicOperatorFlags[0]==0:
	# 			value=fuzzyAnd(Inv(oldValue[inputOrder[0]],inputOrderInvert[0],),Inv(oldValue[inputOrder[1]],inputOrderInvert[1],))
	# 		else:
	# 			value=fuzzyOr(Inv(oldValue[inputOrder[0]],inputOrderInvert[0],),Inv(oldValue[inputOrder[1]],inputOrderInvert[1],))
	# 		for i in range(2,len(inputOrder)):
	# 			if logicOperatorFlags[i-1]==0:
	# 				value=fuzzyAnd(value,Inv(oldValue[inputOrder[i]],inputOrderInvert[i],))
	# 			else:
	# 				value=fuzzyOr(value,Inv(oldValue[inputOrder[i]],inputOrderInvert[i],))
	# 	return value
	# else:
	# 	# print('origvalue')
	# 	# print(oldValue[currentNode])
	# 	return oldValue[currentNode]						
def writeFuzzyNode(currentNode,individual, individualParse, nodeList,inputOrderList, inputOrderInvertList, possibilityNumList):
	#write out evaluation instructions in BooleanNet format. This follows the exact same code as fuzzyUpdate, but writes a string instead of actually updating the values of the nodes
	triple=individualParse[currentNode]
	inputOrder=inputOrderList[currentNode]
	inputOrderInvert=inputOrderInvertList[currentNode]
	writenode=''+nodeList[currentNode]+'='
	if possibilityNumList[currentNode]>0:
		logicOperatorFlags=individual[triple[1]:triple[2]]
		inputOrder=inputOrder[bit2int(individual[triple[0]:triple[1]])%possibilityNumList[currentNode]]
		inputOrderInvert=inputOrderInvert[bit2int(individual[triple[0]:triple[1]])%possibilityNumList[currentNode]]
		if len(inputOrder)==1:
			#print(inputOrder[0])
			if individual[triple[0]]==1:
				if inputOrderInvert[0]:
					writenode=writenode+' not '
				writenode=writenode+nodeList[inputOrder[0]]+' '
			else:
				writenode=writenode+' '+ nodeList[currentNode]
		else:
			counter =0
			if inputOrderInvert[0]:
				writenode=writenode+' not '
			writenode=writenode+nodeList[inputOrder[0]]+' '
			if logicOperatorFlags[0]==0:
				writenode=writenode+' and '
			else:
				writenode=writenode+' or '
			if inputOrderInvert[1]:
				writenode=writenode+' not '
			writenode=writenode+nodeList[inputOrder[1]]+' '
			for i in range(2,len(inputOrder)):
				if logicOperatorFlags[i-1]==0:
					writenode=writenode+' and '
				else:
					writenode=writenode+' or '
				if inputOrderInvert[i]:
					writenode=writenode+' not '
				writenode=writenode+nodeList[inputOrder[i]]+' '
	else:
		writenode=writenode+' '+ nodeList[currentNode]
	return writenode					
	
def writeFuzzyModel(individual, individualParse, nodeList,inputOrderList, inputOrderInvertList, possibilityNumList):
	#iterate over nodes to generate a BooleanNet representation for the entire model
	addString=''
	for i in range(0,len(nodeList)):
		addString=addString+writeFuzzyNode(i,individual, individualParse, nodeList,inputOrderList, inputOrderInvertList, possibilityNumList)
		addString=addString+'\n'
	return addString

def runFuzzyModel( individual, async, individualParse, nodeList, inputOrderList, inputOrderInvertList, possibilityNumList, initValues, steps,sigmaNode, sigmaNetwork,earlyEvalNodes):
	#do fuzzy simulation. individual is just a long bitstring... need to seperate it out.  
	#list of triples which gives the total list of nodes. 
	counter=0;
	newValue=list(initValues)
	simData=[]
	simData.append(list(newValue))
	seq=range(0,len(individualParse))
	earlyEvalNodes=[]
	for node in range(0,len(earlyEvalNodes)):
		newValue[earlyEvalNodes[node]]=propUpdateNet(earlyEvalNodes[node],newValue,individual, individualParse, nodeList,inputOrderList, inputOrderInvertList, possibilityNumList)
	for step in range(0,steps):
		oldValue=list(newValue)
		if async:
			shuffle(seq)
		for i in range(0,len(individualParse)):	
			if async:
				temp=propUpdateNet(seq[i],newValue,individual, individualParse[seq[i]], nodeList,inputOrderList, inputOrderInvertList, possibilityNumList)
			else:
				temp=propUpdateNet(seq[i],oldValue,individual, individualParse[seq[i]], nodeList,inputOrderList, inputOrderInvertList, possibilityNumList)

			if sigmaNode==0:
				newValue[seq[i]]=temp
			else:
				newValue[seq[i]]=temp+sigmaNode*(1-temp/(.1+temp)-temp/10)*random()
		simData.append(list(newValue))
	avg= [0 for x in range(0,len(newValue))]
	stdev= [0 for x in range(0,len(newValue))]
	for step in range(0,len(simData)):
		for element in range(0,len(avg)):
			simData[step][element]=simData[step][element]+sigmaNetwork*random()
	for step in range(len(simData)-10,len(simData)):
		for element in range(0,len(avg)):
			avg[element]=avg[element]+simData[step][element]
	for element in range(0,len(avg)):
		avg[element]=avg[element]/10
	#for step in range(len(simData)-10,len(simData)):
		#for element in range(0,len(avg)):
			#stdev[element]=stdev[element]+(simData[step][element]-avg[element])^2
	#return (avg, stdev)
	return avg

def iterateFuzzyModel(async, iters, individual, individualParse, nodeList, inputOrderList, inputOrderInvertList, possibilityNumList, initValues, steps ,sigmaNode, sigmaNetwork,earlyEvalNodes):
	sum=[0 for x in range(0,len(initValues))]
	for i in range(0,iters):
		avg,stdev=runFuzzyModel( individual, async, individualParse, nodeList, inputOrderList, inputOrderInvertList, possibilityNumList, initValues, steps, sigmaNode, sigmaNetwork, earlyEvalNodes)
		for j in range(0,len(sum)):
			sum[j]=sum[j]+avg[j]
	# print(sum)
	avgs=list(sum)
	for i in range(0,len(sum)):
		avgs[i]=sum[i]/float(iters)
	return avgs
		
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

def simplifyNetwork(graph, ss):
	#network simplification algorithm. 
	# # 1. remove nodes which are neither part of input nor have input to them...
	# # 2. remove straigth paths. 
	# # 3. 


	#collapse complexes of nodes already in the list
	for node in graph.nodes():
		predlist=graph.predecessors(node)
		for pred in predlist:
			if '-' in pred:
				# # print(pred)
				genes=pred.split('-')
				flag=True
				for gene in genes:
					if not gene in predlist:
						flag=False
				if flag:
					graph.remove_edge(pred,node)
	flag=True
	
	# remove nodes with no predecessors or value in given steady state data
	while(flag):
		flag=False	
		# # print(len(graph.nodes()))
		newNodes = [x for x in graph.nodes() if not (len(graph.predecessors(x))==0 and (not (x in ss.keys())))]
		if(len(newNodes)<len(graph.nodes())):
			flag=True
		graph=graph.subgraph(newNodes)
		
		print(len(graph.nodes()))
			
	#  collapse straight lines
	removeNodeList= [x for x in graph.nodes() if (len(graph.predecessors(x))==1 and (len(graph.successors(x))==1))]
	for rm in removeNodeList:
		before=graph.predecessors(rm)[0]
		after=graph.successors(rm)[0]
		edge1=graph.get_edge_data(before,rm)['signal']
		edge2=graph.get_edge_data(rm,after)['signal']
		inhCount=0
		if edge1=='i':
			inhCount=inhCount+1
		if edge2=='i':
			inhCount=inhCount+1
		if inhCount==1:
			graph.add_edge(before,after,signal='i')
		else:
			graph.add_edge(before,after,signal='a')
		graph.remove_node(rm)
	flag=True

	#rewire nodes that have only one upstream node
	print(len(graph.nodes()))
	removeNodeList= [x for x in graph.nodes() if (len(graph.predecessors(x))==1) ]
	for rm in removeNodeList:
		for after in graph.successors(rm):
			before=graph.predecessors(rm)[0]
			edge1=graph.get_edge_data(before,rm)['signal']
			edge2=graph.get_edge_data(rm,after)['signal']
			inhCount=0
			if edge1=='i':
				inhCount=inhCount+1
			if edge2=='i':
				inhCount=inhCount+1
			if inhCount==1:
				graph.add_edge(before,after,signal='i')
			else:
				graph.add_edge(before,after,signal='a')
		graph.remove_node(rm)
	flag=True

	#rewire nodes that have only one downstream node
	print(len(graph.nodes()))
	removeNodeList= [x for x in graph.nodes() if (len(graph.successors(x))==1) ]
	for rm in removeNodeList:
		for start in graph.predecessors(rm):
			finish=graph.successors(rm)[0]
			edge1=graph.get_edge_data(start,rm)['signal']
			edge2=graph.get_edge_data(rm,finish)['signal']
			inhCount=0
			if edge1=='i':
				inhCount=inhCount+1
			if edge2=='i':
				inhCount=inhCount+1
			if inhCount==1:
				graph.add_edge(start,finish,signal='i')
			else:
				graph.add_edge(start,finish,signal='a')
		graph.remove_node(rm)
	flag=True
	return graph
	
def buildFatimaNetwork(): #build a network from Fatimas data and from IL1 pathways in KEGG
	dataFileName='C:/Users/Rohith/Desktop/Rohith_data/ruleMaker_GA/Data/fatima_fpkms.csv'
	#dataFileName='/home/rohith/Documents/fatima_fpkms.csv'
	#two dicts for the models
	nodeUpdateDict={}
	data=loadFpkms(dataFileName)
	sss=sortFpkms(data)
	gc.collect()
	#print(ss.keys())
	graph = nx.DiGraph()
	dict={}
	aliasDict={}
	KEGGdict=constructor.parseKEGGdict('ko00001.keg', aliasDict, dict)	
	
	currentfile='IL1b_pathways.txt'
	inputfile = open(currentfile, 'r')
	line = inputfile.read()
	codelist=re.findall('ko\d\d\d\d\d',line)	
	print(codelist)
	constructor.uploadKEGGcodes(codelist, graph, KEGGdict)
	for node in graph.nodes():
		if node in graph.successors(node):
			graph.remove_edge(node,node)
	nodeList=graph.nodes()
	
	simplifyNetwork(graph,sss[0])
	return graph, sss

			
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
			
def runFatimaGA():
	#very simple function to just run GA from fatima data with no frills. just see output from a single call to eaMuCommaLambda
	graph, ss = buildFatimaNetwork()
	iters, simSteps, generations, popSize, bitFlipProb, crossoverProb,mutationProb = returnSimParams()
	individualLength, individualParse, nodeList, inputOrderList, initValueList, evaluateNodes, possibilityNumList, inputOrderInvertList= setupGAparams(graph, ss)
	toolbox, hof, stats=buildToolbox(iters,individualLength, bitFlipProb, ss, evaluateNodes, individualParse, nodeList, inputOrderList, inputOrderInvertList, possibilityNumList, initValueList,simSteps, )
	toolbox.register("evaluate", evaluate, iters=iters,ss=ss,  evaluateNodes=evaluateNodes,individualParse= individualParse, nodeList=nodeList, inputOrderList=inputOrderList, inputOrderInvertList=inputOrderInvertList, possibilityNumList=possibilityNumList, initValueList=initValues, steps=simSteps, h=h,p=p,hillOn=hillOn)
	population=toolbox.population(n=popSize)
	algo.eaMuCommaLambda(population, toolbox, mu=100, lambda_=200, stats=stats, cxpb=crossoverProb, mutpb=mutationProb, ngen=generations, verbose=True, halloffame=hof)


def synthesizeInputs(graph,samples): # generates synthetic random inputs for steady state simulations
	sss=[]
	for i in range(0,samples):
		sss.append({})
	for node in graph.nodes():
		for i in range(0,samples):
			sss[i][node]=random()
	return sss

def runFatimaSim(networkNoises, nodeNoises, trials):
	#simulate data on Fatima network. Output how well it goes. 
	samples=20
	graph, sss=buildFatimaNetwork()
	sss1=synthesizeInputs(graph,samples)
	async, iters, simSteps, generations, popSize, bitFlipProb, crossoverProb ,mutationProb, mu,lambd= returnSimParams()
	individualLength, individualParse, nodeList, inputOrderList, initValueList, evaluateNodes, possibilityNumList, inputOrderInvertList, earlyEvalNodes=setupGAparams(graph, sss1)
	print("individualLength")
	print(individualLength)
	#open file to save, use noises in name
	#save graph seperately
	
	toolbox, stats=buildToolbox(individualLength, bitFlipProb, samples)

	for networkNoise in networkNoises:
		print('networkNoise')
		print(networkNoise)
		for nodeNoise in nodeNoises:
			print('nodeNoise')
			print(nodeNoise)
			fuzzySimulator(async, iters, toolbox, stats, graph, trials, sss1,individualLength,evaluateNodes, individualParse,nodeList,inputOrderList,inputOrderInvertList,possibilityNumList,initValueList,simSteps ,nodeNoise,networkNoise,popSize ,crossoverProb,mutationProb,generations,earlyEvalNodes, mu,lambd)
		
def fuzzySimulator(async, iters, toolbox, stats, graph, trials, sss,individualLength,evaluateNodes, individualParse,nodeList,inputOrderList,inputOrderInvertList,possibilityNumList,initValueList, steps ,nodeNoise,networkNoise,popSize,crossoverProb,mutationProb,generations,earlyEvalNodes, mu,lambd, complexPenalty, genSteps):
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


def setupBootstrapParams(graph):
	repeat=True
	while(repeat):
		repeat=False
		for node in graph.nodes():
			if node in graph.successors(node):
				graph.remove_edge(node,node)
				repeat=True
	
	evaluateNodes=[] 
	individualParse=[] 
	inputOrderList=[] 
	inputOrderInvertList=[]
	possibilityNumList=[]
	nodeList=graph.nodes()
	nodeDict={}
	for i in range(0,len(nodeList)):
		nodeDict[nodeList[i]]=i
	steps = 100
	counter=0
	for i in range(0,len(nodeList)):
		preds=graph.predecessors(nodeList[i])
		if len(preds)>15:
			preds=preds[1:15]
		for j in range(0,len(preds)):
			preds[j]=nodeDict[preds[j]]
		withNones = zip(preds, itertool.repeat('empty'))
		
		possibilities=list(itertool.product(*withNones))
		activities=[]
		for j in range(0,len(possibilities)):
			possibilities[j]=list(possibilities[j])
			while 'empty' in possibilities[j]:
				possibilities[j].remove('empty')
			while [] in possibilities[j]:
				possibilities[j].remove([])
		while [] in possibilities:
			possibilities.remove([])
		for sequence in possibilities:
			activity=[]
			for node in sequence:
			#print(nodeList[node[0]])
			#print(nodeList[i])
			#print(graph.edges())
				
				if graph.edge[nodeList[node]][nodeList[i]]['signal']=='a':
					activity.append(False)
					#print('a')
				else:
					activity.append(True)
					#print('i')
			activities.append(activity)
		inputOrderList.append(possibilities)
		inputOrderInvertList.append(activities)
		possibilityNumList.append(len(possibilities))
		if len(possibilities)==0:
			logNum=0
		else:
			logNum=ceil(log(len(possibilities))/log(2))
		individualParse.append([int(counter),int(counter+logNum),int(counter+logNum+len(possibilities)-1)])
		counter=counter+logNum+len(possibilities)
		# print(counter)
		# print(individualParse)
	return counter, individualParse, nodeList, inputOrderList, evaluateNodes, possibilityNumList, inputOrderInvertList

def testBootstrap(reSampleNum,sampleSize, graph,sss):

	# graph=nx.DiGraph()
	# graph.add_edge('1','3',signal='a')
	# graph.add_edge('2','3',signal='a')

	individualLength, individualParse, nodeList, inputOrderList, initValueList, evaluateNodes, possibilityNumList, inputOrderInvertList, earlyEvalNodes= setupGAparams(graph, sss)
	nodeNoise=0
	networkNoise=0
	async, iters, simSteps, generations, popSize, bitFlipProb, crossoverProb,mutationProb, mu,lambd, complexPenalty, genSteps=returnSimParams()
	#for i in range(0,trials):
	individual=genRandBits(individualLength)
	# print(individualLength)
	# print(individual)
	
	truths=[]
	#print(writeFuzzyModel(individual, individualParse, nodeList,inputOrderList, inputOrderInvertList, possibilityNumList))
	for i in range(0,len(nodeList)):
		truths.append(writeFuzzyNode(i,individual, individualParse, nodeList,inputOrderList, inputOrderInvertList, possibilityNumList))		
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
		#print(evaluate(individual, async, iters, newSSS, evaluateNodes, individualParse, nodeList, inputOrderList, inputOrderInvertList, possibilityNumList, initValueList, steps, ,earlyEvalNodes,complexPenalty))
	
	initValueList=[]
	for j in range(0,len(newSSS)):
		initValueList.append([])
	for i in range(0,len(nodeList)):
		for j in range(0,len(newSSS)):
			ss=newSSS[j]
			if  nodeList[i] in newSSS[0].keys():
				initValueList[j].append(ss[nodeList[i]])
			else:
				initValueList[j].append(0.5)
				
	# print(individual)
	for i in range(0,len(nodeList)):
		# print('new Node')
		# print(individual)
		# print(individualParse[i])
		deviation=0
		# for steadyStateNum in range(0,len(newSSS)):
		# 	derivedVal=propUpdate(i,initValueList[steadyStateNum],individual[individualParse[i][0]:individualParse[i][2]], individualParse, nodeList,inputOrderList, inputOrderInvertList, possibilityNumList, )
		# 	deviation=deviation+(derivedVal-newSSS[steadyStateNum][nodeList[i]])**2
		# print(deviation)

	ruleList=bootstrapRules(graph,newSSS,reSampleNum,sampleSize, individualLength, individualParse, nodeList, inputOrderList, evaluateNodes, possibilityNumList, inputOrderInvertList)
	wrongCount=0
	correctCount=0
	for item in ruleList:
		for i in range(0,len(nodeList)):
			if possibilityNumList[i]>1:
				if truths[i]==writeFuzzyNode(i,item, individualParse, nodeList,inputOrderList, inputOrderInvertList, possibilityNumList):
					correctCount=correctCount+1
				else:
					wrongCount=wrongCount+1
					print(truths[i])
					print(writeFuzzyNode(i,item, individualParse, nodeList,inputOrderList, inputOrderInvertList, possibilityNumList))

	return wrongCount, correctCount
		# print(item)
		# print(writeFuzzyModel(item, individualParse, nodeList,inputOrderList, inputOrderInvertList, possibilityNumList))


def bitList(n):
    return [1 if digit=='1' else 0 for digit in bin(n)[2:]]

def bootstrapRules(graph,sss,reSampleNum,sampleSize, counter, individualParse, nodeList, inputOrderList, evaluateNodes, possibilityNumList, inputOrderInvertList):
	ruleList=[]
	for i in range(0,reSampleNum):
		sss2=numpy.random.choice(sss, size=sampleSize, replace=True, p=None)
		ruleList.append(propBoolModel( individualParse, nodeList, inputOrderList, evaluateNodes, possibilityNumList, inputOrderInvertList,graph,sss2))
	ruleDict=[]

	# for i in range(0,len(nodeList)):
	# 	ruleDict.append({})
	# for trialNum in range(0,len(ruleList)):
	# 	for nodeNum in range(0,len(ruleList[trialNum])):
	# 		ruleDict[nodeNum][ruleList[trialNum][nodeNum]]=0
	# for trialNum in range(0,len(ruleList)):
	# 	for nodeNum in range(0,len(ruleList[trialNum])):
	# 		ruleDict[nodeNum][ruleList[trialNum][nodeNum]]=ruleDict[nodeNum][ruleList[trialNum][nodeNum]]+1
	return ruleList

def propBoolModel(individualParse, nodeList, inputOrderList, evaluateNodes, possibilityNumList, inputOrderInvertList,graph,sss1):
	async, iters, simSteps, generations, popSize, bitFlipProb, crossoverProb,mutationProb, mu,lambd, complexPenalty, genSteps=returnSimParams()
	bestList=[]
	initValueList=[]
	for j in range(0,len(sss1)):
		initValueList.append([])
	for i in range(0,len(nodeList)):
		for j in range(0,len(sss1)):
			ss=sss1[j]
			if  nodeList[i] in sss1[0].keys():
				initValueList[j].append(ss[nodeList[i]])
			else:
				initValueList[j].append(0.5)
	for i in range(0,len(nodeList)):
		print(nodeList[i])
		print(individualParse[i])
		currentDev=10*len(sss1)
		best=[]
		if possibilityNumList[i]>0:
			for j in range(0,int(2**int(individualParse[i][2]-individualParse[i][0]+1))):
				#if possibilityNumList[i]==1:
					# print('single possibility')
					# print(j)
				bits=[]
				bits=bitList(j)
				while len(bits)<(individualParse[i][2]-individualParse[i][0]+1):
					bits.insert(0,0)
				deviation=0
				for steadyStateNum in range(0,len(sss1)):
					derivedVal=propUpdateNode(i,initValueList[steadyStateNum],bits, individualParse, nodeList,inputOrderList, inputOrderInvertList, possibilityNumList, )
					deviation=deviation+(derivedVal-sss1[steadyStateNum][nodeList[i]])**2
				print(bits)
				print(deviation)
				if(deviation<currentDev):
					print("best")
					best=bits
					currentDev=deviation	
		# print(i)
		# print(nodeList[i])
		bestList.append(best)
		# print(best)
	return [item for sublist in bestList for item in sublist]

if __name__ == '__main__':
	# ss={}
	#data=loadFpkms('C:/Users/Rohith/Desktop/ruleMaker_GA/Data/fatima_fpkms.csv')
	#sortFpkms(data,ss)
	# print(ss)
	#teseter.LiuNetwork2SimTest()
	#runFatimaSim([.01,.05],[.01,.02],1)
	#runFatimaSim([0],[0.001],10)
	#tester.LiuNetwork2SimTest(20,1)
	inputNum=50
	counters=[]
	graphStart=tester.LiuNetwork2Builder()
	sss=synthesizeInputs(graphStart,inputNum)
	graph=simplifyNetwork(graphStart,sss[0])
	nx.write_graphml(graph,'Liu2.graphml')
	for i in range(0,1):
		sss=synthesizeInputs(graph,inputNum)
		counters.append(testBootstrap(5,10,graph,sss))
	print(counters)