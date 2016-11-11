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
#from scipy.stats import logistic
import re
import urllib2 
import fuzzyNetworkConstructor as constructor
import csv
import itertools as itertool
import copy as copy
import cPickle as pickle
import gc as gc
#from scoop import futures

#import pyximport; pyximport.install()
#import evaluator
def returnSimParams():
	iters=100
	simSteps=100
	generations=10
	popSize=10
	mu= 10#individuals selected
	lambd= 20#children produced
	p=.5
	h= 3
	hillOn=False
	bitFlipProb=.1
	crossoverProb=.2
	mutationProb=.2
	async=False
	complexPenalty=False
	genSteps=1000
	return async, iters, simSteps, generations, popSize, h, p, hillOn, bitFlipProb, crossoverProb,mutationProb, mu,lambd, complexPenalty, genSteps

def setupGAparams(graph, sss):
	repeat=True
	while(repeat):
		repeat=False
		for node in graph.nodes():
			if node in graph.successors(node):
				graph.remove_edge(node,node)
				repeat=True
	
	evaluateNodes=[] 
	individualParse=[] 
	nodeOrderList=[] 
	nodeOrderInvertList=[]
	possibilityNumList=[]
	earlyEvalNodes=[]
	nodeList=graph.nodes()
	nodeDict={}
	for i in range(0,len(nodeList)):
		nodeDict[nodeList[i]]=i
	steps = 100
	counter=0
	initValueList=[]
	for j in range(0,len(sss)):
			initValueList.append([])
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
		nodeOrderList.append(possibilities)
		nodeOrderInvertList.append(activities)
		possibilityNumList.append(len(possibilities))
		if len(possibilities)==0:
			logNum=0
		else:
			logNum=ceil(log(len(possibilities))/log(2))
		individualParse.append([int(counter),int(counter+logNum),int(counter+logNum+len(possibilities)-1)])
		counter=counter+logNum+len(possibilities)
		# print(counter)
		# print(individualParse)
		for j in range(0,len(sss)):
			ss=sss[j]
			if  nodeList[i] in sss[0].keys():
				initValueList[j].append(ss[nodeList[i]])
				evaluateNodes.append(i)
			else:
				initValueList[j].append(0.5)
				earlyEvalNodes.append(i)
	#print(initValueList)
	return counter, individualParse, nodeList, nodeOrderList, initValueList, evaluateNodes, possibilityNumList, nodeOrderInvertList, earlyEvalNodes
	
	
def genRandBits(individualLength):
	arr = numpy.random.randint(2, size=(int(individualLength),))
	return list(arr)
	
def loadFatimaData(filename,tempDict):
	inputfile = open(filename, 'r')
	lines = inputfile.readlines()
	for line in lines:
		if line[0]!='#':
			kline=line.split(' ')
			tempDict[str.upper(kline[0][3:])]=logistic.cdf(float(kline[1]))

			
def loadFpkms(filename):
	with open(filename) as csvfile:
		data=[]
		reader = csv.reader(csvfile, delimiter='\t')
		for row in reader:
			data.append(row)
		return data

def sortFpkms(data):
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
	
def hill(x,h,p,hillOn):
	if hillOn:
		return ((1+h**p)*x**p)/(h**p+x**p)		
	else:
		return x
		
def hillInv(x, inverter,h,p,hillOn):
	if inverter:
		return hill(1-x, h,p,hillOn)
	else:
		return hill(x, h,p,hillOn)
		
def fuzzyAnd(num1, num2):
	return min(num1,num2)
			
def fuzzyOr(num1, num2):
	return max(num1, num2)

def bit2int(bitlist):
	out = 0
	for bit in bitlist:
		out = (out << 1) | bit
	return out

def fuzzyUpdate(currentNode,oldValue,individual, individualParse, nodeList,nodeOrderList, nodeOrderInvertList, possibilityNumList, h,p,hillOn):
	triple=individualParse[currentNode]
	nodeOrder=nodeOrderList[currentNode]
	nodeOrderInvert=nodeOrderInvertList[currentNode]
	if possibilityNumList[currentNode]>0:
		logicOperatorFlags=individual[triple[1]:triple[2]]
		nodeOrder=nodeOrder[bit2int(individual[triple[0]:triple[1]])%possibilityNumList[currentNode]]
		nodeOrderInvert=nodeOrderInvert[bit2int(individual[triple[0]:triple[1]])%possibilityNumList[currentNode]]
		if len(nodeOrder)==0:
			value=oldValue[currentNode]
		elif len(nodeOrder)==1:
			#print(nodeOrder[0])
			value=hillInv(oldValue[nodeOrder[0]],nodeOrderInvert[0],h,p,hillOn)
		else:
			counter =0
			if logicOperatorFlags[0]==0:
				value=fuzzyAnd(hillInv(oldValue[nodeOrder[0]],nodeOrderInvert[0],h,p,hillOn),hillInv(oldValue[nodeOrder[1]],nodeOrderInvert[1],h,p,hillOn))
			else:
				value=fuzzyOr(hillInv(oldValue[nodeOrder[0]],nodeOrderInvert[0],h,p,hillOn),hillInv(oldValue[nodeOrder[1]],nodeOrderInvert[1],h,p,hillOn))
			for i in range(2,len(nodeOrder)):
				if logicOperatorFlags[i]==0:
					value=fuzzyAnd(value,hillInv(oldValue[nodeOrder[i]],nodeOrderInvert[i],h,p,hillOn))
				else:
					value=fuzzyOr(value,hillInv(oldValue[nodeOrder[i]],nodeOrderInvert[i],h,p,hillOn))
		return value
	else:
		# print('origvalue')
		# print(oldValue[currentNode])
		return oldValue[currentNode]						

def writeFuzzyNode(currentNode,individual, individualParse, nodeList,nodeOrderList, nodeOrderInvertList, possibilityNumList):
	triple=individualParse[currentNode]
	nodeOrder=nodeOrderList[currentNode]
	nodeOrderInvert=nodeOrderInvertList[currentNode]
	writenode=''+nodeList[currentNode]+'='
	if possibilityNumList[currentNode]>0:
		logicOperatorFlags=individual[triple[1]:triple[2]]
		nodeOrder=nodeOrder[bit2int(individual[triple[0]:triple[1]])%possibilityNumList[currentNode]]
		nodeOrderInvert=nodeOrderInvert[bit2int(individual[triple[0]:triple[1]])%possibilityNumList[currentNode]]
		if len(nodeOrder)==1:
			#print(nodeOrder[0])
			if nodeOrderInvert[0]:
				writenode=writenode+' not '
			writenode=writenode+nodeList[nodeOrder[0]]+' '
		else:
			# print(nodeOrder)
			# print(triple)
			# print(individual)
			# print(logicOperatorFlags)
			counter =0
			if nodeOrderInvert[0]:
				writenode=writenode+' not '
			writenode=writenode+nodeList[nodeOrder[0]]+' '
			if logicOperatorFlags[0]==0:
				writenode=writenode+' and '
			else:
				writenode=writenode+' or '
			if nodeOrderInvert[1]:
				writenode=writenode+' not '
			writenode=writenode+nodeList[nodeOrder[1]]+' '
			for i in range(2,len(nodeOrder)):
				if logicOperatorFlags[i-1]==0:
					writenode=writenode+' and '
				else:
					writenode=writenode+' or '
				if nodeOrderInvert[i]:
					writenode=writenode+' not '
				writenode=writenode+nodeList[nodeOrder[i]]+' '
		
	else:
		# print('origvalue')
		# print(oldValue[currentNode])
		writenode=writenode+' '+ nodeList[currentNode]
	return writenode					
	


def writeFuzzyModel(individual, individualParse, nodeList,nodeOrderList, nodeOrderInvertList, possibilityNumList):
	addString=''
	for i in range(0,len(nodeList)):
		addString=addString+writeFuzzyNode(i,individual, individualParse, nodeList,nodeOrderList, nodeOrderInvertList, possibilityNumList)
		addString=addString+'\n'
	return addString

def runFuzzyModel( individual, async, individualParse, nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValues, steps, h,p,hillOn,sigmaNode, sigmaNetwork,earlyEvalNodes):
	#do fuzzy simulation. individual is just a long bitstring... need to seperate it out.  
	#list of triples which gives the total list of nodes. 
	counter=0;
	newValue=list(initValues)
	simData=[]
	simData.append(list(newValue))
	seq=range(0,len(individualParse))
	earlyEvalNodes=[]
	for node in range(0,len(earlyEvalNodes)):
		newValue[earlyEvalNodes[node]]=fuzzyUpdate(earlyEvalNodes[node],newValue,individual, individualParse, nodeList,nodeOrderList, nodeOrderInvertList, possibilityNumList, h,p,hillOn)
	for step in range(0,steps):
		if async:
			shuffle(seq)
		for i in range(0,len(individualParse)):	
			#print(fuzzyUpdate(seq[i],newValue,individual, individualParse, nodeList,nodeOrderList, nodeOrderInvertList, possibilityNumList, h,p,hillOn))
			temp=fuzzyUpdate(seq[i],newValue,individual, individualParse, nodeList,nodeOrderList, nodeOrderInvertList, possibilityNumList, h,p,hillOn)
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

def iterateFuzzyModel(async, iters, individual, individualParse, nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValues, steps, h,p,hillOn,sigmaNode, sigmaNetwork,earlyEvalNodes):
	sum=[0 for x in range(0,len(initValues))]
	for i in range(0,iters):
		avg,stdev=runFuzzyModel( individual, async, individualParse, nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValues, steps, h,p,hillOn,sigmaNode, sigmaNetwork, earlyEvalNodes)
		for j in range(0,len(sum)):
			sum[j]=sum[j]+avg[j]
	# print(sum)
	avgs=list(sum)
	for i in range(0,len(sum)):
		avgs[i]=sum[i]/float(iters)
	return avgs
		
#this is the evaluation function. Need to create a graph then recreate everything with it.
def evaluate(individual, async, iters, sss, evaluateNodes, individualParse, nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList, steps, h,p,hillOn,earlyEvalNodes,complexPenalty):
	SSEs=[]
	for j in range(0,len(sss)):
		ss=sss[j]
		initValues=initValueList[j]
		SSE=0
		if async:
			boolValues=iterateFuzzyModel(async, iters, individual, individualParse, nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList[j], steps, h,p,hillOn,0,0,earlyEvalNodes)	
		else:
			boolValues=runFuzzyModel(individual, async, individualParse, nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList,initValueList[j], steps, h,p,hillOn,0,0,earlyEvalNodes)	
		#if len(evaluateNodes)<10:
		#	print("so few nodes")
		for i in range(0, len(evaluateNodes)):
			SSE=SSE+(boolValues[evaluateNodes[i]]-ss[nodeList[evaluateNodes[i]]])**2
		SSEs.append(SSE)
	edgeDegree=0
	if complexPenalty:
		for i in range(0,len(nodeList)):
			if possibilityNumList[i]>0:
				edgeDegree=edgeDegree+len(nodeOrderList[i][bit2int(individual[individualParse[i][0]:individualParse[i][1]])%possibilityNumList[i]])
	SSEs.append(edgeDegree/1.)
	gc.collect()
	return SSEs

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
	print(len(graph.nodes()))
	removeNodeList= [x for x in graph.nodes() if (len(graph.predecessors(x))==1 and (len(graph.successors(x))==1)and (not (x in ss.keys())))]
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
	removeNodeList= [x for x in graph.nodes() if (len(graph.predecessors(x))==1)and (not (x in ss.keys()))]
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
	
def buildFatimaNetwork():
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
	
def printSubgraphs(graph):
	nodeList=graph.nodes()
	subgraphnodes=set()
	for node in nodeList:
		if len(graph.predecessors(node))>15:
			subgraphnodes.add(node)
			for element in graph.predecessors(node):
				subgraphnodes.add(element)
			
			print('new node')
			print(node)
			print(len(graph.predecessors(node)))
			print(graph.predecessors(node))

			
def buildToolbox( individualLength, bitFlipProb, samples):
	# # #setup toolbox
	toolbox = base.Toolbox()

	pset = gp.PrimitiveSet("MAIN", arity=1)
	pset.addPrimitive(operator.add, 2)
	pset.addPrimitive(operator.sub, 2)
	pset.addPrimitive(operator.mul, 2)


	weightTup=(-1.0,-1.0)
	if samples>1:
		for i in range(1,samples):
			weightTup=weightTup+(-1.0,)
	# make a fitness minimization function
	creator.create("FitnessMin", base.Fitness, weights=weightTup)
	# create a class of individuals that are lists from networkx
	creator.create("Individual", list, fitness=creator.FitnessMin)


	# how to create aliases for your individuals... 
	# toolbox.register("attr_float", random.random)
	# need this alias to create new graphs... but we should just be using the first one.... 

	toolbox.register("genRandomBitString", genRandBits, individualLength=individualLength)
	toolbox.register("individual", tools.initIterate, creator.Individual, toolbox.genRandomBitString)
	toolbox.register("population", tools.initRepeat, list , toolbox.individual)
	
	# ind1=toolbox.individual()
	# population=toolbox.population(n=50)

	
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
	#algo.eaMuCommaLambda(population, toolbox, mu=100, lambda_=200, stats=stats, cxpb=.2, mutpb=.2, ngen=1, verbose=True)
	
	
	
	# # # graphy=nx.DiGraph()
	# # # graphy.add_edge(1,2,color='blue',type='typing')	
	return toolbox, stats
			
def runFatimaGA():
	graph, ss = buildFatimaNetwork()
	iters, simSteps, generations, popSize, h, p, hillOn, bitFlipProb, crossoverProb,mutationProb = returnSimParams()
	individualLength, individualParse, nodeList, nodeOrderList, initValueList, evaluateNodes, possibilityNumList, nodeOrderInvertList= setupGAparams(graph, ss)
	toolbox, hof, stats=buildToolbox(iters,individualLength, bitFlipProb, ss, evaluateNodes, individualParse, nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList,simSteps,h,p,hillOn)
	toolbox.register("evaluate", evaluate, iters=iters,ss=ss,  evaluateNodes=evaluateNodes,individualParse= individualParse, nodeList=nodeList, nodeOrderList=nodeOrderList, nodeOrderInvertList=nodeOrderInvertList, possibilityNumList=possibilityNumList, initValueList=initValues, steps=simSteps, h=h,p=p,hillOn=hillOn)
	population=toolbox.population(n=popSize)
	algo.eaMuCommaLambda(population, toolbox, mu=100, lambda_=200, stats=stats, cxpb=crossoverProb, mutpb=mutationProb, ngen=generations, verbose=True, halloffame=hof)

def randomizeInputs(sss,samples):
	newSS=[]
	for i in range(0,samples):
		newSS.append({})
	for key in sss[0].keys():
		for i in range(0,samples):
			newSS[i][key]=random()
	return newSS

def synthesizeInputs(graph,samples):
	sss=[]
	for i in range(0,samples):
		sss.append({})
	for node in graph.nodes():
		for i in range(0,samples):
			sss[i][node]=random()
	return sss

def runFatimaSim(networkNoises, nodeNoises, trials):
	samples=20
	graph, sss=buildFatimaNetwork()
	#sss1=randomizeInputs(sss,samples)
	sss1=synthesizeInputs(graph,samples)
	async, iters, simSteps, generations, popSize, h, p, hillOn, bitFlipProb, crossoverProb ,mutationProb, mu,lambd= returnSimParams()
	individualLength, individualParse, nodeList, nodeOrderList, initValueList, evaluateNodes, possibilityNumList, nodeOrderInvertList, earlyEvalNodes=setupGAparams(graph, sss1)
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
			fuzzySimulator(async, iters, toolbox, stats, graph, trials, sss1,individualLength,evaluateNodes, individualParse,nodeList,nodeOrderList,nodeOrderInvertList,possibilityNumList,initValueList,simSteps, h,p,hillOn,nodeNoise,networkNoise,popSize ,crossoverProb,mutationProb,generations,earlyEvalNodes, mu,lambd)
		
def fuzzySimulator(async, iters, toolbox, stats, graph, trials, sss,individualLength,evaluateNodes, individualParse,nodeList,nodeOrderList,nodeOrderInvertList,possibilityNumList,initValueList, steps, h,p,hillOn,nodeNoise,networkNoise,popSize,crossoverProb,mutationProb,generations,earlyEvalNodes, mu,lambd, complexPenalty, genSteps):
	hofSize=10
	avgSSEs=[]
	hofs=[]
	truthIndividuals=[]
	truthAvgs=[]
	hofScores=[]
	newSSS=[]
	truthCounter=[0 for number in xrange(hofSize)]
	for i in range(0,trials):
		individual=genRandBits(individualLength)
		newSSS=[]
		for k in range(0,len(sss)):
			if async:
				boolValues=iterateFuzzyModel(async, iters,individual, individualParse, nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList[k], genSteps, h,p,hillOn,nodeNoise,networkNoise,earlyEvalNodes)
			else:
				boolValues=runFuzzyModel(individual, async,individualParse, nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList[k], genSteps, h,p,hillOn,nodeNoise,networkNoise,earlyEvalNodes)
			newSS=copy.deepcopy(sss[k])
			for j in range(0,len(evaluateNodes)):
				newSS[nodeList[evaluateNodes[j]]]=boolValues[evaluateNodes[j]]
			newSSS.append(newSS)
		#update eval nodes
		newInitValueList=[]
		for j in range(0,len(sss)):
			newInitValueList.append([])
		for i in range(0,len(nodeList)):
			for j in range(0,len(sss)):
				ss=newSSS[j]
				if  nodeList[i] in sss[0].keys():
					newInitValueList[j].append(ss[nodeList[i]])
				else:
					newInitValueList[j].append(0.5)
		#print(len(sss))
		toolbox.register("evaluate", evaluate, iters=iters, async=async, sss=newSSS, evaluateNodes=evaluateNodes,individualParse= individualParse, nodeList=nodeList, nodeOrderList=nodeOrderList, nodeOrderInvertList=nodeOrderInvertList, possibilityNumList=possibilityNumList, initValueList=newInitValueList, steps=steps, h=h,p=p,hillOn=hillOn,earlyEvalNodes=earlyEvalNodes, complexPenalty=complexPenalty)
		#toolbox.register("evaluate", evaluate, ss=ss, evaluateNodes=evaluateNodes,individualParse= individualParse, nodeList=nodeList, nodeOrderList=nodeOrderList, nodeOrderInvertList=nodeOrderInvertList, possibilityNumList=possibilityNumList, initValues=boolValues, steps=steps, h=h,p=p,hillOn=hillOn)
		population=toolbox.population(n=popSize)
		hof = tools.HallOfFame(hofSize, similar=numpy.array_equal)
		output = tools.Logbook()
		output=algo.eaMuCommaLambda(population, toolbox, mu=mu, lambda_=lambd, stats=stats, cxpb=crossoverProb, mutpb=mutationProb, ngen=generations, verbose=False, halloffame=hof)
		truth=(writeFuzzyModel(individual, individualParse, nodeList,nodeOrderList, nodeOrderInvertList, possibilityNumList))
		trueNodeList=[]
		for i in range(0,10):
			bestRun=(writeFuzzyModel(hof[i], individualParse, nodeList,nodeOrderList, nodeOrderInvertList, possibilityNumList))
			if truth==bestRun:
				truthCounter[i]+=1
				break
			elif i==0:
				truthlines=truth.split('\n')
				newlines=bestRun.split('\n')
				trueNodes=0
				for j in range(0,len(truthlines)):
					if truthlines[j]==newlines[j]:
						trueNodes=trueNodes+1
				trueNodeList.append(trueNodes)
		# else:
			# print('WRONG')
			# print(truth)
			# print(bestRun)
		avgs=[output[1][i]['min'] for i in range(0,len(output[1]))]
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

def LiuNetwork1Builder():
	graph = nx.DiGraph()
	
	graph.add_edge('a','c', signal='a')
	graph.add_edge('b','d', signal='a')
	graph.add_edge('c','f', signal='a')	
	graph.add_edge('c','h', signal='a')
	graph.add_edge('d','f', signal='a')	
	graph.add_edge('d','g', signal='a')
	graph.add_edge('g','k', signal='a')	
	graph.add_edge('h','j', signal='a')
	graph.add_edge('j','c', signal='i')	
	graph.add_edge('f','k', signal='i')
	return graph
def LiuSSS1Builder():
	sss=[]
	ss={}
	ss['a']=1
	ss['b']=0
	ss['c']=0
	ss['d']=0
	ss['f']=0
	ss['g']=0
	ss['h']=0
	ss['j']=0
	ss['k']=0
	sss.append(ss)
	ss={}
	ss['a']=0
	ss['b']=1
	ss['c']=0
	ss['d']=0
	ss['f']=0
	ss['g']=0
	ss['h']=0
	ss['j']=0
	ss['k']=0
	sss.append(ss)
	ss={}
	ss['a']=0
	ss['b']=0
	ss['c']=1
	ss['d']=0
	ss['f']=0
	ss['g']=0
	ss['h']=0
	ss['j']=0
	ss['k']=0
	sss.append(ss)
	ss={}
	ss['a']=0
	ss['b']=0
	ss['c']=0
	ss['d']=1
	ss['f']=0
	ss['g']=0
	ss['h']=0
	ss['j']=0
	ss['k']=0
	sss.append(ss)
	ss={}
	ss['a']=0
	ss['b']=0
	ss['c']=0
	ss['d']=0
	ss['f']=1
	ss['g']=0
	ss['h']=0
	ss['j']=0
	ss['k']=0
	# sss.append(ss)
	# ss={}
	# ss['a']=0
	# ss['b']=0
	# ss['c']=0
	# ss['d']=0
	# ss['f']=0
	# ss['g']=1
	# ss['h']=0
	# ss['j']=0
	# ss['k']=0
	# sss.append(ss)
	# ss={}
	# ss['a']=0
	# ss['b']=0
	# ss['c']=0
	# ss['d']=0
	# ss['f']=0
	# ss['g']=0
	# ss['h']=1
	# ss['j']=0
	# ss['k']=0
	# sss.append(ss)
	# ss={}
	# ss['a']=0
	# ss['b']=0
	# ss['c']=0
	# ss['d']=0
	# ss['f']=0
	# ss['g']=0
	# ss['h']=0
	# ss['j']=1
	# ss['k']=0
	# sss.append(ss)
	# ss={}
	# ss['a']=0
	# ss['b']=0
	# ss['c']=0
	# ss['d']=0
	# ss['f']=0
	# ss['g']=0
	# ss['h']=0
	# ss['j']=0
	# ss['k']=1
	# sss.append(ss)
	# ss={}
	# ss['a']=1
	# ss['b']=1
	# ss['c']=1
	# ss['d']=1
	# ss['f']=1
	# ss['g']=1
	# ss['h']=1
	# ss['j']=1
	# ss['k']=1
	# sss.append(ss)
	# ss={}
	# ss['a']=0
	# ss['b']=0
	# ss['c']=0
	# ss['d']=0
	# ss['f']=0
	# ss['g']=0
	# ss['h']=0
	# ss['j']=0
	# ss['k']=0
	# sss.append(ss)
	return sss
def LiuNetwork3Builder():
	graph = nx.DiGraph()
	
	graph.add_edge('EGFR','RAS', signal='a')
	graph.add_edge('EGFR','PI3K', signal='a')
	graph.add_edge('TNFR','FADD', signal='a')	
	graph.add_edge('TNFR','PI3K', signal='a')
	graph.add_edge('TNFR','RIP1', signal='a')	
	graph.add_edge('DNA_DAMAGE','ATM', signal='a')
	graph.add_edge('RAS','ERK', signal='a')	
	graph.add_edge('RAS','JNK', signal='a')
	graph.add_edge('FADD','CASP8', signal='a')	
	graph.add_edge('PI3K','AKT', signal='a')
	graph.add_edge('RIP1','P38', signal='a')
	graph.add_edge('RIP1','JNK', signal='a')
	graph.add_edge('ATM','CHK', signal='a')	
	graph.add_edge('ATM','P53', signal='a')	

	graph.add_edge('ERK','RAS', signal='a')
	graph.add_edge('ERK','MYC', signal='a')	
	graph.add_edge('ERK','BID', signal='i')
	graph.add_edge('CASP8','BID', signal='a')	
	graph.add_edge('CASP8','CASP3', signal='a')
	graph.add_edge('CASP8','RIP1', signal='i')
	graph.add_edge('JNK','P53', signal='a')	
	graph.add_edge('JNK','FOXO1', signal='a')
	graph.add_edge('AKT','BID', signal='i')	
	graph.add_edge('AKT','CASP8', signal='i')
	graph.add_edge('AKT','BIM', signal='i')
	graph.add_edge('AKT','4EBP1', signal='i')	
	graph.add_edge('AKT','S6K', signal='a')
	graph.add_edge('AKT','P53', signal='i')	
	graph.add_edge('AKT','FOXO1', signal='i')
	graph.add_edge('P38','CDC25', signal='a')	
	graph.add_edge('P38','FOXO1', signal='a')
	graph.add_edge('P38','P53', signal='a')
	graph.add_edge('CHK','CDC25', signal='i')	
	graph.add_edge('CHK','P53', signal='a')	
	graph.add_edge('CDC25','CDK', signal='a')
	graph.add_edge('P53','PUMA', signal='a')
	graph.add_edge('P53','FADD', signal='a')
	graph.add_edge('P53','ATM', signal='a')	
	graph.add_edge('P53','CHK', signal='i')
	graph.add_edge('P53','CYCLIN', signal='i')	
	graph.add_edge('P53','BIM', signal='a')
	graph.add_edge('FOXO1','P27', signal='a')
	graph.add_edge('BIM','BAX', signal='a')	
	graph.add_edge('BID','BAX', signal='a')
	graph.add_edge('MYC','PROLIFERATION', signal='a')	
	graph.add_edge('PUMA','BAX', signal='a')
	graph.add_edge('S6K','S6', signal='a')	
	graph.add_edge('S6K','PI3K', signal='a')
	graph.add_edge('P27','CYCLIN', signal='i')
	graph.add_edge('P27','CELL_CYCLE', signal='a')	





	graph.add_edge('d','g', signal='a')
	graph.add_edge('g','k', signal='a')	
	graph.add_edge('h','j', signal='a')
	graph.add_edge('j','c', signal='a')	
	graph.add_edge('f','k', signal='a')
	graph.add_edge('c','h', signal='a')
	graph.add_edge('d','f', signal='a')	
	graph.add_edge('d','g', signal='a')
	graph.add_edge('g','k', signal='a')	
	graph.add_edge('h','j', signal='a')
	graph.add_edge('j','c', signal='a')	
	graph.add_edge('f','k', signal='a')
	graph.add_edge('c','h', signal='a')
	graph.add_edge('d','f', signal='a')	
	graph.add_edge('d','g', signal='a')
	graph.add_edge('g','k', signal='a')	
	graph.add_edge('h','j', signal='a')
	graph.add_edge('j','c', signal='a')	
	graph.add_edge('f','k', signal='a')
	graph.add_edge('d','f', signal='a')	
	graph.add_edge('d','g', signal='a')
	graph.add_edge('g','k', signal='a')	
	graph.add_edge('h','j', signal='a')
	graph.add_edge('j','c', signal='i')	
	graph.add_edge('f','k', signal='i')
	graph.add_edge('c','h', signal='a')
	graph.add_edge('d','f', signal='a')	
	graph.add_edge('d','g', signal='a')
	graph.add_edge('g','k', signal='a')	
	graph.add_edge('h','j', signal='a')
	graph.add_edge('j','c', signal='a')	
	graph.add_edge('f','k', signal='a')
	graph.add_edge('c','h', signal='a')
	graph.add_edge('d','f', signal='a')	
	graph.add_edge('d','g', signal='a')
	graph.add_edge('g','k', signal='a')	
	graph.add_edge('h','j', signal='a')
	graph.add_edge('j','c', signal='a')	
	graph.add_edge('f','k', signal='a')
	graph.add_edge('c','h', signal='a')
	graph.add_edge('d','f', signal='a')	
	graph.add_edge('d','g', signal='a')
	graph.add_edge('g','k', signal='a')	
	graph.add_edge('h','j', signal='a')
	graph.add_edge('j','c', signal='a')	
	graph.add_edge('f','k', signal='a')
	graph.add_edge('c','h', signal='a')
	graph.add_edge('d','f', signal='a')	
	graph.add_edge('d','g', signal='a')
	graph.add_edge('g','k', signal='a')	
	graph.add_edge('h','j', signal='a')
	graph.add_edge('j','c', signal='a')	
	graph.add_edge('f','k', signal='a')

	return graph
def LiuNetwork1SimTest():
	#two dicts for the models
	#sss=LiuSSS1Builder()
	graph = LiuNetwork1Builder()
	#runNoiseSimTest(graph, sss, [.001,.002,.005,.01,.02,.05,.1,.2,.5], [0], 10)
	sss1=synthesizeInputs(graph,5)
	runNoiseSimTest(graph, sss1, [.001], [0], 500)

def LiuNetwork2SimTest():
	#two dicts for the models
	sss, graph=LiuNetwork2Builder()
	runNoiseSimTest(graph, sss, [.01,.02,.05,.1,.2], [.01,.02,.05,.1,.2], 25)

def LiuNetwork2Builder():
	graph = nx.DiGraph()
	graph.add_edge('igf1','ras', signal='a')
	graph.add_edge('tgfa','ras', signal='a')
	graph.add_edge('igf1','pi3k', signal='a')
	graph.add_edge('tgfa','pi3k', signal='a')
	graph.add_edge('ras','pi3k', signal='a')
	graph.add_edge('ras','map3k1', signal='a')
	graph.add_edge('ras','mek12', signal='a')
	graph.add_edge('tnfa','pi3k', signal='a')
	graph.add_edge('tnfa','jnk12', signal='a')
	graph.add_edge('tnfa','map3k1', signal='a')
	graph.add_edge('tnfa','map3k7', signal='a')
	graph.add_edge('tnfa','mkk4', signal='a')
	graph.add_edge('il1a','map3k7', signal='a')
	graph.add_edge('il1a','map3k1', signal='a')
	graph.add_edge('map3k7','ikk', signal='a')
	graph.add_edge('map3k7','mkk4', signal='a')
	graph.add_edge('map3k7','p38', signal='a')
	graph.add_edge('map3k7','hsp27', signal='a')
	graph.add_edge('map3k1','ikk', signal='a')
	graph.add_edge('map3k1','jnk12', signal='a')
	graph.add_edge('map3k1','mkk4', signal='a')
	graph.add_edge('pi3k','map3k1', signal='a')
	graph.add_edge('pi3k','akt', signal='a')
	graph.add_edge('pi3k','mek12', signal='a')
	graph.add_edge('akt','mek12', signal='a')
	graph.add_edge('akt','ikk', signal='a')
	graph.add_edge('mek12','erk12', signal='a')
	graph.add_edge('ikk','ikb', signal='a')
	graph.add_edge('mkk4','jnk12', signal='a')
	graph.add_edge('mkk4','p38', signal='a')
	graph.add_edge('erk12','hsp27', signal='a')
	graph.add_edge('p38','hsp27', signal='a')

	listy=['igf1','tgfa','ras','tnfa','il1a','map3k7','pi3k','map3k1','akt','mek12','ikk','mkk4','erk12','jk12','ikb','p38','hsp27']
	sss=[]
	for i in range(0,10):
		sss.append({})
	for element in listy:
		for i in range(0,10):
			sss[i][element]=random()
	return sss, graph

def runNoiseSimTest(graph, sss, networkNoises, nodeNoises, trials):
	async, iters, steps, generations, popSize, h, p, hillOn, bitFlipProb, crossoverProb,mutationProb,mu,lambd, complexPenalty, genSteps= returnSimParams()
	individualLength, individualParse, nodeList, nodeOrderList, initValueList, evaluateNodes, possibilityNumList, nodeOrderInvertList, earlyEvalNodes=setupGAparams(graph, sss)
	#open file to save, use noises in name
	#save graph seperately
	toolbox, stats=buildToolbox(individualLength, bitFlipProb, len(sss))
	for networkNoise in networkNoises:
		print('networkNoise')
		print(networkNoise)
		for nodeNoise in nodeNoises:
			print('nodeNoise')
			print(nodeNoise)
			fuzzySimulator(async, iters, toolbox, stats, graph, trials, sss,individualLength,evaluateNodes, individualParse,nodeList,nodeOrderList,nodeOrderInvertList,possibilityNumList,initValueList, steps, h,p,hillOn,nodeNoise,networkNoise,popSize,crossoverProb,mutationProb,generations,earlyEvalNodes, mu,lambd, complexPenalty, genSteps)




if __name__ == '__main__':
	# ss={}
	#data=loadFpkms('C:/Users/Rohith/Desktop/ruleMaker_GA/Data/fatima_fpkms.csv')
	#sortFpkms(data,ss)
	# print(ss)
	#LiuNetwork2SimTest()
	#runFatimaSim([.01,.05],[.01,.02],1)
	#runFatimaSim([0],[0.001],10)
	LiuNetwork1SimTest()