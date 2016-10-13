#import stuff
from random import random
from random import shuffle
import numpy as numpy
from math import floor, ceil, log
import operator
from deap import base
from deap import creator
from deap import gp
from deap import tools
from deap import algorithms as algo
import networkx as nx
from scipy.stats import logistic
import re
import urllib2 
import fuzzyNetworkConstructor as constructor
import csv
import itertools as itertool

def setupGAparams(graph):
	global individualParse #keep a triplet of how to parse out each individual
	global steps
	global nodeList
	global nodeOrderList
	global initValues
	global ss
	global evaluateNodes #list of indices in nodeList to be evaluated from nodelist
	global possibilityNumList
	global h
	global p
	global hillOn
	global nodeOrderInvertList
	
	repeat=True
	while(repeat):
		repeat=False
		for node in graph.nodes():
			if node in graph.successors(node):
				graph.remove_edge(node,node)
				repeat=True
	hillOn = False
	p=.5
	h= 3
	evaluateNodes=[] 
	individualParse=[] 
	nodeOrderList=[] 
	nodeOrderInvertList=[]
	possibilityNumList=[]
	nodeList=graph.nodes()
	nodeDict={}
	for i in range(0,len(nodeList)):
		nodeDict[nodeList[i]]=i
	steps = 5
	counter=0
	initValues=[]
	for i in range(0,len(nodeList)):
		preds=graph.predecessors(nodeList[i])
		if len(preds)>15:
			preds=preds[1:15]
			#print(preds)
		for j in range(0,len(preds)):
			preds[j]=nodeDict[preds[j]]
		# print(preds)
		# print(nodeList[i])
		withNones = zip(preds, itertool.repeat('empty'))
		
		possibilities=list(itertool.product(*withNones))
		activities=[]
		for j in range(0,len(possibilities)):
			possibilities[j]=list(possibilities[j])
			while 'empty' in possibilities[j]:
				possibilities[j].remove('empty')
			while [] in possibilities[j]:
				possibilities[j].remove([])
			# print('possible')
			# print possibilities[j]
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
		if  nodeList[i] in ss.keys():
			initValues.append(ss[nodeList[i]])
			evaluateNodes.append(i)
		else:
			initValues.append(0.5)
	# print(initValues)
	# print(ss)
	# print(len(ss.keys()))
	global individualLength
	individualLength=counter
	return counter
	
	
def genRandBits():
	global individualLength
	arr = numpy.random.randint(2, size=(individualLength,))
	#print(arr)
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

def sortFpkms(data,ss):
	for i in range(1,len(data)):
		maxdata=0
		for j in range(1,len(data[i])):
			if float(data[i][j])>maxdata:
				maxdata=float(data[i][j])
		if maxdata==0:
			maxdata=1
		ss[str.upper(data[i][0])]=float(data[i][1])/maxdata
	
def hill(x):
	global h
	global p
	global hillOn
	if hillOn:
		return ((1+h**p)*x**p)/(h**p+x**p)		
	else:
		return x
	
	
def hillInv(x, inverter):
	if inverter:
		return hill(1-x)
	else:
		return hill(x)

# def hill(x, h, p):
	# return ((1+h**p)*x**p)/(h**p+x**p)
			
def fuzzyAnd(num1, num2):
	return min(num1,num2)
			
def fuzzyOr(num1, num2):
	return max(num1, num2)

	
	
def bit2int(bitlist):
	out = 0
	for bit in bitlist:
		out = (out << 1) | bit
	return out


def fuzzyUpdate(currentNode,oldValue,individual):
	global individualParse
	global nodeList
	global nodeOrderList
	global nodeOrderInvertList
	global possibilityNumList
	# print(nodeList[currentNode])
	triple=individualParse[currentNode]
	nodeOrder=nodeOrderList[currentNode]
	nodeOrderInvert=nodeOrderInvertList[currentNode]
	#print(nodeOrderList)
	if possibilityNumList[currentNode]>0:
		#print(bit2int(individual[triple[0]:triple[1]])%possibilityNumList[currentNode])
		logicOperatorFlags=individual[triple[1]:triple[2]]
		nodeOrder=nodeOrder[bit2int(individual[triple[0]:triple[1]])%possibilityNumList[currentNode]]
		nodeOrderInvert=nodeOrderInvert[bit2int(individual[triple[0]:triple[1]])%possibilityNumList[currentNode]]
		# print(nodeList[currentNode])
		# print(nodeOrderInvert)
		# print(nodeOrder)
		if len(nodeOrder)==0:
			value=oldValue[currentNode]
		elif len(nodeOrder)==1:
			#print(nodeOrder[0])
			value=hillInv(oldValue[nodeOrder[0]],nodeOrderInvert[0])
		else:
			counter =0
			if logicOperatorFlags[0]==0:
				value=fuzzyAnd(hillInv(oldValue[nodeOrder[0]],nodeOrderInvert[0]),hillInv(oldValue[nodeOrder[1]],nodeOrderInvert[1]))
			else:
				value=fuzzyOr(hillInv(oldValue[nodeOrder[0]],nodeOrderInvert[0]),hillInv(oldValue[nodeOrder[1]],nodeOrderInvert[1]))
			for i in range(2,len(nodeOrder)):
				if logicOperatorFlags[i]==0:
					value=fuzzyAnd(value,hillInv(oldValue[nodeOrder[i]],nodeOrderInvert[i]))
				else:
					value=fuzzyOr(value,hillInv(oldValue[nodeOrder[i]],nodeOrderInvert[i]))
		# print('value')
		# print(value)
		return value
	else:
		# print('origvalue')
		# print(oldValue[currentNode])
		return oldValue[currentNode]						
				
def runFuzzySim(individual):
	#do fuzzy simulation. individual is just a long bitstring... need to seperate it out. 
	global individualParse
	global nodeList
	global initValues
	#list of triples which gives the total list of nodes. 
	counter=0;
	newValue=list(initValues)
	simData=[]
	simData.append(list(newValue))
	global steps
	seq=range(0,len(individualParse))
	for step in range(0,steps):
		shuffle(seq)
		for i in range(0,len(individualParse)):	
			newValue[seq[i]]=fuzzyUpdate(seq[i],newValue,individual)
		simData.append(list(newValue))
	array= [0 for x in range(0,len(newValue))]
	for step in range(len(simData)-5,len(simData)):
		#print(simData[step])
		for element in range(0,len(array)):
			array[element]=array[element]+simData[step][element]
	for element in range(0,len(array)):
		array[element]=array[element]/5
	return array
	
	
			
#this is the evaluation function. Need to create a graph then recreate everything with it.
def evaluate(individual):
	global ss
	global evaluateNodes
	#print evaluateNodes
	RME=1.0
	#for i in range(0,len(individual)):
		# print(individual)
		# print(i)
		# print(individual[int(i)])
	boolValues=runFuzzySim(individual)	
	#print(len(boolValues))
	for i in range(0, len(evaluateNodes)):
		#print(evaluateNodes[i])
		RME=RME+(boolValues[evaluateNodes[i]]-ss[nodeList[evaluateNodes[i]]])**2
		#print(RME)
	return RME,

def simplifyNetwork(graph, ss):
#network simplification algorithm. 
# # 1. remove nodes which are neither part of input nor have input to them...
# # 2. remove straigth paths. 
# # 3. 

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
		
		
# remove nodes with no predecessors or value in given steady state data
	print(len(graph.nodes()))
	newNodes = [x for x in graph.nodes() if not (len(graph.predecessors(x))==0 and (not (x in ss.keys())))]
	graph=graph.subgraph(newNodes)
	
	print(len(graph.nodes()))

	
	
def mutate(individual):
	global nodeList
	
	return individual,

	
	
def buildFatimaNetwork():
	dataFileName='C:/Users/Rohith/Desktop/ruleMaker_GA/Data/fatima_fpkms.csv'
	
	#two dicts for the models
	nodeUpdateDict={}
	global ss
	ss={} 
	data=loadFpkms(dataFileName)
	sortFpkms(data,ss)
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
	
	simplifyNetwork(graph,ss)
	return graph
	
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
			# for predy in graph.predecessors(node):
				# if node in graph.predecessors(predy):
					# print(predy)
	# graph1=graph.subgraph(subgraphnodes)
	# constructor.drawGraph2(graph1,'dense_subgraph.png')
		# if len(graph.predecessors(node))==0:
			# counter=counter+1
			# if node in ss.keys():
				# counter2=counter2+1
	# print('nodes with no incoming edges')
	# print(counter)
	# print('nodes with no incoming edges that are part of transcriptomics data set')
	# print(counter2)
	# print('number of datapoints in Fatima data')
	# print(len(ss.keys()))
	# print('nodes in overlap')
	# print(len(evaluateNodes))



def runFatimaSim():
	graph=buildFatimaNetwork()
	setupGAparams(graph)
	
	# # #setup toolbox
	
	toolbox = base.Toolbox()

	pset = gp.PrimitiveSet("MAIN", arity=1)
	pset.addPrimitive(operator.add, 2)
	pset.addPrimitive(operator.sub, 2)
	pset.addPrimitive(operator.mul, 2)

	# make a fitness minimization function
	creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
	# create a class of individuals that are lists from networkx
	creator.create("Individual", list, fitness=creator.FitnessMin)


	# how to create aliases for your individuals... 
	# toolbox.register("attr_float", random.random)
	# need this alias to create new graphs... but we should just be using the first one.... 

	toolbox.register("genRandomBitString", genRandBits)
	toolbox.register("individual", tools.initIterate, creator.Individual, toolbox.genRandomBitString)
	toolbox.register("population", tools.initRepeat, list , toolbox.individual)
	
	ind1=toolbox.individual()
	population=toolbox.population(n=100)

	
	stats = tools.Statistics(key=lambda ind: ind.fitness.values)
	stats.register("avg", numpy.mean)
	stats.register("std", numpy.std)
	stats.register("min", numpy.min)
	stats.register("max", numpy.max)
	hof = tools.HallOfFame(1, similar=numpy.array_equal)
	
	# finish registering the toolbox functions
	toolbox.register("mate", tools.cxTwoPoint)
	toolbox.register("mutate", tools.mutFlipBit, indpb=.1)
	toolbox.register("select", tools.selNSGA2)
	toolbox.register("evaluate", evaluate)
	toolbox.register("similar", numpy.array_equal)
	algo.eaMuCommaLambda(population, toolbox, mu=100, lambda_=200, stats=stats, cxpb=.2, mutpb=.2, ngen=100, verbose=True)
	
	
	
	# # # graphy=nx.DiGraph()
	# # # graphy.add_edge(1,2,color='blue',type='typing')	

	# # # ind2=population[2]
	# # # child1, child2 = [toolbox.clone(ind) for ind in (ind1, ind2)]
	# # # tools.cxBlend(child1, child2, 0.5)
	# # # del child1.fitness.values
	# # # del child2.fitness.values
	
	
	
def testFuzzySim():
	filename='inputDataFatima.txt'
	KEGGfileName='ko04060.xml'
		
	#two dicts for the models
	nodeUpdateDict={}
	global ss
	ss={}
	for i in range(1,5):
		ss[str(i)]=0
	ss['zero']=1
	ss['one']=0
	ss['two']=0
	print(ss.keys())
	graph = nx.DiGraph()
	dict={}
	aliasDict={}
	KEGGdict=constructor.parseKEGGdict('ko00001.keg', aliasDict, dict)
	graph.add_edge('zero','one', signal='a')
	graph.add_edge('zero','two', signal='a')
	graph.add_edge('one','two', signal='i')	

	global individualLength
	individualLength=setupGAparams(graph)
	# global nodeOrderInvertList
	# global nodeOrderList
	# global nodeList
	# print(nodeList)
	# print(nodeOrderList)
	# print(nodeOrderInvertList)
	# print(individualLength)
	global nodeList
	print(nodeList)
	print(graph.edges())
	individual=[0,0,0,0,0,0]
	print(individual)
	print(runFuzzySim(individual))
	
	individual=[0,1,0,0,0,0]
	print(individual)
	print(runFuzzySim(individual))
	
	individual=[1,0,0,0,0,0]
	print(individual)
	print(runFuzzySim(individual))
	
	
	
	graph = nx.DiGraph()
	dict={}
	aliasDict={}
	KEGGdict=constructor.parseKEGGdict('ko00001.keg', aliasDict, dict)
	graph.add_edge('zero','one', signal='a')
	graph.add_edge('zero','two', signal='i')
	graph.add_edge('one','two', signal='i')	
	individualLength=setupGAparams(graph)
	print(nodeList)
	print(graph.edges())
	individual=[0,0,0,0,0,0]
	print(individual)
	print(runFuzzySim(individual))
	
	individual=[0,1,0,0,0,0]
	print(individual)
	print(runFuzzySim(individual))
	
	individual=[1,0,0,0,0,0]
	print(individual)
	print(runFuzzySim(individual))
	
	

	
	
	ss['zero']=1
	ss['one']=1
	ss['two']=1
	individualLength=setupGAparams(graph)
	print(nodeList)
	print(graph.edges())
	individual=[0,0,0,0,0,0]
	print(individual)
	print(runFuzzySim(individual))
	
	individual=[0,1,0,0,0,0]
	print(individual)
	print(runFuzzySim(individual))
	
	individual=[1,0,0,0,0,0]
	print(individual)
	print(runFuzzySim(individual))
	




if __name__ == '__main__':
	# ss={}
	# data=loadFpkms('C:/Users/Rohith/Desktop/ruleMaker_GA/Data/fatima_fpkms.csv')
	# sortFpkms(data,ss)
	# print(ss)
	runFatimaSim()